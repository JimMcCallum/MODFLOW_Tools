# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 09:30:11 2024

@author: 00098687
"""

import numpy as np
import scipy as sp
import pylab as plt
import flopy
from modflowapi import ModflowApi

class Model:
    def __init__(self,
                 simname,
                 modelname,
                 tmodelname,
                 t,
                 z,
                 Kx, #probably irrelevent
                 Kz,
                 alphaL,
                 alphaT,
                 alphaTv,
                 Ss,
                 Sy,
                 theta_sat,
                 Water_Table,
                 theta_r,
                 eps,
                 AEP,
                 theta_i,
                 finf,
                 pet, 
                 extdp, 
                 extwc, 
                 hroot, 
                 rootact):
    
        
        #initiale simultation
        self.sim = flopy.mf6.MFSimulation(sim_name=simname, 
                             exe_name="mf6",
                             version="mf6",) 
        
        #set up tims disctretisation
        self.simname = simname
        self.fmod = modelname
        self.tmod = tmodelname
        dt = t[1:] - t[:-1]
        nper = len(dt)
        self.nper = nper
        perioddata = []
        for i in range(nper):
            perioddata.append((dt[i],1,1.0))
        self.tdis = flopy.mf6.ModflowTdis(self.sim, time_units='DAYS', nper=nper, perioddata=perioddata)
        #initiale GWF
        #modelname="gwf-tidal"
        self.gwf = flopy.mf6.ModflowGwf(self.sim, modelname = modelname, model_nam_file = f'{modelname}.nam',
                                   save_flows = True,
                                   newtonoptions='UNDER_RELAXATION') 
        
        #initiate GWT
        if type(tmodelname) != type(None):
            self.gwt = flopy.mf6.ModflowGwt(self.sim, modelname=tmodelname)
         
        #Flow solver
        self.ims_gwf = ims= flopy.mf6.ModflowIms(self.sim, pname="ims",complexity="moderate",
                          filename=f"{simname}.ims",
                          print_option = "all",
                          outer_maximum=(15000)
                          )
        self.sim.register_ims_package(self.ims_gwf, [self.gwf.name])
        
        #Transport Solver
        
        if type(tmodelname) != type(None):
            self.ims_gwt = flopy.mf6.ModflowIms(self.sim,
                                   complexity = 'complex',
                                   linear_acceleration='BICGSTAB',
                                   print_option = "all",
                                   filename="{}.ims".format(self.gwt.name))
            self.sim.register_ims_package(self.ims_gwt, [self.gwt.name])
            
        #Dis packages
        vertices = [[0,0.,1.],
                    [1,1.,1.],
                    [2,0.,0.],
                    [3,1.,0.]]
        cell2d = [[0,0.5,0.5,4,0,1,3,2]]
        top = z[0]
        self.top = top
        botm = z[1:]
        ncpl =len(cell2d)
        nvert = len(vertices)
        nlay = len(botm)
        self.disv = flopy.mf6.ModflowGwfdisv(self.gwf, 
                                length_units='meters', 
                                nlay=nlay, 
                                ncpl=ncpl, 
                                nvert=nvert, 
                                top=top, 
                                botm=botm, 
                                idomain=1, 
                                vertices=vertices, 
                                cell2d=cell2d,
                                filename=f"{modelname}.disv")
        if type(tmodelname) != type(None):
            self.disvt = flopy.mf6.ModflowGwfdisv(self.gwt, 
                                    length_units='meters', 
                                    nlay=nlay, 
                                    ncpl=ncpl, 
                                    nvert=nvert, 
                                    top=top, 
                                    botm=botm, 
                                    idomain=1, 
                                    vertices=vertices, 
                                    cell2d=cell2d,
                                    filename=f"{modelname}_t.disv")


        #NPF, advection and dispersion Packages
        self.npf = flopy.mf6.ModflowGwfnpf(self.gwf, 
                              icelltype=1, 
                              k= Kx, #horizontal conductivity 
                              k33=Kz,#vertical conductivity
                              save_flows=True,
                              save_specific_discharge = True,)
        if type(tmodelname) != type(None):
            self.adv = flopy.mf6.ModflowGwtadv(self.gwt, scheme='UPSTREAM')
        
            ####Dispersion #####
            self.dsp = flopy.mf6.ModflowGwtdsp(self.gwt,diffc=0.0001,alh=alphaL,alv=alphaTv,ath1 = alphaT)
        
        #STO packages
        self.sto = flopy.mf6.ModflowGwfsto(self.gwf,
                                           ss = Ss, #Specific Storage
                                           sy = Sy, #Specific Yield
                                           iconvert = 1)
        if type(tmodelname) != type(None):
            self.tmst = flopy.mf6.ModflowGwtmst(self.gwt,porosity=theta_sat)

        #Drain
        # The top, which will be modified in API
        tide_GHB_SPD = [((0,0),top,10000.)]
        self.tide_DRN = flopy.mf6.ModflowGwfdrn(self.gwf,  
                                                stress_period_data = {0 :tide_GHB_SPD},
                                                pname="DRN_1",
                                                save_flows = True)
        
        #GHB base...
        #then the one for the water table
        WT_GHB_SPD = [((nlay-1,0),Water_Table,1.)]
        self.WT_GHB = flopy.mf6.ModflowGwfghb(self.gwf,  
                                                stress_period_data = {0 : WT_GHB_SPD},
                                                pname="GHB_0",
                                                save_flows = True)
        
        #UZF packages
        uztpackagedata = []
        packagedata = []
        """n = np.argmin(abs(botm-Water_Table))
        if botm[n] > Water_Table:
            n+=1"""
        n = nlay-1
        uzfno = 0
        icpl = 0
        for i in range(n+1):
            lflag = 0
            surfdep = 0.01
            if i < n:
                ivertcon = uzfno+1
            else:
                ivertcon = -1
            if i == 0:
                lflag = 1
            row = [uzfno,
                   (i,icpl),
                   lflag,
                   ivertcon,
                   surfdep,
                   Kz,
                   theta_r,
                   theta_sat,
                   theta_i,
                   eps]
             
            packagedata.append(row)
            uztpackagedata.append((uzfno,35.))
            uzfno +=1
        
        pd0 = [(0, finf, pet, extdp, extwc, AEP, hroot, rootact)]
        self.uzf = flopy.mf6.ModflowGwfuzf(self.gwf, 
                                          packagedata = packagedata,
                                          pname="uzf",
                                          simulate_et = True,
                                          simulate_gwseep=False,
                                          unsat_etwc = True,
                                          print_input=True,
                                          print_flows=True,
                                          save_flows=True,
                                          mover=False,
                                          perioddata={0: pd0},
                                          ntrailwaves=7, 
                                          nwavesets=1000,
                                          nuzfcells=len(packagedata),
                                          budget_filerecord=f"{modelname}.uzf.bud",
                                          wc_filerecord = f"{modelname}.wcf" ,
                                          filename=f"{modelname}.uzf")        
        #initial conditions
        if type(tmodelname) != type(None):
            pd0 = [(0,"infiltration",35.),
                   (0,"uzet",0.)]
            uztperioddata = {0: pd0}
            
            self.uzt = flopy.mf6.modflow.ModflowGwtuzt(self.gwt,
                                                       save_flows=True,
                                                       print_input=True,
                                                       print_flows=True,
                                                       print_concentration=True,
                                                       concentration_filerecord= tmodelname + ".uzt.bin",
                                                       budget_filerecord= tmodelname + ".uzt.bud",
                                                       packagedata=uztpackagedata,
                                                       uztperioddata=uztperioddata,
                                                       pname="uzf",
                                                       filename=f"{tmodelname}.uzt")
            
        
        
        #Initial conditions
        self.icf=flopy.mf6.ModflowGwfic(self.gwf, pname="ic", strt=Water_Table)
        if type(tmodelname) != type(None):
            self.ict=flopy.mf6.ModflowGwtic(self.gwt, pname="ic", strt=35.)
        
            #Flow transport coupling
            
            #sourcerecarray = [("GHB_1", "AUX", "SALT"),]
            self.ssm = flopy.mf6.ModflowGwtssm(self.gwt)#,sources = sourcerecarray)"""
            
            self.gwfgwt = flopy.mf6.ModflowGwfgwt(self.sim, 
                                                 exgtype='GWF6-GWT6',
                                                 exgmnamea=self.gwf.name, 
                                                 exgmnameb=self.gwt.name)
            
            #Oservations
        
        #output control packages
        head_filerecord = f"{simname}.hds"
        budget_filerecord = f"{simname}.cbc"
        self.foc =flopy.mf6.ModflowGwfoc(self.gwf,
                                         head_filerecord=head_filerecord,
                                         budget_filerecord=budget_filerecord,
                                         saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")]
                                         )
        if type(tmodelname) != type(None):
            self.toc = flopy.mf6.ModflowGwtoc(self.gwt,
                                              concentration_filerecord='SALT.ucn',
                                              saverecord=[('CONCENTRATION', 'last')])    
            
    def add_burrow(self,):
    
        print('to do')
        
    def model_make(self):
        self.sim.write_simulation()
    
    
    def model_run(self,t,h,c = None):
        mf6_config_file =  './mfsim.nam'
        mf6 = ModflowApi('libmf6.dll')
        mf6.initialize()  
        #td = np.loadtxt('../tides.csv', delimiter = ',')

        for i in range(self.nper):
            
            t_c = mf6.get_current_time()
            n = np.argmin((t_c - t)**2.)
            sw = h[n]
            
            
            address = ["sinf", self.fmod,"uzf"]
            tag = mf6.get_var_address(*address)
            bound = mf6.get_value_ptr(tag)
            if sw > self.top:
                bound[0] = 1000.0      
            else:
                bound[0] = 0.
            mf6.set_value(tag,bound)
            #print(sw,self.top,bound)
            
            address = ["elev",self.fmod,"drn_1"]
            tag = mf6.get_var_address(*address)
            if sw > self.top:
                bound[0] = sw
                #bound[0,1] = 10000.
            else:
                bound[0] = self.top
                #bound[0,1] = 10000.                
            mf6.set_value(tag,bound)
            
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            #Flow solution
            address = ["MXITER", "SLN_1"]
            mxittag = mf6.get_var_address(*address)
            mxit = mf6.get_value_ptr(mxittag)
            kiter = 0
            mf6.prepare_solve(1)

            while kiter < mxit:
                has_converged = mf6.solve(1)
                kiter += 1
                if has_converged:
                    break
            if not has_converged:
                return None
            mf6.finalize_solve(1)
            if type(self.tmod) != type(None):
                #Transport solution
                address = ["MXITER", "SLN_2"]
                mxittag = mf6.get_var_address(*address)
                mxit = mf6.get_value_ptr(mxittag)
                kiter = 0
                mf6.prepare_solve(2)
    
                while kiter < mxit:
                    has_converged = mf6.solve(2)
                    kiter += 1
                    if has_converged:
                        break
                if not has_converged:
                    return None
                mf6.finalize_solve(2)
            
            mf6.finalize_time_step()
        
        
        mf6.finalize()
        hds = self.gwf.output.head()
        self.h= hds.get_alldata()
        cbc = self.gwf.output.budget()
        ghb = cbc.get_data(text="GHB")
        if type(self.tmod) != type(None):
            conc = self.gwt.output.concentration()
            self.c = conc.get_alldata()