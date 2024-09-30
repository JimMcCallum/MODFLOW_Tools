#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Functions for LOOP2Flopy
#------------------
# You will also need following notebooks. Need to clone this repo and change path to point towards the notebooks. 
# https://github.com/JimMcCallum/MODFLOW_Tools
#----------------
# modelling_routines.ipynb
# geomodel_routines.ipynb
# meshing_routines.ipynb
# disv2disu.py


# In[2]:


class Project:
    
    def __init__(self, 
                 projectname, # string
                 boundingbox, # for geo model: tuple(x0, x1, y0, y1, z0, z1)
                 ):
        
        self.projectname = projectname
        self.boundingbox = boundingbox
        self.x0, self.x1 = boundingbox[0], boundingbox[1]
        self.y0, self.y1 = boundingbox[2], boundingbox[3]
        self.z0, self.z1 = boundingbox[4], boundingbox[5]
        self.Lx = self.x1 - self.x0
        self.Ly = self.y1 - self.y0
        self.Lz = self.z1 - self.z0


# In[1]:


class Model:
    
    def __init__(self,
                 modelname,
                 P,                # Project parameters
                 plan,             # plan (Car - cartesian, Tri - triangular, Vor - voronoi), 
                 transect,         # transect (Vox - voxel, Con - conformable)
                 ):     
    
        #---------------- CREATE PLAN DISCRETISATION ARRAYS ------------------------------------#
        
        self.modelname = modelname             # e.g. 'REF', 'SS', 'SU'...
        self.plan = plan             # Car - cartesian, Tri - triangular, Vor - voronoi
        self.transect = transect     # Vox - voxel, Con - conformable 
        
        # Retrieve these details from meshing
        if plan == 'car':
            self.cell2d   = P.cell2dcar
            self.ncpl     = len(self.cell2d)
            self.vertices = P.verticescar
            self.xcyc     = P.xcyccar
            
        if plan == 'tri':
            self.cell2d   = P.cell2dtri
            self.ncpl     = len(self.cell2d)
            self.vertices = P.verticestri
            self.xcyc     = P.xcyctri
            
        if plan == 'vor':
            self.cell2d   = P.cell2dvor
            self.ncpl     = len(self.cell2d)
            self.vertices = P.verticesvor
            self.xcyc     = P.xcycvor
            
#---------- FUNCTION TO EVALUATE GEO MODEL AND POPULATE HYDRAULIC PARAMETERS ------#
# LOOP2MF take a flexible grid in plan and regular grid in transect and evaluates geo model    

    def create_lith_dis_arrays(self, P): # Takes the project parameters and model class. 
        
        print('   Creating lithology and discretisation arrays for ', self.modelname, ' ...')
        t0 = datetime.now()
        
        def reshape_loop2mf(array):
            array = array.reshape((self.nlay, self.ncpl))
            array = np.flip(array, 0)
            return(array)

#---------- VOX - DIS ARRAY ------#

        if self.transect == 'vox':
            
            self.nlay = P.nlv
            self.dz = (P.z1 - P.z0) / P.nlv
            self.ncell3d = self.ncpl * self.nlay
            self.idomain = np.ones((self.nlay, self.ncpl)) 
            self.top = P.z1 * np.ones((self.ncpl), dtype=float)
            
            self.zc = np.arange(P.z0 + self.dz / 2, P.z1, self.dz)  # Cell centres
            self.zbot = np.arange(P.z1 - self.dz, P.z0 - self.dz, -self.dz)
            
            self.botm = np.zeros((self.nlay, self.ncpl)) 
            for lay in range(self.nlay):
                self.botm[lay,:] = self.zbot[lay]

            #----- VOX - LITH AND VF ------#

            xyz = []                         
            for k in range(self.nlay):
                z = self.zc[k]
                for i in range(self.ncpl):    
                    x, y = self.xcyc[i][0], self.xcyc[i][1]
                    xyz.append([x,y,z])
                    #xyz=np.array(xyz)
            
            litho = P.geomodel.evaluate_model(xyz)  # generates an array indicating lithology for every cell
            vf = P.geomodel.evaluate_model_gradient(xyz) # generates an array indicating gradient for every cell
            
            # Reshape to lay, ncpl   
            litho = np.asarray(litho)
            litho = litho.reshape((self.nlay, self.ncpl))
            litho = np.flip(litho, 0)
            self.lith = litho
            
            angle1, angle2 = [], []
            for i in range(len(vf)):  
                angle1.append(find_angle1(vf[i]))
                angle2.append(find_angle2(vf[i]))
            self.angle1  = reshape_loop2mf(np.asarray(angle1))
            self.angle2  = reshape_loop2mf(np.asarray(angle2))
            
#---------- CON - EVALUATE LITH ON FINE REGULAR LAYERS ------#

        if self.transect == 'con': # CREATING DIS AND NPF ARRAYS
            
            res = 2 # 2m resolution
            nlay = int((P.z1 - P.z0)/res)
            dz = (P.z1 - P.z0)/nlay # actual resolution
            zc = np.arange(P.z0 + dz / 2, P.z1, dz)  # Cell centres
            
            xyz = []  
            for k in range(nlay):
                z = zc[k]
                for i in range(self.ncpl):    
                    x, y = self.xcyc[i][0], self.xcyc[i][1]
                    xyz.append([x,y,z])
                    #xyz=np.array(xyz)
            
            litho = P.geomodel.evaluate_model(xyz)  # generates an array indicating lithology for every cell
            litho = np.asarray(litho)
            litho = litho.reshape((nlay, self.ncpl)) # Reshape to lay, ncpl
            litho = np.flip(litho, 0)
            
            #----- CON - CREATE LITH, BOTM AND IDOMAIN ARRAYS (PILLAR METHOD, PICKS UP PINCHED OUT LAYERS) ------#

            def start_stop_arr(initial_list): # Function to look down pillar and pick geo bottoms
                a = np.asarray(initial_list)
                mask = np.concatenate(([True], a[1:] != a[:-1], [True]))
                idx = np.flatnonzero(mask)
                l = np.diff(idx)
                start = np.repeat(idx[:-1], l)
                stop = np.repeat(idx[1:]-1, l)
                return(start, stop)
            
            # Arrays for geo arrays
            top         = P.z1 * np.ones((self.ncpl), dtype=float)
            botm_geo    = np.zeros((P.nlg, self.ncpl), dtype=float) # bottom elevation of each geological layer
            thick_geo   = np.zeros((P.nlg, self.ncpl), dtype=float) # geo layer thickness
            sublays     = np.zeros((P.nlg, self.ncpl), dtype=float) # number of sublayers
            dz_sublays  = np.zeros((P.nlg, self.ncpl), dtype=float) # geo layer thickness
            idomain_geo = np.zeros((nlay, self.ncpl), dtype=float)      # idomain array for each lithology
            
                             
            for icpl in range(self.ncpl): 
                #get strat column
                strat_log = litho[:,icpl]
                present = np.unique(strat_log)
                start, stop =  start_stop_arr(strat_log)
                start = np.unique(start)
                stop = np.unique(stop)
                for i, lith in enumerate(present):           
                    idomain_geo[lith, icpl] = 1
                    botm_geo[lith, icpl] = P.z1 - (stop[i]+1) * dz
                for lay_geo in range(P.nlg):
                    if idomain_geo[lay_geo, icpl] == 0: # if pinched out geological layer...
                        botm_geo[lay_geo, icpl] = botm_geo[lay_geo-1, icpl]  
                            
            for lay_geo in range(P.nlg):
                    if lay_geo == 0:
                        thick_geo[lay_geo, :] = top - botm_geo[lay_geo,:]
                    else:
                        thick_geo[lay_geo, :] = botm_geo[lay_geo-1,:] - botm_geo[lay_geo,:]
                        
            ####### AMENDED FOR MAXIMUM CELL THICKNESSES ################################

            for lay_geo in range(P.nlg):
                for icpl in range(self.ncpl):
                    min_lay_thick = P.min_modlay_thick[lay_geo]
                    if thick_geo[lay_geo, icpl]/2 > min_lay_thick:
                        sublays[lay_geo, icpl] = math.ceil(thick_geo[lay_geo, icpl]/ min_lay_thick) # geo layer has a minimum of 2 model layers per geo layer
                    else: 
                        sublays[lay_geo, icpl]= 2 # geo layer has a minimum of 2 model layers per geo layer
                    dz_sublays[lay_geo, icpl] = thick_geo[lay_geo, icpl] / sublays[lay_geo, icpl]
                        
            max_sublays = np.ones((P.nlg),  dtype=int)
            for lay_geo in range(P.nlg):
                max_sublays[lay_geo] = sublays[lay_geo, :].max()
            nlay = max_sublays.sum()
            #print('number of DISV model layers = ', nlay)       
            
            # Arrays for flow model
            botm        = np.zeros((nlay, self.ncpl), dtype=float) # bottom elevation of each model layer
            lith        = np.zeros((nlay, self.ncpl), dtype=float) # bottom elevation of each model layer
            idomain     = np.ones((nlay, self.ncpl), dtype=int)    # idomain for each model layer
            
            # Here we make bottom arrays - pinched out cells have the same bottom as layer above
            for icpl in range(self.ncpl):
                lay = 0 # model layer
                for lay_geo in range(P.nlg):
                    #if icpl == 500: print('GEO LAY = ', lay_geo)
                    nsublay    = sublays[lay_geo, icpl]
                    dz         = thick_geo[lay_geo, icpl] / nsublay
                    max_sublay = max_sublays[lay_geo]
                    for s in range(max_sublay): # marches through each sublayer of geo layer
                        if s < nsublay: # active cell
                            if lay == 0:
                                #if icpl == 500: print('Top layer, lay = ', lay)
                                #if icpl == 500: print(top[icpl] - dz)
                                botm[lay, icpl] = top[icpl] - dz
                                lith[lay, icpl] = lay_geo
                            else:
                                #if icpl == 500: print('Not top layer, lay = ', lay)
                                #if icpl == 500: print(botm[lay-1, icpl] - dz)
                                botm[lay, icpl] = botm[lay-1, icpl] - dz
                                lith[lay, icpl] = lay_geo
 
                        else:  # pinched out cell
                            #if icpl == 500: print('PINCHOUT, lay = ', lay)
                            botm[lay, icpl] = botm[lay-1, icpl] # use the bottom before it
                            lith[lay, icpl] = lay_geo
                        lay += 1 # increase mode layer by 1
                        
            # Now we make idomain so that pinched out cells have an idomain of -1
            for icpl in range(self.ncpl):
                if botm[0, icpl] == top[icpl]:
                    idomain[lay, icpl] = -1
                for lay in range(1,nlay):
                    if botm[lay, icpl] == botm[lay-1, icpl]:
                        idomain[lay, icpl] = -1
                  
            self.thick_geo = thick_geo
            self.botm = botm
            self.top = top
            self.idomain = idomain
            self.lith = lith
            self.nlay = nlay
            self.botm_geo = botm_geo
            
#            print('\nmax_sublays =',max_sublays)
#            
#            print('\nicpl = 500')
#            print('top = ', top[500])
#            print('botm_geo = ', botm_geo[:,500])
#            print('thick_geo = ', thick_geo[:,500])
#            print('botm = ', botm[:,500])
#            print('lith = ', lith[:,500])
#            print('idomain = ', idomain[:,500])
#            print('sublays = ', sublays[:, 500])
#            print('nlay (disv) = ', nlay)
                                
            # Get gradients by reevaluationg vector field at each cell
            xyz = []                         
            for lay in range(self.nlay-1, -1, -1):
                for icpl in range(self.ncpl):  
                    x, y, z = self.xcyc[icpl][0], self.xcyc[icpl][1], self.botm[lay, icpl] 
                    xyz.append([x,y,z])
                    #xyz=np.array(xyz)
            vf = P.geomodel.evaluate_model_gradient(xyz) # generates an array indicating gradient for every cell
            
            angle1, angle2 = [], []
            for i in range(len(vf)):  
                angle1.append(find_angle1(vf[i]))
                angle2.append(find_angle2(vf[i]))
            
            self.angle1  = reshape_loop2mf(np.asarray(angle1))
            self.angle2  = reshape_loop2mf(np.asarray(angle2))
            
        self.nnodes_div = len(self.botm.flatten())   
        t1 = datetime.now()
        run_time = t1 - t0
        #print('Time taken = ', run_time.total_seconds())
        
# ------------------ MAKE MODEL GRID CLASS-----------------#
        self.vgrid = flopy.discretization.VertexGrid(vertices=self.vertices,
                               cell2d=self.cell2d,
                               top = self.top,
                               idomain=self.idomain,
                               botm=self.botm,
                               nlay=self.nlay, #number of layers in model
                               ncpl=self.ncpl,)
        
        print(self.angle2.shape)
        # looking into steep dips
        extdip = np.zeros((self.ncpl)) # 2d array to identify pillars with steep dips
        extdip_points = []             # list to append xcyc of points in high dip area
        for icpl in range(self.ncpl):
            for lay in range(nlay):
                if abs(self.angle2[lay, icpl]) >15:
                    extdip[icpl] = 1
                    extdip_points.append((Point(self.xcyc[icpl][0], self.xcyc[icpl][1])))
        extdip_gs = gpd.GeoSeries(extdip_points, crs=P.crs)
        P.extdip = extdip
        P.extdip_points = extdip_points
        P.extdip_gs = extdip_gs
    
    ################## PROP ARRAYS TO BE SAVED IN DISU FORMAT ##################        
    def create_prop_arrays(self, P): # Uses lithology codes to populate arrays 
        
        print('   Creating property arrays for ', self.modelname, ' ...')
        t0 = datetime.now()
        
        # First create an array for cellids in layered version  (before we pop cells that are absent)
        self.cellid_disv = np.empty_like(self.lith, dtype = int)
        self.cellid_disu = -1 * np.ones_like(self.lith, dtype = int)
        i = 0
        for lay in range(self.nlay):
            for icpl in range(self.ncpl):
                self.cellid_disv[lay, icpl] = lay * self.ncpl + icpl
                if self.idomain[lay, icpl] != -1:
                    self.cellid_disu[lay, icpl] = i
                    i += 1
        self.ncell_disv = len(self.cellid_disv.flatten())
        self.ncell_disu = np.count_nonzero(self.cellid_disu != -1)
        
#---------- PROP ARRAYS (VOX and CON) -----   
        self.k11    = np.empty_like(self.lith, dtype = float)
        self.k22    = np.empty_like(self.lith, dtype = float)
        self.k33    = np.empty_like(self.lith, dtype = float)
        self.ss     = np.empty_like(self.lith, dtype = float)
        self.sy     = np.empty_like(self.lith, dtype = float)

        for n in range(P.nlg):  # replace lithologies with parameters
            self.k11[self.lith==n] = P.hk[n] 
            self.k22[self.lith==n] = P.hk[n] 
            self.k33[self.lith==n] = P.vk[n] 
            self.ss[self.lith==n]  = P.ss[n]
            self.sy[self.lith==n]  = P.sy[n]
                   
        # Force all K tensor angles in fault zone to 0 (Loop can't calculate angles in faulted area properly yet!)
        if 'P.fault_poly' in globals(): #if hassattr(P,"fault_poly"):
            for icpl in range(self.ncpl):
                point = Point(self.xcyc[icpl])
                if P.fault_poly.contains(point):
                    for lay in range(self.nlay):
                        self.angle1[lay,icpl] = 0  
                        self.angle2[lay,icpl] = 0   
        ######################################
        
        self.k11    = self.k11[self.cellid_disu != -1].flatten()
        self.k22    = self.k22[self.cellid_disu != -1].flatten()
        self.k33    = self.k33[self.cellid_disu != -1].flatten()
        self.ss     = self.ss[self.cellid_disu != -1].flatten()
        self.sy     = self.sy[self.cellid_disu != -1].flatten()
        self.angle1 = self.angle1[self.cellid_disu != -1].flatten()
        self.angle2 = self.angle2[self.cellid_disu != -1].flatten()
        self.angle3 = np.zeros_like(self.angle1, dtype = float)  # Angle 3 always at 0
        #print('len(k11) =', len(self.k11))
        #print('ncell_disv = ', self.ncell_disv)
        #print('ncell_disu = ', self.ncell_disu)
        #print(np.unique(self.k11))
        self.logk11    = logfunc(self.k11)
        self.logk22    = logfunc(self.k22)
        self.logk33    = logfunc(self.k33)
        
        t1 = datetime.now()
        run_time = t1 - t0
        #print('Time taken = ', run_time.total_seconds())

################## FLOW PACKAGES TO BE SAVED IN DISU FORMAT ##################    
    def create_flow_package_arrays(self, P, rch = True, chd = True, obs = True, wel = True): # e.g. SS.add_flow_packages(recharge = True, chd = True)        
            
        t0 = datetime.now()
        Car = flopy.discretization.VertexGrid(vertices = self.vertices, cell2d = self.cell2d, top = self.top)
        
        print('   Adding flow packages to ', self.modelname, ' ...')
        
        self.strt = P.strt * np.ones_like(self.k11, dtype = float)
        
        if rch:
            self.rch_rec = [] # RECHARGE
            for icpl in range(self.ncpl): 
                self.rch_rec.append([(icpl,), P.rch])
                
        if chd:
            self.vgrid = flopy.discretization.VertexGrid(ncpl = self.ncpl, vertices = self.vertices, cell2d = self.cell2d,
                                                     nlay = self.nlay, botm = self.botm, top = self.top)
            self.gi = flopy.utils.GridIntersect(self.vgrid)
            
            west_bd = LineString([(P.x0, P.y0), (P.x0, P.y1)]) # Western edge
            west_cells = self.gi.intersects(west_bd, shapetype="linestring")
            west_cells = west_cells.cellids.tolist()
                     
            east_bd = LineString([(P.x1, P.y0), (P.x1, P.y1)]) # Eastern edge
            east_cells = self.gi.intersects(east_bd, shapetype="linestring")
            east_cells = east_cells.cellids.tolist()
            
            self.chd_rec = [] # HEAD BOUNDARY
            for icpl in west_cells:
                x,y = self.xcyc[icpl][0], self.xcyc[icpl][1]
                for lay in range(self.nlay):
                    z = lay_to_z(self.botm, self.top, lay, icpl=icpl)
                    zb = self.botm[lay,icpl]
                    if zb < P.chfunc(x,z):                           
                        cell_disv = icpl + lay*(self.ncpl) #8/5
                        cell_disu = self.cellid_disu.flatten()[cell_disv]
                        if cell_disu != -1:
                            self.chd_rec.append([cell_disu, P.chfunc(x,z)]) 
                            
            for icpl in east_cells:
                x,y = self.xcyc[icpl][0], self.xcyc[icpl][1]
                for lay in range(self.nlay):
                    z = lay_to_z(self.botm, self.top, lay, icpl=icpl)
                    zb = self.botm[lay,icpl]
                    if zb < P.chfunc(x,z):                           
                        cell_disv = icpl + lay*(self.ncpl)
                        cell_disu = self.cellid_disu.flatten()[cell_disv]
                        if cell_disu != -1:
                            self.chd_rec.append([cell_disu, P.chfunc(x,z)]) 
        if obs:        
            self.obslist = []
            
            def find_lay_pillar(pillar,z):
                lay = 0  # find the layer z is in          
                while lay < self.nlay-1:
                    if pillar[lay] > z:
                        lay += 1
                        continue
                    else:
                         break
                return (lay)   
            
            for i in range(P.nobs): # Plan
                for j in range(P.nzobs):    # Depth
                    x, y = P.xyobsbores[i][0], P.xyobsbores[i][1] # Get coords of obs point
                    z = P.zobs[i][j]
                    
                    point = Point(x,y) 
                    self.gi = flopy.utils.GridIntersect(self.vgrid)
                    icpl = self.gi.intersect(point)["cellids"]
                    icpl = np.array(list(icpl)) 
                    icpl = icpl[0]

                    pillar = self.botm[:,icpl]
                    lay = find_lay_pillar(pillar,z)  
                    cell_disv = icpl + lay*self.ncpl 
                    cell_disu = self.cellid_disu.flatten()[cell_disv]
                    self.obslist.append([str(P.idobsbores[i] + '_' + str(z)),'head',(cell_disu+1)])

        if wel:
            self.spd_wel_past, self.spd_wel_future = [], [] 

            for n in range(len(P.xypumpbores)):
                x, y = P.xypumpbores[n][0], P.xypumpbores[n][1] # Get coords of obs point
                point = Point(x,y) 
                self.gi = flopy.utils.GridIntersect(self.vgrid)
                icpl = self.gi.intersect(point)["cellids"]
                icpl = np.array(list(icpl)) 
                icpl = icpl[0]
                #print('############ ', icpl)
                #print(self.botm_geo[:,icpl])
                #print(self.botm[:,icpl])
                
                #wel_icpl, wel_icplcoords = find_cell_disv(wel_coords[0], wel_coords[1], self.xcyc)
                wel_top, wel_bot = P.wel_screens[n][0], P.wel_screens[n][1]     
                #print('wel_top and wel_bot = ', wel_top, wel_bot)
                
                if self.transect == 'vox':
                    nwell_cells = int((wel_top - wel_bot)/self.dz)
                    for lay in range(int((0-wel_top)/self.dz), int((0-wel_top)/self.dz) + nwell_cells):   
                        cell_disv = icpl + lay*self.ncpl
                        cell_disu = self.cellid_disu.flatten()[cell_disv]
                        self.spd_wel_past.append([cell_disu, P.qwell_past/nwell_cells])
                        self.spd_wel_future.append([cell_disu, P.qwell_future/nwell_cells])

                if self.transect == 'con':       
                    lay = 0
                    well_layers = []
                    nwell_cells = 0
                    
                    while self.botm[lay, icpl] >= wel_top-0.1: # above top of screen
                        lay += 1
                    while self.botm[lay, icpl] > wel_bot: # above bottom of screen
                        if self.idomain[lay, icpl] != -1: # skips pinched out cells
                            nwell_cells += 1
                            well_layers.append(lay)
                        lay += 1
                    #print(nwell_cells)
                   # print(well_layers)
                    
                    for lay in well_layers:
                        cell_disv = icpl + lay*self.ncpl
                        cell_disu = self.cellid_disu.flatten()[cell_disv]
                        #print('------', icpl, cell_disv, cell_disu)
                        self.spd_wel_past.append([cell_disu, P.qwell_past/nwell_cells])
                        self.spd_wel_future.append([cell_disu, P.qwell_future/nwell_cells])
                        
        t1 = datetime.now()
        run_time = t1 - t0
        print('   Time taken = ', run_time.total_seconds())
                
    def write_run_model(self, P, period, ic_array, staggered = True, complexity='Complex', outer_dvclose=1e-2, inner_dvclose=1e-3, 
                                   outer_maximum=200, newtonoptions = ['UNDER_RELAXATION']):
        
        print('   Writing simulation and gwf for ', self.modelname, ' ...')
        #print(self.modelname)
        t0 = datetime.now()
        
        # -------------- SIM -------------------------
        sim = flopy.mf6.MFSimulation(sim_name='sim', version='mf6',exe_name=P.mfexe_name, sim_ws=P.workspace)

        # -------------- TDIS -------------------------
        if period == 'Steady': tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim)           
        if period == 'Past':   tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, nper=len(P.tdis_past), 
                                                                           perioddata=P.tdis_past)
        if period == 'Future': tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, nper=len(P.tdis_future), 
                                                                           perioddata=P.tdis_future)
        
        # -------------- IMS -------------------------
        # Make linear solver (inner) an order of magnitude tighter than non-linear solver (outer)
        if period == 'Steady': 
            ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', 
                                       complexity    = 'Moderate',
                                       outer_dvclose = 1e-2, 
                                       inner_dvclose = 1e-3, 
                                       outer_maximum = 400, 
                                       linear_acceleration = "BICGSTAB",
                                       preconditioner_levels=5, #1 to 5... PLAY WITH THIS FOR SPEED UP!
                                       preconditioner_drop_tolerance=0.01, # ...if fill 7-18 (hard), DT 1e-2 (7) to 1e-5 (18)
                                       number_orthogonalizations=2,
                                      )
        if period == 'Past' or period == 'Future': 
            ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', 
                                       complexity    = 'Moderate',
                                       outer_dvclose = 1e-2, 
                                       inner_dvclose = 1e-3, 
                                       outer_maximum = 60, 
                                       linear_acceleration = "BICGSTAB",
                                       preconditioner_levels=5, #1 to 5... PLAY WITH THIS FOR SPEED UP!
                                       preconditioner_drop_tolerance=0.01, # ...if fill 7-18 (hard), DT 1e-2 (7) to 1e-5 (18)
                                       number_orthogonalizations=2, # NORTH - increase if hard!
                                       )

        # -------------- GWF -------------------------
        gwf = flopy.mf6.ModflowGwf(sim, modelname=self.modelname, save_flows=True, newtonoptions = newtonoptions,) 

        # -------------- DIS -------------------------       

        import disv2disu
        from importlib import reload
        reload(disv2disu)
        Disv2Disu = disv2disu.Disv2Disu           
        #print(self.idomain[:,500])    
        #print(self.botm[:,500])
        dv2d = Disv2Disu(self.vertices, self.cell2d, self.top, self.botm, staggered=staggered, disv_idomain = self.idomain)
        disu_gridprops = dv2d.get_gridprops_disu6()
        disu = flopy.mf6.ModflowGwfdisu(gwf, **disu_gridprops) # This is the flow package
        
        # -------------- NPF -------------------------
               
        npf = flopy.mf6.modflow.mfgwfnpf.ModflowGwfnpf(gwf, xt3doptions=P.xt3d, k=self.k11, k22=self.k22, k33=self.k33, 
                                                       angle1 = self.angle1, angle2 = self.angle2, angle3 = self.angle3, 
                                                       #angle1 = 0., angle2 = 0., angle3 = 0.,
                                                       icelltype = 1,
                                                       save_flows=False, save_specific_discharge=False,)
                                                       #dev_minimum_saturated_thickness = 1)# try 0.1 then 0.001... no more than 1m!
        
        # -------------- IC / STO / WEL-------------------------
        ic = flopy.mf6.ModflowGwfic(gwf, strt = ic_array)
        if period == 'Steady': 
            csv_file = self.modelname + "_steady.csv"   #To write observation to  
        if period == 'Past':  
            sto = flopy.mf6.modflow.mfgwfsto.ModflowGwfsto(gwf, storagecoefficient=None, iconvert=1, 
                                                           ss=self.ss, sy = self.sy)
            #print(self.spd_wel_past)
            wel = flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(gwf, print_input=True, print_flows=True, 
                                                           stress_period_data = self.spd_wel_past, 
                                                           save_flows=True,)
            csv_file = self.modelname + "_past.csv"   #To write observation to      
        if period == 'Future':
            sto = flopy.mf6.modflow.mfgwfsto.ModflowGwfsto(gwf, storagecoefficient=None, iconvert=1, 
                                                           ss=self.ss, sy = self.sy)
            wel = flopy.mf6.modflow.mfgwfwel.ModflowGwfwel(gwf, print_input=True, print_flows=True, 
                                                           stress_period_data = self.spd_wel_future, 
                                                           save_flows=True,)
            csv_file = self.modelname + "_future.csv" # To write observation to
            
        # -------------- CH / RCH -------------------------
        ch = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, maxbound=len(self.chd_rec),stress_period_data=self.chd_rec,)
        rch = flopy.mf6.modflow.mfgwfrch.ModflowGwfrch(gwf, maxbound=len(self.rch_rec),stress_period_data=self.rch_rec,)
               
        # -------------- OBS / OC -------------------------
        obsdict = {csv_file: self.obslist}
        obs = flopy.mf6.ModflowUtlobs(gwf, filename=self.modelname, print_input=True, continuous=obsdict,)

        oc = flopy.mf6.ModflowGwfoc(gwf, budget_filerecord='{}.bud'.format(self.modelname), 
                                    head_filerecord='{}.hds'.format(self.modelname),
                                    saverecord=[('HEAD', 'LAST'),], printrecord=None,)
        
        # -------------- WRITE AND RUN SIMULATION -------------------------
        sim.write_simulation(silent = True)   

        #print('Running simulation for ', self.modelname, ' ...')
        success, buff = sim.run_simulation(silent = True)   
        print('Period = ', period, '\n   Model success = ', success)

        if success:
            '''# SAVE MFSIM.LST ------------------
            from pathlib import Path
            def copy_and_rename_pathlib(src_path, dest_path, new_name):
                # Create Path objects
                src_path = Path(src_path)
                dest_path = Path(dest_path)
                new_path = dest_path / new_name # Copy and rename the file
                src_path.rename(new_path)
                
            source_file = "../modelfiles/mfsim.lst"
            destination_folder = "../modelfiles/"
            
            if period == 'Steady': new_file_name = "mfsim_success_steady.txt"
            if period == 'Past':   new_file_name = "mfsim_success_past.txt"
            if period == 'Future': new_file_name = "mfsim_success_future.txt"
                    
            path = os.path.join(destination_folder, new_file_name) # delete old file if exists
            if os.path.isfile(path): os.remove(path)            
            copy_and_rename_pathlib(source_file, destination_folder, new_file_name)
            # ------------------------------'''
            
            fname = '{}.hds'.format(self.modelname)
            hds = flopy.utils.binaryfile.HeadFile(os.path.join(P.workspace, fname))  
            times = hds.get_times()
            head = hds.get_data(totim=times[-1]) 
            #bdobj = gwf.output.budget()
            #spdis = bdobj.get_data(text="DATA-SPDIS")[0]
            obs_data = gwf.obs
            t1 = datetime.now()
            run_time = t1 - t0
            print('   run_time = ', run_time.total_seconds())
            return(gwf, head, obs_data, run_time.total_seconds()) 
        
        else:
            '''# SAVE MFSIM.LST ------------------
            from pathlib import Path
            def copy_and_rename_pathlib(src_path, dest_path, new_name):
                # Create Path objects
                src_path = Path(src_path)
                dest_path = Path(dest_path)
                new_path = dest_path / new_name # Copy and rename the file
                src_path.rename(new_path)
                
            source_file = "../modelfiles/mfsim.lst"
            destination_folder = "../modelfiles/"
            
            if period == 'Steady': new_file_name = "mfsim_fail_steady.txt"
            if period == 'Past':   new_file_name = "mfsim_fail_past.txt"
            if period == 'Future': new_file_name = "mfsim_fail_future.txt"
                    
            path = os.path.join(destination_folder, new_file_name) # delete old file if exists
            if os.path.isfile(path): os.remove(path) 
            copy_and_rename_pathlib(source_file, destination_folder, new_file_name)
            # ------------------------------'''
            
            print('   Re-writing IMS - Take 2')
            sim.remove_package(package_name='ims')
            
            ims = flopy.mf6.ModflowIms(sim, print_option='ALL', 
                                       complexity    = 'Complex',
                                       outer_dvclose = 1e-5, 
                                       inner_dvclose = 1e-6, 
                                       inner_maximum = 200,
                                       outer_maximum = 60, 
                                       #no_ptcrecord = ['ALL'],
                                       #reordering_method=['RCM'],
                                       # BACKTRACKING (STOPS SOLUTION RUNNING AWAY WITH NEWTON RAPHSON!) 
                                       backtracking_number = 10,
                                       backtracking_tolerance = 10, #1.01 (aggressive) to 10000
                                       backtracking_reduction_factor = 0.5, # 0.1-0.3, or 0.9 when non-linear convergence HARD 
                                       # LINEAR SOLVER
                                       linear_acceleration = "BICGSTAB",
                                       preconditioner_levels = 7, #1 to 18, usually 7-8 for optimal
                                       preconditioner_drop_tolerance=0.0001, # ...if fill 7-18 (hard), DT 1e-2 (7) to 1e-5 (18)
                                       number_orthogonalizations=10, # NORTH (4-10) - increase if hard!)
                                       relaxation_factor = 0.97,
                                       # DBD USED FOR OSCILLATING NON-LINEAR
                                       under_relaxation = 'DBD',
                                       under_relaxation_kappa = 0.1, #0.05 (aggressive) to 0.3
                                       under_relaxation_theta = 0.7, # 0.5 - 0.9
                                       under_relaxation_gamma = 0.1, # 0-0.2 doesnt make big difference
                                       under_relaxation_momentum = 0.001, #0-0.001 doesn't make big difference
                                      )

            sim.ims.write()
            success2, buff = sim.run_simulation(silent = True)   
            print('Model success2 = ', success2)
            
            if success2:
                '''# SAVE MFSIM.LST ------------------
                from pathlib import Path
                def copy_and_rename_pathlib(src_path, dest_path, new_name):
                    # Create Path objects
                    src_path = Path(src_path)
                    dest_path = Path(dest_path)
                    new_path = dest_path / new_name # Copy and rename the file
                    src_path.rename(new_path)

                source_file = "../modelfiles/mfsim.lst"
                destination_folder = "../modelfiles/"

                if period == 'Steady': new_file_name = "mfsim_success2_steady.txt"
                if period == 'Past':   new_file_name = "mfsim_success2_past.txt"
                if period == 'Future': new_file_name = "mfsim_success2_future.txt"

                path = os.path.join(destination_folder, new_file_name) # delete old file if exists
                if os.path.isfile(path): os.remove(path)            
                copy_and_rename_pathlib(source_file, destination_folder, new_file_name)
                # ------------------------------'''
            
                fname = '{}.hds'.format(self.modelname)
                hds = flopy.utils.binaryfile.HeadFile(os.path.join(P.workspace, fname))  
                times = hds.get_times()
                head = hds.get_data(totim=times[-1]) 
                obs_data = gwf.obs
                t1 = datetime.now()
                run_time = t1 - t0
                print('   run_time = ', run_time.total_seconds())
                return(gwf, head, obs_data, run_time.total_seconds())      
            '''else:
                # SAVE MFSIM.LST ------------------
                from pathlib import Path
                def copy_and_rename_pathlib(src_path, dest_path, new_name):
                    # Create Path objects
                    src_path = Path(src_path)
                    dest_path = Path(dest_path)
                    new_path = dest_path / new_name # Copy and rename the file
                    src_path.rename(new_path)

                source_file = "../modelfiles/mfsim.lst"
                destination_folder = "../modelfiles/"

                if period == 'Steady': new_file_name = "mfsim_fail2_steady.txt"
                if period == 'Past':   new_file_name = "mfsim_fail2_past.txt"
                if period == 'Future': new_file_name = "mfsim_fail2_future.txt"

                path = os.path.join(destination_folder, new_file_name) # delete old file if exists
                if os.path.isfile(path): os.remove(path) 
                copy_and_rename_pathlib(source_file, destination_folder, new_file_name)
                # ------------------------------'''
            
            '''else:
                print('   Re-writing IMS - Take 3')
                
                if period == 'Past':   # Increase number of timesteps to help convergence
                    future_years = 5
                    nts_future = future_years * 12
                    tdis_future = [(future_years * 365, nts_future, 1.2)] # period length, number of timesteps, tsmult
                    sim.remove_package(package_name='tdis')
                    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim, nper=len(tdis_future), perioddata=tdis_future)
                
                # More aggressive solver settings
                sim.remove_package(package_name='ims')
                
                ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', 
                            complexity    = 'Complex',
                            outer_dvclose = 1e-2, 
                            inner_dvclose = 1e-6, 
                            outer_maximum = 60,
                            inner_maximum = 200,
                            linear_acceleration = "BICGSTAB",
                            reordering_method=['RCM'],
                            no_ptcrecord = ['ALL'],
                            under_relaxation = 'DBD',
                            under_relaxation_kappa = 0.3, #0.05 (aggressive) to 0.3
                            under_relaxation_theta = 0.7, # 0.5 - 0.9
                            under_relaxation_gamma = 0.1, # 0-0.2 doesnt make big difference
                            under_relaxation_momentum = 0.001, #0-0.001 doesn't make big difference
                            backtracking_number = 10,
                            backtracking_tolerance = 10, #1.01 (aggressive) to 10000
                            backtracking_reduction_factor = 0.7, # 0.1-0.3, or 0.9 when non-linear convergence HARD 
                            preconditioner_levels=5, #1 to 5... PLAY WITH THIS FOR SPEED UP!
                            preconditioner_drop_tolerance=0.01, # ...if fill 7-18 (hard), DT 1e-2 (7) to 1e-5 (18)
                            number_orthogonalizations=2,)
                sim.ims.write()
                success3, buff = sim.run_simulation(silent = True)   
                print('Model success3 = ', success3)
                
                if success3:
                    fname = '{}.hds'.format(self.modelname)
                    hds = flopy.utils.binaryfile.HeadFile(os.path.join(P.workspace, fname))  
                    times = hds.get_times()
                    head = hds.get_data(totim=times[-1]) 
                    obs_data = gwf.obs
                    t1 = datetime.now()
                    run_time = t1 - t0
                    print('   run_time = ', run_time.total_seconds())
                    return(gwf, head, obs_data, run_time.total_seconds())   '''
        


# In[2]:


# angle 1 (DIP DIRECTION) rotates around z axis counterclockwise looking from +ve z.
def find_angle1(nv): # nv = normal vector to surface
    # The dot product of perpencicular vectors = 0
    # A vector perpendicular to nv would be [a,b,c]
    if nv[2] == 0:
        angle1 = 0.
    else:
        a = nv[0]
        b = nv[1]
        c = -(a*nv[0]+b*nv[1])/nv[2]
        v = [a,b,c]
        if np.isnan(v[0]) == True or np.isnan(v[1]) == True: 
            angle1 = 0.
        if v[0] == 0.:
            if v[1] > 0:
                angle1 = 90
            else:
                angle1 = -90
        else:             
            tantheta = v[1]/v[0] 
            angle1 = np.degrees(math.atan(tantheta))
    return(angle1)

# angle 2 (DIP) rotates around y axis clockwise looking from +ve y.
def find_angle2(nv): # nv = normal vector to surface
    # The dot product of perpencicular vectors = 0
    # A vector perpendicular to nv would be [a,b,c]
    if nv[2] == 0:
        angle2 = 0.
    else:
        a = nv[0]
        b = nv[1]
        c = -(a*nv[0]+b*nv[1])/nv[2]
        v = [a,b,c]
        if np.isnan(v[0]) == True or np.isnan(v[1]) == True or np.isnan(v[2]) == True:
            angle2 = 0.
        else:
            v_mag = (v[0]**2 + v[1]**2 + v[2]**2)**0.5 
            costheta = v[2]/v_mag
            angle2 = 90-np.degrees(math.acos(costheta)) 
    return(angle2)


# In[5]:


print('loop2flopy routines loaded!')

