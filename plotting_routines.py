#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def multiplot_vgrid_transect(P, models, array, X0, Y0, X1, Y1, vmin = None, vmax = None): # array needs to be a string of a property eg. 'k11', 'angle2'
    nmodels = len(models)
    if nmodels > 1: fig = plt.figure(figsize = (10,2*nmodels))
    if nmodels ==1: fig = plt.figure(figsize = (10,2.25))
    fig.suptitle("TRANSECT - " + array)
    for i in range(nmodels):
        M = models[i]
        a = getattr(M, array)
        ax = plt.subplot(nmodels, 1, i+1)
        ax.set_title(M.modelname, size = 10) 
        xsect = flopy.plot.PlotCrossSection(modelgrid=M.vgrid, line={"line": [(X0, Y0),(X1, Y1)]}, 
                                            extent = [P.x0,P.x1,P.z0,P.z1], geographic_coords=True)
        csa = xsect.plot_array(a = a, cmap = 'Spectral', alpha=0.8, vmin = vmin, vmax = vmax)
        if i == nmodels-1: ax.set_xlabel('x (m)', size = 10)
        if i == int(nmodels/2): ax.set_ylabel('z (m)', size = 10)
        if nmodels>1: linecollection = xsect.plot_grid(lw = 0.1, color = 'black') # Don't plot grid for reference
        plt.colorbar(csa, shrink = 0.7)
    plt.tight_layout()  
    plt.show()   


# In[ ]:


def plot_flow_features(P, models):
    
    fig = plt.figure(figsize = (12,5))
    fig.suptitle("Flow Features")
    for i, M in enumerate(models):
        ax = plt.subplot(1,3,i+1)
        
        M.ibd = np.zeros(M.vgrid.ncpl, dtype=int) # for plotting special cells! 0-all, 1-wel, 2-obs, 3-chd
                
        for coords in P.xypumpbores:
            point = Point(coords)
            cells = M.gi.intersect(point)["cellids"]
            cells = np.array(list(cells))
            M.ibd[cells] = 1

        for coords in P.xyobsbores:
            point = Point(coords)
            cells = M.gi.intersect(point)["cellids"]
            cells = np.array(list(cells))
            M.ibd[cells] = 2
                      
        west_bd = LineString([(P.x0, P.y0), (P.x0, P.y1)]) # Western edge
        cells = M.gi.intersects(west_bd, shapetype="linestring")
        cells = cells.cellids.tolist()
        M.ibd[cells] = 3
        
        east_bd = LineString([(P.x1, P.y0), (P.x1, P.y1)]) # Eastern edge
        cells = M.gi.intersects(east_bd, shapetype="linestring")
        cells = cells.cellids.tolist()
        M.ibd[cells] = 3
        
        pmv = flopy.plot.PlotMapView(modelgrid=M.vgrid)
        pmv.plot_array(M.ibd, alpha = 0.6)
        if M.plan == 'car': P.sg.plot(ax=ax, edgecolor='black', lw = 0.2)
        if M.plan == 'tri': P.tri.plot(ax=ax, edgecolor='black', lw = 0.2)
        if M.plan == 'vor': P.vor.plot(ax=ax, edgecolor='black', lw = 0.2)
            
        ax.set_title(M.modelname, size = 10) 
        ax.set_xlabel('x (m)', size = 10)
        if i == 0: ax.set_ylabel('y (m)', size = 10)
        if i == 1 or i ==2: ax.set_yticks([])
    plt.savefig('Flow_features.jpg')


# In[ ]:


## PLOTTING HYDRAULIC PROPERTIES

def multiplot_prop_plan(P, models, array, layer, vmin = None, vmax = None):   # array needs to be a string of a property eg. 'k11', 'logk11'  
    fig = plt.figure(figsize = (10,12))
    fig.suptitle("PLAN - " + array)
    nmodels = len(models)   
    for i in range(nmodels):
        ax = plt.subplot(3,2,i+1)
        M = models[i]
        model = M.gwf
        a = getattr(M, array)
                 
        ax.set_title(M.modelname, size = 10)
        mapview = flopy.plot.PlotMapView(model=model, layer = layer)
        plan = mapview.plot_array(a, cmap='Spectral', alpha=0.8, vmin = vmin, vmax = vmax)
        linecollection = mapview.plot_grid(lw = 0.1, color = 'black')
        if i == 4 or i == 5: ax.set_xlabel('x (m)', size = 10)
        if i == 0 or i == 2 or i == 4: ax.set_ylabel('y (m)', size = 10)
        #if transient == True:
        #    plt.plot(wel_cell_coords[0], wel_cell_coords[1], 'o', c = 'black', ms = 5)
        plt.colorbar(plan, shrink = 0.4)
    plt.tight_layout()  
    
def multiplot_prop_transect(P, models, array, X0, Y0, X1, Y1, vmin = None, vmax = None): # array needs to be a string of a property eg. 'k11', 'angle2'
    nmodels = len(models)
    if nmodels > 1: fig = plt.figure(figsize = (10,2*nmodels))
    if nmodels ==1: fig = plt.figure(figsize = (10,2.5))
    fig.suptitle("TRANSECT - " + array)
    for i in range(nmodels):
        M = models[i]
        model = M.gwf
        a = getattr(M, array)
        
        ax = plt.subplot(nmodels, 1, i+1)
        ax.set_title(M.modelname, size = 10) 
        xsect = flopy.plot.PlotCrossSection(model=model, line={"line": [(X0, Y0),(X1, Y1)]}, 
                                            extent = [P.x0,P.x1,P.z0,P.z1], geographic_coords=True)
        csa = xsect.plot_array(a = a, cmap = 'Spectral', alpha=0.8, vmin = vmin, vmax = vmax)
        linecollection = xsect.plot_grid(lw = 0.1, color = 'black')
        if i == nmodels-1: ax.set_xlabel('x (m)', size = 10)
        if i == int(nmodels/2): ax.set_ylabel('z (m)', size = 10)
        if nmodels>1: linecollection = xsect.plot_grid(lw = 0.1, color = 'black') # Don't plot grid for reference
        plt.colorbar(csa, shrink = 0.7)
    plt.tight_layout()  
    plt.show()    
    


# In[8]:


### PLOTTING HEADS

def multiplot_watertable(P, models, period): 
    nmodels = len(models)
    fig = plt.figure(figsize = (10,12))
    #contours = np.arange(0, 60, 5)
    from flopy.plot import styles
    
    fig.suptitle("PLAN")
    for i in range(nmodels):
        M = models[i]
        model = M.gwf
        ax = plt.subplot(3, 2, i+1)
        ax.set_title(M.modelname, size = 10) 
        if period == 'Steady' : water_table = flopy.utils.postprocessing.get_water_table(M.head_ss, hdry=-1e30)  
        if period == 'Past'   : water_table = flopy.utils.postprocessing.get_water_table(M.head_present, hdry=-1e30)  
        if period == 'Future' : water_table = flopy.utils.postprocessing.get_water_table(M.head_future, hdry=-1e30)  
            
        m = flopy.modflow.Modflow.load(str(M.modelname + '.nam'), model_ws=P.workspace)
        pmv = flopy.plot.PlotMapView(modelgrid = m.modelgrid, ax=ax)
        #pmv.plot_vector(M.ss_spdis["qx"], M.spdis["qy"], alpha=0.5)
        linecollection = pmv.plot_grid(lw = 0.1)
        h = pmv.plot_array(water_table, cmap='Spectral')#, vmin=hmin, vmax=hmax, )    
        #if period == 'Steady' : water_table = M.head_ss
        #if period == 'Past'   : water_table = M.head_present
        #if period == 'Future' : water_table = M.head_future  
        #hmin, hmax = -10,60 #water_table.min(), water_table.max()
        
        for j in range(len(P.xyobsbores)):
            ax.plot(P.xyobsbores[j][0], P.xyobsbores[j][1],'o', ms = '4', c = 'black')
            ax.annotate(j, (P.xyobsbores[j][0], P.xyobsbores[j][1]+60), c = 'black', size = 12) #, weight = 'bold')
        
        for j in range(len(P.xypumpbores)):
            ax.plot(P.xypumpbores[j][0], P.xypumpbores[j][1],'o', ms = '4', c = 'red')
            ax.annotate(j, (P.xypumpbores[j][0], P.xypumpbores[j][1]+60), c = 'red', size = 12) #, weight = 'bold')
            
        if i == 2 or i == 3: ax.set_xlabel('x (m)', size = 10)
        if i == 0 or i == 2: ax.set_ylabel('y (m)', size = 10)
        plt.plot([P.fx1, P.fx2],[P.fy1, P.fy2], c = 'black', lw = 0.5)
        
        '''with styles.USGSMap():
            pmv = flopy.plot.PlotMapView(modelgrid = model.modelgrid, ax=ax)
            #pmv.plot_vector(M.ss_spdis["qx"], M.spdis["qy"], alpha=0.5)
            linecollection = pmv.plot_grid(lw = 0.1)
            h = pmv.plot_array(water_table, cmap='Spectral')#, vmin=hmin, vmax=hmax, )
            #c = pmv.contour_array(water_table, levels=contours, colors="black", linewidths=0.75, linestyles=":", )
            #plt.clabel(c, fontsize=8)
            pmv.plot_inactive()
            plt.colorbar(h, ax=ax, shrink=0.5)'''
        
        #pmv.plot_vector(M.ss_spdis["qx"], M.spdis["qy"], alpha=0.5)
        linecollection = pmv.plot_grid(lw = 0.1)
        
    plt.tight_layout()  

def multiplot_heads_plan(P, layer, models, period):#, obs_points):#, vmin, vmax):   
    nmodels = len(models)
    fig = plt.figure(figsize = (10,6))
    fig.suptitle("PLAN")
    for i in range(nmodels):
        ax = plt.subplot(2, int(nmodels/2)+1,i+1)
        M = models[i]
        model = M.gwf
        if period == 'Steady': array = M.head_ss
        if period == 'Past': array = M.head_present
        if period == 'Future': array = M.head_future
        ax.set_title(M.modelname, size = 10)
        mapview = flopy.plot.PlotMapView(model=model, layer = layer)
        plan = mapview.plot_array(array, cmap='Spectral', alpha=0.8)#, vmin = vmin, vmax = vmax)
        #mapview.plot_vector(M.ss_spdis["qx"], M.ss_spdis["qy"], alpha=0.5)
        linecollection = mapview.plot_grid(lw = 0.1)
        #if i == 2 or i == 3: ax.set_xlabel('x (m)', size = 10)
        #if i == 0 or i == 2: ax.set_ylabel('y (m)', size = 10)
        plt.plot([P.fx1, P.fx2],[P.fy1,P.fy2], c = 'black', lw = 0.5)
        #if transient == True:
        #    plt.plot(wel_cell_coords[0], wel_cell_coords[1], 'o', c = 'black', ms = 5)
        #plt.colorbar(plan, shrink = 0.4)
           
        for j in range(len(P.xyobsbores)):
            ax.plot(P.xyobsbores[j][0], P.xyobsbores[j][1],'o', ms = '4', c = 'black')
            ax.annotate(j, (P.xyobsbores[j][0], P.xyobsbores[j][1]+60), c = 'black', size = 12) #, weight = 'bold')
    plt.tight_layout()  

def multiplot_heads_transect(period, X0, Y0, X1, Y1):#, vmin, vmax):    
    nmodels = len(models)
    fig = plt.figure(figsize = (10, 2*nmodels))
    fig.suptitle("TRANSECT")
    for i in range(nmodels):
        ax = plt.subplot(nmodels,1,i+1)
        M = models[i]
        model = M.gwf
        if period == 'Steady': array = M.ss_head
        if period == 'Past': array = M.head_present
        if period == 'Future': array = M.head_future
        ax.set_title(M.modelname, size = 10) 
        xsect = flopy.plot.PlotCrossSection(model=model, line={"line": [(X0, Y0),(X1, Y1)]}, 
                                            extent = [x0,x1,z0,z1], geographic_coords=True)
        csa = xsect.plot_array(a = array, cmap = 'Spectral', alpha=0.8)#, vmin = vmin, vmax = vmax)
        if i == nmodels-1: ax.set_xlabel('x (m)', size = 10)
        if i == int(nmodels/2): ax.set_ylabel('z (m)', size = 10)
        linecollection = xsect.plot_grid(lw = 0.1, color = 'black')
        plt.colorbar(csa, shrink = 0.7)
    plt.tight_layout()  


# In[9]:


def vtk_make(self, result = False, pfx = None):
    f = open(self.modelname + '.vtk','w')
    f.write('# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n')
    NE = self.disu_gridprops["nodes"]
    NP = NE*6
    f.write('POINTS %i float\n' % NP)
    for k in range(NE):
        for j in range(4,7,1):
            f.write('%.10g %.10g %.10g\n' % (self.disu_gridprops['vertices'][self.disu_gridprops['cell2d'][k][j]][1],
                                    self.disu_gridprops['vertices'][self.disu_gridprops['cell2d'][k][j]][2],
                                    self.disu_gridprops['top'][k]))
        for j in range(4,7,1):
            f.write('%.10g %.10g %.10g\n' % (self.disu_gridprops['vertices'][self.disu_gridprops['cell2d'][k][j]][1],
                                    self.disu_gridprops['vertices'][self.disu_gridprops['cell2d'][k][j]][2],
                                    self.disu_gridprops['bot'][k]))

    f.write('CELLS %i %i\n' % (NE,NE*7))
    ct = 0
    for i in range(NE):
        f.write('%i %i %i %i %i %i %i \n' % (6,ct,ct+1,ct+2,ct+3,ct+4,ct+5))
        ct+=6
    f.write('CELL_TYPES %i\n' % NE)
    for i in range(NE):
        f.write('%i \n' % 13)

    f.write('CELL_DATA %i\nSCALARS Layer int 1 \nLOOKUP_TABLE default\n' % NE )
    for i in range(NE):
        if i == 0:
            f.write('%i\n' % self.disu_gridprops['ihc'][0])
        else:
            f.write('%i\n' % self.disu_gridprops['ihc'][np.cumsum(self.disu_gridprops['iac'])[i-1]])
    f.close()

    if result:
        nstp = np.shape(self.head)[0]
        for i in range(nstp):
            shutil.copy(self.modelname + '.vtk', self.modelname + pfx + '_' + str(i) + '.vtk')
            f = open(self.modelname + pfx + '_' + str(i) + '.vtk','a')
            f.write('SCALARS head float 1 \nLOOKUP_TABLE default\n')
            for j in range(NE):
                f.write('%g\n' % self.head[i,0,0,j])
            f.close()
            
def make_vtk(P, nam_file): # from https://flopy.readthedocs.io/en/latest/Notebooks/export_vtk_tutorial.html
    from flopy.export import vtk
    from pathlib import Path
    from tempfile import TemporaryDirectory
    workspace = P.workspace
    
    ml = flopy.modflow.Modflow.load(nam_file, model_ws=workspace, check=False)

    tempdir = TemporaryDirectory()
    workspace = Path(tempdir.name)

    output_dir = P.workspace / "arrays_test"
    output_dir.mkdir(exist_ok=True)
    
    ml.dis.top.export(output_dir / "TOP", fmt="vtk")
    ml.dis.botm.export(model_bottom_dir = output_dir / "BOTM", fmt="vtk")
    ml.rch.rech.export(output_dir / "RECH", fmt="vtk", pvd=True)
    ml.upw.hk.export(model_hk_dir = output_dir / "HK", smooth=True, fmt="vtk", name="HK", point_scalars=True)
    
    # set up package export folder
    output_dir = workspace / "package_output_test"
    output_dir.mkdir(exist_ok=True)

    # export
    ml.dis.export(output_dir / "DIS", fmt="vtk")
    ml.upw.export(output_dir / "UPW", fmt="vtk", point_scalars=True, xml=True)
    ml.export(workspace / "model_output_test", fmt="vtk")
    
    # create a binary XML VTK object and enable PVD file writing
    vtkobj = vtk.Vtk(ml, xml=True, pvd=True, vertical_exageration=10)
    vtkobj = vtk.Vtk(ml, vertical_exageration=10) # Create a vtk object

    ## create some random array data
    r_array = np.random.random(ml.modelgrid.nnodes) * 100
    vtkobj.add_array(r_array, "random_data") ## add random data to the VTK object
    vtkobj.add_array(ml.dis.botm.array, "botm") ## add the model botom data to the VTK object
    vtkobj.write(output_dir / "Array_example" / "model.vtu") ## write the vtk object to file
    vtkobj = vtk.Vtk(ml, xml=True, pvd=True, vertical_exageration=10) # create a vtk object

    recharge = ml.rch.rech.transient_2ds ## add recharge to the VTK object
    vtkobj.add_transient_array(recharge, "recharge", masked_values=[0,],)
    vtkobj.write(output_dir / "tr_array_example" / "recharge.vtu") ## write vtk files
    vtkobj = vtk.Vtk(ml, xml=True, pvd=True, vertical_exageration=10) # create the vtk object

    spd = ml.wel.stress_period_data ## add well fluxes to the VTK object
    vtkobj.add_transient_list(spd, masked_values=[0,],)
    vtkobj.write(output_dir / "tr_list_example" / "wel_flux.vtu") ## write vtk files


# In[ ]:


def plot_runtime_complexity():   
    titles = ['Steady', 'Past', 'Future']
    fig = plt.figure(figsize = (10,4))
    fig.suptitle('Model run times')
        
    for i in range(3): # each time period
        ax = plt.subplot(1, 3, i+1)
        ax.set_title(titles[i], size = 10)
        
        for m in range(4): 
            for n in range(nruns):
                for c in range(len(complex_options)): 
                    ax.plot(m, run_time_results[m, i, 0, n],'o', ms = '4', alpha = 0.6, c = 'green') # moderate
                    ax.plot(m, run_time_results[m, i, 1, n],'o', ms = '4', alpha = 0.6, c = 'blue')  # complex
        ax.set_ylim(0, 70)
        if i ==0: ax.set_ylabel('run_time (s)', size = 10)
        ax.set_xticks([0,1,2,3])
        ax.set_xticklabels(['SS', 'US', 'SU', 'UU'])
    plt.legend(['Moderate', 'Complex'])
    plt.tight_layout()  
    fig.savefig('../figures/complexity_runtime.tif', dpi=300)


# In[ ]:


def plot_bylayer(P, models, layer, vmin = None, vmax = None):
    
    fig = plt.figure(figsize=(12, 8))
    nmodels = len(models)
    for i in range(nmodels):
        ax = plt.subplot(2,3,i+1, aspect="equal")
        M = models[i]

        model = M.gwf
        water_table = flopy.utils.postprocessing.get_water_table(M.gwf.output.head().get_data())
        M.heads_disv = -1e30 * np.ones_like(M.idomain, dtype=float) 
        for i, h in enumerate(water_table):
            if math.isnan(h) == False: 
                M.heads_disv[M.cellid_disu==i] = h        
        pmv = flopy.plot.PlotMapView(modelgrid=M.vgrid)
        H = pmv.plot_array(M.heads_disv[layer], vmin = vmin, vmax = vmax, cmap = 'Spectral', alpha = 0.6)
        for j in range(len(P.xyobsbores)):
            ax.plot(P.xyobsbores[j][0], P.xyobsbores[j][1],'o', ms = '4', c = 'black')
            ax.annotate(P.idobsbores[j], (P.xyobsbores[j][0], P.xyobsbores[j][1]+100), c='black', size = 12) #, weight = 'bold')
        
        for j in range(len(P.xypumpbores)):
            ax.plot(P.xypumpbores[j][0], P.xypumpbores[j][1],'o', ms = '4', c = 'red')
            ax.annotate(P.idpumpbores[j], (P.xypumpbores[j][0], P.xypumpbores[j][1]+100), c='red', size = 12) #, weight = 'bold')
            
        if M.plan == 'car': P.sg.plot(ax=ax, edgecolor='black', lw = 0.2)
        if M.plan == 'tri': P.tri.plot(ax=ax, edgecolor='black', lw = 0.2)
        if M.plan == 'vor': P.vor.plot(ax=ax, edgecolor='black', lw = 0.2)
        ax.set_title(M.modelname, size = 10)
        plt.colorbar(H, shrink = 0.4)


# In[ ]:


print('Plotting routines loaded!')

