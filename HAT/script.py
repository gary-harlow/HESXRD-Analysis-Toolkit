# -*- coding: utf-8 -*-
'''This how HAT supports rudimentary scripting, the run script button will run the script_main() func
actions can be scripted be calling functions within hat, some documentation and example
scripts are at end of file'''
import glob
from http.client import REQUEST_URI_TOO_LONG
import numpy as np
import math
import time
from xrayhat.fileio import *
import mayavi.mlab as mlab



def script_main(hat):
    import numpy as np
    grid=hat.showData
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(grid),vmin=0, vmax=100)
    #mlab.outline()
    #mlab.show()
    np.save("d:/2022/papers/hat/20_35_35_3d_axis.npy", hat.zgrid)
    return
    
    print(hat.binBounds)
    return

    import plotly.graph_objects as go
    import numpy as np



    
    X, Y, Z = np.ogrid[-1*hat.qmax:hat.qmax:200j, -1*hat.qmax:hat.qmax:200j, -1*hat.qmax:hat.qmax:200j]
    values = hat.showData
    
    twenty = 0.2*np.max(hat.showData)
    eighty = 0.8*np.max(hat.showData)


    #export values
    np.save("d:/test.npy", values)

    #obj = contour3d(values, contours=4, transparent=True)
    #mlab.contour3d(values, contours=4, transparent=True)
    mlab.pipeline.volume(mlab.pipeline.scalar_field(values),vmin=140, vmax=250)

    #mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
    #                        plane_orientation='x_axes',
    #                        slice_index=10,
    #                    )
    #mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
    #                        plane_orientation='y_axes',
    #                        slice_index=10,
    #                   )
    mlab.outline()
    mlab.show()
    print("done")


    return
    print("We have made the mesh grid")

    fig = go.Figure(data=go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        isomin=0.1,
        isomax=0.8,
        opacity=0.1, # needs to be small to see through all surfaces
        surface_count=17, # needs to be a large number for good volume rendering
        ))
    fig.show()

    return
    hat.restore() 
    #hat.restoreImageStack() 
    #hat.load()
    #hat.loadMasks("E:/data/Dec2021/inplane.msk")
    hat.setProjection(0)
    #hat.makeProjection()

    h_list= [0,0, 2,0,1,1,1,2,2,3]
    k_list= [1,2,-1,3,1,0,2,0,1,0]
    path = "E:/data/Dec2021/pt111_HClO4_3/519"
    
    file_list1 = glob.glob(path+"/4/"+"/*.cbf")
    file_list2 = glob.glob(path+"/3/"+"/*.cbf")
    hat.image_stack.select_images(file_list1, file_list2)
    hat.image_stack.select_dark(path+"/d4_dark/WAXS_00001.dark.cbf",path+"/d3_dark/WAXS_00001.dark.cbf")
    hat.load()
    hat.loadMasks("E:/data/Dec2021/inplane.msk")
    hat.makeProjection()
    return
    
    for i in range(len(h_list)):
        h = h_list[i]
        k = k_list[i]

        ctr = str(h)+str(k)
        print(ctr)

        hat.clearMasks()
        hat.loadMasks("E:/data/Dec2021/masks/"+ctr+".msk")
        hat.setProjection(1)
        hat.makeProjection()

        hat.enableProfile(True) 

        hat.loadroi("E:/data/Dec2021/rois/"+ctr+"_ctr.roi")
        hat.saveprofile("E:/data/Dec2021/analysis/"+ctr+".csv") 

        hat.loadroi("E:/data/Dec2021/rois/"+ctr+"_bak1.roi")
        hat.saveprofile("E:/data/Dec2021/analysis/"+ctr+"_bak1.csv") 
    
        hat.loadroi("E:/data/Dec2021/rois/"+ctr+"_bak2.roi")
        hat.saveprofile("E:/data/Dec2021/analysis/"+ctr+"_bak2.csv") 

        hat.enableProfile(False)
    








"""    #Is this normalised?
    background = "D:/au_111_naoh_2_00023/au_111_naoh_2_00023-00380.tif"
    bknorm  = 1.1300110772777499
    bright = "D:/p07/bright.tif"

    folder_list = glob.glob("D:/CV_NaOH_3/*")
    file_display_list =  [None] * 80
    direction = True
    hat.loadMasks("D:/p07/quickmask.msk")
    hat.setUnits(0) #q units   
    hat.setColorMap("viridis") 

    for i,folder in enumerate(folder_list):

        #get list of files in folder if the name contains dark ignore it
        file_list = glob.glob(folder+"/*.tif")
        #we should sort by number here
        files_no_dark = [x for x in file_list if "dark" not in x]

        k=0        
        if i == 1:                  
            hat.p.param('Experiment', 'Beamline Preset').setValue(5) 
            hat.p.param('Experiment', 'Manual Start Angle').setValue(-47.01417)
            hat.p.param('Experiment', 'Manual End Angle').setValue(-43.09816)

        for k,file in enumerate(files_no_dark): 
            expt,monitor,omega,time = hat.image_stack.get_p07_attributes(file)
            mon = monitor *expt
            if i != 0:
                if direction == False:
                    #file_display_list[k]=file
                    hat.image_stack.image_data[k]=hat.image_stack.pubload_image(file)/mon-hat.image_stack.subtract_image/mon
                    #if i > 0 and k < 79:
                    #    hat.image_stack.image_data[k+1]=hat.image_stack.pubload_image(bright)-hat.image_stack.subtract_image
                        #file_display_list[k+1]=bright
                else:
                    #file_display_list[-k]=file
                    hat.image_stack.image_data[-k]=hat.image_stack.pubload_image(file)/mon-hat.image_stack.subtract_image/mon
                    #if i > 0 and k < 79:
                    #    hat.image_stack.image_data[-1*(k+1)]=hat.image_stack.pubload_image(bright)-hat.image_stack.subtract_image
                        #file_display_list[-1*(k+1)]=bright

            #load files 
            if i != 0:                   
                print(i*70+k+25*i+1)
                if (i*70+k+25*i+1) != 176:
                    continue  

                hat.single_thread == True
                print(i,k)
                #hat.image_stack.select_images(show)  
                #hat.load()                 
                hat.makeProjection()
                hat.setLevels(0, 200) 
                index = str(i*70+k+25*i+1).zfill(6)
                hat.makeFigure("D:/output/"+index)   
                #next we want to keep the meta data
                meta = [omega, time,hat.qmax]                          
                np.savez("D:/output/"+index+"_meta.npz",metadata=meta)                

        #for temp in range(1,26):    
        #    hat.makeFigure("D:/output/"+str(i*70+k+temp+1).zfill(4))      
        
        direction^=True        
        #return()"""


'''Documented functions
=====================
Useful actions for scripting

views
------
hat.makeProjection() - Create in-plane projection
hat.detectorTransform() - Out-of-plane map
hat.detectorMode() - just show the max/sum of the images using pixel coordinates, i.e. Detector View
hat.setProjection() - set the projection mode 0: hk, 1: h-l, 2: h-k, 3: sqrt(h^2+k^2)-l 
hat.makeFigure(filename) - save a figure using plotting.py, filename will be the output file name
hat.setLevels(cmin, cmax) - set the histogram minium and maxium colour levels
hat.setAngleRange(start, end) - set the angle region, can only do in pixel view

hat.setColorMap("viridis") - viridis is one example color map you can set others
                             (not so useful as you have to set your own colourmap in plotting.py)
hat.setUnits(value) - value is either 0 for Q or 1 for HKL

masks
------
hat.clearMasks()
hat.loadMasks(filename) - load masks from a file called fileName

files
------
hat.restore() - load saved parameters 
hat.restoreImageStack() - reselect last saved files

hat.image_stack.select_images(files, files2) - select the images to load, files is a list of image names (with path) to load
                                                  files2 is an optional 2nd list if there is a 2nd detector
                                                  if hdf5 file files should be a string and not a list
hat.image_stack.select_dark(file1, files2)  - similar syntax to above, dark is subtracted from each image
hat.image_stack.select_background(filename, filename2) - similar syntax to above, background is also subtracted from each image
hat.image_stack.select_background_dark(filename, filename2) - similar syntax to above, dark is subtracted  background image
hat.image_stack.read_log_file(filename)  - read log file, only works for some predefined beamlines
hat.load() - load the files defined in image_stack

Profile tools
-----------
hat.setProfileLog(value)     - Pass either true or false to change if the profile should have a log y axis
hat.enableProfile(value)  - Pass either true or false to enable/disable the box profile
hat.loadroi(filename) - load a ROI from a file called fileName
hat.saveprofile(filename) - save roi as xy ascii, filename is output filename
hat.saverocks(foldername) - save rocking scans form region of interest, pass a folder name
                            each file is a csv with angle and intensity for each verticle line along y                        
hat.setProfileAxis(value) - change the axis of interest for the profile 
                            "H (Qx)" = 1, "K (Qy)" = 2, "HK (Qr)" = 3,"L (qz)" =  4
                            Default is 4.

'''

