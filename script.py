# -*- coding: utf-8 -*-
'''This how HAT supports rudimentary scripting, the run script button will run the script_main() func
actions can be scripted be calling functions within hat, some documentation and example
scripts are at end of file'''
import glob
import numpy as np
import math

def script_main(hat):
    hat.restore()
    folder_list = glob.glob("D:/CV_NaOH_2/*")
    for i,folder in enumerate(folder_list):

        #get list of file names in directory ending with tif
        file_list = glob.glob(folder+"/*.tif")

        #filter out any dark images
        files_no_dark = [x for x in file_list if "dark" not in x]

        #load files 
        hat.image_stack.select_images(files_no_dark)      
        hat.load()    
        hat.setUnits(0) #q units

        #make hkl map        
        hat.loadMasks("quickrotate.msk") 
        hat.setProjection(0) 
        hat.makeProjection()
        hat.setColorMap("viridis")
        hat.setLevels(20, 70) 
        hat.makeFigure("D:/output/"+str(i))


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
hat.extractrois(filename) - save roi as xy ascii, filename is output filename
hat.saverocks(foldername) - save rocking scans form region of interest, pass a folder name
                            each file is a csv with angle and intensity for each verticle line along y                        
hat.setProfileAxis(value) - change the axis of interest for the profile 
                            "H (Qx)" = 1, "K (Qy)" = 2, "HK (Qr)" = 3,"L (qz)" =  4
                            Default is 4.

'''

