#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import sys,os
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QCheckBox, QProgressBar
from PyQt5.QtGui import QIcon, QMatrix2x3, QPixmap
from PyQt5 import QtGui
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QLabel, QGridLayout, QWidget,QPushButton, QMenuBar, QToolBar, QAction,QDialog
from PyQt5.QtCore import QSize, Qt, QObject, QThread, pyqtSignal, QRunnable  
from pyqtgraph.widgets.ValueLabel import ValueLabel    
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtCore import QUrl
from scipy import interpolate
import pyqtgraph.parametertree.parameterTypes as pTypes
import pyqtgraph as pg
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import matplotlib.pyplot as plt
import numpy as np
import fabio
import pickle
import fileio
import importlib
import locale
locale.setlocale(locale.LC_ALL, '')
from locale import atof

import plotting
import script
from crystal import *
from fileio import *


import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot
import matplotlib.ticker as ticker
from scipy import interpolate


class Ram(QWidget):
    image_stack = None
    def __init__(self):
        super().__init__()
        self.title = 'HESXRD Analysis Toolkit (HAT)'

        self.ImgState=0 #0 normal, 1 hkl, 2 inlplane, 5 combined background and dark images
        self.showData = None
        self.pixMaps=[]
        self.roixdata=None
        self.roiydata=None
        self.grid_h = None
        self.grid_k = None
        self.select_l_slice = True

        self.crystal=CrystallographyParameters(name='Crystal')
        self.experiment=ExperimentParameter(name='Experiment')
        self.initUI()
        

    def imageHoverEvent(self,event):
        
        """Show the position, pixel, and value under the mouse cursor.
        """
        binning=self.image_stack.binning
        xSize=self.image_stack.x_size
        ySize=self.image_stack.y_size
        
        if event.isExit():
            self.p3.setTitle("")
            return
        pos = event.pos()
        i, j = int(pos.y()), int(pos.x())
        ppos = self.imgLeft.mapToParent(pos)
        x, y = ppos.x(), ppos.y()        
        if self.ImgState==5 or self.ImgState==0: 
            #showing background image, or in pixel coordinates
            self.p3.setTitle(" (%0.3f, %0.3f) / (%d, %d), Binned / Est. Real I: %0.3f / %0.3f" % (x, ySize-y, x*binning, binning*(ySize-y),self.imgLeft.image[j,i],self.imgLeft.image[j,i]/binning**2))
        elif self.ImgState==2:
            #in-plane map
            self.p3.setTitle("h: %0.3f, k: %0.3f, dist. from origin: %0.3f Intensity: %0.3f " % (x, y,np.sqrt(x*x+y*y),self.imgLeft.image[j,i]))  
        elif self.ImgState==1:
            #reciprocal space (out-of-plane)
            self.p3.setTitle("hk: %0.3f, l: %0.3f, dist. from origin: %0.3f Intensity: %0.3f" % (x, y,np.sqrt(x*x+y*y),self.imgLeft.image[j,i]))

    def detectorMode(self):
        """Go back to pixel coordiante mode"""
        self.ImgState=0
        self.select_l_slice = True
        self.angleRegion.setMovable(True)  
        #self.hist.setImageItem(self.imgLeft)        
        #self.imgLeft.show()        
        self.imgLeft.resetTransform()
        self.updateRegion()
        self.ROICalc

        self.p3.getAxis('bottom').setScale(None)
        self.p3.getAxis('left').setScale(None)
        self.p3.getAxis('bottom').setGrid(0)
        self.p3.getAxis('left').setGrid(0)
        self.p3.getAxis('bottom').setLabel("PIXELS")       
        self.p3.getAxis('left').setLabel("PIXELS")        
        self.statusLabel.setText('Viewing in pixel coordinates')
        self.p3.getViewBox().setAspectLocked(False)
        self.p3.autoRange()

    def setupProjectionAxis(self):
        projection = self.p.param('Data Processing', 'Select Projection').value()

        if self.crystal.param('Reciprocal Units').value() == 0:
            if projection == 0:
                self.p3.getAxis('bottom').setLabel("Qx(Å⁻¹)")    
                self.p3.getAxis('left').setLabel("Qy (Å⁻¹)")  
            else:
                self.p3.getAxis('left').setLabel("Qz (Å⁻¹)")        
                if projection == 1:
                    self.p3.getAxis('bottom').setLabel("Qx(Å⁻¹)")  
                if projection == 2:
                    self.p3.getAxis('bottom').setLabel("Qy(Å⁻¹)")  
                if projection == 3:
                    self.p3.getAxis('bottom').setLabel("Sqrt(Qx^2 + Qy^2) (Å⁻¹)") 
        else:
            if projection == 0:
                self.p3.getAxis('left').setLabel("K (RLU)")
                self.p3.getAxis('bottom').setLabel("H (RLU)") 
            else:
                self.p3.getAxis('left').setLabel("L (RLU)")        
                if projection == 1:
                    self.p3.getAxis('bottom').setLabel("H (RLU)")  
                if projection == 2:
                    self.p3.getAxis('bottom').setLabel("K (RLU)") 
                if projection == 3:
                    self.p3.getAxis('bottom').setLabel("sqrt(H^2 + k^2) (RLU)") 

    def __projection_progress(self, n):        
        self.progressBar.setValue(n)

    def makeProjection(self):    
        
        """ Fnction for creating an in-plane map by generating coordinates for each pixel and place them in the correct histogram bin"""
        #first we check if it makes sense to do the calculation
        if self.busy == True:
            print('Error: Hat is busy and can not perform requested binning. This is probably because you are loading files')
            return            
        if self.image_stack.angle_mode == False:
            self.statusLabel.setText('Error: No angle information you probably need to read the log file or use manual mode')   
            return
        #we need some masks to make a projects
        if not self.mask_list:
            self.statusLabel.setText('Error: You must first select a mask')   
            return
        else:   
            projection = self.p.param('Data Processing', 'Select Projection').value()
            self.binBounds = []          
            self.convertMasks()
            if projection != 0 :
                binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
                second_det = self.p.param('Data Processing', 'Use 2nd detector').value()  
                self.pixMaps=self.experiment.dector_frame_to_hkl_frame(self,binning,second_det,0)

            if self.single_thread == False:
                worker = fileio.Worker(self.__makeProjection)
                worker.signals.finished.connect(self.__projection_complete)
                worker.signals.progress.connect(self.__projection_progress)
                self.busy = True
                self.threadpool.start(worker) 
            else:
                empty_callback = 0
                self.statusLabel.setText("SCRIPT MODE: Loading Files! (GUI DISABLED)") 
                self.__makeProjection(empty_callback)
                self.__projection_complete()
        
    def __projection_complete(self):        
        self.imgLeft.resetTransform()
        self.hist.setImageItem(self.imgLeft)
        self.p3.getViewBox().setAspectLocked(True)          
        self.ctrROI.hide()
        self.p4.hide()
        self.ImgState=2
        self.angleRegion.setMovable(False)              
        self.setupProjectionAxis()      
        grid_res = self.p.param('Data Processing', 'Grid Size').value()

        self.imgLeft.scale(2*self.qmax/grid_res,2*self.qmax/grid_res)
        self.imgLeft.setPos(-1*self.qmax,-1*self.qmax)    

        self.qRegion.hide()
        self.statusLabel.setText('Vewing in Reciprocal lattice units, angluar integration is fixed, press detector view (CTRL+1) to go back')   
        self.busy = False
        self.progressBar.setValue(self.progressBar.maximum())
        self.hideMasks()
        self.p3.autoRange()
        self.update()   
        

    def __makeProjection(self,progress_callback):
        """Intenral function for creating an in-plane map by generating coordinates for each pixel and place them in the correct histogram bin"""
        #now we should launch this in a seperate thread
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()  
        
        #angle information from the slider
        start_angle = self.angleRegion.getRegion()[0]
        end_angle = self.angleRegion.getRegion()[1]
        number_of_images =self.image_stack.number_of_images_in_range(start_angle,end_angle)        
        step_size = self.image_stack.get_step()

        self.progressBar.setMaximum(number_of_images) 
        #image numbers for our angular range
        from_image = int(np.floor(self.image_stack.angle2image(self.angleRegion.getRegion()[0])))
        to_image = int(np.ceil(self.image_stack.angle2image(self.angleRegion.getRegion()[1])))      

        #images we are interested in binning
        image_data=self.image_stack.image_data[from_image:to_image,::]   

        #interpolate the angles in our range
        angles=np.linspace(start_angle-step_size,end_angle+step_size,number_of_images)

        if self.p.param('Data Processing', 'Select Projection').value() == 0:
            qmax = self.qrmax+0.1
            if self.qzmax > self.qmax:
                qmax = self.qzmax+0.1  
        else:
            qmax = np.max(self.pixMaps[2])+0.1 #qz     

        self.qmax = qmax       

        #grid res is the number of bins in each direciton
        grid_res = self.p.param('Data Processing', 'Grid Size').value()
        gridstep = (2*qmax/grid_res) # x2 for negative values
        print("GRID STEP:",gridstep, grid_res)
        #these grids are coordinate grids that are used when selecting box profiles
        self.grid_qy = np.tile((np.arange(-1*qmax,qmax,gridstep)),(math.ceil(2*qmax/gridstep),1))
        print(np.shape(self.grid_qy))
        self.grid_qx = np.transpose(self.grid_qy)
        self.grid_qr = np.sqrt(self.grid_qx**2 + self.grid_qy**2)                
        
        #call the binning function
        grid_i = self.experiment.projection_binner(self,binning,second_det,angles,image_data,qmax,progress_callback)        
        
        #pass the grid intensities back to showData
        self.showData=grid_i
        self.imgLeft.setImage(grid_i,autoLevels=False)

     
    def detectorTransform(self):
        """Convert pixel coordinates to Reciprocal coordinates"""
        import time
        self.statusLabel.setText('Calculating reciprocal lattice coordinates (Please wait)..')
        self.update()
        self.ImgState=1
        self.imgLeft.resetTransform()
        self.qRegion.hide()
        self.angleRegion.setMovable(False) 
        self.select_l_slice = True
        self.showData=None

        self.p3.getAxis('bottom').setGrid(255)
        self.p3.getAxis('left').setGrid(255)     
        self.p3.getAxis('bottom').setZValue(1)        
        self.p3.getAxis('left').setZValue(1)

        if self.crystal.param('Reciprocal Units').value() == 0:
            self.p3.getAxis('bottom').setLabel("Qᵣ(Å⁻¹)")    
            self.p3.getAxis('left').setLabel("Qz (Å⁻¹)")         
        else:
            self.p3.getAxis('left').setLabel("L (RLU)")
            if self.crystal.param('β₃').value() == 90:
                   self.p3.getAxis('bottom').setLabel("sqrt((Hphi/b1)^2 + (Kphi)^2) (RLU)")   
            if self.crystal.param('b₁').value() != self.crystal.param('b₂').value():
                self.p3.getAxis('bottom').setLabel("sqrt(H^2 + (rK)^2) (RLU)")    
            else:
                self.p3.getAxis('bottom').setLabel("sqrt(H^2 + K^2) (RLU)")   
                      
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()  
        res=int(self.p.param('Data Processing', 'hkl Resolution').value())      

        #calculate the transformation[0]==h, [1]==k, [2]==l,[3]==sqrt(h^2 +k^2),[4]= intensity correciton
        self.pixMaps=self.experiment.dector_frame_to_hkl_frame(self,binning,second_det,0)
        #h_values = self.pixMaps[0].ravel()
        #k_values = self.pixMaps[1].ravel()
        l_values = self.pixMaps[2].ravel()
        hk_values = self.pixMaps[3].ravel()
        c_values = self.pixMaps[4]
      
        self.showData = self.image_stack.get_image(self.angleRegion.getRegion()[0], self.angleRegion.getRegion()[1])
   
        #create grid space for interpolation 
        self.grid_hk, self.grid_l=np.mgrid[np.min(hk_values):np.max(hk_values):res*1j,np.min(l_values):np.max(l_values):res*1j]        
        grid_i = interpolate.griddata((hk_values,l_values),np.ravel(np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='nearest')

        #do 2d interpolation, and apply intensity corrections if needed
        if self.p.param('Data Processing', 'Apply Intensity Corrections').value():
            #we unravel each pixmap to get lists of coordinates (hk, l, correction)
            grid_i = interpolate.griddata((hk_values,l_values),np.ravel(c_values*np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='nearest')
        else:
            grid_i = interpolate.griddata((hk_values,l_values),np.ravel(np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='nearest')           
            

        #grid_i[grid_i<=0]=np.nan
        #replace showData with the 
        self.showData=grid_i
        self.gridl = self.grid_l

        #set the image and scale the axis
        self.imgLeft.setImage(grid_i,autoLevels=False)
        #self.hist.setImageItem(self.imgLeft)
        self.imgLeft.scale((np.max(hk_values)-np.min(hk_values))/res,(np.max(l_values)-np.min(l_values))/res)
        self.imgLeft.setPos(np.min(hk_values),np.min(l_values))        
        self.statusLabel.setText('Viewing in reciprocal lattice units, angluar integration is fixed, press Detector View to go back')  
        self.showMasks()  

        #depending on the on the log intensity state, set ImgLeft to a logged showData
        if self.p.param('Data Processing', 'Log Intensity').value():
            self.logIntensity()  
        self.p3.autoRange()
        self.update()


    def updateRegion(self):
        """When the part of the image stack to show is changed we need to update the image displayed, this however is a little complicated
        since the images are discreate and not continous. Also angles need to be handled differently to image numbers. 
        The selection can be continous providing that it is not smaller than one image. 
       """       
 
        if self.image_stack.images_read:
            from_region=float(self.angleRegion.getRegion()[0])
            to_region=float(self.angleRegion.getRegion()[1])
            if self.image_stack.is_region_ok(from_region,to_region):
                self.showData = self.image_stack.get_image(from_region,to_region)                
                self.imgLeft.setImage(self.showData,autoLevels=False)  
                self.logIntensity() # check if we need to log intensity                   
            else:
                if self.image_stack.angle_mode:
                    if from_region == self.image_stack.start_angle:    
                        self.angleRegion.setRegion((from_region,to_region+self.image_stack.get_step()))
                    elif to_region == self.image_stack.end_angle:
                        self.angleRegion.setRegion((from_region-self.image_stack.get_step(),to_region))   
                    else:
                        #since we don't know if we can move the left of right handle we check which half we are in
                        if to_region < (self.image_stack.end_angle - self.image_stack.start_angle)/2:
                            self.angleRegion.setRegion((from_region,to_region+self.image_stack.get_step()/2))
                        else:
                            self.angleRegion.setRegion((from_region-self.image_stack.get_step()/2,to_region))  
                #so the box is less that one image in width
                else:
                    if from_region == 0:
                        self.angleRegion.setRegion((from_region,to_region+self.image_stack.get_step()/2))                                                
                    elif from_region == self.image_stack.number_of_images:
                        self.angleRegion.setRegion((from_region-self.image_stack.get_step()/2,to_region))
                    else:
                        #since we don't know if we can move the left of right handle we check which half we are in
                        if to_region < int(self.image_stack.number_of_images/2):
                            #left half so increase the right handle one half step
                            self.angleRegion.setRegion((from_region,to_region+self.image_stack.get_step()/2))  
                        else:
                            #right half so decrease the left handle one half step
                            self.angleRegion.setRegion((from_region-self.image_stack.get_step()/2,to_region)) 

                from_region=float(self.angleRegion.getRegion()[0])
                to_region=float(self.angleRegion.getRegion()[1])
                self.showData = self.image_stack.get_image(from_region,to_region)
                self.imgLeft.setImage(self.showData,autoLevels=False)   
                self.logIntensity() # check if we need to log intensity                

            self.p3.getAxis('bottom').setScale(None)
            self.p3.getAxis('left').setScale(None)
            self.p3.getAxis('bottom').setGrid(0)
            self.p3.getAxis('left').setGrid(0)                        
            self.ROICalc()
            if self.ImgState==0:           
                self.p3.getAxis('bottom').setLabel("PIXELS")       
                self.p3.getAxis('left').setLabel("PIXELS") 
                self.updateAngle()

    def setProfileLog(self, value):
        """Simple function to change file to make scripting cleaner, not used internally"""
        self.p.param('Profile tools', 'Log profile').setValue(value)

    def setBinning(self, value):
        """Simple function to change file to make scripting cleaner, not used internally"""
        self.p.param('Data Processing', 'Binning (X by X)').setValue(value)

    def enableProfile(self, value):
        """Simple function to change file to make scripting cleaner, not used internally"""
        self.p.param('Profile tools', 'Enable Box Profile').setValue(value)

    def setProfileAxis(self, value):
        """Simple function to change profile axis to make scripting cleaner, not used internally"""
        self.p.param('Profile tools', 'Axis of interest').setValue(value)
    
    def setProjection(self, value):
        """Simple function to change profile axis to make scripting cleaner, not used internally"""
        self.p.param('Data Processing', 'Select Projection').setValue(value)

    def setLevels(self, cmin, cmax):
        self.hist.setLevels(min=cmin, max=cmax)   

    def setColorMap(self, colormap):
        self.hist.gradient.loadPreset(colormap)

    def setAngleRange(self, start, end):
        if self.ImgState==0:
            self.angleRegion.setRegion((start,end))
        else:            
            self.statusLabel.setText('You can only change the region in detector mode!')     
            self.update()
    
    def setUnits(self, unit):
        self.crystal.param('Reciprocal Units').setValue(unit)
        
    def ROICalc(self):
        """Function to update the box profile display"""
        res=int(self.p.param('Data Processing', 'hkl Resolution').value())
        
        if self.p.param('Profile tools', 'Log profile').value():
            self.p4.setLogMode(x=None,y=True)
        else:
             self.p4.setLogMode(x=None,y=False)

        if not self.p.param('Profile tools', 'Enable Box Profile').value():
            return
        #hkl view
        if self.ImgState==1:
            if self.p.param('Profile tools', 'Axis of interest').value() == 1:
                grid_h, tmp=np.mgrid[np.min(self.pixMaps[0]):np.max(self.pixMaps[0]):res*1j,np.min(self.pixMaps[2]):np.max(self.pixMaps[2]):res*1j]
                xRoiData=self.ctrROI.getArrayRegion(grid_h, self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 2:
                grid_k, tmp=np.mgrid[np.min(self.pixMaps[1]):np.max(self.pixMaps[1]):res*1j,np.min(self.pixMaps[2]):np.max(self.pixMaps[2]):res*1j]
                xRoiData=self.ctrROI.getArrayRegion(grid_k, self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 3:
                xRoiData=self.ctrROI.getArrayRegion(self.grid_hk, self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 4:
                xRoiData=self.ctrROI.getArrayRegion(self.grid_l, self.imgLeft)

            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft)  
            xdata=np.mean(xRoiData,axis=1)
            self.roixdata = xdata
            #we take the average across the box rather than just the single max pixel
            #this allows for background subtraction.             
            self.roiydata = np.mean(ROIData,axis=1)                
            self.CTRCurve.setData(x=xdata,y=self.roiydata)

        #projectionview 
        elif self.ImgState==2:
            projection = self.p.param('Data Processing', 'Select Projection').value()
            
            #axis = H
            if self.p.param('Profile tools', 'Axis of interest').value() == 1:
                #h-k
                if projection == 0:                
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #h-l
                elif projection == 1:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #k-l    
                elif projection == 2:
                    self.statusLabel.setText('This projection does not calculate that axis, using qy/k instead')     
                    self.update()
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    self.statusLabel.setText('This projection does not calculate that axis, using qr/(sqrt(h^2+k^2)) instead')     
                    self.update()
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)

            #axis = K     
            elif self.p.param('Profile tools', 'Axis of interest').value() == 2:
                #h-k
                if projection == 0:                                 
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qy, self.imgLeft)
                #h-l
                elif projection == 1:
                    self.statusLabel.setText("This projection does not calculate that axis, using qx/h instead") 
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #k-l    
                elif projection == 2:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    self.statusLabel.setText('This projection does not calculate that axis, using qr/(sqrt(h^2+k^2)) instead')     
                    self.update()                    
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)

            #axis = qr
            elif self.p.param('Profile tools', 'Axis of interest').value() == 3:
                #h-k
                if projection == 0:                          
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qr, self.imgLeft)
                #h-l
                elif projection == 1: 
                    self.statusLabel.setText('This axis could be misleading showing h/qx instead')     
                    self.update()    
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #k-l    
                elif projection == 2: 
                    self.statusLabel.setText('This axis could be misleading showing k/qx instead')     
                    self.update()    
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)          
            #axis = L
            elif self.p.param('Profile tools', 'Axis of interest').value() == 4:
                #h-k
                if projection == 0:                       
                    self.statusLabel.setText('This projection does not calculate that axis, using qx/h instead')     
                    self.update()                         
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qx, self.imgLeft)
                #h-l
                elif projection == 1:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qy, self.imgLeft)
                #k-l    
                elif projection == 2:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qy, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qy, self.imgLeft)   

    
            
            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft)  
            xdata=np.mean(xRoiData,axis=1)
            self.roixdata = xdata
            self.roiydata = np.mean(ROIData,axis=1)
            self.CTRCurve.setData(x=xdata,y=self.roiydata)#-np.min(ROIData[0],axis=1))
            
        #pixel view
        else:      
            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft)     
            self.roiydata = np.mean(ROIData,axis=1)                
            self.roixdata = range(0,len(np.mean(ROIData,axis=1)))  
            self.CTRCurve.setData(self.roiydata)



    def data_format_changed(self):
        """We update the interface if a beamline preset is chosen"""
        #ID31 Beamline
        if self.experiment.param('Beamline Preset').value() == 4:
            #self.p.param('File', 'Select Dark Image').hide()
            #self.p.param('File', 'Select Background Dark Image').hide()
            #self.p.param('File', 'Select Data Images').setName("Select HDF5 file")
            #self.p.param('File', 'Select Log File').hide()
            self.p.param('Data Processing', 'Use 2nd detector').setValue(False)
            self.p.param('Data Processing', 'Use 2nd detector').hide()           
        #P07 Beamline
        if self.experiment.param('Beamline Preset').value() == 3:
            #self.p.param('File', 'Select Log File').hide()
            self.p.param('Experiment', 'Y Pixels').setValue(2048)
            self.p.param('Experiment', 'X Pixels').setValue(2048)
            self.p.param('Experiment', 'Pixel Size').setValue(200e-6)
            self.p.param('Experiment', 'Energy').setValue(73700)
            self.p.param('Experiment', 'Sample-Detector Dist.').setValue(1.6)
            #There was no 2nd detector and dark images subtracted normaly 
            self.p.param('Data Processing', 'Use 2nd detector').hide()
            #self.p.param('File', 'Select Dark Image').hide()
            #self.p.param('File', 'Select Background Dark Image').hide()   

        if self.experiment.param('Beamline Preset').value() == 6:
            #self.p.param('File', 'Select Log File').hide()
            self.p.param('Experiment', 'Y Pixels').setValue(2880)
            self.p.param('Experiment', 'X Pixels').setValue(2880)
            self.p.param('Experiment', 'Pixel Size').setValue(150e-6)
            self.p.param('Experiment', 'Energy').setValue(73700)
            self.p.param('Experiment', 'Sample-Detector Dist.').setValue(2.1)
            #There was no 2nd detector and dark images subtracted normaly 
            self.p.param('Data Processing', 'Use 2nd detector').hide()
            #self.p.param('File', 'Select Dark Image').hide()
            #self.p.param('File', 'Select Background Dark Image').hide()   

        #P21.2 Beamline
        if self.experiment.param('Beamline Preset').value() == 1 or self.experiment.param('Beamline Preset').value() == 2:  
            self.p.param('Experiment', 'Y Pixels').setValue(2880)
            self.p.param('Experiment', 'X Pixels').setValue(2880)
            self.p.param('Experiment', 'Pixel Size').setValue(150e-6)

        if self.experiment.param('Beamline Preset').value() == 5:
            self.experiment.param('Manual Start Angle').show()
            self.experiment.param('Manual End Angle').show()
            #self.p.param('File', 'Select Dark Image').show()
            #self.p.param('Data Processing', 'Use 2nd detector').show()  
            #self.p.param('File', 'Select Background Dark Image').show()
            #self.p.param('File', 'Select Log File').hide()
        else:
            self.experiment.param('Manual Start Angle').hide()
            self.experiment.param('Manual End Angle').hide()


    
    def sample_preset(self):
        """Apply presets for sample if changed, eventually this should be read from a text file"""
        #Au 111 (surface units)
        if self.p.param('Crystal','Preset').value() == 1:
            self.p.param('Crystal', 'a₁').setValue(2.885)
            self.p.param('Crystal', 'a₂').setValue(2.885)
            self.p.param('Crystal', 'a₃').setValue(7.064)

            self.p.param('Crystal', 'α₁').setValue(90)
            self.p.param('Crystal', 'α₂').setValue(90)
            self.p.param('Crystal', 'α₃').setValue(120)
        #Au 100 (surface units)
        if self.p.param('Crystal','Preset').value() == 2:
            self.p.param('Crystal', 'a₁').setValue(2.88)
            self.p.param('Crystal', 'a₂').setValue(2.88)
            self.p.param('Crystal', 'a₃').setValue(4.08)

            self.p.param('Crystal', 'α₁').setValue(90)
            self.p.param('Crystal', 'α₂').setValue(90)
            self.p.param('Crystal', 'α₃').setValue(90)
        #TiO2
        if self.p.param('Crystal','Preset').value() == 3:
            self.p.param('Crystal', 'a₁').setValue(6.496)
            self.p.param('Crystal', 'a₂').setValue(2.959)
            self.p.param('Crystal', 'a₃').setValue(6.496)

            self.p.param('Crystal', 'α₁').setValue(90)
            self.p.param('Crystal', 'α₂').setValue(90)
            self.p.param('Crystal', 'α₃').setValue(90)


    def saveprofile(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","csv (*.csv);;All Files (*)", options=options)  
        data = np.asarray([self.roixdata,self.roiydata])
        np.savetxt(filename,np.transpose(data),fmt='%10.5f', delimiter=',',newline='\n')       

    def saverocks(self, paramHandle, folderName=None):
        if type(paramHandle) == str:
            folderName =  paramHandle
        if not folderName:
            folderName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        #get the image indices to use
        if self.image_stack.angle_mode == False:
            from_image=int(np.round(self.angleRegion.getRegion()[0]))
            to_image=int(np.round(self.angleRegion.getRegion()[1]))
        else:
            from_image = int(np.floor(self.image_stack.angle2image(self.angleRegion.getRegion()[0])))
            to_image = int(np.ceil(self.image_stack.angle2image(self.angleRegion.getRegion()[1])))

        #images in our angular range
        tmp = self.image_stack.image_data[from_image:to_image+1,:,:]
        profiles = []
        angles = []        

        for i,image in enumerate(tmp):
            ROIData=self.ctrROI.getArrayRegion(image, self.imgLeft) 
            profiles.append(np.sum(ROIData,axis=1))
            angles.append(self.image_stack.start_angle+i*self.image_stack.get_step())

        profiles = np.transpose(np.asarray(profiles))
        res=int(self.p.param('Data Processing', 'hkl Resolution').value())
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())      
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()                     
        self.pixMaps=self.experiment.dector_frame_to_hkl_frame(self,binning,second_det,0)
        

        xRoiData=self.ctrROI.getArrayRegion(np.rot90(self.pixMaps[2],3), self.imgLeft)
        xRoiData2=self.ctrROI.getArrayRegion(np.rot90(self.pixMaps[4],3), self.imgLeft)
        qz = np.mean(xRoiData,axis=1)  
        corr = np.mean(xRoiData2,axis=1) 
        np.savetxt(folderName+'/axis.csv',[qz,corr], delimiter=',',newline='\n')         
        for i in range(len(profiles)):
            fileName = str(folderName)+'/'+str(i)+'.csv'
            tmp2 = profiles[i]
            data = np.transpose([angles,tmp2])  
            np.savetxt(fileName,data, delimiter=',',newline='\n') 

        text = str(i)+" rocking scans saved in folder: "+ folderName
        self.statusLabel.setText(text)     
        self.update()               
         
    def saveroi(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","roi (*.roi);;All Files (*)", options=options)     
        if filename:
            state = [self.ctrROI.saveState(),self.angleRegion.getRegion()]
            pickle.dump(state, open(filename, "wb" ))   

    def loadroi(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Select a ROI file too use", "","ROI File (*.roi);;All Files (*)", options=options)
        
        if filename:
            state = pickle.load(open(filename, "rb" ))
            self.ctrROI.setState(state[0])
            #self.angleRegion.setRegion(state[1])    

    def load(self):
        if self.image_stack.images_read:
            self.image_stack.load_aux_files()
            self.image_stack.load_images() 

    def loadUpdate(self, n):
        self.angleRegion.setMovable(0)
        self.showData = self.image_stack.get_image(0,n)
        self.imgLeft.setImage(self.showData, autoLevels=False)
        self.imgLeft.show() 

    def loadComplete(self):        
        if self.image_stack.images_read:    
            if self.image_stack.angle_mode:
                self.angleRegion.setMovable(1)  
                self.angleRegion.setRegion((self.image_stack.start_angle,self.image_stack.end_angle))
                self.angleRegion.setBounds((self.image_stack.start_angle,self.image_stack.end_angle))                
                self.showData = self.image_stack.get_image(self.image_stack.start_angle,self.image_stack.end_angle)

                self.imgLeft.setImage(self.showData, autoLevels=False)
                self.imgLeft.show()                    
                #self.hist.setImageItem(self.imgLeft)     
                self.p2.setLabel('bottom',text='Angle (degrees)')     
                
            else:
                self.angleRegion.setMovable(1)  
                self.angleRegion.setRegion((1,self.image_stack.number_of_images))
                self.angleRegion.setBounds((1,self.image_stack.number_of_images))
                self.showData = self.image_stack.get_image(0,self.image_stack.number_of_images)
                self.imgLeft.setImage(self.showData, autoLevels=False)
                self.imgLeft.show()                    
                self.hist.setImageItem(self.imgLeft)       
            
        else:  
            self.statusLabel.setText('Background that will be subtracted. Select images to load')
            self.image_stack.load_aux_files()
            self.imgLeft.setImage(self.image_stack.subtract_image, autoLevels=False)
            self.imgLeft.show()    
            self.hist.setImageItem(self.imgLeft,)   
            self.ImgState = 5

    def makeFigure(self, paramHandle, filename=None):
        projection = self.p.param('Data Processing', 'Select Projection').value()
        if type(paramHandle) == str:
            filename =  paramHandle
        """This function reloads plotting.py and runs the correct plotting function to
        generate a figure with matplotlib"""
        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","png (*.png);;All Files (*)", options=options)     

        if filename:
            importlib.reload(plotting)
            if self.ImgState == 2:
                cmin,cmax=self.imgLeft.getLevels()
                self.grid_qy 
                if projection == 0:
                    plotting.plot_projection_hk(self, self.grid_qx,self.grid_qy,self.showData,cmin,cmax, filename)
                    text = "In-plane qx/qy map saved as: "+str(filename)
                    self.statusLabel.setText(text)
                    
                if projection == 1:
                    plotting.plot_projection_hl(self, self.grid_qx,self.grid_qy,self.showData,cmin,cmax, filename)
                    text = "qx/qz map saved as: "+ str(filename)
                    self.statusLabel.setText(text)
                if projection == 2:
                    plotting.plot_projection_kl(self, self.grid_qx,self.grid_qy,self.showData,cmin,cmax, filename)
                    text = "qy/qz map saved as: "+str(filename)
                    self.statusLabel.setText(text)
                if projection == 3:
                    plotting.plot_projection_qrl(self, self.grid_qx,self.grid_qy,self.showData,cmin,cmax, filename)
                    text="qr/qz map saved as: "+str(filename)
                    self.statusLabel.setText(text)
            if self.ImgState == 1:
                cmin,cmax=self.imgLeft.getLevels()
                plotting.plot_transformed_detector(self.grid_hk,self.grid_l,self.showData,cmin,cmax, filename)
                text = "Transfomred Detector View saved as: "+str(filename)
                self.statusLabel.setText(text)
            if self.ImgState == 0:
                cmin,cmax=self.imgLeft.getLevels()
                plotting.plot_out_of_plane(self.showData,cmin,cmax, filename)
                text = "Image saved as: "+str(filename)
                self.statusLabel.setText(text)

               
    def runScript(self):
        importlib.reload(script)
        self.single_thread = True
        text = "Running Script"
        self.update()
        script.script_main(self)
        text = "Script Finished"
        self.statusLabel.setText(text) 
        self.single_thread = False

    def addMask(self): 
        #when the button is pressed add a new ROI  
        i = len(self.mask_list)     
        self.mask_list.append(pg.RectROI([0, 5], [1, 1], pen=(i,9)))
        self.mask_list[i].show()
        self.p3.addItem(self.mask_list[i]) 

    def convertMasks(self):
        #this converts the ROIs into bounds for the binning algorithm
        self.qrmax = 0
        self.qzmax = 0
        for mask in self.mask_list:
            xmin = mask.pos()[0]
            xmax = mask.size()[0] + xmin
            ymin = mask.pos()[1]
            ymax = mask.size()[1] + ymin
            self.binBounds.append([xmin,xmax,ymin,ymax])  
            if xmax > self.qrmax:
                self.qrmax = xmax
            if abs(xmin) > self.qrmax:
                self.qrmax = abs(xmin)
            if ymax > self.qzmax:
                self.qzmax = ymax
            if abs(ymin) > self.qzmax:
                self.qzmax = abs(ymax)

    def clearMasks(self):
        self.binBounds = []   
        for mask in self.mask_list:
            self.p3.removeItem(mask)
        self.mask_list = []
    
    def hideMasks(self): 
        for mask in self.mask_list:
            self.p3.removeItem(mask)
    
    def showMasks(self): 
        for mask in self.mask_list:
            self.p3.addItem(mask)

    def saveMasks(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","msk (*.msk);;All Files (*)", options=options)     
        if filename:
            mask_states = []
            for mask in self.mask_list:
                mask_states.append(mask.saveState())
            pickle.dump(mask_states, open( filename, "wb" ))  
            
             
    def loadMasks(self, paramHandle, filename=None): 
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:       
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Select a Mask file too use", "","Mask File (*.msk);;All Files (*)", options=options)
        if filename:
            mask_states = pickle.load(open(filename, "rb" ))
            for mask_state in mask_states:
                self.mask_list.append(pg.RectROI([0, 5], [1, 1], pen=(len(mask_states),9)))
                self.mask_list[-1].setState(mask_state)
                self.mask_list[-1].show()
                self.p3.addItem(self.mask_list[-1])
        self.convertMasks()            
      


    def save(self):
        state = self.p.saveState()
        pickle.dump(state, open( 'params.par', "wb" ))      
        pickle.dump(self.image_stack.save_state(),open('imagestack.bin', "wb" ) )           
        self.statusLabel.setText("Parameters and file selection saved")    

    def restore(self):
        self.p.restoreState(pickle.load(open( 'params.par', "rb" )),removeChildren=False)
        self.statusLabel.setText("Saved parameters restored") 

    def importParameters(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:       
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Select a parameter file to load", "","Binary Parameter *.par);;All Files (*)", options=options)
        if filename:
            self.p.restoreState(pickle.load(open(filename, "rb" )),removeChildren=False)
            self.statusLabel.setText("Saved parameters restored") 

    def importImageStack(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:       
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getOpenFileName(self,"Select an Image Stack to Load", "","Binary File (*.bin);;All Files (*)", options=options)
        if filename:        
            self.image_stack.restore_state(pickle.load(open(filename, "rb" )))   
            self.image_stack.load_aux_files()       
            text = str(self.image_stack.number_of_images)+" images selected with flags"+str(self.image_stack.flags)+"  Press Load (CTRL+L)"            
            self.statusLabel.setText(text) 

    def exportParameters(self, paramHandle, filename=None):
        if type(paramHandle) == str:
                filename =  paramHandle
        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","par (*.par);;All Files (*)", options=options)     
        if filename:
            state = self.p.saveState()
            pickle.dump(state, open(filename, "wb" ))
            self.statusLabel.setText("Parameters exported to "+filename)             

    def exportImageStack(self, paramHandle, filename=None):
        if type(paramHandle) == str:
                filename =  paramHandle
        if not filename:
            options = QFileDialog.Options()
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","bin (*.bin);;All Files (*)", options=options)     
        if filename:
            pickle.dump(self.image_stack.save_state(),open(filename, "wb" ) ) 
            self.statusLabel.setText("File selection exported to "+filename)             

    def restoreImageStack(self):        
        self.image_stack.restore_state(pickle.load(open( 'imagestack.bin', "rb" )))   
        self.image_stack.load_aux_files()           
        text = str(self.image_stack.number_of_images)+" images selected with flags"+str(self.image_stack.flags)+"  Press Load (CTRL+L)"            
        self.statusLabel.setText(text) 

    def reloadWarning(self):
        if self.image_stack.image_file_names_1 is not None:
            self.statusLabel.setText("For this change to take affect you need to reaload the files (CTRL+L)") 
            self.warnDialog("Manually type the whole number. For this change to take affect you need to reaload the files (CTRL+L)")
    
    def updateCrange(self):
        cmin, cmax = self.hist.getLevels()
        self.cminbox.setText(str(cmin))
        self.cmaxbox.setText(str(cmax))

    def updateAngle(self):
        angleMin, AngleMax = self.angleRegion.getRegion()
        self.startangle.setText(str(angleMin))
        self.endangle.setText(str(AngleMax ))

    def updateAngleManual(self):
        start = atof(self.startangle.text())
        end = atof(self.endangle.text())  
        self.angleRegion.setRegion((start,end))

    def updateCrangeManual(self):
        cmin = atof(self.cminbox.text())
        cmax = atof(self.cmaxbox.text())      
        self.setLevels(cmin,cmax)     
        self.updateRegion()

    def logIntensity(self):
        '''Function to show the logged intensity'''         
        if self.p.param('Data Processing', 'Log Intensity').value(): 
            if self.log_active == False:
                self.log_active = True
                #if there are any negative numbers we should shift the intensities
                min_value = np.min(self.showData) 
                cmin, cmax = self.hist.getLevels()     
      
                if min_value <= 0:
                    self.statusLabel.setText("Negative value in log detected shifting intensities")
                    self.cshift = abs(min_value)+1
                    logdata = np.log10(self.showData+self.cshift)
                    self.setLevels(np.log10(cmin+self.cshift),np.log10(cmax+self.cshift)) 

                else:
                    logdata = np.log10(self.showData)
                    self.setLevels(np.log10(cmin),np.log10(cmax)) 
            else:
                logdata = np.log10(self.showData+self.cshift)

            self.imgLeft.setImage(logdata, autoLevels=False)
        else:
            if self.log_active:
                cmin, cmax = self.hist.getLevels()
                self.imgLeft.setImage(self.showData, autoLevels=False)
                self.setLevels(10**(cmin)-self.cshift,10**(cmax)-self.cshift) 
                self.cshift = 0
                self.log_active = False
                self.imgLeft.show()   

    def toggleBg(self):            
            if self.white_background:
                self.win.setBackground('k')
                self.white_background = False
            else:
                self.win.setBackground('w')   
                self.white_background = True    

    def warnDialog(self,text):
        d = QDialog()
        d.setWindowTitle("Warning")
        d.setFixedSize(550,50)
        t = QLabel(d)
        t.setText(text)
        t.move(20,2)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()
                
    def aboutHelp(self):
        """Function to display and about help dialog"""
        d = QDialog()
        d.setWindowTitle("About")
        d.setFixedSize(500,160)

        t1 = QLabel("HESXRD Analysis Toolkit.\n\n(c) 2021 Gary Harlow. \n\nFree to use under MIT License\n\nCreated by Gary Harlow and Sebastian Pfaff at Lund University, Sweden ",d)        
        t1.setWordWrap(True)

        t2 = QLabel(d)             
        t2.setText('<a href="https://github.com/gary-harlow/HESXRD-Analysis-Toolkit">Goto GitHuB Page</a>')
        t2.setOpenExternalLinks(True)

        t3 = QLabel(d)
        t3.setText('Eventually a link to a readthedocs.io page will come here')

        t1.move(20,2)
        t2.move(20,100)
        t3.move(20,130)
        d.setWindowModality(Qt.ApplicationModal)
        d.exec_()

    def exithat(self):
        print("Goodbye. Have a nice day!")
        quit()

                              
    def initUI(self):
        self.image_stack=ImageStack(self)
        self.setWindowTitle(self.title)
        gridLayout = QGridLayout(self)     
        self.setLayout(gridLayout)   
        self.single_thread = False    
        self.cshift = 0 #for a potential color shift when taking log
        self.log_active = False
        self.qmax =0
        self.white_background = False
        
            
        params = [
            self.experiment,
            self.crystal,
            {'name': 'Data Processing', 'type': 'group', 'children': [
                {'name': 'Binning (X by X)', 'type': 'int', 'value': 4, 'step': 1},
                {'name': 'hkl Resolution', 'type': 'int', 'value': 1000, 'step': 1},
                {'name': 'Grid Size', 'type': 'int', 'value': 800, 'step': 1},
                {'name': 'Select Projection', 'type': 'list', 'values': {"h-k": 0, "h-l": 1, "k-l": 2,"h^2+k^2 - l": 3},'value': 0},               
                {'name': 'Mean Images Instead of Max', 'type': 'bool','valie':False},
                {'name': 'Acceleration', 'type': 'list', 'values': {"none": 0, "numba (cpu)": 1, "cuda (gpu)": 2},'value': 1},                    
                {'name': 'Bin From Full Images', 'type': 'bool', 'value': False},             
                {'name': 'Apply Intensity Corrections', 'type': 'bool', 'value': False},
                {'name': 'Correct for Refraction', 'type': 'bool', 'value': False},
                {'name': 'Intensity Offset', 'type': 'float', 'value': 0, 'step':10},
                {'name': 'Multiply intensity by', 'type': 'float', 'value': 1, 'step':0.1},
                {'name': 'Perform Background Subtraction', 'type': 'bool', 'value': False},
                {'name': 'Log Intensity', 'type': 'bool', 'value': False},
                {'name': 'Image Rotation', 'type': 'list', 'values': {"0 degrees": 0, "90 degrees": 1, "180 degrees": 2,"270 degrees": 3},'value': 0},
                {'name': 'Image Flip U/D', 'type': 'list', 'values': {"True": True,"False": False},'value': False},
                {'name': 'Image Flip L/R', 'type': 'list', 'values': {"True": True,"False": False},'value': False},
                {'name': 'Use 2nd detector', 'type': 'bool', 'value': False},      
            ]},
            {'name': 'Profile tools', 'type': 'group', 'children': [
                {'name': 'Enable Box Profile', 'type': 'bool',},
                {'name': 'Log profile', 'type': 'bool', 'value': True},
                {'name': 'Extract profile', 'type': 'action'},
                {'name': 'Save Rocking Scans', 'type': 'action'},              
                {'name': 'Axis of interest', 'type': 'list',  'values': {"H (Qx)": 1, "K (Qy)": 2, "HK (Qr)": 3,"L (qz)": 4}, 'value': 3} ,
                {'name': 'Save ROI', 'type': 'action'},
                {'name': 'Load ROI', 'type': 'action'},
            ]}
        ]


        #p.param('Save/Restore functionality', 'Save State').sigActivated.connect(save)


        self.p = Parameter.create(name='params', type='group', children=params)

        
        #Menu Bar
        menubar = QMenuBar()       
        gridLayout.addWidget(menubar, 0, 0)
        actionFile = menubar.addMenu("&File")

        loadAct = QAction('Load', self)
        loadAct.setShortcut('Ctrl+L')
        loadAct.setStatusTip('Load currently selected files')
        loadAct.triggered.connect(self.load)
        actionFile.addAction(loadAct)

        qsaveAct = QAction('Quick Save', self)
        qsaveAct.setShortcut('Ctrl+S')
        qsaveAct.setStatusTip('Quick Save Parameters and File')
        qsaveAct.triggered.connect(self.save)
        actionFile.addAction(qsaveAct)

        qparmAct = QAction('Quick Restore Parameters', self)
        qparmAct.setShortcut('Ctrl+P')
        qparmAct.setStatusTip('Quick Restore Parameters previously quick saved')
        qparmAct.triggered.connect(self.restore)
        actionFile.addAction(qparmAct)

        qfileAct = QAction('Quick Restore Image Stack', self)
        qfileAct.setShortcut('Ctrl+R')
        qfileAct.setStatusTip('Quick Restore Files previously quick saved')
        qfileAct.triggered.connect(self.restoreImageStack)
        actionFile.addAction(qfileAct)

        actionFile.addSeparator()

        selectImagesAct = QAction('Select Data Images', self)
        selectImagesAct.setShortcut('Ctrl+I')
        selectImagesAct.setStatusTip('Select Image Files to be Loaded')
        selectImagesAct.triggered.connect(self.image_stack.select_images)
        actionFile.addAction(selectImagesAct)

        selectBkImagesAct = QAction('Select Background Image', self)
        selectBkImagesAct.setShortcut('Ctrl+B')
        selectBkImagesAct.setStatusTip('Select background image files to be subtracted from each image (optional)')
        selectBkImagesAct.triggered.connect(self.image_stack.select_background)
        actionFile.addAction(selectBkImagesAct)

        selectLogAct = QAction('Select Log File', self)
        #selectBkImagesAct.setShortcut('Ctrl+B')
        selectLogAct.setStatusTip('Select a supported log file containing angles and monitor values')
        selectLogAct.triggered.connect(self.image_stack.read_log_file)
        actionFile.addAction(selectLogAct)

        actionFile.addSeparator()

        selectDarkAct = QAction('Select Dark Image', self)
        selectDarkAct.setShortcut('Ctrl+D')
        selectDarkAct.setStatusTip('Chose a dark image to be subtracted (optional)')
        selectDarkAct.triggered.connect(self.image_stack.select_dark)
        actionFile.addAction(selectDarkAct)

        selectBkDarkAct = QAction('Select Background Dark Image', self)
        selectBkDarkAct.setStatusTip('Chose a dark image to be subtracted from any background image(optional)')
        selectBkDarkAct.triggered.connect(self.image_stack.select_background_dark)
        actionFile.addAction(selectBkDarkAct)

        actionFile.addSeparator()

        exportParmAct = QAction('Export Current Parameters', self)
        exportParmAct.setStatusTip('Save current parameters to a file')
        exportParmAct.triggered.connect(self.exportParameters)
        actionFile.addAction(exportParmAct)

        exportFileAct = QAction('Export File Selection', self)
        exportFileAct.setStatusTip('Save list of currently selected dark images,files, and background to a file')
        exportFileAct.triggered.connect(self.exportImageStack)
        actionFile.addAction(exportFileAct)

        importParmAct = QAction('Export Current Parameters', self)
        importParmAct.setStatusTip('Save current parameters to a file')
        importParmAct.triggered.connect(self.importParameters)
        actionFile.addAction(importParmAct)

        importFileAct = QAction('Export File Selection', self)
        importFileAct.setStatusTip('Save list of currently selected dark images,files, and background to a file')
        importFileAct.triggered.connect(self.importImageStack)
        actionFile.addAction(importFileAct)

        actionFile.addSeparator()
        exitAct = QAction('Quit', self)
        exitAct.setShortcut('Ctrl+Q')
        exitAct.setStatusTip('Exit application')
        exitAct.triggered.connect(self.exithat)
        actionFile.addAction(exitAct)
        
        actionView = menubar.addMenu("&View")

        detViewAct = QAction('Detector View', self)
        detViewAct.setShortcut('Ctrl+1')
        detViewAct.setStatusTip('View the image stack in detector mode')
        detViewAct.triggered.connect(self.detectorMode)
        actionView.addAction(detViewAct)

        transViewAct = QAction('Transformed Detector View', self)
        transViewAct.setShortcut('Ctrl+2')
        transViewAct.setStatusTip('Summed or averaged image stack is corrected for curvature of the Ewald Sphere and given coordinates')
        transViewAct.triggered.connect(self.detectorTransform)
        actionView.addAction(transViewAct)
        
        binViewAct = QAction('Binned Projection', self)
        binViewAct.setShortcut('Ctrl+3')
        binViewAct.setStatusTip('Use all pixels across the selected angular range to create a projection')
        binViewAct.triggered.connect(self.makeProjection)
        actionView.addAction(binViewAct)

        actionView .addSeparator()
        toggleBackgroundColour = QAction('Toggle Background Color (white/black)', self)
        toggleBackgroundColour.triggered.connect(self.toggleBg)
        actionView.addAction(toggleBackgroundColour)


        
        actionRun = menubar.addMenu("&Script")

        scriptAct = QAction('Run script.py', self)
        scriptAct.setShortcut('F5')
        scriptAct.setStatusTip('Run the file called script.py')
        scriptAct.triggered.connect(self.runScript)
        actionRun.addAction(scriptAct)

        plotAct = QAction('Run plotting.py', self)
        plotAct.setShortcut('F6')
        plotAct.setStatusTip('Use matplotlib to plot the current view as defined in plotting.py')
        plotAct.triggered.connect(self.makeFigure)
        actionRun.addAction(plotAct)

        actionHelp = menubar.addMenu("&Help")
        
        actionAbout = QAction('About', self)
        actionAbout.setShortcut('F1')
        actionAbout.triggered.connect(self.aboutHelp)
        actionHelp.addAction(actionAbout)
        
        #MASK Toolbar
        toolbar = QToolBar()
        gridLayout.addWidget(toolbar, 0, 1)
        text = QLabel("Mask Tools")
        toolbar.addWidget(text)   
        
        addMask = toolbar.addAction("Add Mask")
        addMask.setIcon(QIcon("icons/addmask.svg"))
        addMask.triggered.connect(self.addMask)
        
        saveMask = toolbar.addAction("Save Masks")
        saveMask.setIcon(QIcon("icons/savemask.svg"))
        saveMask.triggered.connect(self.saveMasks)
        
        loadMask = toolbar.addAction("Load Masks")
        loadMask.setIcon(QIcon("icons/loadmask.svg"))
        loadMask.triggered.connect(self.loadMasks)
        
        clearMask = toolbar.addAction("Clear Masks")
        clearMask.setIcon(QIcon("icons/clearmask.svg"))
        clearMask.triggered.connect(self.clearMasks)
        
        toolbar1 = QToolBar()
        gridLayout.addWidget(toolbar1, 2, 0)
        onlydub = QtGui.QDoubleValidator()
        
        text3 = QLabel("  Start Angle ")
        toolbar1.addWidget(text3)        
        self.startangle = QLineEdit()
        self.startangle.setValidator(onlydub)
        startanglewidget = toolbar1.addWidget(self.startangle)
        self.startangle.editingFinished.connect(self.updateAngleManual)
        
        text4 = QLabel("  End Angle ")
        toolbar1.addWidget(text4)        
        self.endangle = QLineEdit()
        self.endangle.setValidator(onlydub)
        endanglewidget = toolbar1.addWidget(self.endangle)
        self.endangle.editingFinished.connect(self.updateAngleManual)        
        
        text1 = QLabel("  Cmin ")
        toolbar1.addWidget(text1)        
        self.cminbox = QLineEdit()
        self.cminbox.setValidator(onlydub)
        cminwidget = toolbar1.addWidget(self.cminbox)
        self.cminbox.editingFinished.connect(self.updateCrangeManual)
        
        
        text2 = QLabel("  Cmax ")
        toolbar1.addWidget(text2)        
        self.cmaxbox = QLineEdit()
        self.cmaxbox.setValidator(onlydub)
        cmaxwidget = toolbar1.addWidget(self.cmaxbox)   
        self.cmaxbox.editingFinished.connect(self.updateCrangeManual)            
        
        self.experiment.param('Manual Start Angle').sigValueChanged.connect(self.reloadWarning)
        self.experiment.param('Manual End Angle').sigValueChanged.connect(self.reloadWarning)
                
        self.p.param('Data Processing', 'Log Intensity').sigValueChanged.connect(self.logIntensity)
        self.p.param('Data Processing', 'Perform Background Subtraction').sigValueChanged.connect(self.updateRegion)
        self.p.param('Data Processing', 'Mean Images Instead of Max').sigValueChanged.connect(self.updateRegion)   
        self.p.param('Data Processing', 'Multiply intensity by').sigValueChanged.connect(self.updateRegion) 
        self.p.param('Data Processing', 'Intensity Offset').sigValueChanged.connect(self.updateRegion) 

        self.p.param('Profile tools', 'Extract profile').sigActivated.connect(self.saveprofile)
        self.p.param('Profile tools', 'Save Rocking Scans').sigActivated.connect(self.saverocks)
        self.p.param('Profile tools', 'Log profile').sigValueChanged.connect(self.ROICalc)  
        self.p.param('Profile tools', 'Axis of interest').sigValueChanged.connect(self.ROICalc)  
        self.p.param('Profile tools', 'Save ROI').sigActivated.connect(self.saveroi)
        self.p.param('Profile tools', 'Load ROI').sigActivated.connect(self.loadroi)     

        self.experiment.param('Beamline Preset').sigValueChanged.connect(self.data_format_changed)
        self.crystal.param('Preset').sigValueChanged.connect(self.sample_preset)

        t = ParameterTree()
        t.setMaximumWidth(320)
        t.setParameters(self.p, showTop=False)

        gridLayout.addWidget(t,1,1)
        win = pg.GraphicsLayoutWidget()
        self.win = win
        gridLayout.addWidget(win,1,0)

        self.progressBar=QProgressBar()
        gridLayout.addWidget(self.progressBar,3,0)

        self.statusLabel=QLabel('Select Images ... (Press Load Saved Files to load the files you had when you last pressed \'Save Parameters\')')
        gridLayout.addWidget(self.statusLabel,3,0)

        win.nextRow()
        self.p2 = win.addPlot(row=4,colspan=3)
        self.angleRegion=pg.LinearRegionItem()
        self.angleRegion.sigRegionChanged.connect(self.updateRegion)
     

        self.p2.setMaximumHeight(80)
        self.p2.addItem(self.angleRegion)
        #p2.setXRange(float(scanRange[0]),float(scanRange[1]))
        self.p2.showAxis('left',show=False)

        self.qRegion=pg.LinearRegionItem(orientation='horizontal')
        self.p3 = win.addPlot(row=1,col=0)
        self.imgLeft = pg.ImageItem()
        #self.p3.getViewBox().setAspectLocked(True)
        self.p3.addItem(self.imgLeft)
        self.p3.addItem(self.qRegion)


        self.inPlaneImage=pg.ImageItem()
        self.qRegion.hide()
        self.p4 = win.addPlot(row=3,col=0)
        self.p4.hide()
        self.ctrROI=pg.LineROI([0, 5], [0, 0], width=2, pen=(1,9))

        def ctrROIClick():
            #if self.ImgState==2:
            #    return
            binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
            xSize=int(self.p.param('Experiment', 'X Pixels').value()/binning)
            if self.p.param('Profile tools', 'Enable Box Profile').value():
                self.ctrROI.show()
                self.p4.show()
                if self.ImgState == 0:
                    self.ctrROI.setSize(100)
                    self.ctrROI.setPos([xSize/2,xSize/2])
                else:
                    self.ctrROI.setSize(1)
                    self.ctrROI.setPos([0,0])
                                      

            else:
                self.ctrROI.hide()
                self.p4.hide()
            return
        self.p4.setMaximumHeight(300)
        self.p4.setLogMode(x=None,y=True)
        self.CTRCurve=self.p4.plot()
        self.p.param('Profile tools', 'Enable Box Profile').sigValueChanged.connect(ctrROIClick)
        self.ctrROI.hide()
        self.ctrROI.sigRegionChanged.connect(self.ROICalc)
        self.p3.addItem(self.ctrROI)          


        self.hist = pg.HistogramLUTItem()        
        self.hist.setImageItem(self.imgLeft)        
        self.hist.setMaximumWidth(200)
        self.hist.sigLevelChangeFinished.connect(self.updateCrange)
        win.addItem(self.hist,row=1,col=2)
        self.imgLeft.hoverEvent = self.imageHoverEvent   
        self.mask_list = []
        self.binBounds = []
        self.data_format_changed()
        self.threadpool = QThreadPool()
        self.busy = False
        #self.setFixedSize(1024,768)
        self.show() 
