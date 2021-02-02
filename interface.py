#!/usr/bin/env python3

# -*- coding: utf-8 -*-
import sys,os
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QCheckBox, QProgressBar
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QLabel, QGridLayout, QWidget,QPushButton
from PyQt5.QtCore import QSize, Qt    
from scipy import interpolate
import pyqtgraph.parametertree.parameterTypes as pTypes
import pyqtgraph as pg
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import matplotlib.pyplot as plt
import numpy as np
import fabio
import pickle
import fileio

from crystal import *
from fileio import *



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
            self.p3.setTitle("pos: (%0.3f, %0.3f)  real pixel: (%d, %d), Binned intensity: %0.3f" % (x, y, x*binning, (ySize-y*binning),self.imgLeft.image[j,i]))
        elif self.ImgState==2:
            #in-plane map
            self.p3.setTitle("h: %0.3f, k: %0.3f, dist. from origin: %0.3f" % (x, y,np.sqrt(x*x+y*y)))  
        elif self.ImgState==1:
            #reciprocal space (out-of-plane)
            self.p3.setTitle("hk: %0.3f, l: %0.3f, dist. from origin: %0.3f" % (x, y,np.sqrt(x*x+y*y)))


    def makeInPlane(self):
        """Create an in-plane map based on the selected region"""
        if self.select_l_slice:
            self.qRegion.show()
            self.select_l_slice = False
            self.statusLabel.setText('First you must choose a slice to compress along L')
        else:
            if self.image_stack.angle_mode == False:
                self.statusLabel.setText('Error: You probably need to read the log file in')   
                return

            self.ctrROI.hide()
            self.p4.hide()
            self.ImgState=2
            self.imgLeft.resetTransform()
            self.angleRegion.setMovable(False)  

            binning=int(self.p.param('Data Processing', 'Binning').value())
            res=int(self.p.param('Data Processing', 'hkl Resolution').value())
    

            start_angle = self.angleRegion.getRegion()[0]
            end_angle = self.angleRegion.getRegion()[1]
            number_of_images =self.image_stack.number_of_images_in_range(start_angle,end_angle)
            step_size = self.image_stack.get_step()

            #we add extra angles either side so to avoid problems with interpolation 
            #if angluar range is small
            angles=np.linspace(start_angle-step_size,end_angle+step_size,number_of_images+2)
            
            #caluclate momentum transfer as this does not depend on the sample rotation
            qMaps=np.array(self.experiment.dector_frame_to_lab_frame(self,binning))

            #We only need one row since we assume 90 degrees in plane        
            qCentreRow=qMaps[:,int(self.p.param('Experiment', 'Center Pixel Y').value()/binning)]

            #we now return list of hlk coorinates for each angle, along the one row  
            pixMaps=self.experiment.get_coords(angles,qCentreRow,self)

            #We average along the 'L' direction, so we just get one row        
            #qSliceTmp=np.mean(self.ImageData[fromImage:toImage,:,int(self.qRegion.getRegion()[0]):int(self.qRegion.getRegion()[1])],axis=2)
            from_image = int(np.floor(self.image_stack.angle2image(self.angleRegion.getRegion()[0])))
            to_image = int(np.ceil(self.image_stack.angle2image(self.angleRegion.getRegion()[1])))
            print(from_image,to_image)
            qSliceTmp=np.mean(self.image_stack.image_data[from_image:to_image,:,int(self.qRegion.getRegion()[0]):int(self.qRegion.getRegion()[1])],axis=2)

            #add zeros to the ends so we don't try and interperlate past the angles
            emptyr = np.zeros_like(qCentreRow[0])        
            qSlice = np.vstack((emptyr, qSliceTmp,emptyr))      
                
            #create coordinate grids, and interpolate data
            self.grid_h, self.grid_k=np.mgrid[np.min(pixMaps[0]):np.max(pixMaps[0]):res*1j,np.min(pixMaps[1]):np.max(pixMaps[1]):res*1j]
            grid_i = interpolate.griddata((pixMaps[0].ravel(),pixMaps[1].ravel()),np.ravel(qSlice),(self.grid_h, self.grid_k),method='linear')

            self.p3.getAxis('bottom').setGrid(255)
            self.p3.getAxis('left').setGrid(255)

            self.p3.getAxis('bottom').setZValue(1)
            self.p3.getAxis('bottom').setLabel("H (RLU)")

            
            self.p3.getAxis('left').setZValue(1)
            self.p3.getAxis('left').setLabel("K (RLU)")

            
            #grid_i[grid_i<=0]=np.nan

            hScale=(np.max(pixMaps[0])-np.min(pixMaps[0]))/res
            kScale=(np.max(pixMaps[1])-np.min(pixMaps[1]))/res
            self.p3.getViewBox().setAspectLocked(True)
            self.showData=grid_i
            grid_i[grid_i<=0]=np.nan
            self.lShift=np.min(pixMaps[2])
            self.hist.setImageItem(self.imgLeft)
            self.imgLeft.setImage(grid_i)
            self.imgLeft.scale(hScale,kScale)
            
            self.imgLeft.setPos(np.min(self.grid_h),np.min(self.grid_k))
            self.qRegion.hide()
            self.statusLabel.setText('Vewing in recirpocal lattice units, angluar integration is fixed, press pixel view to go back')   
            self.imgLeft.setLevels([0.99*np.mean(self.showData), 1.01*np.mean(self.showData)]) 
            self.update()
            

    def makehkl(self):
        """Convert pixel coordinates to reciprical coordinates"""
        import time
        self.statusLabel.setText('Calculating reciprocal lattice coordinates (Please wait)..')
        self.update()
        self.ImgState=1
        self.imgLeft.resetTransform()
        self.qRegion.hide()
        self.angleRegion.setMovable(False) 
        self.showData=None

        self.p3.getAxis('bottom').setGrid(255)
        self.p3.getAxis('left').setGrid(255)     
        self.p3.getAxis('bottom').setZValue(1)
        self.p3.getAxis('bottom').setLabel("sqrt(H^2 + K^2) (RLU)")    
        self.p3.getAxis('left').setZValue(1)
        self.p3.getAxis('left').setLabel("L (RLU)")     
                      
        binning=int(self.p.param('Data Processing', 'Binning').value())
        numba = False
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()  
        res=int(self.p.param('Data Processing', 'hkl Resolution').value())      

        #calculate the transformation[0]==h, [1]==k, [2]==l,[3]==sqrt(h^2 +k^2),[4]= intensity correciton

        self.pixMaps=self.experiment.dector_frame_to_hkl_frame(self,binning,second_det,0)
        pixMaps = self.pixMaps #just don't want to type self all the time       
        self.showData = self.image_stack.get_image(self.angleRegion.getRegion()[0], self.angleRegion.getRegion()[1])
   
        #create grid space for interpolation 
        self.grid_hk, self.grid_l=np.mgrid[np.min(pixMaps[3]):np.max(pixMaps[3]):res*1j,np.min(pixMaps[2]):np.max(pixMaps[2]):res*1j]

  
        grid_i = interpolate.griddata((pixMaps[3],pixMaps[2]),np.ravel(np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='nearest')

        """#do 2d interpolation, and apply intensity corrections if needed
        if self.p.param('Data Processing', 'Apply intensity corrections').value():
            #we unravel each pixmap to get lists of coordinates (hk, l, correction)
            grid_i = interpolate.griddata((pixMaps[3].ravel(),pixMaps[2].ravel()),np.ravel(pixMaps[4]*np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='linear')
        else:
            grid_i = interpolate.griddata((pixMaps[3].ravel(),pixMaps[2].ravel()),np.ravel(np.rot90(self.showData,k=1)),(self.grid_hk, self.grid_l),method='linear')"""
            

        grid_i[grid_i<=0]=np.nan
        #replace showData with the 
        self.showData=grid_i
        self.gridl = self.grid_l

        #set the image and scale the axis
        self.imgLeft.setImage(grid_i,autoRange=True, autoLevels=True)
        self.hist.setImageItem(self.imgLeft)
        self.imgLeft.scale((np.max(pixMaps[3])-np.min(pixMaps[3]))/res,(np.max(pixMaps[2])-np.min(pixMaps[2]))/res)
        self.imgLeft.setPos(np.min(pixMaps[3]),np.min(pixMaps[2]))        
        self.statusLabel.setText('Vewing in recirpocal lattice units, angluar integration is fixed, press Pixel Coordinates to go back')        
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
                self.imgLeft.setImage(self.image_stack.get_image(from_region,to_region),autoLevels=False)
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
                self.imgLeft.setImage(self.image_stack.get_image(from_region,to_region),autoLevels=False)

            self.p3.getAxis('bottom').setScale(None)
            self.p3.getAxis('left').setScale(None)
            self.p3.getAxis('bottom').setGrid(0)
            self.p3.getAxis('left').setGrid(0)
            self.ctrROICalc()
            if self.ImgState==0:           
                self.p3.getAxis('bottom').setLabel("PIXELS")       
                self.p3.getAxis('left').setLabel("PIXELS")     

    def ctrROICalc(self):
        res=int(self.p.param('Data Processing', 'hkl Resolution').value())
        
        if self.p.param('Profile tools', 'Log profile').value():
            self.p4.setLogMode(x=None,y=True)
        else:
             self.p4.setLogMode(x=None,y=False)

        if not self.p.param('Profile tools', 'Enable Box Profile').value():
            return
        #pixel view
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
        #hkl view
        elif self.ImgState==2:
            if self.p.param('Profile tools', 'Axis of interest').value() == 1:
                xRoiData=self.ctrROI.getArrayRegion(self.grid_h, self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 2:
                xRoiData=self.ctrROI.getArrayRegion(self.grid_k, self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 3:
                xRoiData=self.ctrROI.getArrayRegion(np.sqrt(self.grid_h**2+self.grid_k**2), self.imgLeft)
            elif self.p.param('Profile tools', 'Axis of interest').value() == 4:                
                xRoiData=self.ctrROI.getArrayRegion(self.grid_h, self.imgLeft)

            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft)  
            xdata=np.mean(xRoiData,axis=1)
            self.roixdata = xdata
            self.roiydata = np.mean(ROIData,axis=1)
            self.CTRCurve.setData(x=xdata,y=self.roiydata)#-np.min(ROIData[0],axis=1))
        #inplane view
        else:      
            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft)     
            self.roiydata = np.mean(ROIData,axis=1)                
            self.roixdata = range(0,len(np.mean(ROIData,axis=1)))  
            self.CTRCurve.setData(self.roiydata)

    def roipop(self):       
        """This is a quite testing function to create a matplotlib windows, however
        it needs integrating into pyqt as plt.show() currently generates warining due
        to Qt conflicts"""

        #Are we in hkl mode?   
        if self.ImgState == 1:  
            res=int(self.p.param('Data Processing', 'hkl Resolution').value())
            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft,returnMappedCoords=False)

            #h axis
            if self.p.param('Profile tools', 'Axis of interest').value() == 1:
                #we make a mesh grid using the min/max of H coorindate array (pixMaps[0]) and the min/max of L coordinate array pixMaps[2]
                grid_h, tmp=np.mgrid[np.min(self.pixMaps[0]):np.max(self.pixMaps[0]):res*1j,np.min(self.pixMaps[2]):np.max(self.pixMaps[2]):res*1j]
                #we then use the region of interest to take h values for each position in the roi
                xRoiData=self.ctrROI.getArrayRegion(grid_h, self.imgLeft)
            #k
            elif self.p.param('Profile tools', 'Axis of interest').value() == 2:
                grid_k, tmp=np.mgrid[np.min(self.pixMaps[1]):np.max(self.pixMaps[1]):res*1j,np.min(self.pixMaps[2]):np.max(self.pixMaps[2]):res*1j]
                xRoiData=self.ctrROI.getArrayRegion(grid_k, self.imgLeft)
            #hk
            elif self.p.param('Profile tools', 'Axis of interest').value() == 3:
                xRoiData=self.ctrROI.getArrayRegion(self.grid_hk, self.imgLeft)

            #L axis doesn't make too much sense
            elif self.p.param('Profile tools', 'Axis of interest').value() == 4:
                print("L is already the vertical axis please select a diffrent axis")
                return

            #we already have a mesh grid of l values, so we exctract l values using that
            yRoiData=self.ctrROI.getArrayRegion(self.gridl, self.imgLeft)           

            #Here we plot
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cmin,cmax=self.imgLeft.getLevels()
            ax.pcolormesh(xRoiData,yRoiData,ROIData,vmin=cmin, vmax=cmax,cmap='viridis',shading='auto')
            ax.grid(True,color='w',linestyle='-', linewidth=0.04)          
            plt.show()  #This is bad practice!! will cuase error:
            #QCoreApplication::exec: The event loop is already running
       
        #if we are view in pixel coordinates we don't need to do anything too fancy
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cmin,cmax=self.imgLeft.getLevels()
            ROIData=self.ctrROI.getArrayRegion(self.showData, self.imgLeft,vmin=cmin, vmax=cmax,cmap='viridis',shading='auto')
            ax.imshow(ROIData)
            ax.grid(True,color='w',linestyle='-', linewidth=0.04)
            ax.set_aspect('equal')    
            plt.show()          

    def rawrefresh(self):
        """Go back to pixel coordiante mode"""
        self.ImgState=0
        self.select_l_slice = True
        self.angleRegion.setMovable(True)  
        self.hist.setImageItem(self.imgLeft)        
        self.imgLeft.show()        
        self.imgLeft.resetTransform()
        self.updateRegion()
        self.ctrROICalc

        self.p3.getAxis('bottom').setScale(None)
        self.p3.getAxis('left').setScale(None)
        self.p3.getAxis('bottom').setGrid(0)
        self.p3.getAxis('left').setGrid(0)
        self.p3.getAxis('bottom').setLabel("PIXELS")       
        self.p3.getAxis('left').setLabel("PIXELS")        
        self.statusLabel.setText('Viewing in pixel coordinates')
        self.p3.getViewBox().setAspectLocked(False)

    def data_format_changed(self):
        """We update the interface if a beamline preset is chosen"""
        #ID31 Beamline
        if self.experiment.param('Beamline Preset').value() == 4:
            self.p.param('File', 'Select Dark Image').hide()
            self.p.param('File', 'Select Background Dark Image').hide()
            self.p.param('File', 'Select Data Images').setName("Select HDF5 file")
            self.p.param('File', 'Select Log File').hide()
            self.p.param('Data Processing', 'Use 2nd detector').setValue(False)
            self.p.param('Data Processing', 'Use 2nd detector').hide()           
        #P07 Beamline
        if self.experiment.param('Beamline Preset').value() == 3:
            self.p.param('File', 'Select Log File').hide()
            self.p.param('Experiment', 'Y Pixels').setValue(2048)
            self.p.param('Experiment', 'X Pixels').setValue(2048)
            self.p.param('Experiment', 'Pixel Size').setValue(200e-6)
            self.p.param('Experiment', 'Energy').setValue(73700)
            self.p.param('Experiment', 'Sample-Detector Dist.').setValue(1.6)
            #There was no 2nd detector and dark images subtracted normaly 
            self.p.param('Data Processing', 'Use 2nd detector').hide()
            self.p.param('File', 'Select Dark Image').hide()
            self.p.param('File', 'Select Background Dark Image').hide()
    
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


    def saveprofile(self):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","csv (*.csv);;All Files (*)", options=options)      
        data = np.asarray([self.roixdata,self.roiydata])
        np.savetxt(fileName,np.transpose(data),fmt='%10.5f', delimiter=',',newline='\n')       

    def saverocks(self):
        options = QFileDialog.Options()
        #options |= QFileDialog.DontUseNativeDialog
        folderName = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        #get the image indices to use
        if self.LogFile is None:
            fromImage=int(np.round(self.angleRegion.getRegion()[0]))
            toImage=int(np.round(self.angleRegion.getRegion()[1]))
        else:
            angle2image=interpolate.interp1d((self.scanRange[0],self.scanRange[1]),(0,len(self.Images)))
            fromImage=int(np.round(angle2image(self.angleRegion.getRegion()[0])))
            toImage=int(np.round(angle2image(self.angleRegion.getRegion()[1])))

        #images in our angular range
        tmp = self.ImageData[fromImage:toImage+1,:,:]
        profiles = []
        angles = []


        for i,image in enumerate(tmp):
            ROIData=self.ctrROI.getArrayRegion(image, self.imgLeft) 
            profiles.append(np.sum(ROIData,axis=1))
            if self.LogFile is None:
                angles.append(fromImage+i)
            else:
                angles.append(self.allangles[fromImage+i+1])

        #next we loop through each line in
        profiles = np.asarray(profiles)
        plt.plot(angles,profiles[:,0])
        plt.plot(angles,profiles[:,-1])
        plt.show()
        #for i in range(len(profiles))
        #data = np.asarray([self.roixdata,self.roiydata])
        print(angles)     
        #np.savetxt(fileName,np.transpose(data),fmt='%10.5f', delimiter=',',newline='\n')    
              

    def extractrois(self):
        options = QFileDialog.Options()
        roifiles, _ = QFileDialog.getOpenFileNames(self,"Select ROIs to use", "","ROI Files (*.roi);;All Files (*)", options=options)
        bkgroifiles, _ = QFileDialog.getOpenFileNames(self,"Select Background ROIs to use", "","ROI Files (*.roi);;All Files (*)", options=options)      
        l_values =  np.array([])
        i_values = np.array([])
        bkg_values =np.array([])
        if roifiles:
            self.progressBar.setMaximum(len(roifiles))
            for i in range(len(roifiles)):
                self.progressBar.setValue(i)
                self.update()
                print(i+1,"/",str(),len(roifiles))
                self.rawrefresh()
                #load the roi
                state = pickle.load(open(roifiles[i], "rb" ))
                self.ctrROI.setState(state[0])
                self.angleRegion.setRegion(state[1])
                #calculate hkl and the roi again
                self.makehkl()
                self.ctrROICalc()
                self.update()
                #extract the profile                
                i_values = np.append(i_values,self.roiydata)
                l_values = np.append(l_values,self.roixdata)              
                #now for the background
                self.rawrefresh()
                state = pickle.load(open(bkgroifiles[i], "rb" ))
                self.ctrROI.setState(state[0])
                self.angleRegion.setRegion(state[1])
                self.makehkl()
                self.ctrROICalc()
                self.update()                
                #keep the background subtracted values
                bkg_values = np.append(bkg_values,self.roiydata)  
        print(l_values)     

        y = i_values - bkg_values
        x = l_values
        fig = plt.figure(figsize=(5,5),dpi=300)        
        plt.plot(x,y)
        plt.show()
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","csv (*.csv);;All Files (*)", options=options)  
        data = np.vstack((x,y))
        np.savetxt(fileName,np.transpose(data),fmt='%10.5f', delimiter=',',newline='\n')
         
    def saveroi(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","roi (*.roi);;All Files (*)", options=options)     
        if fileName:
            state = [self.ctrROI.saveState(),self.angleRegion.getRegion()]
            pickle.dump(state, open( fileName, "wb" ))   

    def loadroi(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self,"Select a ROI file too use", "","ROI File (*.roi);;All Files (*)", options=options)
        if fileName:
            state = pickle.load(open(fileName, "rb" ))
            self.ctrROI.setState(state[0])
            self.angleRegion.setRegion(state[1])
    
    def roictr(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileNames(self,"Select ROI files too use", "","ROI File (*.roi);;All Files (*)", options=options)
        if fileName:
            state = pickle.load(open(fileName, "rb" ))
            self.ctrROI.setState(state[0])
            self.angleRegion.setRegion(state[1])

    def load(self):
        if self.image_stack.images_read:
            self.image_stack.load_aux_files()
            self.image_stack.load_images()  
            
            if self.image_stack.angle_mode:
                self.angleRegion.setRegion((self.image_stack.start_angle,self.image_stack.end_angle))
                self.angleRegion.setBounds((self.image_stack.start_angle,self.image_stack.end_angle))
                self.showData = self.image_stack.get_image(self.image_stack.start_angle,self.image_stack.end_angle)
                self.imgLeft.setImage(self.showData)
                self.imgLeft.show()                    
                self.hist.setImageItem(self.imgLeft)     
                self.p2.setLabel('bottom',text='Angle (degrees)')       
            else:
                self.angleRegion.setRegion((1,self.image_stack.number_of_images))
                self.angleRegion.setBounds((1,self.image_stack.number_of_images))
                self.showData = self.image_stack.get_image(0,self.image_stack.number_of_images)
                self.imgLeft.setImage(self.showData)
                self.imgLeft.show()                    
                self.hist.setImageItem(self.imgLeft)       
            
        else:  
            self.statusLabel.setText('Background that will be subtracted. Select images to load')
            self.image_stack.load_aux_files()
            self.imgLeft.setImage(self.image_stack.subtract_image)
            self.imgLeft.show()    
            self.hist.setImageItem(self.imgLeft)   
            self.ImgState = 5

    def save(self):
        state = self.p.saveState()
        pickle.dump(state, open( 'params.bin', "wb" ))      
        pickle.dump(self.image_stack.save_state(),open('imagestack.bin', "wb" ) ) 
        print("Parameters and file selection saved\n")

    def restore(self):
        self.p.restoreState(pickle.load(open( 'params.bin', "rb" )),removeChildren=False)
        print("Saved parameters restored\n")

    def restoreImageStack(self):
        print("Restoring file selection:")
        self.image_stack.restore_state(pickle.load(open( 'imagestack.bin', "rb" )))   
        self.image_stack.load_aux_files()           
        print(self.image_stack.flags)
        print(self.image_stack.number_of_images, "images selected")
        print("Complete (Press Load)\n")
        
        
    def initUI(self):
        self.setWindowTitle(self.title)
        gridLayout = QGridLayout(self)     
        self.setLayout(gridLayout)
        
        params = [
            self.experiment,
            self.crystal,
            {'name': 'Data Processing', 'type': 'group', 'children': [
                {'name': 'Binning', 'type': 'int', 'value': 4, 'step': 1},
                {'name': 'hkl Resolution', 'type': 'int', 'value': 1000, 'step': 1},
                {'name': 'White background', 'type': 'bool', 'value': False},
                {'name': 'Mean Images instead of Max', 'type': 'bool',},
                {'name': 'Acceleration', 'type': 'list', 'values': {"none": 0, "numba (cpu)": 1, "cuda (gpu)": 2},'value': 0},
                 {'name': 'Use 2nd detector', 'type': 'bool', 'value': False},                
                {'name': 'Apply intensity corrections', 'type': 'bool', 'value': False},
                {'name': 'Correct for refraction', 'type': 'bool', 'value': False},
                {'name': 'Intensity offset', 'type': 'float', 'value': 0, 'step':10},
                {'name': 'Multiply intensity by', 'type': 'float', 'value': 1, 'step':0.1},
                {'name': 'Image Rotation', 'type': 'list', 'values': {"0 degrees": 0, "90 degrees": 1, "180 degress": 2,"270 degress": 3},'value': 0},
                {'name': 'Image Flip U/D', 'type': 'list', 'values': {"True": True,"False": False},'value': False},
                {'name': 'Image Flip L/R', 'type': 'list', 'values': {"True": True,"False": False},'value': False},
            ]},
            {'name': 'File', 'type': 'group', 'children': [  
                {'name': 'Load', 'type': 'action'},
                {'name': 'Restore File Selection', 'type': 'action'},
                {'name': 'Select Dark Image', 'type': 'action'},
                {'name': 'Select Background Image', 'type': 'action'},
                {'name': 'Select Background Dark Image', 'type': 'action'},
                {'name': 'Select Data Images', 'type': 'action'},
                {'name': 'Select Log File', 'type': 'action'},
                {'name': 'Save Parameters and Files', 'type': 'action'},
                {'name': 'Load Parameters', 'type': 'action'}, 
            ]},

            {'name': 'View Mode', 'type': 'group', 'children': [             
                {'name': 'Pixel Coordinates', 'type': 'action'},
                {'name': 'Out-of-plane Map', 'type': 'action'},                
                {'name': 'In-plane Map', 'type': 'action'},

            ]},
            {'name': 'Profile tools', 'type': 'group', 'children': [
                {'name': 'Enable Box Profile', 'type': 'bool',},
                {'name': 'Log profile', 'type': 'bool', 'value': True},
                {'name': 'ROI Popout', 'type': 'action'},
                #{'name': 'Export ROI as TIF', 'type': 'action'},
                {'name': 'Extract profile', 'type': 'action'},
                {'name': 'Save Rocking Scans', 'type': 'action'},
                {'name': 'Extract CTR ROIs', 'type': 'action'},
                #{'name': 'Batch Extract', 'type': 'action'},                
                {'name': 'Axis of interest', 'type': 'list',  'values': {"H": 1, "K": 2, "HK": 3,"L": 4}, 'value': 4} ,
                {'name': 'Save ROI', 'type': 'action'},
                {'name': 'Load ROI', 'type': 'action'},
            ]}
        ]


        #p.param('Save/Restore functionality', 'Save State').sigActivated.connect(save)
        def toggleBg():
            if self.p.param('Data Processing', 'White background').value():
                win.setBackground('w')
            else:
                win.setBackground('k')       

        self.p = Parameter.create(name='params', type='group', children=params)

        self.image_stack=ImageStack(self)

        self.p.param('File', 'Select Dark Image').sigActivated.connect(self.image_stack.select_dark)
        self.p.param('File', 'Select Background Image').sigActivated.connect(self.image_stack.select_dark)
        self.p.param('File', 'Select Background Dark Image').sigActivated.connect(self.image_stack.select_background_dark)
        self.p.param('File', 'Select Data Images').sigActivated.connect(self.image_stack.select_images)
        self.p.param('File', 'Select Log File').sigActivated.connect(self.image_stack.read_log_file)
        self.p.param('File', 'Load').sigActivated.connect(self.load)
        self.p.param('File', 'Save Parameters and Files').sigActivated.connect(self.save)
        self.p.param('File', 'Load Parameters').sigActivated.connect(self.restore)
        self.p.param('File', 'Restore File Selection').sigActivated.connect(self.restoreImageStack)

        self.p.param('View Mode', 'Pixel Coordinates').sigActivated.connect(self.rawrefresh)
        self.p.param('View Mode', 'Out-of-plane Map').sigActivated.connect(self.makehkl)        
        self.p.param('View Mode', 'In-plane Map').sigActivated.connect(self.makeInPlane)


        self.p.param('Data Processing', 'White background').sigValueChanged.connect(toggleBg)
        self.p.param('Data Processing', 'Mean Images instead of Max').sigValueChanged.connect(self.updateRegion)   
        
        self.p.param('Profile tools', 'ROI Popout').sigActivated.connect(self.roipop)
        self.p.param('Profile tools', 'Extract profile').sigActivated.connect(self.saveprofile)
        self.p.param('Profile tools', 'Save Rocking Scans').sigActivated.connect(self.saverocks)
        self.p.param('Profile tools', 'Log profile').sigValueChanged.connect(self.ctrROICalc)  
        self.p.param('Profile tools', 'Axis of interest').sigValueChanged.connect(self.ctrROICalc)  
        self.p.param('Profile tools', 'Extract CTR ROIs').sigActivated.connect(self.extractrois)
        self.p.param('Profile tools', 'Save ROI').sigActivated.connect(self.saveroi)
        self.p.param('Profile tools', 'Load ROI').sigActivated.connect(self.loadroi) 

        self.experiment.param('Beamline Preset').sigValueChanged.connect(self.data_format_changed)
        self.crystal.param('Preset').sigValueChanged.connect(self.sample_preset)

        t = ParameterTree()
        t.setMaximumWidth(320)
        t.setParameters(self.p, showTop=False)

        gridLayout.addWidget(t,0,2)
        win = pg.GraphicsLayoutWidget()
        gridLayout.addWidget(win,0,1)

        self.progressBar=QProgressBar()
        gridLayout.addWidget(self.progressBar,1,1)

        self.statusLabel=QLabel('Select Images ... (Press Load Saved Files to load the files you had when you last pressed \'Save Parameters\')')
        gridLayout.addWidget(self.statusLabel,1,1)

        win.nextRow()
        self.p2 = win.addPlot(row=3,colspan=3)
        self.angleRegion=pg.LinearRegionItem()
        self.angleRegion.sigRegionChanged.connect(self.updateRegion)
     

        self.p2.setMaximumHeight(80)
        self.p2.addItem(self.angleRegion)
        #p2.setXRange(float(scanRange[0]),float(scanRange[1]))
        self.p2.showAxis('left',show=False)

        self.qRegion=pg.LinearRegionItem(orientation='horizontal')
        self.p3 = win.addPlot(row=0,col=0)
        self.imgLeft = pg.ImageItem()
        #self.p3.getViewBox().setAspectLocked(True)
        self.p3.addItem(self.imgLeft)
        self.p3.addItem(self.qRegion)


        self.inPlaneImage=pg.ImageItem()
        self.qRegion.hide()
        self.p4 = win.addPlot(row=2,col=0)
        self.p4.hide()
        self.ctrROI=pg.LineROI([0, 5], [0, 0], width=2, pen=(1,9))

        def ctrROIClick():
            #if self.ImgState==2:
            #    return
            binning=int(self.p.param('Data Processing', 'Binning').value())
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
        self.ctrROI.sigRegionChanged.connect(self.ctrROICalc)
        self.p3.addItem(self.ctrROI)          


        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.imgLeft)
        self.hist.setMaximumWidth(200)
        win.addItem(self.hist,row=0,col=2)
        self.imgLeft.hoverEvent = self.imageHoverEvent   
        


        win.show()
        
        self.setFixedSize(1920,1080)
        self.show() 
