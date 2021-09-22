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
import math
import cmath


class CrystallographyParameters(pTypes.GroupParameter):

    def calcBInv(self):
        B = np.matrix([[self.param('b₁').value(), self.param('b₂').value()*np.cos(np.deg2rad(self.param('β₃').value())), self.param('b₃').value()*np.cos(np.deg2rad(self.param('β₂').value()))],
        [0, self.param('b₂').value()*np.sin(np.deg2rad(self.param('β₃').value())), -1*self.param('b₃').value()*np.sin(np.deg2rad(self.param('β₂').value()))*np.cos(np.deg2rad(self.param('α₁').value()))],
        [0, 0, 2*np.pi/self.param('a₃').value()]])

        #we don't really know this so it's just the identity
        U = np.matrix( [[1, 0,0],
                        [0, 1, 0],
                        [0, 0, 1]])
            
        UB = U*B
        Binv = np.linalg.inv(UB)
        return Binv

    def __init__(self, **opts):
        opts['type'] = 'bool'
        opts['value'] = True
        pTypes.GroupParameter.__init__(self, **opts)  
        self.addChild({'name': 'Preset', 'type': 'list', 'values': {"Manual": 0, "Au (111) Surf.": 1, "Au (001) Surf.": 2, "TiO2": 3},'value': 0})
        self.addChild({'name': 'a₁', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.001})
        self.addChild({'name': 'a₂', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.001})
        self.addChild({'name': 'a₃', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.001})

        self.addChild({'name': 'α₁', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'α₂', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'α₃', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})

        self.addChild({'name': 'b₁', 'type': 'float', 'value': 2.52,'suffix': 'Å⁻¹', 'step': 0.001})
        self.addChild({'name': 'b₂', 'type': 'float', 'value': 2.53,'suffix': 'Å⁻¹', 'step': 0.001})
        self.addChild({'name': 'b₃', 'type': 'float', 'value': 0.89,'suffix': 'Å⁻¹', 'step': 0.001})

        self.addChild({'name': 'β₁', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'β₂', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'β₃', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'Reciprical Units', 'type': 'list', 'values': {"Q": 0, "HKL": 1}, 'value':1})

        self.a1=self.param('a₁')
        self.a2=self.param('a₂')
        self.a3=self.param('a₃')

        self.A1=self.param('α₁')
        self.A2=self.param('α₂')
        self.A3=self.param('α₃')

        self.b1=self.param('b₁')
        self.b2=self.param('b₂')
        self.b3=self.param('b₃')

        self.B1=self.param('β₁')
        self.B2=self.param('β₂')
        self.B3=self.param('β₃')

        self.a1.sigValueChanged.connect(self.getReciprocal)
        self.a2.sigValueChanged.connect(self.getReciprocal)
        self.a3.sigValueChanged.connect(self.getReciprocal)

        self.A1.sigValueChanged.connect(self.getReciprocal)
        self.A2.sigValueChanged.connect(self.getReciprocal)
        self.A3.sigValueChanged.connect(self.getReciprocal)

        self.b1.sigValueChanged.connect(self.getReal)
        self.b2.sigValueChanged.connect(self.getReal)
        self.b3.sigValueChanged.connect(self.getReal)

        self.B1.sigValueChanged.connect(self.getReal)
        self.B2.sigValueChanged.connect(self.getReal)
        self.B3.sigValueChanged.connect(self.getReal)

        self.getReciprocal()

    def getReciprocal(self,event=None):
        b1,b2,b3,B1,B2,B3=self.invertLattice(self.a1.value(),self.a2.value(),self.a3.value(),np.deg2rad(self.A1.value()),\
            np.deg2rad(self.A2.value()),np.deg2rad(self.A3.value()))

        self.b1.setValue(b1*2*np.pi,blockSignal=self.getReal)
        self.b2.setValue(b2*2*np.pi,blockSignal=self.getReal)
        self.b3.setValue(b3*2*np.pi,blockSignal=self.getReal)
        self.B1.setValue(np.rad2deg(B1),blockSignal=self.getReal)
        self.B2.setValue(np.rad2deg(B2),blockSignal=self.getReal)
        self.B3.setValue(np.rad2deg(B3),blockSignal=self.getReal)

    def getReal(self):
        a1,a2,a3,A1,A2,A3=self.invertLattice(self.b1.value(),self.b2.value(),self.b3.value(),np.deg2rad(self.B1.value()),\
            np.deg2rad(self.B2.value()),np.deg2rad(self.B3.value()))
        self.a1.setValue(a1/2/np,blockSignal=self.getReciprocal)
        self.a2.setValue(a2/2/np,blockSignal=self.getReciprocal)
        self.a3.setValue(a3/2/np,blockSignal=self.getReciprocal)
        self.A1.setValue(np.rad2deg(A1),blockSignal=self.getReciprocal)
        self.A2.setValue(np.rad2deg(A2),blockSignal=self.getReciprocal)
        self.A3.setValue(np.rad2deg(A3),blockSignal=self.getReciprocal)
    
    def invertLattice(self,a,b,c,A,B,C):
        V=a*b*c*np.sqrt(1-np.cos(A)**2-np.cos(B)**2-np.cos(C)**2\
            +2*np.cos(A)*np.cos(B)*np.cos(C))

        b1=b*c*np.sin(A)/V
        b2=a*c*np.sin(B)/V
        b3=a*b*np.sin(C)/V

        B1=np.arccos(np.cos(B)*np.cos(C)-np.cos(A)/(np.sin(B)*np.sin(C)))
        B2=np.arccos(np.cos(A)*np.cos(C)-np.cos(B)/(np.sin(A)*np.sin(C)))
        B3=np.arccos(np.cos(A)*np.cos(B)-np.cos(C)/(np.sin(A)*np.sin(B)))

        return(b1,b2,b3,B1,B2,B3)

        self.addChild({'name': 'Real Centre Angle', 'type': 'float', 'value': 0,'suffix': '°','suffixGap': ''})

class ExperimentParameter(pTypes.GroupParameter):
    
    def dector_frame_to_hkl_frame(self,window,binning,twodetectors,midangle):
        """Function to generate coordinate grid for detector image"""
        angi=np.deg2rad(self.param('Angle of Incidence').value())
        omega_rad=np.deg2rad(self.param('Angle offset').value())
        rotDirection=(self.param('Axis directon').value())
        pixel_count_x=int(self.param('X Pixels').value()/binning)
        pixel_count_y=int(self.param('Y Pixels').value()/binning)
        acceleration=window.p.param('Data Processing', 'Acceleration').value() 
        correct_refraction=window.p.param('Data Processing', 'Correct for Refraction').value() 
        correct=window.p.param('Data Processing', 'Apply Intensity Corrections').value()
        angi_c = np.deg2rad(self.param('Critical Angle').value()*binning)
        if twodetectors:
            pixel_count_x=int((2*self.param('X Pixels').value()+self.param('Detector Gap (pixels)').value())/binning)
        pixel_x=self.param('Pixel Size').value()*binning
        pixel_y=self.param('Pixel Size').value()*binning
        
        SDD=self.param('Sample-Detector Dist.').value()
        invSDD=1/self.param('Sample-Detector Dist.').value()   
        q0=np.array([self.param('Center Pixel X').value()/binning,self.param('Center Pixel Y').value()/binning]).copy()
        k0=(np.pi*2)/self.param('Wavelength').value()
        p_h=self.param('Horizontal Polarisation').value()
        Binv=window.crystal.calcBInv()
        beta3 = window.crystal.param('β₃').value()
        b1=window.crystal.param('b₁').value()
        b2=window.crystal.param('b₂').value()
        recp_units=window.crystal.param('Reciprical Units').value()

        if  acceleration==2:
            if "cuda" not in sys.modules:
                from numba import cuda
            @cuda.jit(fastmath=True)
            def detector_to_hkl_kernel(h_glob,k_glob,l_glob,hk_glob,c_glob,midang):
                #get the current thread position
                row,col = cuda.grid(2)

                if row < hk_glob.shape[0] and col < hk_glob.shape[1]:
                    delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                    delR = math.sqrt(delta_x**2 + delta_z**2)            
                    dist = math.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                    #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                    #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007   

                    gam_pix  = math.atan(delta_x*invSDD)
                    del_pix =  math.atan(delta_z/math.sqrt(delta_x**2 + SDD**2))-angi*math.cos(gam_pix)                                    
                                     
                    qx = k0*(math.cos(del_pix)*math.cos(gam_pix)-math.cos(angi))                    
                    qy = k0*(math.cos(del_pix)*math.sin(gam_pix)) 
                    qz = k0*(math.sin(del_pix)+math.sin(angi))

                    #refraction correction
                    # Busch et al., J. Appl. Cryst. 39, 433-442 (2006)
                    if correct_refraction:
                        kzi = k0*math.sqrt(abs(math.sin(angi)**2 - math.sin(angi_c)**2))
                        kzf = k0*math.sqrt(abs(math.sin(del_pix)**2 - math.sin(angi_c)**2))
                        qz = kzf - kzi       

                    if correct:
                        Phor = 1- (math.cos(del_pix)*math.sin(gam_pix))**2
                        Pver = 1 - math.sin(del_pix)**2
                        pol = p_h*Phor + (1-p_h)*Pver
                        #L= (math.sin(del_pix)*math.cos(angi))*math.copysign(1,del_pix)    
                        L=math.sin(gam_pix)*math.cos(del_pix)*math.cos(angi)*math.copysign(1,gam_pix)                            
                        Cd = dist**2/SDD**2    
                        Ci = 1/math.cos(math.atan(delR/SDD))                            
                        Crod = (math.sin(angi-gam_pix)**2 * math.cos(del_pix) + math.cos(angi)*math.cos(angi-gam_pix))/math.cos(angi)         
                        c_glob[row,col] = L*Crod*Cd*pol*Ci
                    else:
                        c_glob[row,col] = 1                 
                
                    """apply sample rotations to frame of reference
                    rotations are phi and omega_h which is angle of incidence
                    this is eqn 44. in ref. 1"""
                    so = math.sin(rotDirection*(omega_rad+midang))
                    co = math.cos(rotDirection*(omega_rad+midang))
                    ci = 1 #angle of incidence set ealier q calc
                    si = 0 #angle of incidence set ealier q calc

                    hphi_1 = so*qy+co*qx
                    hphi_2 = co*qy-so*qx
                    hphi_3 = ci*qz-si*qy

                    if recp_units == 1:
                        h=  Binv[0][0]*hphi_1+Binv[0][1]*hphi_2+Binv[0][2]*hphi_3
                        k = Binv[1][0]*hphi_1+Binv[1][1]*hphi_2+Binv[1][2]*hphi_3
                        l = Binv[2][0]*hphi_1+Binv[2][1]*hphi_2+Binv[2][2]*hphi_3
                    else:
                        h = hphi_1
                        k = hphi_2
                        l = hphi_3                           

                    h_glob[row,col] = h
                    k_glob[row,col] = k
                    l_glob[row,col] = l

                    if recp_units == 1 and b1 != b2:
                        hk_glob[row,col] = math.copysign(1,gam_pix)*math.sqrt((hphi_1/b1)**2 + (hphi_2/b1)**2)
                    elif recp_units == 1 and abs(90 - beta3) > 0.1:
                        hk_glob[row,col] = math.copysign(1,gam_pix)*math.sqrt((hphi_1/b1)**2 + (hphi_2/b2)**2)
                    else:
                        hk_glob[row,col] = math.copysign(1,gam_pix)*math.sqrt((h)**2 + (k)**2)                      
                    
                
            h_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            k_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            l_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x))) 
            hk_global_mem = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            c_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x))) #array of pixel intensity corrections           

            # Configure the blocks
            threadsperblock = (32, 8)
            blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)
     
            detector_to_hkl_kernel[blockspergrid, threadsperblock](h_global_mem,k_global_mem,l_global_mem,hk_global_mem,c_global_mem,midangle)        
            return [h_global_mem.copy_to_host(),k_global_mem.copy_to_host(),l_global_mem.copy_to_host(),hk_global_mem.copy_to_host(),c_global_mem.copy_to_host()]  

        #this is the none gpu based function 
        def _dector_frame_to_hkl_frame():
                h_img = np.zeros((pixel_count_y,pixel_count_x)) 
                k_img = np.zeros((pixel_count_y,pixel_count_x))
                l_img = np.zeros((pixel_count_y,pixel_count_x)) 
                hk_img = np.zeros((pixel_count_y,pixel_count_x))
                c_img = np.zeros((pixel_count_y,pixel_count_x)) #array of pixel intensity corrections    
             
                """Function calculate the pixel positon in lab frame"""
                for row in range(pixel_count_y):    
                    for col in range(pixel_count_x):
                        delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                        #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007        
                        gam_pix  = np.arctan(delta_x/ SDD)
                        del_pix = np.arctan(delta_z/np.sqrt(delta_x**2 + SDD**2))-angi*np.cos(gam_pix)   

                        qx = k0*(np.cos(del_pix)*np.cos(gam_pix)-np.cos(angi))
                        qy = k0*(np.cos(del_pix)*np.sin(gam_pix)) 
                        qz = k0*(np.sin(del_pix)+np.sin(angi))

                        if correct_refraction:
                            kzi = k0*np.sqrt(abs(np.sin(angi)**2 - np.sin(angi_c)**2))
                            kzf = k0*np.sqrt(abs(np.sin(del_pix)**2 - np.sin(angi_c)**2))
                            qz = kzf - kzi
                        if correct:
                            #tth = np.arccos(np.cos(del_pix)*np.cos(gam_pix))
                        
                            #polarization correction
                            Phor = 1- (np.cos(del_pix)*np.sin(gam_pix))**2
                            Pver = 1 - np.sin(del_pix)**2
                            pol = p_h*Phor + (1-p_h)*Pver   
                            L=math.sin(gam_pix)*math.cos(del_pix)*math.cos(angi)*math.copysign(1,gam_pix)      
                            Cd = dist**2/SDD**2
                            Ci = 1/np.cos(math.atan(delR/SDD))                            
                            Crod = (np.sin(angi-gam_pix)**2 * np.cos(del_pix) + np.cos(angi)*np.cos(angi-gam_pix))/np.cos(angi)              
                    
                            c_img[row][col] = L*Crod*Cd*pol*Ci
                        else:
                            c_img[row][col] = 1

                        so = np.sin(rotDirection*(omega_rad+midangle))
                        co = np.cos(rotDirection*(omega_rad+midangle))
                        ci = 1
                        si = 0  
                        """apply sample rotations to frame of reference
                            rotations are phi and omega_h which is angle of incidence
                            this is eqn 44. in ref. 1"""

                        hphi_1 = so*(ci*qy+si*qz)+co*qx
                        hphi_2 = co*(ci*qy+si*qz)-so*qx
                        hphi_3 = ci*qz-si*qy  

                        if recp_units == 1:
                            h=  Binv[0,0]*hphi_1+Binv[0,1]*hphi_2+Binv[0,2]*hphi_3
                            k = Binv[1,0]*hphi_1+Binv[1,1]*hphi_2+Binv[1,2]*hphi_3
                            l = Binv[2,0]*hphi_1+Binv[2,1]*hphi_2+Binv[2,2]*hphi_3
                        else:
                            h = hphi_1
                            k = hphi_2
                            l = hphi_3    

                        h_img[row,col] = h
                        k_img[row,col] = k
                        l_img[row,col] = l    

                        if recp_units == 1 and b1 != b2:
                            hk_img[row,col] = math.copysign(1,gam_pix)*math.sqrt((hphi_1/b1)**2 + (hphi_2/b1)**2)
                        elif recp_units == 1 and abs(90 - beta3) > 0.1:
                            hk_img[row,col] = math.copysign(1,gam_pix)*math.sqrt((hphi_1/b1)**2 + (hphi_2/b2)**2)
                        else:
                            hk_img[row,col] = math.copysign(1,gam_pix)*math.sqrt((h)**2 + (k)**2)         
                return [h_img,k_img,l_img,hk_img,c_img]

        #use numba on cpu instead
        if acceleration == 1:
            if "jit" not in sys.modules:                
                from numba import jit
            detector_frame_to_hkl=jit(_dector_frame_to_hkl_frame,nopython=True)
        else:
            detector_frame_to_hkl=_dector_frame_to_hkl_frame

        tmp = np.array(detector_frame_to_hkl())
        return (tmp)  

    def projection_binner(self,window,binning,twodetectors,angles,imagedata,qxymax,qzmax):
        angi=np.deg2rad(self.param('Angle of Incidence').value())
        omega_rad=np.deg2rad(self.param('Angle offset').value())
        rotDirection=(self.param('Axis directon').value())
        use_raw_files= window.p.param('Data Processing', 'Bin From Full Images').value()   

        #some parameters will be different with we don't use the binned data in memory
        if use_raw_files:
            pixel_count_x=int(self.param('X Pixels').value())
            pixel_count_y=int(self.param('Y Pixels').value())
            if twodetectors:
                pixel_count_x=int((2*self.param('X Pixels').value()+self.param('Detector Gap (pixels)').value()))
            q0=np.array([self.param('Center Pixel X').value(),self.param('Center Pixel Y').value()]).copy()
            pixel_x=self.param('Pixel Size').value()
            pixel_y=self.param('Pixel Size').value()
            
        else:
            pixel_count_x=int(self.param('X Pixels').value()/binning)
            pixel_count_y=int(self.param('Y Pixels').value()/binning)
            if twodetectors:
                pixel_count_x=int((2*self.param('X Pixels').value()+self.param('Detector Gap (pixels)').value())/binning)
            q0=np.array([self.param('Center Pixel X').value()/binning,self.param('Center Pixel Y').value()/binning]).copy()

            pixel_x=self.param('Pixel Size').value()*binning
            pixel_y=self.param('Pixel Size').value()*binning

        acceleration=window.p.param('Data Processing', 'Acceleration').value() 
        correct_refraction=window.p.param('Data Processing', 'Correct for Refraction').value()
        correct=window.p.param('Data Processing', 'Apply Intensity Corrections').value()
        bounds = np.asarray(window.binBounds,dtype=np.single)
        projection=window.p.param('View Mode', 'Select Projection').value()
        norm_hist=window.p.param('Data Processing', 'Divide Bins By Frequency').value() 
        
        angi_c = np.deg2rad(self.param('Critical Angle').value())           
        from_image = int(np.floor(window.image_stack.angle2image(window.angleRegion.getRegion()[0])))

        SDD=self.param('Sample-Detector Dist.').value()
        invSDD=1/self.param('Sample-Detector Dist.').value()        
        k0=(np.pi*2)/self.param('Wavelength').value()
        p_h=self.param('Horizontal Polarisation').value()
        Binv=window.crystal.calcBInv()
        recp_units=window.crystal.param('Reciprical Units').value()
        
        gridsize=int(window.p.param('Data Processing', 'Grid Size').value())
        index_offset = int(math.ceil(gridsize/2))
        #currently it is just a square grid so it needs to be big enough
        if qzmax > qxymax:
            qxymax = qzmax
        invgridstepx = 1/(2*qxymax/gridsize)       
        invgridstepy = 1/(2*qxymax/gridsize)   
 
        if  acceleration==2:
            if "cuda" not in sys.modules:
                from numba import cuda  

            @cuda.jit(fastmath=True)
            def  _bin_kernel(image, midang, histogram_result,histogram_weights,bin_bounds):
                #get the current thread position
                row,col = cuda.grid(2)

                if row < image.shape[0] and col < image.shape[1]:
                    delta_z= (q0[1]-row)*pixel_y #real-space dinstance from centre pixel y
                    delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                    delR = math.sqrt(delta_x**2 + delta_z**2)            
                    dist = math.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                    #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                    #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007      

                    gam_pix  = math.atan(delta_x*invSDD)
                    del_pix =  math.atan(delta_z/math.sqrt(delta_x**2 + SDD**2))-angi*math.cos(gam_pix)                                    
                                     
                    qx = k0*(math.cos(del_pix)*math.cos(gam_pix)-math.cos(angi))                    
                    qy = k0*(math.cos(del_pix)*math.sin(gam_pix)) 
                    qz = k0*(math.sin(del_pix)+math.sin(angi))

                    #refraction correction
                    # Busch et al., J. Appl. Cryst. 39, 433-442 (2006)
                    if correct_refraction:
                        kzi = k0*math.sqrt(abs(math.sin(angi)**2 - math.sin(angi_c)**2))
                        kzf = k0*math.sqrt(abs(math.sin(del_pix)**2 - math.sin(angi_c)**2))
                        qz = kzf - kzi
                    
                    if correct:
                        Phor = 1- (math.cos(del_pix)*math.sin(gam_pix))**2
                        Pver = 1 - math.sin(del_pix)**2
                        pol = p_h*Phor + (1-p_h)*Pver
                        L= 1#(math.sin(del_pix)*math.cos(angi))*math.copysign(1,del_pix)  
                        #L=math.sin(gam_pix)*math.cos(del_pix)*math.cos(angi)                        
                        Cd = dist**2/SDD**2    
                        Ci = 1/math.cos(math.atan(delR/SDD))                            
                        Crod = (math.sin(angi-gam_pix)**2 * math.cos(del_pix) + math.cos(angi)*math.cos(angi-gam_pix))/math.cos(angi)         
                        cfac = L*Crod*Cd*pol*Ci
                    else:
                        cfac = 1
                
                    """apply sample rotations to frame of reference
                    rotations are phi and omega_h which is angle of incidence
                    this is eqn 44. in ref. 1"""
                    so = math.sin(rotDirection*(omega_rad+midang))
                    co = math.cos(rotDirection*(omega_rad+midang))

                    hphi_1 = so*qy+co*qx
                    hphi_2 = co*qy-so*qx
                    hphi_3 = qz

                    if recp_units == 1:
                        h=  Binv[0][0]*hphi_1+Binv[0][1]*hphi_2+Binv[0][2]*hphi_3
                        k = Binv[1][0]*hphi_1+Binv[1][1]*hphi_2+Binv[1][2]*hphi_3
                        l = Binv[2][0]*hphi_1+Binv[2][1]*hphi_2+Binv[2][2]*hphi_3
                    else:
                        h = hphi_1
                        k = hphi_2
                        l = hphi_3                  
 
                    #Now we we have the coordinates of our pixel we need to decide which bin to put it
                    #Projection 0: hk - we assume the bounds correspond qr and l 
                    if projection == 0:                                           
                        xbound = math.sqrt(h**2 + k**2)*math.copysign(1,gam_pix)
                        ybound = l
                        xcord = h
                        ycord = k
                    #Projection 1: h-l, this is selected from an in-plane map so should be h k, and output to h and l
                    elif projection == 1:                        
                        xbound = h
                        ybound = k
                        xcord = h
                        ycord = l
                    
                    #Projection 2: k-l, this is selected from an in-plane map so should be h k, and output to k and l
                    elif projection == 2:
                        xbound = h
                        ybound = k
                        xcord = k
                        ycord = l
                    #Projection 3: qr-l, this is selected from an in-plane map so should be h k, and output to qr and l
                    else:
                        xbound = h
                        ybound = k
                        xcord = math.sqrt(h**2 + k**2)
                        ycord = l                       


                    #qr = math.sqrt(h**2 + k**2)*math.copysign(1,gam_pix)
                    for bound_index in range(len(bin_bounds)):                            
                        if xbound > bin_bounds[bound_index][0] and xbound < bin_bounds[bound_index][1]:   
                            if ybound > bin_bounds[bound_index][2] and ybound < bin_bounds[bound_index][3]:                         
                                histx = int(round(xcord*invgridstepx))+index_offset
                                histy = int(round(ycord*invgridstepy))+index_offset  
                                intensity = image[row,col]*cfac
                                cuda.atomic.add(histogram_result,(histx,histy),intensity)
                                cuda.atomic.add(histogram_weights,(histx,histy),1)

            # Configure the blocks
            threadsperblock = (32, 8)
            blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)            
            histogram_result  = cuda.to_device(np.zeros((gridsize,gridsize),dtype=np.single)) 
            histogram_weights = cuda.to_device(np.zeros((gridsize,gridsize),dtype=np.single))            

            for i, angle in enumerate(angles):                 
                if use_raw_files:              
                    #we can also use every pixel but we have to load each image again one at a time
                    #however we should not keep them all in memory.       
                    tmp = np.ascontiguousarray(np.rot90(window.image_stack.get_image_unbinned(i+from_image)))
                else:
                    tmp = np.ascontiguousarray(np.rot90(imagedata[i]))
                imagei = cuda.to_device(tmp)  
                print("Calculating angle ",i, " = ", angle)  
                _bin_kernel[blockspergrid, threadsperblock](imagei, np.deg2rad(angle), histogram_result, histogram_weights,bounds)          
            #normalise histogram
            #cuda.synchronize()    
            hist = np.asarray(histogram_result.copy_to_host())     
            weight = np.asarray(histogram_weights.copy_to_host())     
            if norm_hist:
                hist2 = np.divide(hist, weight, where=weight!=0)                
            else:
                hist2=hist
            return hist2

        else:     
            def _binner(image, angle):
                hist = np.zeros((gridsize,gridsize),dtype=np.single)
                weights = np.zeros((gridsize,gridsize),dtype=np.single) 
                """Function calculate the pixel positon in lab frame"""
                for row in range(pixel_count_y):    
                    for col in range(pixel_count_x):
                        delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                        #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007        
                        gam_pix  = np.arctan(delta_x/ SDD)
                        del_pix = np.arctan(delta_z/np.sqrt(delta_x**2 + SDD**2))-angi*np.cos(gam_pix)   

                        qx = k0*(np.cos(del_pix)*np.cos(gam_pix)-np.cos(angi))
                        qy = k0*(np.cos(del_pix)*np.sin(gam_pix)) 
                        qz = k0*(np.sin(del_pix)+np.sin(angi))

                        if correct_refraction:
                            kzi = k0*np.sqrt(abs(np.sin(angi)**2 - np.sin(angi_c)**2))
                            kzf = k0*np.sqrt(abs(np.sin(del_pix)**2 - np.sin(angi_c)**2))
                            qz = kzf - kzi
                        if correct:
                            Phor = 1- (np.cos(del_pix)*np.sin(gam_pix))**2
                            Pver = 1 - np.sin(del_pix)**2
                            pol = p_h*Phor + (1-p_h)*Pver   
                            #L= (np.sin(del_pix)*np.cos(angi))*np.copysign(1,del_pix)
                            L=math.sin(gam_pix)*math.cos(del_pix)*math.cos(angi)
                            Cd = dist**2/SDD**2
                            Ci = 1/np.cos(math.atan(delR/SDD))                            
                            Crod = (np.sin(angi-gam_pix)**2 * np.cos(del_pix) + np.cos(angi)*np.cos(angi-gam_pix))/np.cos(angi)              
                    
                            cfac = L*Crod*Cd*pol*Ci
                        else:
                            cfac = 1

                        so = np.sin(rotDirection*(omega_rad+angle))
                        co = np.cos(rotDirection*(omega_rad+angle))
                        ci = 1
                        si = 0  
                        """apply sample rotations to frame of reference
                            rotations are phi and omega_h which is angle of incidence
                            this is eqn 44. in ref. 1"""

                        hphi_1 = so*(ci*qy+si*qz)+co*qx
                        hphi_2 = co*(ci*qy+si*qz)-so*qx
                        hphi_3 = ci*qz-si*qy  

                        if recp_units == 1:
                            h=  Binv[0,0]*hphi_1+Binv[0,1]*hphi_2+Binv[0,2]*hphi_3
                            k = Binv[1,0]*hphi_1+Binv[1,1]*hphi_2+Binv[1,2]*hphi_3
                            l = Binv[2,0]*hphi_1+Binv[2,1]*hphi_2+Binv[2,2]*hphi_3
                        else:
                            h = hphi_1
                            k = hphi_2
                            l = hphi_3

                        #Now we we have the coordinates of our pixel we need to decide which bin to put it
                        #Projection 0: hk - we assume the bounds correspond qr and l 
                        if projection == 0:                                           
                            xbound = math.sqrt(h**2 + k**2)*math.copysign(1,gam_pix)
                            ybound = l
                            xcord = h
                            ycord = k
                        #Projection 1: h-l, this is selected from an in-plane map so should be h k, and output to h and l
                        elif projection == 1:                        
                            xbound = h
                            ybound = k
                            xcord = h
                            ycord = l
                        
                        #Projection 2: k-l, this is selected from an in-plane map so should be h k, and output to k and l
                        elif projection == 2:
                            xbound = h
                            ybound = k
                            xcord = k
                            ycord = l
                        #Projection 3: qr-l, this is selected from an in-plane map so should be h k, and output to qr and l
                        else:
                            xbound = h
                            ybound = k
                            xcord = math.sqrt(h**2 + k**2)
                            ycord = l

                            
                        #qr = math.sqrt(h**2 + k**2)*math.copysign(1,gam_pix)
                        for bound_index in range(len(bounds)):                            
                            if xbound > bounds[bound_index][0] and xbound < bounds[bound_index][1]:   
                                if ybound > bounds[bound_index][2] and ybound < bounds[bound_index][3]:                         
                                    histx = int(round(xcord*invgridstepx))+index_offset
                                    histy = int(round(ycord*invgridstepy))+index_offset  
                                    intensity = image[row,col]*cfac
                                    hist[histx][histy] += intensity
                                    weights[histx][histy] += 1

                return hist,weights  

            #use numba on cpu instead
            if acceleration == 1:
                if "jit" not in sys.modules:                
                    from numba import jit
                binner = jit(_binner,nopython=True)
            else:
                print("Not even Numba ...  reeallly sloooww")
                binner =_binner

            final_hist = np.zeros((gridsize,gridsize),dtype=np.single)
            final_weights = np.zeros((gridsize,gridsize),dtype=np.single) 
            
            for i, angle in enumerate(angles): 
                print("Calculating angle ",i, " = ", angle)
                if use_raw_files:           
                    image = np.rot90(window.image_stack.get_image_unbinned(i+from_image))   

                else:
                    image = np.rot90(imagedata[i])

                image_bin,image_weights = binner(image,np.deg2rad(angle))

                final_hist= final_hist + image_bin
                final_weights = final_weights + image_weights

            if norm_hist:
                hist2 = np.divide(final_hist, final_weights, where=final_weights!=0)                
            else:
                hist2=final_hist
            return hist2    

    def getLambda(self,energy):
        h=6.62607004e-34
        e=1.60217662e-19
        c=299792458
        return h*c/(energy*e)*1e10

    def geteV(self,lambdas):
        h=6.62607004e-34
        e=1.60217662e-19
        c=299792458
        return (h*c/lambdas*1e10)/e
    
    def __init__(self, **opts):
        opts['type'] = 'bool'
        opts['value'] = True
        pTypes.GroupParameter.__init__(self, **opts)
        
        self.addChild({'name': 'Energy', 'type': 'float', 'value': 67000,'siPrefix': True, 'suffix': 'eV', 'step': 100})
        self.addChild({'name': 'Wavelength', 'type': 'float', 'value': 1., 'suffix': 'Å', 'step': 0.001})
        self.addChild({'name': 'Sample-Detector Dist.', 'type': 'float', 'value': 1.26, 'suffix': 'm', 'siPrefix': True, 'step': 0.0001,'dec': True})
        self.addChild({'name': 'Horizontal Polarisation', 'type': 'float', 'value': 0.98})
        
        self.addChild({'name': 'Center Pixel X', 'type': 'int', 'value': 1424, 'step': 1})
        self.addChild({'name': 'Center Pixel Y', 'type': 'int', 'value': 2730, 'step': 1})
        self.addChild({'name': 'Pixel Size', 'type': 'float', 'value': 150e-6, 'step': 1e-6, 'suffix': 'm', 'siPrefix': True})
        self.addChild({'name': 'X Pixels', 'type': 'int', 'value': 2880, 'step': 1, })
        self.addChild({'name': 'Y Pixels', 'type': 'int', 'value': 2880, 'step': 1, })
        self.addChild({'name': 'Detector Gap (pixels)', 'type': 'int', 'value': 1, 'step': 1, })
        self.addChild({'name': 'Critical Angle', 'type': 'float', 'value': 0.08, 'step': 0.01, })

        self.addChild({'name': 'Angle of Incidence', 'type': 'float', 'value': 0.07, 'suffix': '°', 'step': 0.001,'suffixGap': ''})
        self.addChild({'name': 'Axis directon', 'type': 'list', 'values': {"Positive": 1, "Negative": -1}, 'value': -1})
        self.addChild({'name': 'Angle offset', 'type': 'float', 'value': 0,'suffix': '°','suffixGap': '', 'step': 0.01})
        self.addChild({'name': 'Beamline Preset', 'type': 'list', 'values': {"P21.2 1": 1, "P21.2 2": 2, "P07 2019": 3,"ID31 HDF5": 4},'value': 1})

        self.energy = self.param('Energy')
        self.wavelength = self.param('Wavelength')
        self.energy.sigValueChanged.connect(self.energyChanged)
        self.wavelength.sigValueChanged.connect(self.wavelengthChanged)
        self.wavelength.setValue(self.getLambda(self.energy.value()))

    def energyChanged(self):
        self.wavelength.setValue(self.getLambda(self.energy.value()), blockSignal=self.wavelengthChanged)

    def wavelengthChanged(self):
        self.energy.setValue(self.geteV(self.wavelength.value()), blockSignal=self.energyChanged)
