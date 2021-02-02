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
        self.addChild({'name': 'a₁', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.01})
        self.addChild({'name': 'a₂', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.01})
        self.addChild({'name': 'a₃', 'type': 'float', 'value': 4.081,'suffix': 'Å', 'step': 0.01})

        self.addChild({'name': 'α₁', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'α₂', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'α₃', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})

        self.addChild({'name': 'b₁', 'type': 'float', 'value': 2.52,'suffix': 'Å⁻¹', 'step': 0.01})
        self.addChild({'name': 'b₂', 'type': 'float', 'value': 2.53,'suffix': 'Å⁻¹', 'step': 0.01})
        self.addChild({'name': 'b₃', 'type': 'float', 'value': 0.89,'suffix': 'Å⁻¹', 'step': 0.01})

        self.addChild({'name': 'β₁', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'β₂', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})
        self.addChild({'name': 'β₃', 'type': 'float', 'value': 90,'suffix': '°','suffixGap': ''})

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
        angi=np.deg2rad(self.param('Angle of Incidence').value())
        omega_rad=np.deg2rad(self.param('Angle offset').value()+midangle)
        rotDirection=(self.param('Axis directon').value())
        pixel_count_x=int(self.param('X Pixels').value()/binning)
        pixel_count_y=int(self.param('Y Pixels').value()/binning)
        acceleration=window.p.param('Data Processing', 'Acceleration').value() 
        correct_refraction=window.p.param('Data Processing', 'Correct for refraction').value() 
        angi_c = np.deg2rad(self.param('Critical Angle').value()*binning)
        if twodetectors:
            pixel_count_x=int((2*self.param('X Pixels').value()+self.param('Detector Gap (pixels)').value())/binning)
        pixel_x=self.param('Pixel Size').value()*binning
        pixel_y=self.param('Pixel Size').value()*binning
        
        SDD=self.param('Sample-Detector Dist.').value()
        q0=np.array([self.param('Center Pixel X').value()/binning,self.param('Center Pixel Y').value()/binning]).copy()
        k0=(np.pi*2)/self.param('Wavelength').value()
        p_h=self.param('Horizontal Polarisation').value()
        Binv=window.crystal.calcBInv()
        b1=window.crystal.param('b₁').value()

        if  acceleration==2:
            if "cuda" not in sys.modules:
                from numba import cuda
            @cuda.jit
            def detector_to_hkl_kernel(h_glob,k_glob,l_glob,hk_glob,c_glob):
                #get the current thread position
                j,i = cuda.grid(2)

                if j < h_glob.shape[0] and i < h_glob.shape[1]:
                    delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (i-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                    delR = math.sqrt(delta_x**2 + delta_z**2)            
                    dist = math.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                    #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                    #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007            

                    del_pix  = math.atan(delta_x/ SDD)
                    gam_pix = math.atan(delta_z/math.sqrt(delta_x**2 + SDD**2))-angi*math.cos(del_pix)                                    
                                  
                    tth  =  math.acos(math.cos(del_pix)* math.cos(gam_pix))    

                    qx = k0*(math.cos(gam_pix)*math.cos(del_pix)-math.cos(angi))
                    qy = k0*(math.cos(gam_pix)*math.sin(del_pix)) 
                    qz = k0*(math.sin(gam_pix)+math.sin(angi))

                    #refraction correction
                    # Busch et al., J. Appl. Cryst. 39, 433-442 (2006)
                    if correct_refraction:
                        kzi = k0*math.sqrt(abs(math.sin(angi)**2 - math.sin(angi_c)**2))
                        kzf = k0*math.sqrt(abs(math.sin(gam_pix)**2 - math.sin(angi_c)**2))
                        qz = kzf - kzi

                    #polarization correction
                    P = p_h*(1-math.cos(del_pix)**2*math.sin(gam_pix)**2)+(1-p_h)*(1-math.sin(del_pix)**2)           
                    
                    #Lorentz factor
                    if abs(del_pix) <= 0.005:
                        L = 1/math.sin(2*(angi+del_pix)) #z-axis reflectivity
                    else:
                        L = math.cos(tth/2)/math.sin(tth)
                    L= 1/L                                

                    #correction factor for change in distance 
                    #due to flat detector (open-slit)
                    Cd = dist**2/SDD**2    
                    #correction factor for projected pixel size 
                    #due to beam inclination (open-slit)    
                    Ci = 1/math.cos(math.atan(delR/SDD))
                    
                    #Rod interception z-axis mode
                    Crod = 1/math.cos(abs(gam_pix))           

                    c_glob[j,i] = Ci*Cd*Crod*P*L
                    
                    """apply sample rotations to frame of reference
                    rotations are phi and omega_h which is angle of incidence
                    this is eqn 44. in ref. 1"""
                    so = math.sin(rotDirection*omega_rad)
                    co = math.cos(rotDirection*omega_rad)
                    # we deal with the angle of incidence in the momentum transfer calc
                    ci = 1 #math.cos(angi) 
                    si = 0 #math.sin(angi)  

                    hphi_1 = so*(ci*qy+si*qz)+co*qx
                    hphi_2 = co*(ci*qy+si*qz)-so*qx
                    hphi_3 = ci*qz-si*qy
                    
                    #H= Hphi #should be Binv dot Hphi 
                    # compute the dot product manual since cuda doens't
                    # support np.dot       
                    h_glob[j,i] = Binv[0][0]*hphi_1+Binv[0][1]*hphi_2+Binv[0][2]*hphi_3
                    k_glob[j,i] = Binv[1][0]*hphi_1+Binv[1][1]*hphi_2+Binv[1][2]*hphi_3
                    l_glob[j,i] = Binv[2][0]*hphi_1+Binv[2][1]*hphi_2+Binv[2][2]*hphi_3

                    #this works for non-orthongal coorindates
                    hk_glob[j,i] = math.copysign(1,del_pix)*math.sqrt((hphi_1/b1)**2 + (hphi_2/b1)**2)                     
                    
                
            h_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            k_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            l_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x))) 
            hk_global_mem = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            c_global_mem  = cuda.to_device(np.ones((pixel_count_y,pixel_count_x))) #array of pixel intensity corrections           

            # Configure the blocks
            threadsperblock = (16, 16)
            blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)
            import time
            t0 = time.time()
            detector_to_hkl_kernel[blockspergrid, threadsperblock](h_global_mem,k_global_mem,l_global_mem,hk_global_mem,c_global_mem) 
            t1 = time.time()  
            print("inner",t1-t0)
            return [h_global_mem.ravel(),k_global_mem.ravel(),l_global_mem.ravel(),hk_global_mem.ravel(),c_global_mem.copy_to_host()]  

        #use numba on cpu instead
        if acceleration == 1:
            print("ONLY GPU ACCELERATION IS CURRENTLY FUNCTINAL!!")
            return
            if "jit" not in sys.modules:                
                from numba import jit
            @jit(nopython=True)
            def _dector_frame_to_hkl_frame():
                h_img = np.zeros((pixel_count_y,pixel_count_x)) 
                k_img = np.zeros((pixel_count_y,pixel_count_x))
                l_img = np.zeros((pixel_count_y,pixel_count_x)) 
                hk_img = np.zeros((pixel_count_y,pixel_count_x))
                c_img = np.ones((pixel_count_y,pixel_count_x)) #array of pixel intensity corrections
                
                """Function calculate the pixel positon in lab frame"""
                for j in range(pixel_count_y):    
                    for i in range(pixel_count_x):
                        delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (q0[0]-i)*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        
                        del_pix =  np.arctan(delta_z/SDD)
                        gam_pix = np.arcsin(delta_x/dist) 


                        qx = -1*k0*(np.cos(del_pix)*np.sin(gam_pix))
                        qy = k0*(np.cos(del_pix)*np.cos(gam_pix)-1)    
                        qz = k0*(np.sin(del_pix)) 
                       
                        tth = np.arccos(np.cos(del_pix)*np.cos(gam_pix))
                        
                        #polarization correction
                        P = p_h*(1-np.cos(del_pix)**2*np.sin(gam_pix)**2)+(1-p_h)*(1-np.sin(del_pix)**2)           
                        
                        #Lorentz factor
                        if abs(tth) <= 0.005:
                            L = 1/np.sin(2*(angi+del_pix)) #z-axis reflectivity
                        else:
                            L = np.cos(tth/2)/np.sin(tth)
                        L= 1/L
                                    
                        #correction factor for change in distance 
                        #due to flat detector (open-slit)
                        Cd = dist**2/SDD**2    
                        #correction factor for projected pixel size 
                        #due to beam inclination (open-slit)    
                        Ci = 1/np.cos(np.arctan(delR/SDD))
                        
                        #Rod interception z-axis mode
                        Crod = 1/np.cos(abs(gam_pix))           
                    
                        c_img[j][i] = Ci*Cd*Crod*P*L
                        
                        """apply sample rotations to frame of reference
                        rotations are phi and omega_h which is angle of incidence
                        this is eqn 44. in ref. 1"""
                        so = np.sin(rotDirection*omega_rad)
                        co = np.cos(rotDirection*omega_rad)
                        ci = np.cos(angi)
                        si = np.sin(angi)  
                        
                        Hphi = np.array([[so*(ci*qy+si*qz)+co*qx],
                                        [co*(ci*qy+si*qz)-so*qx],
                                        [ci*qz-si*qy]])
                        
                        H= Binv.dot(Hphi)
                        h_img[j][i] = H[0][0]
                        k_img[j][i] = H[1][0]
                        l_img[j][i] = H[2][0]
                        
                    #this works for non-orthongal coorindates
                        hk_img[j][i] = np.sign(gam_pix)*np.sqrt((Hphi[0][0]/b1)**2 + (Hphi[1][0]/b1)**2)
                return [h_img,k_img,l_img,hk_img,c_img]
        else:
            print("ONLY GPU ACCELERATION IS CURRENTLY FUNCTINAL!!")
            return
            def _dector_frame_to_hkl_frame():
                h_img = np.zeros((pixel_count_y,pixel_count_x)) 
                k_img = np.zeros((pixel_count_y,pixel_count_x))
                l_img = np.zeros((pixel_count_y,pixel_count_x)) 
                hk_img = np.zeros((pixel_count_y,pixel_count_x))
                c_img = np.ones((pixel_count_y,pixel_count_x)) #array of pixel intensity corrections

                """Function calculate the pixel positon in lab frame"""
                for j in range(pixel_count_y):    
                    for i in range(pixel_count_x):
                        delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (q0[0]-i)*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        
                        del_pix =  np.arctan(delta_z/SDD)
                        gam_pix = np.arcsin(delta_x/dist)  
                
                        qx =-1* k0*(np.cos(del_pix)*np.sin(gam_pix))
                        qy = k0*(np.cos(del_pix)*np.cos(gam_pix)-1)    
                        qz = k0*(np.sin(del_pix)) 
                        
                        tth = np.arccos(np.cos(del_pix)*np.cos(gam_pix))
                        
                        #polarization correction
                        P = p_h*(1-np.cos(del_pix)**2*np.sin(gam_pix)**2)+(1-p_h)*(1-np.sin(del_pix)**2)           
                        
                        #Lorentz factor
                        if abs(tth) <= 0.005:
                            L = 1/np.sin(2*(angi+del_pix)) #z-axis reflectivity
                        else:
                            L = np.cos(tth/2)/np.sin(tth)
                        L= 1/L
                                    
                        #correction factor for change in distance 
                        #due to flat detector (open-slit)
                        Cd = dist**2/SDD**2    
                        #correction factor for projected pixel size 
                        #due to beam inclination (open-slit)    
                        Ci = 1/np.cos(np.arctan(delR/SDD))
                        
                        #Rod interception z-axis mode
                        Crod = 1/np.cos(abs(gam_pix))           
                    
                        c_img[j][i] = Ci*Cd*Crod*P*L
                        
                        """apply sample rotations to frame of reference
                        rotations are phi and omega_h which is angle of incidence
                        this is eqn 44. in ref. 1"""
                        so = np.sin(rotDirection*omega_rad)
                        co = np.cos(rotDirection*omega_rad)
                        ci = np.cos(angi)
                        si = np.sin(angi)  
                        
                        Hphi = np.array([[so*(ci*qy+si*qz)+co*qx],
                                        [co*(ci*qy+si*qz)-so*qx],
                                        [ci*qz-si*qy]])
                        
                        H= Binv.dot(Hphi)
                        h_img[j][i] = H[0][0]
                        k_img[j][i] = H[1][0]
                        l_img[j][i] = H[2][0]
                        
                    #this works for non-orthongal coorindates
                        hk_img[j][i] = np.sign(gam_pix)*np.sqrt((Hphi[0][0]/b1)**2 + (Hphi[1][0]/b1)**2)                      
    
                return [h_img,k_img,l_img,hk_img,c_img]
        import time

        t0 = time.time()
        tmp = np.array(_dector_frame_to_hkl_frame())
        t1 = time.time()
        print(t1-t0)
        return (tmp)

    def dector_frame_to_lab_frame(self,window,binning):
        angi=np.deg2rad(self.param('Angle of Incidence').value())
        omega_rad=np.deg2rad(self.param('Angle offset').value())
        rotDirection=(self.param('Axis directon').value())
        pixel_count_x=int(self.param('X Pixels').value()/binning)
        pixel_count_y=int(self.param('Y Pixels').value()/binning)
        acceleration=window.p.param('Data Processing', 'Acceleration').value() 
        correct_refraction=window.p.param('Data Processing', 'Correct for refraction').value() 
        angi_c = np.deg2rad(self.param('Critical Angle').value()*binning)
        if window.p.param('Data Processing', 'Use 2nd detector').value() :
            pixel_count_x=int((2*self.param('X Pixels').value()+self.param('Detector Gap (pixels)').value())/binning)
        pixel_x=self.param('Pixel Size').value()*binning
        pixel_y=self.param('Pixel Size').value()*binning
        
        SDD=self.param('Sample-Detector Dist.').value()
        q0=np.array([self.param('Center Pixel X').value()/binning,self.param('Center Pixel Y').value()/binning]).copy()
        k0=(np.pi*2)/self.param('Wavelength').value()
        p_h=self.param('Horizontal Polarisation').value()
        Binv=window.crystal.calcBInv()
        b1=window.crystal.param('b₁').value()


        if  acceleration==2:
            if "cuda" not in sys.modules:
                from numba import cuda
            @cuda.jit
            def detector_to_hkl_kernel(qx_glob,qy_glob,qz_glob):
                #get the current thread position
                j,i = cuda.grid(2)

                if j < qx_glob.shape[0] and i < qx_glob.shape[1]:
                    delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (i-q0[0])*pixel_x  #real-space dinstance from centre pixel x
                    #delR = math.sqrt(delta_x**2 + delta_z**2)            
                    #dist = math.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                    #https://www.classe.cornell.edu/~dms79/D-lineNotes/GISAXS-at-D-line/GIWAXS@D1_Conversions.html
                    #i.e Smilgies & Blasini, J. Appl. Cryst. 40, 716-718 (2007            

                    del_pix  = math.atan(delta_x/ SDD)
                    gam_pix = math.atan(delta_z/math.sqrt(delta_x**2 + SDD**2))-angi*math.cos(del_pix)                                    
                                  
                    tth  =  math.acos(math.cos(del_pix)* math.cos(gam_pix))    

                    qx_glob[j,i] = k0*(math.cos(gam_pix)*math.cos(del_pix)-math.cos(angi))
                    qy_glob[j,i] = k0*(math.cos(gam_pix)*math.sin(del_pix)) 
                    qz_glob[j,i] = k0*(math.sin(gam_pix)+math.sin(angi))

          
            qx_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            qy_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
            qz_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))   

            # Configure the blocks
            threadsperblock = (16, 16)
            blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)

            detector_to_hkl_kernel[blockspergrid, threadsperblock](qx_global_mem,qy_global_mem,qz_global_mem) 
            return [qx_global_mem.copy_to_host(),qy_global_mem.copy_to_host(),qz_global_mem.copy_to_host()]  


        if  acceleration==1:
            print("ONLY GPU ACCELERATION IS CURRENTLY FUNCTIONAL")
            return
            from numba import jit

            @jit(nopython=True)    
            def _dector_frame_to_lab_frame():
                qx_img = np.zeros((pixel_count_x,pixel_count_y)) 
                qy_img = np.zeros((pixel_count_x,pixel_count_y))
                qz_img = np.zeros((pixel_count_x,pixel_count_y)) 
                
                """Function calculate the pixel positon in lab frame"""
                for j in range(pixel_count_y):    
                    for i in range(pixel_count_x):
                        delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (q0[0]-i)*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        
                        del_pix =  np.arctan(delta_z/SDD)
                        gam_pix = np.arcsin(delta_x/dist)  
                
                        qx =-1* k0*(np.cos(del_pix)*np.sin(gam_pix))
                        qy = k0*(np.cos(del_pix)*np.cos(gam_pix)-1)    
                        qz = k0*(np.sin(del_pix))  

                        qx_img[j][i] = qx
                        qy_img[j][i] = qy
                        qz_img[j][i] = qz

                return [qx_img,qy_img,qz_img]
        else:
            print("ONLY GPU ACCELERATION IS CURRENTLY FUNCTIONAL")
            return
            def _dector_frame_to_lab_frame():
                qx_img = np.zeros((pixel_count_x,pixel_count_y)) 
                qy_img = np.zeros((pixel_count_x,pixel_count_y))
                qz_img = np.zeros((pixel_count_x,pixel_count_y)) 
                
                """Function calculate the pixel positon in lab frame"""
                for j in range(pixel_count_y):    
                    for i in range(pixel_count_x):
                        delta_z= (q0[1]-j)*pixel_y  #real-space dinstance from centre pixel y
                        delta_x = (q0[0]-i)*pixel_x  #real-space dinstance from centre pixel x
                        delR = np.sqrt(delta_x**2 + delta_z**2)            
                        dist = np.sqrt(delta_x**2+SDD**2 + delta_z**2) #distance to pixel
                        
                        del_pix =  np.arctan(delta_z/SDD)
                        gam_pix = np.arcsin(delta_x/dist)  
                
                        qx =-1* k0*(np.cos(del_pix)*np.sin(gam_pix))
                        qy = k0*(np.cos(del_pix)*np.cos(gam_pix)-1)    
                        qz = k0*(np.sin(del_pix))  

                        qx_img[j][i] = qx
                        qy_img[j][i] = qy
                        qz_img[j][i] = qz

                return [qx_img,qy_img,qz_img]
        return(np.array(_dector_frame_to_lab_frame()))

    def get_coords(self,corrected_ang,q,window):
        Binv=window.crystal.calcBInv()
        qx_list=q[0]
        qy_list=q[1]
        qz_list=q[2]
        angi=np.deg2rad(self.param('Angle of Incidence').value())
        angleoff=np.deg2rad(self.param('Angle offset').value())
        rotDirection=(self.param('Axis directon').value())
        corrected_angles=corrected_ang
        acceleration=window.p.param('Data Processing', 'Acceleration').value() 

        if  acceleration==2 or acceleration==1:
            from numba import jit
            @jit(nopython=True) #save time and compile
            def _get_coords():
                points_h = []
                points_k = []
                points_l = []
                #angle of incidence now down in previous step
                ci = 1 #np.cos(angi)
                si = 0 #np.sin(angi)  
                
                for i, angle in enumerate(corrected_angles):
                    omega_rad = np.deg2rad(angle) + angleoff
                    so = np.sin(rotDirection*omega_rad) #the sign is important
                    co = np.cos(rotDirection*omega_rad)
                    
                    for i in range(len(qx_list)):
                        qx = qx_list[i]
                        qy = qy_list[i]
                        qz = qz_list[i]           
                        Hphi = np.array([[so*(ci*qy+si*qz)+co*qx],
                                        [co*(ci*qy+si*qz)-so*qx],
                                        [ci*qz-si*qy]])            
                        H= Binv.dot(Hphi)
                        points_h.append(H[0][0])
                        points_k.append(H[1][0])
                        points_l.append(H[2][0])  
                return np.array([points_h,points_k,points_l])
            return(_get_coords())
        else:
            def _get_coords():
                points_h = []
                points_k = []
                points_l = []
                ci = np.cos(angi)
                si = np.sin(angi)  
                
                for i, angle in enumerate(corrected_angles):
                    omega_rad = np.deg2rad(angle) + angleoff
                    so = np.sin(rotDirection*omega_rad) #the sign is important
                    co = np.cos(rotDirection*omega_rad)
                    
                    for i in range(len(qx_list)):
                        qx = qx_list[i]
                        qy = qy_list[i]
                        qz = qz_list[i]           
                        Hphi = np.array([[so*(ci*qy+si*qz)+co*qx],
                                        [co*(ci*qy+si*qz)-so*qx],
                                        [ci*qz-si*qy]])            
                        H= Binv.dot(Hphi)
                        points_h.append(H[0][0])
                        points_k.append(H[1][0])
                        points_l.append(H[2][0])  
                return np.array([points_h,points_k,points_l])
            return(_get_coords())


    def get_coords_np(self,corrected_ang,q,crystal):
        Binv=crystal.calcBInv()
        angi=-np.deg2rad(self.param('Angle of Incidence').value())
        angleoff=np.deg2rad(self.param('Angle offset').value())
        rotDirection=(self.param('Axis directon').value())
        corrected_angles=corrected_ang
        def _get_coords():
            ci = np.cos(angi)
            si = np.sin(angi)
            omega_rad = np.deg2rad(corrected_ang) + angleoff
            so = np.sin(rotDirection*omega_rad) #the sign is important
            co = np.cos(rotDirection*omega_rad)
            Hs=np.zeros((len(omega_rad),3))
            for i in range(len(omega_rad)):
                Hphi = np.array([[so[i]*(ci*q[1]+si*q[2])+co[i]*q[0]],
                                 [co[i]*(ci*q[1]+si*q[2])-so[i]*q[0]],
                                 [ci*q[2]-si*q[1]]])            
                H= Binv.dot(Hphi)
                Hs[i,:]=H[:,0]
            return (Hs)
        return(_get_coords())

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
        self.addChild({'name': 'Wavelength', 'type': 'float', 'value': 1., 'suffix': 'Å', 'step': 0.01})
        self.addChild({'name': 'Sample-Detector Dist.', 'type': 'float', 'value': 1.26, 'suffix': 'm', 'siPrefix': True, 'step': 0.01,'dec': True})
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
        self.addChild({'name': 'Angle offset', 'type': 'float', 'value': 0,'suffix': '°','suffixGap': ''})
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
