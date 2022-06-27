#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math 
import sys
import numpy as np

def dector_frame_to_theta2chi(window,binning,twodetectors):
    """Function to transform detector image into either q or RLU units"""
    experiment = window.experiment
    #collect some variables together
    alpha=np.deg2rad(experiment.param('Angle of Incidence').value())
    pixel_count_x=int(experiment.param('X Pixels').value()/binning)
    pixel_count_y=int(experiment.param('Y Pixels').value()/binning)
    acceleration=window.p.param('Data Processing', 'Acceleration').value() 
    correct_refraction=window.p.param('Data Processing', 'Correct for Refraction').value() 
    correct=window.p.param('Data Processing', 'Apply Intensity Corrections').value()
    critical_angle = np.deg2rad(experiment.param('Critical Angle').value())
    if twodetectors:
        pixel_count_x=int((2*experiment.param('X Pixels').value()+experiment.param('Detector Gap (pixels)').value())/binning)
    pixel_x=experiment.param('Pixel Size').value()*binning
    pixel_y=experiment.param('Pixel Size').value()*binning
    
    dist=experiment.param('Sample-Detector Dist.').value()
    q0=np.array([experiment.param('Center Pixel X').value()/binning,experiment.param('Center Pixel Y').value()/binning]).copy()
    k0=(np.pi*2)/experiment.param('Wavelength').value()
    b1=window.crystal.param('b₁').value()
    b2=window.crystal.param('b₂').value()
    b3=window.crystal.param('b₃').value()
    recp_units=window.crystal.param('Reciprocal Units').value()

    #if we are to use the GPU
    if  acceleration==2:
        if "cuda" not in sys.modules:
            #check if we need to load the library
            from numba import cuda
            
        @cuda.jit(fastmath=True)
        def detector_to_hkl_kernel(twotheta_glob,chi_glob,c_glob):
            #get the current thread position
            row,col = cuda.grid(2)

            if row < pixel_count_y and col < pixel_count_x:
                delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x 
                delta_x_2 = delta_x*delta_x
                dist_2 =dist*dist            
                
                #lab frame
                gam = math.atan(delta_x/dist)
                delta = math.atan(delta_z/math.sqrt(dist_2+delta_x_2))               
                #surface frame
                psi = gam
                beta = delta - alpha*math.cos(psi)
                
                sin_psi = math.sin(psi)
                cos_psi = math.cos(psi)
                sin_beta = math.sin(beta)
                cos_beta = math.cos(beta)
                cos_alpha = math.cos(alpha)
                sin_alpha = math.sin(alpha)
                cos_gam = math.cos(gam)
                cos_del = math.cos(delta)

                rad2deg = 180/3.1415926535
                theta2 = rad2deg*math.acos(cos_gam*cos_del)
                if math.tan(gam):
                    chi = rad2deg*(math.atan(math.tan(delta)/math.tan(gam)))
                            
                if correct:   
                    delta_z_2 = delta_z*delta_z                
                    delR = math.sqrt(delta_x_2 + delta_z_2)            
                    distR = math.sqrt(delta_x_2+dist_2 + delta_z_2) #distance to pixel
                    #Here we use the geometry independed correction factors from:
                    #doi.org/10.1063/1.1461876
                    Cpol = 1/(cos_beta*cos_beta*cos_psi*cos_psi+sin_beta*sin_beta)
                    Crod = 1/cos_beta
                    Cd = distR**2/dist**2    
                    Ci = 1./math.cos(math.atan(delR/dist))
                    #For integration around phi
                    Clor = cos_alpha*sin_psi*cos_beta 
                    C = Cpol*Crod*Cd*Ci*Clor
                else:
                    C = 1  

                chi_glob[row,col] = chi
                twotheta_glob[row,col] = theta2
                c_glob[row,col] = C
         
            
        chi_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
        twotheta_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x)))
        c_global_mem  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x))) #array of pixel intensity corrections           

        # Configure the blocks
        threadsperblock = (32, 8)
        blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
        blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
        blockspergrid = (blockspergrid_x, blockspergrid_y)
    
        detector_to_hkl_kernel[blockspergrid, threadsperblock](twotheta_global_mem,chi_global_mem,c_global_mem)        
        return [twotheta_global_mem.copy_to_host(),chi_global_mem.copy_to_host(),c_global_mem.copy_to_host()]  

    #Here we define the none gpu based function 
    def _dector_frame_to_hkl_frame():
            chi_img = np.zeros((pixel_count_y,pixel_count_x)) 
            theta2_img = np.zeros((pixel_count_y,pixel_count_x))            
            c_img = np.zeros((pixel_count_y,pixel_count_x)) #array of pixel intensity corrections    
            
            """Function calculate the pixel positon in lab frame"""
            for row in range(pixel_count_y):    
                for col in range(pixel_count_x):
                    delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x 
                    delta_x_2 = delta_x*delta_x
                    dist_2 =dist*dist            
                    
                    #lab frame
                    gam = math.atan(delta_x/dist)
                    delta = math.atan(delta_z/math.sqrt(dist_2+delta_x_2))               
                    #surface frame
                    psi = gam
                    beta = delta - alpha*math.cos(psi)
                    
                    sin_psi = math.sin(psi)
                    cos_psi = math.cos(psi)
                    sin_beta = math.sin(beta)
                    cos_beta = math.cos(beta)
                    cos_alpha = math.cos(alpha)
                    sin_alpha = math.sin(alpha)
                    cos_gam = math.cos(gam)
                    cos_del = math.cos(delta)

                    rad2deg = 180/3.1415926535
                    theta2 = rad2deg*math.acos(cos_gam*cos_del)
                    if math.tan(gam):
                        chi = rad2deg*(math.atan(math.tan(delta)/math.tan(gam)))
                                   
                    if correct:   
                        delta_z_2 = delta_z*delta_z                
                        delR = math.sqrt(delta_x_2 + delta_z_2)            
                        distR = math.sqrt(delta_x_2+dist_2 + delta_z_2) #distance to pixel
                        #Here we use the geometry independed correction factors from:
                        #doi.org/10.1063/1.1461876
                        Cpol = 1/(cos_beta*cos_beta*cos_psi*cos_psi+sin_beta*sin_beta)
                        Crod = 1/cos_beta
                        Cd = distR**2/dist**2    
                        Ci = 1./math.cos(math.atan(delR/dist))
                        #For integration around phi
                        Clor = cos_alpha*sin_psi*cos_beta 
                        C = Cpol*Crod*Cd*Ci*Clor
                    else:
                        C = 1  
 
                    chi_img[row,col] = chi
                    theta2_img[row,col] = theta2
                    c_img[row,col] = C

            return [theta2_img,chi_img,c_img]

    #use numba on cpu instead
    if acceleration == 1:
        if "jit" not in sys.modules: 
            #if numba not loaded load it               
            from numba import jit
        detector_frame_to_hkl=jit(_dector_frame_to_hkl_frame,nopython=True)
    else:
        detector_frame_to_hkl=_dector_frame_to_hkl_frame

    tmp = np.array(detector_frame_to_hkl())
    
    return (tmp)