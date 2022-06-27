#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from concurrent.futures.process import _MAX_WINDOWS_WORKERS
import math 
import sys
import numpy as np

def projection_binner(window,angles,imagedata,progress_callback):

    #First we start by collecting variables and handles from the main windows
    crystal = window.crystal
    experiment = window.experiment
    
    alpha=np.deg2rad(experiment.param('Angle of Incidence').value())
    omega_rad=np.deg2rad(experiment.param('Angle offset').value())+np.deg2rad(experiment.param('Aux Angle offset').value())
    
    rotDirection=(experiment.param('Axis directon').value())
    use_raw_files= window.p.param('Data Processing', 'Bin From Full Images').value() 
    binning= window.p.param('Data Processing', 'Binning (X by X)').value() 
    twodetectors=window.p.param('Data Processing', 'Use 2nd detector').value()
    imgState = window.ImgState

    #some parameters will be different if we don't use the binned data in memory
    #in this the binner re-reads the file from disk
    if use_raw_files:
        print("RAW")
        pixel_count_x=int(experiment.param('X Pixels').value())
        pixel_count_y=int(experiment.param('Y Pixels').value())
        if twodetectors:
            pixel_count_x=int((2*experiment.param('X Pixels').value()+experiment.param('Detector Gap (pixels)').value()))
        q0=np.array([experiment.param('Center Pixel X').value(),experiment.param('Center Pixel Y').value()]).copy()
        pixel_x=experiment.param('Pixel Size').value()
        pixel_y=experiment.param('Pixel Size').value()
        
    else:
        pixel_count_x=int(experiment.param('X Pixels').value()/binning)
        pixel_count_y=int(experiment.param('Y Pixels').value()/binning)
        if twodetectors:
            pixel_count_x=int((2*experiment.param('X Pixels').value()+experiment.param('Detector Gap (pixels)').value())/binning)
        q0=np.array([experiment.param('Center Pixel X').value()/binning,experiment.param('Center Pixel Y').value()/binning]).copy()

        pixel_x=experiment.param('Pixel Size').value()*binning
        pixel_y=experiment.param('Pixel Size').value()*binning

    acceleration=window.p.param('Data Processing', 'Acceleration').value() 
    correct_refraction=window.p.param('Data Processing', 'Correct for Refraction').value()
    correct=window.p.param('Data Processing', 'Apply Intensity Corrections').value()
    bounds = np.asarray(window.binBounds,dtype=np.single)
    projection=window.p.param('Data Processing', 'Select Projection').value()
    
    critical_angle = np.deg2rad(experiment.param('Critical Angle').value())  
    from_image = int(np.floor(window.image_stack.angle2image(window.angleRegion.getRegion()[0])))

    dist=experiment.param('Sample-Detector Dist.').value()     
    invDist = 1/dist
    dist_2 = dist*dist
    k0=(np.pi*2)/experiment.param('Wavelength').value()
    Binv=crystal.calcBInv()
    recp_units=crystal.param('Reciprocal Units').value()
    b1=window.crystal.param('b₁').value()
    b2=window.crystal.param('b₂').value()
    b3=window.crystal.param('b₃').value()
    
    gridsizex=int(window.p.param('Data Processing', 'Grid Size Qx').value())
    gridsizey=int(window.p.param('Data Processing', 'Grid Size Qy').value())
    gridsizez=int(window.p.param('Data Processing', 'Grid Size Qz').value())

    xmax = window.xmax
    xmin = window.xmin
    ymax = window.ymax
    ymin = window.ymin
    zmax = window.zmax
    zmin = window.zmin

    #The transformed detector view has an xaxis of Qr so we need to be prepared for rotation
    #this is for transfomred detector view      
    #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
    if imgState == 1:
        xmax = 10
        ymax = 10
        xmin = -10
        ymin = -1


    elif imgState == 3:         
        qrs = []
        mask_list = window.mask_list
        #We loop throught he masks and find the maxium 
        for mask in mask_list:
            if mask[1] == 1:
                qrs.append(abs(mask[0].pos()[0]))
                qrs.append(abs(mask[0].size()[0] + mask[0].pos()[0]))

        #but we need to be careful here since we dont know the orinatation so it could be negative
        qrmax = np.max(qrs)
        xmax = qrmax
        ymax = qrmax
        xmin = -1*qrmax
        ymin = -1*qrmax
        window.qrmax = qrmax #we need this later


    gridstepx = abs((xmax-xmin)/gridsizex)
    gridstepy = abs((ymax-ymin)/gridsizey)
    gridstepz = abs((zmax-zmin)/gridsizez)

    #We add some extra padding incase any rounding erros cause overflows
    #This could be simplifed mathematically but we go for more readable here
    xmax = xmax + gridstepx
    xmin = xmin - gridstepx
    ymax = ymax + gridstepy
    ymin = ymin - gridstepy
    zmax = zmax + gridstepz
    zmin = zmin - gridstepz

    gridstepx = abs((xmax-xmin)/gridsizex)
    gridstepy = abs((ymax-ymin)/gridsizey)
    gridstepz = abs((zmax-zmin)/gridsizez)
   
    #the inverses of the the stepsizes to avoid division later on       
    invgridstepx = 1/gridstepx     
    invgridstepy = 1/gridstepy
    invgridstepz = 1/gridstepz 
    
    #we need the grid size for a given projection

    #Transformed Detector
    if imgState== 1:
        grid_nx = gridsizex
        grid_ny = gridsizez

    #HK
    if imgState== 3:
        grid_nx = gridsizex
        grid_ny = gridsizey
    #HL
    elif imgState == 4:
        grid_nx = gridsizex
        grid_ny = gridsizez
    #KL
    elif imgState == 5:  
        grid_nx = gridsizey
        grid_ny = gridsizez
    #3d
    elif imgState== 6:  
        grid_nx = gridsizex
        grid_ny = gridsizey
        
    grid_nz = gridsizez #this is always true    

    if  acceleration==2:
            if "cuda" not in sys.modules:
                from numba import cuda  

            @cuda.jit(fastmath=True)
            def  _lab_frame(qx_vals,qy_vals,qz_vals,correction,twotheta_vals,chi_vals):                
                """Calcuation the lab/surface frame, it is constant so only need to calc once"""
                #get the current thread position
                row,col = cuda.grid(2)
                if row < pixel_count_y and col < pixel_count_x:
                    delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x 
                    delta_x_2 = delta_x*delta_x         
                    
                    #lab frame
                    gam = math.atan(delta_x*invDist)
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
                    
                    rad2deg = 180/3.1415926535
                    theta2 = rad2deg*math.acos(math.cos(gam)*math.cos(delta))
                    chi = 0
                    if math.tan(gam):
                        chi = rad2deg*(math.atan(math.tan(delta)/math.tan(gam)))                                       
                
                    qx = k0*(cos_beta*cos_psi - cos_alpha)                    
                    qy = k0*(cos_beta*sin_psi) 
                    qz = k0*(sin_beta+sin_alpha) 
        
                    if correct_refraction:
                        kzi = k0*math.sqrt(abs(sin_alpha**2 - math.sin(critical_angle)**2))
                        kzf = k0*math.sqrt(abs(sin_beta**2 - math.sin(critical_angle)**2))
                        qz = kzf - kzi
                        
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
                        #Clor = cos_alpha*sin_psi*cos_beta 
                        Clor = 1 # the integration is in carteasion coordinates for binning
                        C = Cpol*Crod*Cd*Ci*Clor
                    else:
                        C = 1                   
                    
                    qx_vals[row,col] = qx
                    qy_vals[row,col] = qy
                    qz_vals[row,col] = qz
                    correction[row,col] = abs(C)
                    twotheta_vals[row,col] = theta2
                    chi_vals[row,col] = chi                


            @cuda.jit(fastmath=True)
            def  _bin_kernel(image, midang, histogram_result,histogram_weights,bounds, qx_vals, qy_vals, qz_vals,corrections,twotheta_vals,chi_vals):               
                """This function rotates the surface frame based on the rotation of the sample around 
                its surface normal and then bins the intensity in a grid. It is mostly slowed down by
                transfering the detector image to the graphics memory"""
                row,col = cuda.grid(2)
                if row < image.shape[0] and col < image.shape[1]:

                    so = math.sin(rotDirection*(omega_rad+midang))
                    co = math.cos(rotDirection*(omega_rad+midang))
                    ci = float(1)
                    si = 0  

                    qx = qx_vals[row,col]
                    qy = qy_vals[row,col]
                    qz = qz_vals[row,col]
                    theta2 = twotheta_vals[row,col]
                    chi = chi_vals[row,col]
                    C = corrections[row,col]
                    delta_x = (col-q0[0])*pixel_x
                    sign = math.copysign(1,delta_x) #negative or positive qr      
                            
                    hphi_1 = so*(ci*qy+si*qz)+co*qx
                    hphi_2 = co*(ci*qy+si*qz)-so*qx
                    hphi_3 = ci*qz-si*qy  

                    if recp_units == 1:
                        h=  Binv[0,0]*hphi_1+Binv[0,1]*hphi_2+Binv[0,2]*hphi_3
                        k = Binv[1,0]*hphi_1+Binv[1,1]*hphi_2+Binv[1,2]*hphi_3
                        l = Binv[2,0]*hphi_1+Binv[2,1]*hphi_2+Binv[2,2]*hphi_3
                        r = b1/b2
                        qr = np.floor(sign*math.sqrt((hphi_1/b1)**2 + (r*hphi_2/b2)**2))
                    else:
                        h = hphi_1
                        k = hphi_2
                        l = hphi_3                    
                        qr =sign*math.sqrt(h**2 + k**2)

                    include_pixel_in_hk = False
                    include_pixel_in_h = False
                    include_pixel_in_l = False
                    include_pixel_in_pix = True
                    include_pixel_in_2theta = True
                
                    #Next we check if our pixel is within a mask
                    for bound_index in range(len(bounds)): 
                        bound_xmin = bounds[bound_index][0]
                        bound_xmax = bounds[bound_index][1]
                        bound_ymin = bounds[bound_index][2]
                        bound_ymax = bounds[bound_index][3]
                        #Index 4 corresponds to the img state of the mask
                        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 
                        # 3 = HK, 4=HL, 5=KL, 6=3d
                        
                        # 0 detector view pixels (exclusive)
                        if bounds[bound_index][4] == 0:                            
                            if col > bound_xmin and col < bound_xmax:                       
                                if row > bound_ymin and row < bound_ymax:
                                    #pixel masks are exlusive
                                    include_pixel_in_pix = False

                        #1 Mask is on Transformed detector view (inclusive)
                        elif bounds[bound_index][4] == 1:                            
                            if bound_xmin <= qr <= bound_xmax:                                              
                                if bound_ymin <= l <= bound_ymax:                                                                      
                                    include_pixel_in_hk = True
                                    include_pixel_in_l = True 
                        
                        #2 Mask is on 2 theta/chi (exclusive)
                        elif bounds[bound_index][4] == 2:                            
                            if bound_xmin <= theta2 <= bound_xmax:                                              
                                if bound_ymin <= chi <= bound_ymax:                                                                      
                                    include_pixel_in_2theta = False

                        #3 Mask is on HK view (inclusive)
                        elif bounds[bound_index][4] == 3:                            
                            if h > bound_xmin and h < bound_xmax:                       
                                if k > bound_ymin and k < bound_ymax:
                                    include_pixel_in_h = True
                    
                        #we currently don't use masks on view 4,5 or 6
                
                        if include_pixel_in_pix and include_pixel_in_2theta :
                            intensity = image[row,col]*C
                            #Now we we have the coordinates of our pixel we need to decide which bin to put it
                            histx = int((h-xmin)*invgridstepx)
                            histy = int((k-ymin)*invgridstepy)
                            histz = int((l-zmin)*invgridstepz)                                           
                            
                            #we bin differently depending on the projection direction
                            #HK
                            if imgState == 3:    
                                if include_pixel_in_hk and include_pixel_in_l:       
                                    cuda.atomic.add(histogram_result,(histx,histy),intensity)
                                    cuda.atomic.add(histogram_weights,(histx,histy),1)
                            #HL    
                            elif imgState == 4:  
                                #use both the in-plane mask for h and k and the transformed detector for L
                                if include_pixel_in_h and include_pixel_in_l:       
                                    cuda.atomic.add(histogram_result,(histx,histz),intensity)
                                    cuda.atomic.add(histogram_weights,(histx,histz),1)     
                            #KL
                            elif imgState  == 5:    
                                if include_pixel_in_h and include_pixel_in_l:  
                                    cuda.atomic.add(histogram_result,(histy,histz),intensity)
                                    cuda.atomic.add(histogram_weights,(histy,histz),1)     
                            #3d
                            elif imgState == 6:  
                                if include_pixel_in_h and include_pixel_in_l:            
                                    cuda.atomic.add(histogram_result,(histx,histy,histz),intensity)
                                    cuda.atomic.add(histogram_weights,(histx,histy,histz),1)                           
            # Configure the blocks
            threadsperblock = (16, 8) #could be moved to setting or config file
            blockspergrid_x = int(math.ceil(pixel_count_y / threadsperblock[0]))
            blockspergrid_y = int(math.ceil(pixel_count_x / threadsperblock[1]))
            blockspergrid = (blockspergrid_x, blockspergrid_y)   
            qx_vals  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single)) 
            qy_vals = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single))           
            qz_vals  = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single)) 
            corrections = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single))
            twotheta_vals = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single))
            chi_vals = cuda.to_device(np.zeros((pixel_count_y,pixel_count_x),dtype=np.single))
            #for later
            if imgState == 6:
                histogram_result = cuda.to_device(np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double))                       
                histogram_weights = cuda.to_device(np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double))
            else:
                histogram_result = cuda.to_device(np.zeros((grid_nx,grid_ny),dtype=np.double))                       
                histogram_weights = cuda.to_device(np.zeros((grid_nx,grid_ny),dtype=np.double))
            _lab_frame[blockspergrid, threadsperblock](qx_vals,qy_vals,qz_vals,corrections,twotheta_vals,chi_vals)         

      
        
            for i, angle in enumerate(angles): 
                window.statusLabel.setText("Calculating angle "+str(i)+" = "+str(angle))  
                if use_raw_files:           
                    tmp = np.ascontiguousarray(np.rot90(window.image_stack.get_image_unbinned(i+from_image)))
                else:
                    tmp = np.ascontiguousarray(np.rot90(imagedata[i]))
                imagei = cuda.to_device(tmp)            
                if (i%5==0):
                    window.statusLabel.setText("Calculating angle "+str(i)+" = "+str(angle))  

                _bin_kernel[blockspergrid, threadsperblock](imagei, np.deg2rad(angle), histogram_result,
                                                             histogram_weights,bounds, qx_vals, qy_vals, qz_vals,
                                                             corrections,chi_vals,twotheta_vals) 
                if window.single_thread == False:
                    progress_callback.emit(i)  
                else:
                    window.progressBar.setValue(i)  

                del(imagei)                    

            final_hist = np.asarray(histogram_result.copy_to_host()) 
            final_weights = np.asarray(histogram_weights.copy_to_host())     
            hist2 = np.divide(final_hist, final_weights, where=final_weights!=0)     

            xaxis = np.arange(xmin,xmax,gridstepx)
            yaxis = np.arange(ymin,ymax,gridstepy)
            zaxis = np.arange(zmin,zmax,gridstepz) 

            xgrid = np.tile(xaxis,(grid_nx,1))

            if imgState  == 5: 
                ygrid = np.flipud(np.transpose(np.tile(yaxis,(grid_ny,1))))
            else:
                ygrid = np.tile(yaxis,(grid_ny,1))

            zgrid = np.tile(zaxis,(grid_nz,1))
            #zgrid = np.flipud(np.transpose(np.tile(zaxis,(grid_nz,1))))   
            return xgrid,ygrid,zgrid,hist2  
                     
    else:     
        def _lab_frame():
            qx_img = np.zeros((pixel_count_y,pixel_count_x)) 
            qy_img = np.zeros((pixel_count_y,pixel_count_x))
            qz_img = np.zeros((pixel_count_y,pixel_count_x)) 
            theta2_img = np.zeros((pixel_count_y,pixel_count_x)) 
            chi_img = np.zeros((pixel_count_y,pixel_count_x)) 
            c_img = np.zeros((pixel_count_y,pixel_count_x)) #array of pixel intensity corrections            
            """Function calculate the pixel positon in lab frame"""
            for row in range(pixel_count_y):    
                for col in range(pixel_count_x):
                    delta_z= (q0[1]-row)*pixel_y  #real-space dinstance from centre pixel y
                    delta_x = (col-q0[0])*pixel_x  #real-space dinstance from centre pixel x 
                    delta_x_2 = delta_x*delta_x         
                    
                    #lab frame
                    gam = math.atan(delta_x*invDist)
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
                    
                    rad2deg = 180/3.1415926535
                    theta2 = rad2deg*math.acos(math.cos(gam)*math.cos(delta))
                    chi = 0
                    if math.tan(gam):
                        chi = rad2deg*(math.atan(math.tan(delta)/math.tan(gam)))                                       
                
                    qx = k0*(cos_beta*cos_psi - cos_alpha)                    
                    qy = k0*(cos_beta*sin_psi) 
                    qz = k0*(sin_beta+sin_alpha) 
        
                    if correct_refraction:
                        kzi = k0*math.sqrt(abs(sin_alpha**2 - math.sin(critical_angle)**2))
                        kzf = k0*math.sqrt(abs(sin_beta**2 - math.sin(critical_angle)**2))
                        qz = kzf - kzi
                        
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
                        #Clor = cos_alpha*sin_psi*cos_beta 
                        Clor = 1 # the integration is in carteasion coordinates for binning
                        C = Cpol*Crod*Cd*Ci*Clor
                    else:
                        C = 1                   
                    
                    qx_img[row][col] = qx
                    qy_img[row][col] = qy
                    qz_img[row][col] = qz
                    c_img[row][col] = abs(C)
                    theta2_img[row][col] = theta2
                    chi_img[row][col] = chi
                    
            return qx_img,qy_img,qz_img,c_img,theta2_img,chi_img
                
        def  _binner(image, angle, qx_img, qy_img, qz_img,corrections_img,theta2_img,chi_img):
            if imgState == 6:
                hist = np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double)                       
                weights = np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double)
            else:
                hist = np.zeros((grid_nx,grid_ny),dtype=np.double)                       
                weights = np.zeros((grid_nx,grid_ny),dtype=np.double)

            so = math.sin(rotDirection*(omega_rad+angle))
            co = math.cos(rotDirection*(omega_rad+angle))
            ci = float(1)
            si = 0  

            for row in range(pixel_count_y):    
                for col in range(pixel_count_x): 
                    qx = qx_img[row][col]
                    qy = qy_img[row][col]
                    qz = qz_img[row][col]
                    theta2 = theta2_img[row][col]
                    chi = chi_img[row][col]
                    C = corrections_img[row][col]
                    delta_x = (col-q0[0])*pixel_x
                    sign = math.copysign(1,delta_x) #negative or positive qr      
                              
                    hphi_1 = so*(ci*qy+si*qz)+co*qx
                    hphi_2 = co*(ci*qy+si*qz)-so*qx
                    hphi_3 = ci*qz-si*qy  

                    if recp_units == 1:
                        h=  Binv[0,0]*hphi_1+Binv[0,1]*hphi_2+Binv[0,2]*hphi_3
                        k = Binv[1,0]*hphi_1+Binv[1,1]*hphi_2+Binv[1,2]*hphi_3
                        l = Binv[2,0]*hphi_1+Binv[2,1]*hphi_2+Binv[2,2]*hphi_3
                        r = b1/b2
                        qr = np.floor(sign*math.sqrt((hphi_1/b1)**2 + (r*hphi_2/b2)**2))
                    else:
                        h = hphi_1
                        k = hphi_2
                        l = hphi_3                    
                        qr =sign*math.sqrt(h**2 + k**2)

                    include_pixel_in_hk = False
                    include_pixel_in_h = False
                    include_pixel_in_l = False
                    include_pixel_in_pix = True
                    include_pixel_in_2theta = True
                
                    #Next we check if our pixel is within a mask
                    for bound_index in range(len(bounds)): 
                        bound_xmin = bounds[bound_index][0]
                        bound_xmax = bounds[bound_index][1]
                        bound_ymin = bounds[bound_index][2]
                        bound_ymax = bounds[bound_index][3]
                        #Index 4 corresponds to the img state of the mask
                        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 
                        # 3 = HK, 4=HL, 5=KL, 6=3d
                        
                        # 0 detector view pixels (exclusive)
                        if bounds[bound_index][4] == 0:                            
                            if col > bound_xmin and col < bound_xmax:                       
                                if row > bound_ymin and row < bound_ymax:
                                    #pixel masks are exlusive
                                    include_pixel_in_pix = False

                        #1 Mask is on Transformed detector view (inclusive)
                        elif bounds[bound_index][4] == 1:                            
                            if bound_xmin <= qr <= bound_xmax:                                              
                                if bound_ymin <= l <= bound_ymax:                                                                      
                                    include_pixel_in_hk = True
                                    include_pixel_in_l = True 
                        
                        #2 Mask is on 2 theta/chi (exclusive)
                        elif bounds[bound_index][4] == 2:                            
                            if bound_xmin <= theta2 <= bound_xmax:                                              
                                if bound_ymin <= chi <= bound_ymax:                                                                      
                                    include_pixel_in_2theta = False

                        #3 Mask is on HK view (inclusive)
                        elif bounds[bound_index][4] == 3:                            
                            if h > bound_xmin and h < bound_xmax:                       
                                if k > bound_ymin and k < bound_ymax:
                                    include_pixel_in_h = True
                       
                        #we currently don't use masks on view 4,5 or 6

                   
                    if include_pixel_in_pix and include_pixel_in_2theta :
                        intensity = image[row][col]*C
                        #Now we we have the coordinates of our pixel we need to decide which bin to put it
                        histx = int((h-xmin)*invgridstepx)
                        histy = int((k-ymin)*invgridstepy)
                        histz = int((l-zmin)*invgridstepz)                                           
                        
                        #we bin differently depending on the projection direction
                        #HK
                        if imgState == 3:    
                            if include_pixel_in_hk and include_pixel_in_l:                               
                                hist[histx][histy] += intensity
                                weights[histx][histy] += 1
                        #HL    
                        elif imgState == 4:  
                            #use both the in-plane mask for h and k and the transformed detector for L
                            if include_pixel_in_h and include_pixel_in_l:                                              
                                hist[histx][histz] += intensity
                                weights[histx][histz] += 1
                        #KL
                        elif imgState  == 5:    
                            if include_pixel_in_h and include_pixel_in_l:                                            
                                hist[histy][histz] += intensity
                                weights[histy][histz] += 1
                        #3d
                        elif imgState == 6:  
                            if include_pixel_in_h and include_pixel_in_l:
                                hist[histx][histy][histz] += intensity
                                weights[histx][histy][histz] += +1

            return hist,weights  

        #use numba on cpu instead
        if acceleration == 1:
            if "jit" not in sys.modules:                
                from numba import jit                        
            binner = jit(_binner,nopython=True)
            lab_frame = jit(_lab_frame,nopython=True)
        else:
            print("This is probably really slow, you should try Numba or even better the GPU mode")
            binner =_binner
            lab_frame = _lab_frame

        #We need an extra dimension if 3D binning
        if imgState  == 6:
            final_hist = np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double)
            final_weights = np.zeros((grid_nx,grid_ny,grid_nz),dtype=np.double) 
        else:
            final_hist = np.zeros((grid_nx,grid_ny),dtype=np.double)
            final_weights = np.zeros((grid_nx,grid_ny),dtype=np.double) 

        qx_img,qy_img,qz_img,c_img,theta2_img,chi_img = lab_frame()
        #So to get the transformed detector view we can just bin on one image, the summed or averaged image
        if window.ImgState == 1:
            window.statusLabel.setText("Calculating..")
            image = window.showData
            image_bin,image_weights = binner(image,np.deg2rad(0),qx_img,qy_img,qz_img,c_img,theta2_img,chi_img)           
            hist = np.divide(image_bin, image_weights, where=image_bin!=0)  
            xaxis = np.arange(xmin,xmax,gridstepx)
            yaxis = np.arange(zmin,zmax,gridstepz)
            zaxis = np.arange(zmin,zmax,gridstepz) 
            xgrid = np.tile(xaxis,(grid_ny,1))
            ygrid = np.flipud(np.transpose(np.tile(yaxis,(grid_nx,1))))
            zgrid = np.flipud(np.transpose(np.tile(zaxis,(grid_nx,1))))
            return xgrid,ygrid,zgrid,hist

        
        
        for i, angle in enumerate(angles): 
            window.statusLabel.setText("Calculating angle "+str(i)+" = "+str(angle))  
            if use_raw_files:           
                image = np.rot90(window.image_stack.get_image_unbinned(i+from_image))  
            else:
                image = np.rot90(imagedata[i])                    
            image_bin,image_weights = binner(image,np.deg2rad(angle),qx_img,qy_img,qz_img,c_img,theta2_img,chi_img) 
      
            final_hist= final_hist + image_bin

            final_weights = final_weights + image_weights
            if window.single_thread == False:
                progress_callback.emit(i)  
            else:
                window.progressBar.setValue(i)                        

        hist2 = np.divide(final_hist, final_weights, where=final_weights!=0)  

        xaxis = np.arange(xmin,xmax,gridstepx)
        yaxis = np.arange(ymin,ymax,gridstepy)
        zaxis = np.arange(zmin,zmax,gridstepz) 


        xgrid = np.tile(xaxis,(grid_nx,1))

        if imgState  == 5: 
            ygrid = np.flipud(np.transpose(np.tile(yaxis,(grid_ny,1))))
        else:
            ygrid = np.tile(yaxis,(grid_ny,1))
        zgrid = np.tile(zaxis,(grid_nz,1))
        #zgrid = np.flipud(np.transpose(np.tile(zaxis,(grid_nz,1))))

          
        return xgrid,ygrid,zgrid,hist2    