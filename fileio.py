#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from PyQt5.QtWidgets import QMessageBox
import pyqtgraph.parametertree.parameterTypes as pTypes
import pyqtgraph as pg
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog, QCheckBox, QProgressBar
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QLabel, QGridLayout, QWidget,QPushButton,QInputDialog,QListWidget,QListWidgetItem
from PyQt5.QtCore import QSize, Qt    
import fabio
import h5py
import hdf5plugin
import pickle
import util
from scipy import interpolate


class ImageStack():
    """ This Class handles the image stack and background/dark subtraciton. It can also sitch together two images,
    horizontally, so they can be treated as one detector. 
    This is a rewrite of the orginal file input and output code to make it more object oriented."""

    #Flag dictionary
    flags = {'background dark image':False,
             'dark image':False,
             'background image':False}

    #File name dictionary for back and background images
    aux_file_names = {'background dark image 1': None,
                    'background dark image 2': None,
                    'dark image 1': None, 
                    'dark image 2': None,
                    'background image 1': None,
                    'background image 2': None}

    #image file names
    image_file_names_1 = None
    image_file_names_2 = None

    #image data
    image_data = None
    number_of_images = None
    log_file_name = None
    scanno = None
    scan_name = None
    hdf5_file_name = None

    #If there is an error, e.g. in selecting files
    file_error = False       

    #Normalisation  

    def __init__(self,parameters):
        self.images_read = False     
        self.angle_mode = False
        self.window = parameters 
        self.__update_params

    def save_state(self):
        #we will just select the file names and pickle
        state = [self.flags,
                self.aux_file_names,
                self.image_file_names_1,
                self.image_file_names_2,
                self.log_file_name,
                self.images_read,
                self.number_of_images,
                self.hdf5_file_name,
                self.scan_name,
                self.scanno]
        return state

    def restore_state(self,state):
        self.flags,\
        self.aux_file_names,\
        self.image_file_names_1,\
        self.image_file_names_2,\
        self.log_file_name,\
        self.images_read,\
        self.number_of_images,\
        self.hdf5_file_name,\
        self.scan_name,\
        self.scanno = state          

    def __update_params(self):
        #{"P21.2 1": 1, "P21.2 2": 2, "P07 2019": 3,"ID31 HDF5": 4},'
        self.beamline = self.window.p.param('Experiment', 'Beamline Preset').value() 
        self.rotation = self.window.p.param('Data Processing', 'Image Rotation').value() 
        self.image_flipud = self.window.p.param('Data Processing', 'Image Flip U/D').value() 
        self.image_fliplr = self.window.p.param('Data Processing', 'Image Flip L/R').value() 

        self.offset = self.window.p.param('Data Processing', 'Intensity Offset').value() 
        self.multiply = self.window.p.param('Data Processing', 'Multiply intensity by').value() 
        self.binning=int(self.window.p.param('Data Processing', 'Binning').value())

        #image sizes
        self.xSize=int(self.window.p.param('Experiment', 'X Pixels').value())
        self.ySize=int(self.window.p.param('Experiment', 'Y Pixels').value())

        self.gap = int(self.window.p.param('Experiment','Detector Gap (pixels)').value())        
        self.two_detectors = self.window.p.param('Data Processing', 'Use 2nd detector').value()
        self.mean =  self.window.p.param('Data Processing', 'Mean Images Instead of Max').value()
       
        #these are the binned sizes
        extra_x = util.add_to_divide(self.xSize,self.binning)
        extra_y = util.add_to_divide(self.ySize,self.binning)
        self.x_size=int((self.xSize+extra_x)/self.binning)
        self.y_size=int((self.ySize+extra_y)/self.binning)

        if self.two_detectors:
            self.x_size2 = int(2*self.x_size) + int(self.gap/self.binning)
        else:
            self.x_size2 = int(self.x_size)

    def __flip(self,image_data):
        """Chekc if the images should be flipped"""
        if self.image_flipud:
            image_data = np.fliplr(image_data)
        if self.image_fliplr:
            image_data = np.flipud(image_data)
        return image_data

    def __load_image(self,file_name, nobin=False):
        self.__update_params()
        if self.beamline != 4:           
            print("File loaded: ",file_name)
            #we load the image and bin,rotate,flip, multiply and offset as neeeded
            if nobin:
                return self.multiply*self.__flip(np.rot90(np.asarray(fabio.open(file_name).data), self.rotation))+self.offset
            else:
                return self.multiply*self.__flip(np.rot90(util.rebin(np.asarray(fabio.open(file_name).data),int(self.binning)), self.rotation))+self.offset
        else:
            #for HDF5 FILES
            #not quite sure how to do the same binning trick on 3d datasets
            #also don't really want to load the whole unbinned dataset into memory
            dset = (self.hdf5_file[self.scan_name]) 
            self.number_of_images = dset.len()          
            self.image_data=np.zeros([self.number_of_images,self.x_size2,self.y_size])
            for i in range(self.number_of_images):               
                self.image_data[i] = self.multiply*self.__flip(np.rot90(util.rebin(np.asarray(dset[i]),int(self.binning)), self.rotation)/self.monitor[i])+self.offset
            return self.image_data

    def get_image_unbinned(self,image_index):             
        self.__update_params()
        file_name = self.image_file_names_1[image_index]      
        if self.beamline != 4:          
       
            if self.two_detectors:
                unbinned_x_size = 2*self.xSize + self.gap
            else:
                unbinned_x_size = self.xSize 
            #we load the image and bin,rotate,flip, multiply and offset as neeeded       
            mon = (self.multiply*self.monitor[image_index]+self.offset)  
            img1 = self.__load_image(file_name, nobin=True) / mon
            image_data = np.zeros([unbinned_x_size,self.ySize])            
            
            image_data[0:self.xSize,:] += img1

            if self.two_detectors:                
                file_name_2 = self.image_file_names_2[image_index]
                img2 = img1 = self.__load_image(file_name_2, nobin=True) / mon      
                image_data[(self.xSize + self.gap):unbinned_x_size,:] += img2    
          
            return image_data - self.full_subtract_image*(self.multiply*self.monitor[image_index]+self.offset)           
        else:
            #for HDF5 FILES
            #not quite sure how to do the same binning trick on 3d datasets
            #also don't really want to load the whole unbinned dataset into memory
            dset = (self.hdf5_file[self.scan_name]) 
            print(image_index)
            return self.multiply*self.__flip(np.rot90(np.asarray(dset[image_index]), self.rotation)/self.monitor[image_index])+self.offset

    def __select_aux_image(self,description,dictionary_prefix,filename=None, filename2=None):
        """Selection dialogs for backround and dark images"""
        self.__update_params()
        if self.beamline == 4:  
            print("This needs implementing")
        else:
            self.__update_params()
            if not filename:
                options = QFileDialog.Options()
                self.aux_file_names[dictionary_prefix + ' 1'], _ = QFileDialog.getOpenFileName(self.window, description +' 1',"","TIFF File (*.tif);;All Files (*)", options=options) 
            else:
                self.aux_file_names[dictionary_prefix + ' 1'] = filename

            if self.aux_file_names[dictionary_prefix + ' 1']:
                self.flags[dictionary_prefix] = True

            #For 2nd detector
            if self.two_detectors:
                if not filename2:
                    self.aux_file_names[dictionary_prefix + ' 2'], _ = QFileDialog.getOpenFileName(self.window, description +' 2',"","TIFF File (*.tif);;All Files (*)", options=options)   
                    if self.aux_file_names[dictionary_prefix + ' 2'] == None:
                        QMessageBox.about(self, "You have two detectors selected but only selected one image", "Warning")
                        self.file_error = True  
                else:
                    self.aux_file_names[dictionary_prefix + ' 2'] = filename2 
  
    def load_aux_files(self):   
        """load background and dark images"""    
        self.__update_params()
        self.subtract_image=np.zeros((self.x_size2,self.y_size))  
        if self.two_detectors:
            unbinned_x_size = 2*self.xSize + self.gap
        else:
            unbinned_x_size = self.xSize   

        self.full_subtract_image = np.zeros((unbinned_x_size,self.ySize))   

        if self.two_detectors:
            after_gap = self.x_size + int(self.gap/self.binning)
            if self.flags['dark image']: 
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] += self.__load_image(self.aux_file_names['dark image 2'])  
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] += self.__load_image(self.aux_file_names['dark image 2'], nobin=True) 
            if self.flags['background image']:  
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] += self.__load_image(self.aux_file_names['background image 2']) 
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] += self.__load_image(self.aux_file_names['background image 2'], nobin=True) 
            if self.flags['background dark image']:  
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 2']) 
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] -= self.__load_image(self.aux_file_names['background dark image 2'],nobin=True) 
                
        if self.flags['dark image']: 
            self.subtract_image[0:self.x_size,0:self.y_size] += self.__load_image(self.aux_file_names['dark image 1'])  
            self.full_subtract_image[0:self.xSize,:] += self.__load_image(self.aux_file_names['dark image 1'],nobin=True)  

        if self.flags['background image']:       
            self.subtract_image[0:self.x_size,0:self.y_size] += self.__load_image(self.aux_file_names['background image 1'])
            self.full_subtract_image[0:self.xSize,:] += self.__load_image(self.aux_file_names['background image 1'],nobin=True) 

        if self.flags['background dark image']: 
            self.subtract_image[0:self.x_size,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 1'])
            self.full_subtract_image[0:self.xSize,:] -= self.__load_image(self.aux_file_names['background dark image 1'],nobin=True) 
      
    def select_background(self,paramHandle, filename=None, filename2=None):
        if type(paramHandle) == str:
            filename =  filename
            filename2 =  paramHandle
        """Call this to select a background"""
        self.__select_aux_image("Choose a background image for detector ", 'background image',filename, filename2)

    def select_background_dark(self,paramHandle,filename=None, filename2=None):
        if type(paramHandle) == str:
            filename =  filename
            filename2 =  paramHandle
        """Call this to select a dark image for the background image"""
        self.__select_aux_image("Choose background dark image for detector ", 'background dark image',filename, filename2)

    def select_dark(self, paramHandle, filename=None, filename2=None):
        if type(paramHandle) == str:
            filename =  filename
            filename2 =  paramHandle
        """Call this to select a dark image"""
        self.__select_aux_image("Choose dark image for detector ", 'dark image',filename, filename2)    
    
    def select_images(self, paramHandle,files=None, files2 = None):
        if type(paramHandle) == list:
            files2 =  files
            files =  paramHandle
        """Choose the image data"""
        self.__update_params()
        #This assumes with have TIF file, it will ignore anything containing .dark (i.e. for P07)
        if self.beamline != 4:
            if not files:
                options = QFileDialog.Options()
                files, _ = QFileDialog.getOpenFileNames(self.window,"Select images to use", "","TIFF Files (*.tif);;All Files (*)", options=options)
            if files:
                self.image_file_names_1 = []
                for file_name in files:
                    if ".dark" not in file_name:
                        self.image_file_names_1.append(file_name)
                self.number_of_images = len(self.image_file_names_1)
                self.images_read = True 
                
            if self.two_detectors:
                if not files2:
                    options2 = QFileDialog.Options()
                    files2, _ = QFileDialog.getOpenFileNames(self.window,"Select right detector images to use", "","TIFF Files (*.tif);;All Files (*)", options=options2)  
                if files2:
                    self.image_file_names_2 = []
                    for file_name in files2:
                        if ".dark" not in file_name:
                            self.image_file_names_2.append(file_name)
                else:
                    print("Two detectors selected by only one set of images selected")
                    self.images_read = False
                
        else:
            #in this case we use the HDF5 file format from the ESRF - Nov 2020
            #since dual detectors is not an option here it is not supported but 
            #should be simple to change in the future
            if not files:
                options = QFileDialog.Options()
                self.hdf5_file_name, _ = QFileDialog.getOpenFileName(self.window,"Select HDF5 file to use", "","HDF5 File (*.h5);;All Files (*)", options=options)
            else:
                self.hdf5_file_name = files
            if self.hdf5_file_name:
                with h5py.File(self.hdf5_file_name,'r') as hf:
                   items = list(hf.keys())
                   self.scanno, ok = QInputDialog.getItem(self.window, "Scan number",    "Select scan ending in .1", items, 0, False)
                   if ok:                    
                    self.scan_name = '/'+self.scanno+'/measurement/p3'
                    self.images_read = True                     

    def load_images(self): 
        """Loads files selected into memory"""  
        self.__update_params()   
        #esrf id31 beamline
        if self.beamline == 4:             
            self.hdf5_file = h5py.File(self.hdf5_file_name, 'r')   
            #we want the epoch informaiton for each image
            epoch = self.hdf5_file['/'+self.scanno+'/measurement/epoch_trig']
            #the scanno .2 branch has a shorter measurment frequecy, but first we need the scan no
            short_scan = self.scanno.split('.')[0]+'.2'
            #then we get the less frequent epoch            
            short_epoch = self.hdf5_file['/'+short_scan+'/measurement/epoch']
            #and the ring current which we use to normalise
            srcur = self.hdf5_file['/'+short_scan+'/measurement/srcur']
            #we then want to interplote the srcur reported
            epoch2srcur=interpolate.interp1d(short_epoch,srcur)
            self.monitor = epoch2srcur(epoch)
            self.image_data = self.__load_image(self.image_file_names_1)
            self.number_of_images = len(self.image_data)
            #We just make a titled numpy array to subtract the background form each image
            background = np.tile(self.subtract_image,(self.number_of_images,1,1))
            self.image_data -= np.int32(background)
            
            self.__process_log_file()
            self.window.statusLabel.setText(str(self.number_of_images) + ' images selected')
            self.hdf5_file.close()
        #p07
        if self.beamline == 3:
            #We do this to surpress warnings from TIFF file having a bad header
            import logging
            logger = logging.getLogger()
            logger.setLevel(logging.CRITICAL)
            self.log_file_name = "something"
            self.flags['normalize images'] = True

        if self.beamline != 4:
            self.window.statusLabel.setText(str(self.number_of_images) + ' images selected')
            self.window.progressBar.setMaximum(self.number_of_images)  
            #create empty array to place image data
            self.image_data=np.zeros([self.number_of_images,self.x_size2,self.y_size])
            #loop through the images and put in data, also subtract background if needed. 
            for i,filename in enumerate(self.image_file_names_1):
                self.window.progressBar.setValue(i)
                self.window.update()
                self.image_data[i,0:self.x_size,0:self.y_size]=self.__load_image(filename)-self.subtract_image[0:self.x_size,0:self.y_size]

            #We load them in two loops because this is often quicker than the disk 
            #jumping back and forward all the time, depending on how the files are stored. 
            if self.two_detectors:
                for i,filename in enumerate(self.image_file_names_2):
                    after_gap = self.x_size + int(self.gap/self.binning)
                    self.image_data[i,after_gap:self.x_size2,0:self.y_size]=self.__load_image(filename)-self.subtract_image[after_gap:self.x_size2,0:self.y_size]

            if self.log_file_name:
                self.__process_log_file()

            self.window.progressBar.setValue(i+1)
        self.window.update()

    def __get_p07_attributes(self, imgname):
        '''get_p07 attributes gets some useful data out of the metafile'''
        #time, monitor, angle,exposure
        with open(imgname.split('.')[0]+'.tif.metadata', 'r') as metadata:
            for line in metadata:
                ls = line.split('=')
                linesplit = [x.rstrip() for x in ls]
                if(linesplit[0]=='exposureTime'):
                    expt = float(linesplit[1]) #exposure time
                if(linesplit[0]=='userComment1'):
                    omega = float(linesplit[2].strip('"')) #angle
                if(linesplit[0]=='timeStamp'):
                    time = float(linesplit[1]) #time
                if(linesplit[0]=='extraInputs\\1\\extraInputs'):
                    monitor = abs(float(linesplit[1]))   #monitor              
        return expt,monitor,omega,time        

    def read_log_file(self,paramHandle,filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle    
        self.__update_params()
        #p21.2 - 2020
        if self.beamline == 2 or self.beamline == 1:  
            if not filename:      
                options = QFileDialog.Options()
                self.log_file_name, _ = QFileDialog.getOpenFileName(self.window,"Select log file to use", "","Log Files (*.fio);;All Files (*)", options=options)     
            else:
               self.log_file_name = filename 
        elif self.beamline == 4:
            print("Not available for this experiment metadata is extracted automatically")

    def __process_log_file(self):
        """Calls reads in the log file depending on the beamline"""
        if self.beamline == 2:   
            first_file_name_number=str(sorted(self.image_file_names_1)[0]).split('/')[-1]
            last_file_name_number=str(sorted(self.image_file_names_1)[-1]).split('/')[-1]  
            LogFileLines=None
            with open(self.log_file_name) as f:
                LogFileLines = f.readlines()
                self.LogFile=True
                for i, line in enumerate(LogFileLines):
                    if i < 42:
                        continue
                    if self.number_of_images != 0:
                        if line.split()[3]== first_file_name_number:
                            self.start_angle=float(line.split()[0])
                        if line.split()[3]==last_file_name_number:
                            self.end_angle=float(line.split()[0])
            self.monitor = np.ones(self.number_of_images)
            self.angle2image=interpolate.interp1d((self.start_angle,self.end_angle),(0,self.number_of_images))
            self.angle_mode = True

        #P07, DESY
        if self.beamline == 3: 
            all_angles = []
            self.monitor = []
            for i,imgname in enumerate(self.image_file_names_1):
                #exposure, monitor, angle, time
                expt,monitor,omega,time = self.__get_p07_attributes(imgname)
                all_angles.append(omega)
                self.monitor.append(monitor*expt)
                self.angle_mode = True               
                #not currently using the time
                self.image_data[i] = self.image_data[i]/(monitor*expt)   

            #are the angles backwards?       
            if all_angles[0] < all_angles[-1]:
                self.start_angle=all_angles[0]
                self.end_angle=all_angles[-1]
            else:
                self.start_angle=all_angles[-1]
                self.end_angle=all_angles[0]  
                self.image_data = np.flip(self.image_data, axis=0)  
                self.image_file_names_1 = np.flip(self.image_file_names_1)  
                self.image_file_names_2 = np.flip(self.image_file_names_2)  

        if self.beamline == 4:  
            all_angles = np.asarray(self.hdf5_file['/'+self.scanno+'/measurement/th/'])

            if all_angles[0] < all_angles[-1]:
                self.start_angle=all_angles[0]
                self.end_angle=all_angles[-1]
            else:
                self.start_angle=all_angles[-1]
                self.end_angle=all_angles[0]
                self.image_data = np.flip(self.image_data, axis=0) 

        self.angle2image=interpolate.interp1d((self.start_angle,self.end_angle),(0,self.number_of_images))
        self.angle_mode = True

    def max_example_3d(result, values):
        """
        Find the maximum value in values and store in result[0].
        Both result and values are 3d arrays.
        """
        i, j, k = cuda.grid(3)
        # Atomically store to result[0,1,2] from values[i, j, k]
        cuda.atomic.max(result, (0, 1, 2), values[i, j, k])
   
    def get_image(self,start, end):
        """return a summed or max image between with two angles or image numbers,
           depening on self.angle_mode and the 'Mean Images Instead of Max' parameter"""
        self.__update_params()
        if self.angle_mode:            
            self.from_image2=int(np.floor(self.angle2image(start)))
            self.to_image2=int(np.ceil(self.angle2image(end))) 
        else:
            self.from_image2=int(np.floor(start))
            self.to_image2=int(np.ceil(end))
        if self.mean:
            self.img = np.mean(self.image_data[self.from_image2:self.to_image2,:,:],axis=0)
        else:
            self.img =  np.max(self.image_data[self.from_image2:self.to_image2,:,:],axis=0)
        return self.img

    def is_region_ok(self,start,end):
        """Function to check if a given range is less is less than the image increment"""
        if self.angle_mode:      
            from_image=self.angle2image(start)
            to_image=self.angle2image(end)
        else:
            from_image = start
            to_image = end

        #we do 0.9 so it is actually possible to select one image           
        if to_image - from_image < 0.9:
                return False
        return True  

    def get_step(self):
        if self.angle_mode:            
            return (self.end_angle-self.start_angle)/self.number_of_images
        else:
            return 1
    def number_of_images_in_range(self,start_angle,end_angle):
        from_image=int(np.floor(self.angle2image(start_angle)))
        to_image=int(np.ceil(self.angle2image(end_angle)))  
        return to_image-from_image


            


