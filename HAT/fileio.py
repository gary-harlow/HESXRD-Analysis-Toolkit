#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from PyQt6.QtWidgets import QMessageBox
from PyQt6.QtWidgets import QInputDialog, QFileDialog
from PyQt6.QtCore import QObject, pyqtSignal,QRunnable,pyqtSlot
import fabio
import h5py
import hdf5plugin
from xrayhat import util
from scipy import interpolate
import natsort as ns
import traceback,sys


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)

class Worker(QRunnable):
    '''
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function
    '''

    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress

    @pyqtSlot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''

        # Retrieve args/kwargs here; and fire processing using them
        try:
            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)  # Return the result of the processing
        finally:
            self.signals.finished.emit()  # Done


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
        self.image_data = [None]        
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
        self.binning=int(self.window.p.param('Data Processing', 'Binning (X by X)').value())

        #image sizes
        self.xSize=int(self.window.p.param('Experiment', 'X Pixels').value())
        self.ySize=int(self.window.p.param('Experiment', 'Y Pixels').value())

        self.gap = int(self.window.p.param('Experiment','Detector Gap (pixels)').value())        
        self.two_detectors = self.window.p.param('Data Processing', 'Use 2nd detector').value()
        self.mean =  self.window.p.param('Data Processing', 'Mean Images Instead of Max').value()
        self.subtract_bak_image_box =  self.window.p.param('Data Processing', 'Perform Background Subtraction').value()
      
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
            text = "File loaded: "+str(file_name)
            self.window.statusLabel.setText(text) 
            #we load the image and bin,rotate,flip, multiply and offset as neeeded
            if nobin:
                return self.__flip(np.rot90(np.asarray(fabio.open(file_name).data), self.rotation))
            else:
                return self.__flip(np.rot90(util.rebin(np.asarray(fabio.open(file_name).data),int(self.binning)), self.rotation))
        else:
            #for HDF5 FILES
            #not quite sure how to do the same binning trick on 3d datasets
            #also don't really want to load the whole unbinned dataset into memory
            dset = (self.hdf5_file[self.scan_name]) 
            self.number_of_images = dset.len()   
            #if the image is rotated 90 or 270 degrees then the x and y sizes swap

            self.image_data=np.zeros([self.number_of_images,self.x_size2,self.y_size])
            for i in range(self.number_of_images):               
                self.image_data[i] = self.__flip(np.rot90(util.rebin(np.asarray(dset[i]),int(self.binning)), self.rotation)/self.monitor[i])
            return self.image_data

    def pubload_image(self,file_name, nobin=False):
        self.__update_params()
        if self.beamline != 4:                     
            text = "File loaded: "+str(file_name)
            self.window.statusLabel.setText(text) 
            #we load the image and bin,rotate,flip, multiply and offset as neeeded
            if nobin:
                return self.__flip(np.rot90(np.asarray(fabio.open(file_name).data), self.rotation))
            else:
                return self.__flip(np.rot90(util.rebin(np.asarray(fabio.open(file_name).data),int(self.binning)), self.rotation))
        else:
            #for HDF5 FILES
            #not quite sure how to do the same binning trick on 3d datasets
            #also don't really want to load the whole unbinned dataset into memory
            dset = (self.hdf5_file[self.scan_name]) 
            self.number_of_images = dset.len()   
            #if the image is rotated 90 or 270 degrees then the x and y sizes swap

            self.image_data=np.zeros([self.number_of_images,self.x_size2,self.y_size])
            for i in range(self.number_of_images):               
                self.image_data[i] = self.__flip(np.rot90(util.rebin(np.asarray(dset[i]),int(self.binning)), self.rotation)/self.monitor[i])
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
            mon = (self.monitor[image_index])  
            img1 = self.__load_image(file_name, nobin=True) / mon
            image_data = np.zeros([unbinned_x_size,self.ySize])            
            
            image_data[0:self.xSize,:] += img1

            if self.two_detectors:                
                file_name_2 = self.image_file_names_2[image_index]
                img2 = img1 = self.__load_image(file_name_2, nobin=True) / mon      
                image_data[(self.xSize + self.gap):unbinned_x_size,:] += img2    
          
            return image_data - self.full_subtract_image*(self.monitor[image_index])           
        else:
            #for HDF5 FILES
            #not quite sure how to do the same binning trick on 3d datasets
            #also don't really want to load the whole unbinned dataset into memory
            dset = (self.hdf5_file[self.scan_name]) 
            print(image_index)
            return self.__flip(np.rot90(np.asarray(dset[image_index]), self.rotation)/self.monitor[image_index])

    def __select_aux_image(self,description,dictionary_prefix,filename=None, filename2=None):
        """Selection dialogs for backround and dark images"""
        self.__update_params()
        if self.beamline == 4:  
            text ="Not available for this experiment metadata is extracted automatically"
            self.window.statusLabel.setText(text)             
        else:            
            self.__update_params()
            print(description,dictionary_prefix)
            if not filename:
                self.aux_file_names[dictionary_prefix + ' 1'], _ = QFileDialog.getOpenFileName(self.window, description +' 1',"","Image Files (*.tif *.tiff *.cbf *.edf);;All Files (*)") 
            else:
                self.aux_file_names[dictionary_prefix + ' 1'] = filename

            if self.aux_file_names[dictionary_prefix + ' 1']:
                self.flags[dictionary_prefix] = True

            #For 2nd detector
            if self.two_detectors:
                if not filename2:
                    self.aux_file_names[dictionary_prefix + ' 2'], _ = QFileDialog.getOpenFileName(self.window, description +' 2',"","Image Files (*.tif *.tiff *.cbf *.edf);;All Files (*)")   
                    if self.aux_file_names[dictionary_prefix + ' 2'] == None:
                        QMessageBox.about(self, "You have two detectors selected but only selected one image", "Warning")
                        self.file_error = True  
                else:
                    self.aux_file_names[dictionary_prefix + ' 2'] = filename2 
  
    def load_aux_files(self):   
        """load background and dark images"""    
        self.__update_params()
        self.subtract_image=np.zeros((self.x_size2,self.y_size))  
        self.subtract_bak_image=np.zeros((self.x_size2,self.y_size))  
        if self.two_detectors:
            unbinned_x_size = 2*self.xSize + self.gap
        else:
            unbinned_x_size = self.xSize   

        self.full_subtract_image = np.zeros((unbinned_x_size,self.ySize))   
        self.full_subtract_bak_image = np.zeros((unbinned_x_size,self.ySize)) 

        if self.two_detectors:
            after_gap = self.x_size + int(self.gap/self.binning)
            if self.flags['dark image']: 
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] += self.__load_image(self.aux_file_names['dark image 2'])  
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] += self.__load_image(self.aux_file_names['dark image 2'], nobin=True) 
            if self.flags['background image']:  
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] += self.__load_image(self.aux_file_names['background image 2']) 
                self.subtract_bak_image[after_gap:self.x_size2,0:self.y_size] += self.__load_image(self.aux_file_names['background image 2']) 
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] += self.__load_image(self.aux_file_names['background image 2'], nobin=True) 
            if self.flags['background dark image']:  
                self.subtract_image[after_gap:self.x_size2,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 2'])
                self.subtract_bak_image[after_gap:self.x_size2,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 2']) 
                self.full_subtract_image[(self.xSize + self.gap):unbinned_x_size,:] -= self.__load_image(self.aux_file_names['background dark image 2'],nobin=True) 
                
        if self.flags['dark image']: 
            self.subtract_image[0:self.x_size,0:self.y_size] += self.__load_image(self.aux_file_names['dark image 1'])  
            self.full_subtract_image[0:self.xSize,:] += self.__load_image(self.aux_file_names['dark image 1'],nobin=True)  

        if self.flags['background image']:       
            self.subtract_image[0:self.x_size,0:self.y_size] += self.__load_image(self.aux_file_names['background image 1'])
            self.subtract_bak_image[0:self.x_size,0:self.y_size] += self.__load_image(self.aux_file_names['background image 1'])
            self.full_subtract_image[0:self.xSize,:] += self.__load_image(self.aux_file_names['background image 1'],nobin=True) 

        if self.flags['background dark image']: 
            self.subtract_image[0:self.x_size,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 1'])
            self.subtract_bak_image[0:self.x_size,0:self.y_size] -= self.__load_image(self.aux_file_names['background dark image 1'])
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
            filename2 =  filename
            filename =  paramHandle
        """Call this to select a dark image"""
        self.__select_aux_image("Choose dark image for detector ", 'dark image',filename, filename2)    
    
    def select_images(self, paramHandle,files=None, files2 = None):
        if type(paramHandle) == list:
            files2 =  files
            files =  paramHandle
        print(files)
        print(paramHandle)
        """Choose the image data"""
        self.__update_params()
        #This assumes with have TIF file, it will ignore anything containing .dark (i.e. for P07)
        if self.beamline != 4:
            if not files:
                files, _ = QFileDialog.getOpenFileNames(self.window,"Select images to use", "","Image Files (*.tif *.tiff *.cbf *.edf);;All Files (*)")
            if files:
                self.image_file_names_1 = []
                for file_name in files:
                    if ".dark" not in file_name:
                        self.image_file_names_1.append(file_name)
                self.number_of_images = len(self.image_file_names_1)
                self.images_read = True 
                
            if self.two_detectors:
                if not files2:                    
                    files2, _ = QFileDialog.getOpenFileNames(self.window,"Select right detector images to use", "","Image Files (*.tif *.tiff *.cbf *.edf);;All Files (*)")  
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
                self.hdf5_file_name, _ = QFileDialog.getOpenFileName(self.window,"Select HDF5 file to use", "","HDF5 File (*.h5);;All Files (*)")
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
        #scripting is not super compatible with multithreading if we are running a script
        #we should not load the files in the background, otherwise the next script action will run
        #clear current image stack
        if None in self.image_data:
            del(self.image_data)
        try:
            self.image_data=np.zeros([self.number_of_images,self.x_size2,self.y_size],dtype=np.float32)
        except:
            self.window.statusLabel.setText('Insufficent memory!! Try to increase the binning.')                
            return 

        #we should appropriately sort the filelist
        self.image_file_names_1 = ns.natsorted(self.image_file_names_1,alg=ns.PATH)
        if self.image_file_names_2:
            self.image_file_names_2 = ns.natsorted(self.image_file_names_2,alg=ns.PATH)

        if self.window.single_thread == False:
            worker = Worker(self.__load_images) # Any other args, kwargs are passed to the run function
            #worker.signals.result.connect(self.print_output)
            worker.signals.finished.connect(self.thread_complete)
            worker.signals.progress.connect(self.progress_fn)
            self.window.busy = True
            self.window.threadpool.start(worker)   
        else:
            self.window.statusLabel.setText("SCRIPT MODE: Loading Files! (GUI DISABLED)")             
            self.__load_images(None)    

    def __load_images(self, progress_callback): 
        """Loads files selected into memory"""  
        self.__update_params()   
        self.angle_mode = False
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
        if self.beamline == 3 or self.beamline == 6:
            #We do this to surpress warnings from TIFF file having a bad header
            import logging
            logger = logging.getLogger()
            logger.setLevel(logging.CRITICAL)
            self.log_file_name = "something"
            self.flags['normalize images'] = False

        if self.beamline != 4:
            self.window.statusLabel.setText(str(self.number_of_images) + ' images selected')
            self.window.progressBar.setMaximum(self.number_of_images) 
            if self.two_detectors: 
                after_gap = self.x_size + int(self.gap/self.binning)
            #loop through the images and put in data, also subtract background if needed.           
                             
            for i,filename in enumerate(self.image_file_names_1):                               
                self.image_data[i,0:self.x_size,0:self.y_size]=self.__load_image(filename)-self.subtract_image[0:self.x_size,0:self.y_size]
                if self.two_detectors: 
                    self.image_data[i,after_gap:self.x_size2,0:self.y_size]=self.__load_image(self.image_file_names_2[i])-self.subtract_image[after_gap:self.x_size2,0:self.y_size]
                if self.window.single_thread == False:
                    progress_callback.emit(i)
                else:
                    self.window.progressBar.setValue(i)
            #We load them in two loops because this is often quicker than the disk 
            #jumping back and forward all the time, depending on how the files are stored. 
        
        if self.window.single_thread == True:
            text = str(self.number_of_images)+" Files Loaded"
            self.window.statusLabel.setText(text) 
            self.window.progressBar.setValue(self.number_of_images)
            if self.log_file_name or self.beamline == 5:
                self.__process_log_file()
                if self.beamline == 7:
                    #divide images by monitor 
                    self.image_data = self.image_data#/self.monitor
            
            self.window.loadComplete()          

    def progress_fn(self, n):        
        '''File loading progess'''
        self.window.progressBar.setValue(n)
        if n == 1:
             self.window.loadUpdate(n)
        if n!= 0 and n%100 == 0:
            self.window.loadUpdate(n)

    def thread_complete(self):
        '''Function to be run after files are loaded'''
        text = str(self.number_of_images)+" Files Loaded"
        self.window.statusLabel.setText(text) 
        self.window.progressBar.setValue(self.number_of_images)
        if self.log_file_name or self.beamline == 5:
            self.__process_log_file()
            if self.beamline == 7:
                #divide images by monitor 
                for i, mon in enumerate(self.monitor):
                    self.image_data[i,:,:] = self.image_data[i,:,:]/mon
            
        self.window.loadComplete()
        self.window.busy = False

    def __get_p07_attributes(self, imgname):
        '''Internal funciton same as external function - this is scheduled to be removed'''
        #time, monitor, angle,exposure
        with open(imgname.split('.')[0]+'.tif.metadata', 'r') as metadata:
            for line in metadata:
                ls = line.split('=')
                linesplit = [x.rstrip() for x in ls]
                if(linesplit[0]=='exposureTime'):
                    expt = float(linesplit[1]) #exposure time
                if(linesplit[0]=='userComment1'):
                    omega = float(linesplit[2].strip('=')[0:-2]) #angle
                if(linesplit[0]=='timeStamp'):
                    time = float(linesplit[1]) #time
                if(linesplit[0]=='extraInputs\\1\\extraInputs'):
                    monitor = 1#abs(float(linesplit[1]))   #monitor              
        return expt,monitor,omega,time 

    def get_p07_attributes(self, imgname):
        '''get attributes from p07 meta data files'''
        #time, monitor, angle,exposure
        with open(imgname.split('.')[0]+'.tif.metadata', 'r') as metadata:                        
            for line in metadata:                
                ls = line.split('=')
                linesplit = [x.rstrip() for x in ls]
                if(linesplit[0]=='exposureTime'):
                    expt = float(linesplit[1]) #exposure time
                if(linesplit[0]=='userComment1'):
                    omega = float(linesplit[2].strip('=')[0:-2]) #angle
                if(linesplit[0]=='timeStamp'):
                    time = float(linesplit[1]) #time
                if(linesplit[0]=='extraInputs\\1\\extraInputs'):
                    monitor = abs(float(linesplit[1]))   #monitor    
        return expt,monitor,omega,time 

    def __get_p07_attributes_2021(self, imgname):
        '''get_p07 attributes gets some useful data out of the metafile'''
        #time, monitor, angle,exposure
        with open(imgname.split('.')[0]+'.tif.metadata', 'r') as metadata:
            for line in metadata:
                ls = line.split('=')
                linesplit = [x.rstrip() for x in ls]
                if(linesplit[0]=='exposureTime'):
                    expt = float(linesplit[1]) #exposure time
                if(linesplit[0]=='userComment1'):
                    omega = float(linesplit[2].strip('=')[0:-2]) #angle
                if(linesplit[0]=='timeStamp'):
                    time = float(linesplit[1]) #time
                if(linesplit[0]=='extraInputs\\1\\extraInputs'):
                    monitor = abs(float(linesplit[1]))   #monitor              
        return expt,monitor,omega,time        

    def read_log_file(self,paramHandle,filename=None):
        """Function to read log file"""
        if type(paramHandle) == str:
            filename =  paramHandle    
        self.__update_params()
        #p21.2 - 2020
        if self.beamline == 2 or self.beamline == 1 or self.beamline== 7:  
            if not filename:                 
                self.log_file_name, _ = QFileDialog.getOpenFileName(self.window,"Select log file to use", "","Log Files (*.fio);;All Files (*)")     
            else:
               self.log_file_name = filename 
        elif self.beamline == 4:
            self.window.statusLabel.setText("Not available for this experiment metadata is extracted automatically") 

    def __process_log_file(self):
        """Internal function to proess several different log file formats"""
        self.__update_params()
        """Calls reads in the log file depending on the beamline"""
        #P21.2 Dec 21
        if self.beamline == 7:   
            flip = False
            LogFileLines=None
            with open(self.log_file_name) as f:
                #first we find the data section
                line = f.readline()
                while("%d" not in line):
                    line = f.readline()
                    if not line:
                        print("BAD INPUT FILE")
                        quit()
                        
                #read the col names
                angle_col = 0
                time_col = 0
                file_col = 0
                chan_col = 0

                line = f.readline()
                while("Col" in line):
                    col_name= line.split()[2]
                    #next we want the column number
                    if col_name == "idrz1(encoder)":
                        angle_col = int(line.split()[1])-1
                    if col_name == "filename":
                        file_col = int(line.split()[1])-1
                    if col_name == "unix":
                        time_col = int(line.split()[1])-1
                    if col_name == "channel":
                        chan_col = int(line.split()[1])-1
                    if "clearning" in line:
                        print("Warning a clearning image was included!")
                    line = f.readline()
                data = []        
                #read in columns we care about

                while len(line.split())>2:
                    row = line.split()
                    data.append([row[angle_col],row[time_col],row[file_col],row[chan_col]])
                    line = f.readline()                
     
            #we want the angular range
            first_file_name=str(self.image_file_names_1[0]).split('/')[-1]
            last_file_name=str(self.image_file_names_1[-1]).split('/')[-1]     


            for row in data:
                if first_file_name in row[2]:
                    self.start_angle = float(row[0])
                if last_file_name in row[2]:
                    self.end_angle = float(row[0])

            
            #swap the order if they are backwards
            if self.start_angle > self.end_angle:
                self.start_angle,self.end_angle = self.end_angle, self.start_angle
                flip = True
            
            #if the angle is not in the range discard
            data2 = []
            for i,row in enumerate(data):
                if float(row[0]) < self.start_angle or float(row[0])>self.end_angle:
                    continue
                if int(row[3]) == 1:
                    data2.append(row[1])

            #next we open the slower log file and interplote the monitor values as a function of time                   
            #this is a bit janky and will fail is fio is in the path
            timestamp = []
            eh_diode = []
            with open(self.log_file_name.replace("fio","log")) as f:
                line = f.readline()                
                #skip comments
                while (line[0]=="#"):
                    line = f.readline()   

                while(len(line.split()) >3):
                    row = line.split()
                    timestamp.append(float(row[0]))
                    eh_diode.append(float(row[3]))
                    line = f.readline()

            time_to_mon = interpolate.interp1d(np.asarray(timestamp),np.asarray(eh_diode))
            mon = []       
            for row in data2:
                mon.append(time_to_mon(row))            
            self.monitor = np.asarray(mon)
            if flip:
                #also need to revese monitor list and file list should be done automatically also
                self.monitor = np.flip(self.monitor)              
            self.angle2image=interpolate.interp1d((self.start_angle,self.end_angle),(0,self.number_of_images))
            self.LogFile=True
            self.angle_mode = True

        if self.beamline == 2:   
            first_file_name_number=str(self.image_file_names_1)[0].split('/')[-1]
            last_file_name_number=str(self.image_file_names_1)[-1].split('/')[-1]  
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
        if self.beamline == 3  : 
            all_angles = []
            self.monitor = []
            for i,imgname in enumerate(self.image_file_names_1):
                #exposure, monitor, angle, time
                expt,monitor,omega,time = self.__get_p07_attributes(imgname)             
                all_angles.append(omega)
                self.monitor.append(monitor*expt)
                self.angle_mode = True               
                #not currently using the time
                self.image_data[i] = self.image_data[i] #/(monitor*expt)  
            self.start_angle=all_angles[0]
            self.end_angle=all_angles[-1] 

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

        if self.beamline == 6: 
            all_angles = []
            self.monitor = []
            for i,imgname in enumerate(self.image_file_names_1):
                #exposure, monitor, angle, time
                expt,monitor,omega,time = self.__get_p07_attributes_2021(imgname)
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
        if self.beamline == 5:


            self.start_angle= float(self.window.p.param('Experiment', 'Manual Start Angle').value())
            self.end_angle= float(self.window.p.param('Experiment', 'Manual End Angle').value())
            if self.start_angle > self.end_angle:
                self.start_angle,self.end_angle = self.end_angle, self.start_angle
                flip = True
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
            #if there is background subtraction and we don't want it we should add the subtracted background back
            if self.flags['background image'] == True: 
                if self.subtract_bak_image_box == False:
                    self.img = self.img + self.subtract_bak_image
        else:
            self.img =  np.max(self.image_data[self.from_image2:self.to_image2,:,:],axis=0)
            if self.flags['background image'] == True:
                if self.subtract_bak_image_box == False:
                    self.img = self.img + self.subtract_bak_image*(end-start)
        return self.multiply*self.img+self.offset
        
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
        """Return the step size of one image """
        if self.angle_mode:            
            return (self.end_angle-self.start_angle)/self.number_of_images
        else:
            return 1

    def number_of_images_in_range(self,start_angle,end_angle):
        """Return how many images are in a given angular range"""
        from_image=int(np.floor(self.angle2image(start_angle)))
        to_image=int(np.ceil(self.angle2image(end_angle)))  
        return to_image-from_image


            


