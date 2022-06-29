#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy import interpolate
import pickle
import importlib

import locale
locale.setlocale(locale.LC_ALL, '')
from locale import atof
import pkg_resources

from PyQt6.QtCore import QSize, QThreadPool
from PyQt6.QtGui import QAction, QIcon,QDoubleValidator,QTransform
from PyQt6.QtWidgets import (
    QLabel,
    QMainWindow,
    QStatusBar,
    QToolBar,
    QGridLayout,
    QWidget,
    QDialog,
    QProgressBar,
    QLineEdit,
    QFileDialog    
)
from pyqtgraph.parametertree import Parameter, ParameterTree
import pyqtgraph as pg

#import various modules
from xrayhat import  plotting
from xrayhat import script
from xrayhat.crystal import *
from xrayhat import fileio
import xrayhat.binner as proj
from xrayhat import chi2theta
from xrayhat.transformdet import *
from xrayhat import scene3d


class MainWindow(QMainWindow):
    
    def __init__(self):
        '''This method initializes the window and sets up the widgets etc'''
        super().__init__()
        self.setWindowTitle("HESXRD Analysis Toolkit (HAT)")
        self.setStatusBar(QStatusBar(self)) 
      
        self.__init_globals__()   
        self.__init_menu__()
        self.__init_toolbar__()
        self.__init_parameter_tree__()
        
        layout = QGridLayout()        
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)       
        
        #plotting region
        self.win = pg.GraphicsLayoutWidget()        
        layout.addWidget(self.win,0,0)
        self.win.nextRow()    
                       
        self.p2 = self.win.addPlot(row=4,colspan=3)
      
        self.angleRegion=pg.LinearRegionItem()
        self.p2.setMaximumHeight(80)
        self.p2.addItem(self.angleRegion)
        self.p2.showAxis('left',show=False)
        
        self.qRegion=pg.LinearRegionItem(orientation='horizontal')
        self.p3 = self.win.addPlot(row=1,col=0)
        self.imgLeft = pg.ImageItem()
        self.p3.addItem(self.imgLeft)
        self.p3.addItem(self.qRegion)
        
        #parameter tree
        t = ParameterTree()
        t.setMaximumWidth(360)
        t.setParameters(self.p, showTop=False)
        layout.addWidget(t,0,1)
        
        self.progressBar=QProgressBar()
        layout.addWidget(self.progressBar,2,0)      
        
        self.statusLabel=QLabel('Select Images ... (Press Load Saved Files to load the files you had when you last pressed \'Save Parameters\')')
        layout.addWidget(self.statusLabel,2,0)   
        
        self.inPlaneImage=pg.ImageItem()
        self.qRegion.hide()
        self.p4 = self.win.addPlot(row=3,col=0)
        self.p4.hide()
        self.ctrROI=pg.LineROI([0, 5], [0, 0], width=2, pen=(1,9))
        
        self.p4.setMaximumHeight(300)
        self.p4.setLogMode(x=None,y=True)
        self.CTRCurve=self.p4.plot()
        self.ctrROI.hide()
        self.p3.addItem(self.ctrROI)      
        
  
        self.hist = pg.HistogramLUTItem()        
        self.hist.setImageItem(self.imgLeft)        
        self.hist.setMaximumWidth(200)
        self.win.addItem(self.hist,row=1,col=2)
        self.imgLeft.hoverEvent = self.imageHoverEvent   
        
        self.win.addItem(self.hist,row=1,col=2)
        self.imgLeft.hoverEvent = self.imageHoverEvent   
        
        #text boxes        
        toolbar1 = QToolBar()
        layout.addWidget(toolbar1, 1, 0)
        onlydub = QDoubleValidator()
        
        text3 = QLabel("  Start Angle ")
        toolbar1.addWidget(text3)        
        self.startangle = QLineEdit()
        self.startangle.setValidator(onlydub)
        startanglewidget = toolbar1.addWidget(self.startangle)
        
        text4 = QLabel("  End Angle ")
        toolbar1.addWidget(text4)        
        self.endangle = QLineEdit()
        self.endangle.setValidator(onlydub)
        endanglewidget = toolbar1.addWidget(self.endangle)    
        
        text1 = QLabel("  Cmin ")
        toolbar1.addWidget(text1)        
        self.cminbox = QLineEdit()
        self.cminbox.setValidator(onlydub)
        cminwidget = toolbar1.addWidget(self.cminbox)     
        
        text2 = QLabel("  Cmax ")
        toolbar1.addWidget(text2)        
        self.cmaxbox = QLineEdit()
        self.cmaxbox.setValidator(onlydub)
        cmaxwidget = toolbar1.addWidget(self.cmaxbox)   
        
               
        #signals
        self.hist.sigLevelChangeFinished.connect(self.updateCrange)  
        self.angleRegion.sigRegionChanged.connect(self.updateRegion)
        self.cmaxbox.editingFinished.connect(self.updateCrangeManual)
        self.cminbox.editingFinished.connect(self.updateCrangeManual)   
        self.endangle.editingFinished.connect(self.updateAngleManual)    
        self.startangle.editingFinished.connect(self.updateAngleManual)  
        self.ctrROI.sigRegionChanged.connect(self.ROICalc)   
        
        self.experiment.param('Manual Start Angle').sigValueChanged.connect(self.reloadWarning)
        self.experiment.param('Manual End Angle').sigValueChanged.connect(self.reloadWarning)
                
        self.p.param('Data Processing', 'Log Intensity').sigValueChanged.connect(self.logIntensity)
        self.p.param('Data Processing', 'Perform Background Subtraction').sigValueChanged.connect(self.updateRegion)
        self.p.param('Data Processing', 'Mean Images Instead of Max').sigValueChanged.connect(self.updateRegion)   
        self.p.param('Data Processing', 'Multiply intensity by').sigValueChanged.connect(self.updateRegion) 
        self.p.param('Data Processing', 'Intensity Offset').sigValueChanged.connect(self.updateRegion) 

        self.p.param('Profile tools', 'Log profile').sigValueChanged.connect(self.ROICalc)  
        self.p.param('Profile tools', 'Axis of interest').sigValueChanged.connect(self.ROICalc)  
        self.experiment.param('Beamline Preset').sigValueChanged.connect(self.data_format_changed)
        self.crystal.param('Preset').sigValueChanged.connect(self.sample_preset)       
        
                                                    
                              
    def __init_globals__(self):
        '''Initialize various class globals'''
        self.image_stack = fileio.ImageStack(self)     
        self.crystal=CrystallographyParameters(name='Crystal')
        self.experiment=ExperimentParameter(name='Experiment')           
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
        self.ImgState=0 
        self.showData = None #data in current view
        self.pixMaps=[] #coordinate arrays
        self.roixdata=None #data from roi/line profile
        self.roiydata=None #data from roi/line profile
        self.grid_h = None #coordinate arrays
        self.grid_k = None #coordinate arrays
        self.select_l_slice = True
        self.showROI = False #Toggle Line profile
        self.newImgState = 0 # review
        self.lineProfile = False #Line profile is active? Review wrt showROI
        self.single_thread = False    #Single Threaded mode
        self.cshift = 0 #for a potential color shift when taking log
        self.log_active = False #Log intensity?
        self.qmax =0  #max Q value possible, used for axis - REVIEW
        self.white_background = False #Background color mode
        self.qrmax = 0 #max qr from inplane axis size
        
        self.mask_list = []
        self.binBounds = []
        self.threadpool = QThreadPool()
        self.busy = False
        
        
    def __init_toolbar__(self):
        '''Internal function to initialize the toolbar'''
        #toolbar
        self.toolbar = QToolBar("Toolbar")
        self.toolbar.setIconSize(QSize(32, 32))
        self.addToolBar(self.toolbar)
        
        icons = pkg_resources.resource_filename('xrayhat', 'icons/')
        
        addMask = self.toolbar.addAction("Add Mask")
        addMask.setIcon(QIcon(icons+"addmask.svg"))
        addMask.setStatusTip("Add a new mask to the current view")
        addMask.triggered.connect(self.addMask)
        
        saveMask = self.toolbar.addAction("Save Masks")
        saveMask.setIcon(QIcon(icons+"savemask.svg"))
        saveMask.setStatusTip("Save masks across all views to a file")
        saveMask.triggered.connect(self.saveMasks)
                        
        loadMask = self.toolbar.addAction("Load Masks")
        loadMask.setIcon(QIcon(icons+"loadmask.svg"))
        loadMask.setStatusTip("Load masks (possibly for multiple views) from a file, current masks are not cleared")
        loadMask.triggered.connect(self.loadMasks)
        
        clearMask = self.toolbar.addAction("Clear Masks")
        clearMask.setIcon(QIcon(icons+"clearmask.svg"))
        clearMask.setStatusTip("Clear Masks across all views")
        clearMask.triggered.connect(self.clearMasks)
        
        #toolbar.addSeparator()
        self.toolbar2 = QToolBar("Toolbar")
        self.toolbar2.setIconSize(QSize(32, 32))
        self.addToolBar(self.toolbar2)

        
        tglROI = self.toolbar2.addAction("Toogle ROI")
        tglROI.setIcon(QIcon(icons+"roi.svg"))
        tglROI.triggered.connect(self.setProfile)
        tglROI.setStatusTip("Switch box profile on and off")

        loadROI = self.toolbar2.addAction("Load ROI")
        loadROI.setIcon(QIcon(icons+"loadroi.svg"))
        loadROI.triggered.connect(self.loadroi)
        loadROI.setStatusTip("Load the position and shape of a previous box profile")

        saveROI = self.toolbar2.addAction("Save ROI")
        saveROI.setIcon(QIcon(icons+"saveroi.svg"))
        saveROI.triggered.connect(self.saveroi)
        saveROI.setStatusTip("Save the current box profile position and shape to a file")

        extractROI =self.toolbar2.addAction("Extract ROI")
        extractROI.setIcon(QIcon(icons+"extractroi.svg"))
        extractROI.triggered.connect(self.saveprofile)
        extractROI.setStatusTip("Save the data from the profile to a file")

        extractRock = self.toolbar2.addAction("Extract Rocking Scans")
        extractRock.setIcon(QIcon(icons+"rockingscan.svg"))
        extractRock.triggered.connect(self.saverocks)  
        extractRock.setStatusTip("Save rocking scans (as a function of angle) for each pixel along the profile")
        
        #button_action = QAction(QIcon("bug.png"), "&Your button", self)
        #button_action.setStatusTip("This is your button")
        #button_action.triggered.connect(self.onMyToolBarButtonClick)
        #button_action.setCheckable(True)
        #toolbar.addAction(button_action)
        
    def __init_menu__(self):
        '''Internal Function to initialize'''        
        menubar = self.menuBar()
        actionFile = menubar.addMenu("&File")    
        loadAct = QAction('Load',self)
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

        importParmAct = QAction('Import Current Parameters', self)
        importParmAct.setStatusTip('Save current parameters to a file')
        importParmAct.triggered.connect(self.importParameters)
        actionFile.addAction(importParmAct)

        importFileAct = QAction('Import File Selection', self)
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

        thet2chiViewAct = QAction('2Theta / Chi View', self)
        thet2chiViewAct.setShortcut('Ctrl+3')
        thet2chiViewAct.setStatusTip('Summed or averaged image stack is corrected for curvature of the Ewald Sphere and given coordinates')
        thet2chiViewAct.triggered.connect(self.theta2chiTransform)
        actionView.addAction(thet2chiViewAct)
        
        binViewAct = QAction('Binned Projection', self)
        binViewAct.setShortcut('Ctrl+4')
        binViewAct.setStatusTip('Use all pixels across the selected angular range to create a projection')
        binViewAct.triggered.connect(self.makeProjection)
        actionView.addAction(binViewAct)

        full3dViewAct = QAction('3D Binned Projection', self)
        full3dViewAct.setShortcut('Ctrl+5')
        full3dViewAct.setStatusTip('Use all pixels across the selected angular range to create a projection')
        full3dViewAct.triggered.connect(self.make3dProjection)
        actionView.addAction(full3dViewAct)

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
        
    def __init_parameter_tree__(self):
        params = [
            self.experiment,
            self.crystal,
            {'name': 'Data Processing', 'type': 'group', 'children': [
                {'name': 'Binning (X by X)', 'type': 'int', 'value': 4, 'step': 1},
                {'name': 'View 2 Resolution', 'type': 'int', 'value': 1000, 'step': 1},
                {'name': 'Grid Size Qx', 'type': 'int', 'value': 800, 'step': 1},
                {'name': 'Grid Size Qy', 'type': 'int', 'value': 800, 'step': 1},
                {'name': 'Grid Size Qz', 'type': 'int', 'value': 800, 'step': 1},
                {'name': 'Select Projection', 'type': 'list', 'values': {"h-k": 0, "h-l": 1, "k-l": 2,"3D View": 3},'value': 0},               
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
                {'name': 'Log profile', 'type': 'bool', 'value': True},
                {'name': 'Axis of interest', 'type': 'list',  'values': {"H (Qx)": 1, "K (Qy)": 2, "HK (Qr)": 3,"L (qz)": 4}, 'value': 3} ,
            ]}
        ]
        self.p = Parameter.create(name='params', type='group', children=params)   
        
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
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d 
        if self.ImgState==7 or self.ImgState==0: 
            #showing background image, or in pixel coordinates
            self.p3.setTitle(" (%0.3f, %0.3f) / (%d, %d), Binned / Est. Real I: %0.3f / %0.3f" % (x, ySize-y, x*binning, binning*(ySize-y),self.imgLeft.image[j,i],self.imgLeft.image[j,i]/binning**2))
        elif self.ImgState==2:
            #in-plane map
            self.p3.setTitle("h: %0.3f, k: %0.3f, dist. from origin: %0.3f Intensity: %0.3f " % (x, y,np.sqrt(x*x+y*y),self.imgLeft.image[j,i]))  
        elif self.ImgState==1:
            #reciprocal space (out-of-plane)
            self.p3.setTitle("hk: %0.3f, l: %0.3f, dist. from origin: %0.3f Intensity: %0.3f" % (x, y,np.sqrt(x*x+y*y),self.imgLeft.image[j,i]))        
        
    
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



            
    def aboutHelp(self):
        """Function to display and about help dialog"""
        d = QDialog()
        d.setWindowTitle("About")
        d.setFixedSize(500,210)

        t1 = QLabel("HESXRD Analysis Toolkit.(c) 2021 Gary Harlow. \n\nFree to use under MIT License\n\nCreated by Gary Harlow and Sebastian Pfaff",d)        
        t1.setWordWrap(True)

        t2 = QLabel(d)             
        t2.setText('<a href="https://github.com/gary-harlow/HESXRD-Analysis-Toolkit">Link to GitHub Page</a>')
        t2.setOpenExternalLinks(True)

        t3 = QLabel(d)
        t3.setText('<a href="https://xray-hat.readthedocs.io/en/latest">Link to Documentation</a>')
        t3.setOpenExternalLinks(True)

        t1.move(20,2)
        t2.move(20,150)
        t3.move(20,170)
        
        d.exec()
        
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
            filename, _ = QFileDialog.getOpenFileName(self,"Select a parameter file to load", "","par (*.par);;All Files (*)")
        if filename:
            self.p.restoreState(pickle.load(open(filename, "rb" )),removeChildren=False)
            self.statusLabel.setText("Saved parameters restored") 

    def importImageStack(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:       
            filename, _ = QFileDialog.getOpenFileName(self,"Select an Image Stack to Load", "","bin (*.bin);;All Files (*)")
        if filename:        
            self.image_stack.restore_state(pickle.load(open(filename, "rb" )))   
            self.image_stack.load_aux_files()       
            text = str(self.image_stack.number_of_images)+" images selected with flags"+str(self.image_stack.flags)+"  Press Load (CTRL+L)"            
            self.statusLabel.setText(text) 

    def exportParameters(self, paramHandle, filename=None):
        if type(paramHandle) == str:
                filename =  paramHandle
        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","par (*.par);;All Files (*)")     
        if filename:
            state = self.p.saveState()
            pickle.dump(state, open(filename, "wb" ))
            self.statusLabel.setText("Parameters exported to "+filename)             

    def exportImageStack(self, paramHandle, filename=None):
        if type(paramHandle) == str:
                filename =  paramHandle
        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","bin (*.bin);;All Files (*)")     
        if filename:
            pickle.dump(self.image_stack.save_state(),open(filename, "wb" ) ) 
            self.statusLabel.setText("File selection exported to "+filename)             

    def restoreImageStack(self):        
        self.image_stack.restore_state(pickle.load(open( 'imagestack.bin', "rb" )))   
        self.image_stack.load_aux_files()           
        text = str(self.image_stack.number_of_images)+" images selected with flags"+str(self.image_stack.flags)+"  Press Load (CTRL+L)"            
        self.statusLabel.setText(text) 

    def exithat(self):
        print("Goodbye. Have a nice day!")
        quit()
        
    def toggleBg(self):            
        if self.white_background:
            self.win.setBackground('k')
            self.white_background = False
        else:
            self.win.setBackground('w')   
            self.white_background = True  
    
    def runScript(self):
        importlib.reload(script)
        self.single_thread = True
        text = "Running Script"
        self.update()
        script.script_main(self)
        text = "Script Finished"
        self.statusLabel.setText(text) 
        self.single_thread = False 
    
    def makeFigure(self, paramHandle, filename=None):
        projection = self.p.param('Data Processing', 'Select Projection').value()
        if type(paramHandle) == str:
            filename =  paramHandle
        """This function reloads plotting.py and runs the correct plotting function to
        generate a figure with matplotlib"""
        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","png (*.png);;All Files (*)")     

        if filename:
            importlib.reload(plotting)
            if self.ImgState > 2:
                cmin,cmax=self.imgLeft.getLevels()
                if projection == 0:
                    plotting.plot_projection_hk(self, self.xgrid,self.ygrid,self.showData,cmin,cmax, filename)
                    text = "In-plane qx/qy map saved as: "+str(filename)
                    self.statusLabel.setText(text)
                    
                if projection == 1:
                    plotting.plot_projection_hl(self, self.xgrid,self.ygrid,self.showData,cmin,cmax, filename)
                    text = "qx/qz map saved as: "+ str(filename)
                    self.statusLabel.setText(text)
                if projection == 2:
                    plotting.plot_projection_kl(self, self.xgrid,self.ygrid,self.showData,cmin,cmax, filename)
                    text = "qy/qz map saved as: "+str(filename)
                    self.statusLabel.setText(text)

            if self.ImgState == 2:
                cmin,cmax=self.imgLeft.getLevels()
                plotting.plot_2theta_chi(self.theta2grid,self.chigrid,self.showData,cmin,cmax, filename)
                text = "2Theta/Chi View saved as: "+str(filename)
                self.statusLabel.setText(text)

            if self.ImgState == 1:
                cmin,cmax=self.imgLeft.getLevels()
                plotting.plot_transformed_detector(self.grid_hk,self.grid_l,self.showData,cmin,cmax, filename)
                text = "Transfomred Detector View saved as: "+str(filename)
                self.statusLabel.setText(text)

            if self.ImgState == 0:
                cmin,cmax=self.imgLeft.getLevels()
                plotting.plot_det_view(self.showData,cmin,cmax, filename)
                text = "Image saved as: "+str(filename)
                self.statusLabel.setText(text)        
                
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
            #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
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

    def ctrROIClick(self):
        #if self.ImgState==2:
        #    return
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
        xSize=int(self.p.param('Experiment', 'X Pixels').value()/binning)

        if not self.showROI:
            self.ctrROI.show()
            self.p4.show()
            if self.ImgState == 0:
                self.ctrROI.setSize(100)
                self.ctrROI.setPos([xSize/2,xSize/2])
            else:
                self.ctrROI.setSize(1)
                self.ctrROI.setPos([0,0])
            self.showROI = True
        else:
            self.showROI = False
            self.ctrROI.hide()
            self.p4.hide()
        return

    def setProfile(self,value=None):
        """Function to toogle line profile"""
        if value == None:
            self.lineProfile^=True #toogle flag            
        else: 
            self.lineProfile = value       
        self.ctrROIClick()    
        self.ROICalc()          
      
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
        res=int(self.p.param('Data Processing', 'View 2 Resolution').value())
        
        if self.p.param('Profile tools', 'Log profile').value():
            self.p4.setLogMode(x=None,y=True)
        else:
             self.p4.setLogMode(x=None,y=False)

        if not self.showROI:
            return
        #hkl view
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL
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
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL
        elif self.ImgState > 2:
            projection = self.p.param('Data Processing', 'Select Projection').value()
            
            #axis = qx H
            if self.p.param('Profile tools', 'Axis of interest').value() == 1:
                #h-k
                if projection == 0:                
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)
                #h-l
                elif projection == 1:
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)
                #k-l    
                elif projection == 2:
                    self.statusLabel.setText('This projection does not calculate that axis, using qy/k instead')     
                    self.update()
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)
                #2theta-chi
                elif projection == 3:
                    self.statusLabel.setText('Feature currently not available')     
                    self.update()
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)

            #axis = qy 
            elif self.p.param('Profile tools', 'Axis of interest').value() == 2:
                #h-k
                if projection == 0:                                 
                    xRoiData=self.ctrROI.getArrayRegion(self.ygrid, self.imgLeft)
                #h-l
                elif projection == 1:
                    self.statusLabel.setText("This projection does not calculate that axis, using qx/h instead") 
                    xRoiData=self.ctrROI.getArrayRegion(self.ygrid, self.imgLeft)
                #k-l    
                elif projection == 2:
                    xRoiData=self.ctrROI.getArrayRegion(self.ygrid, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    self.statusLabel.setText('This projection does not calculate that axis, using qr/(sqrt(h^2+k^2)) instead')     
                    self.update()                    
                    xRoiData=self.ctrROI.getArrayRegion(self.ygrid, self.imgLeft)

            #axis = qr
            elif self.p.param('Profile tools', 'Axis of interest').value() == 3:
                #h-k
                if projection == 0:                          
                    xRoiData=self.ctrROI.getArrayRegion(self.grid_qr, self.imgLeft)
                #h-l
                elif projection == 1: 
                    self.statusLabel.setText('This axis could be misleading showing h/qx instead')     
                    self.update()    
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)
                #k-l    
                elif projection == 2: 
                    self.statusLabel.setText('This axis could be misleading showing k/qx instead')     
                    self.update()    
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    xRoiData=self.ctrROI.getArrayRegion(self.xgrid, self.imgLeft)          
            #axis = qz
            elif self.p.param('Profile tools', 'Axis of interest').value() == 4:
                #h-k
                if projection == 0:                       
                    self.statusLabel.setText('This projection does not calculate that axis, using qx/h instead')     
                    self.update()                         
                    xRoiData=self.ctrROI.getArrayRegion(self.zgrid, self.imgLeft)
                #h-l
                elif projection == 1:
                    xRoiData=self.ctrROI.getArrayRegion(self.zgrid, self.imgLeft)
                #k-l    
                elif projection == 2:
                    xRoiData=self.ctrROI.getArrayRegion(self.zgrid, self.imgLeft)
                #h^2 + k ^2
                elif projection == 3:
                    xRoiData=self.ctrROI.getArrayRegion(self.zgrid, self.imgLeft)   
    
            
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
            
    def addMask(self): 
        #when the button is pressed add a new ROI 
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL
        view_mode =  self.ImgState #it's type depends on the view it was created in         
        i = len(self.mask_list)     
        self.mask_list.append([pg.RectROI([0, 5], [1, 1], pen=(i,9)),view_mode])
        self.mask_list[i][0].show()
        self.p3.addItem(self.mask_list[i][0]) 

    def convertMasks(self):
        #this converts the ROIs into bounds for the binning algorithm
        xmaxs = []
        xmins = []
        zmaxs = []
        zmins = []
        ymaxs = []
        ymins = []
        for mask in self.mask_list:   
            view =  mask[1]    
            mask_xmin = mask[0].pos()[0]
            mask_xmax = mask[0].size()[0] + mask[0].pos()[0]
            mask_ymin = mask[0].pos()[1]
            mask_ymax = mask[0].size()[1] + mask[0].pos()[1]                        
            #Transformed detector    
            if view == 1:
                zmins.append(mask_ymin)
                zmaxs.append(mask_ymax)
            #HK
            if view == 3:
                xmins.append(mask_xmin)
                xmaxs.append(mask_xmax)
                ymins.append(mask_ymin)
                ymaxs.append(mask_ymax)
            #generate a list of bounds
            self.binBounds.append([mask_xmin,mask_xmax,mask_ymin,mask_ymax,view])   

        if xmaxs:
            self.xmin = np.min(xmins)
            self.xmax = np.max(xmaxs)
        else:
            self.xmin = 0
            self.xmax = 1
        
        if ymaxs:
            self.ymin = np.min(ymins)
            self.ymax = np.max(ymaxs)
        else:
            self.ymin = 0
            self.ymax = 1

        if zmaxs:
            self.zmin = np.min(zmins)
            self.zmax = np.max(zmaxs)
        else:
            self.zmin = 0
            self.zmax = 1      
        
                    
    def clearMasks(self):
        self.binBounds = []   
        for mask in self.mask_list:
            self.p3.removeItem(mask[0])
        self.mask_list = []
    
    def hideMasks(self): 
        for mask in self.mask_list:
            self.p3.removeItem(mask[0])
    
    def showMasks(self): 
        for mask in self.mask_list:
            #We should only add masks relevant to the current view
            if mask[1]== self.ImgState:
                self.p3.addItem(mask[0])

    def saveMasks(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","msk (*.msk);;All Files (*)")     
        if filename:
            mask_states = []
            for mask in self.mask_list:
                mask_states.append([mask[0].saveState(),mask[1]])
            pickle.dump(mask_states, open( filename, "wb" ))  
            
             
    def loadMasks(self, paramHandle, filename=None): 
        if type(paramHandle) == str:
            filename =  paramHandle
        if not filename:       
            filename, _ = QFileDialog.getOpenFileName(self,"Select a Mask file too use", "","Mask File (*.msk);;All Files (*)")
        if filename:
            mask_states = pickle.load(open(filename, "rb" ))
            for mask_state in mask_states:
                self.mask_list.append([pg.RectROI([0, 5], [1, 1], pen=(len(mask_states),9)),mask_state[1]])
                self.mask_list[-1][0].setState(mask_state[0])
                self.mask_list[-1][0].show()
                self.p3.addItem(self.mask_list[-1][0])
        self.convertMasks()            

    def saveroi(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","roi (*.roi);;All Files (*)")     
        if filename:
            state = [self.ctrROI.saveState(),self.angleRegion.getRegion()]
            pickle.dump(state, open(filename, "wb" ))   

    def loadroi(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            filename, _ = QFileDialog.getOpenFileName(self,"Select a ROI file too use", "","ROI File (*.roi);;All Files (*)")
        
        if filename:
            state = pickle.load(open(filename, "rb" ))
            self.ctrROI.setState(state[0])
            #self.angleRegion.setRegion(state[1])   

    def saveprofile(self, paramHandle, filename=None):
        if type(paramHandle) == str:
            filename =  paramHandle

        if not filename:
            filename, _ = QFileDialog.getSaveFileName(self,"Chose file name", "","csv (*.csv);;All Files (*)")  
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
        res=int(self.p.param('Data Processing', 'View 2 Resolution').value())
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())      
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()                     
        self.pixMaps=dector_frame_to_hkl_frame(self,binning,second_det)
        

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

    """Different Views are defined below"""    
        
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
        self.hideMasks()
        self.showMasks()

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
        if not self.mask_list and self.ImgState > 3:
            self.statusLabel.setText('Error: You must first select a mask')   
            return
        else:   
            projection = self.p.param('Data Processing', 'Select Projection').value()
            #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
            #Some messy variables this should be tidied up
            if self.newImgState == 1:
                self.newImgState = 0
            else: 
                self.ImgState=3+projection
            
            self.binBounds = []          
            self.convertMasks()

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
        projection = self.p.param('Data Processing', 'Select Projection').value()
        self.imgLeft.resetTransform()
                
        self.p3.getViewBox().setAspectLocked(True)          
        self.ctrROI.hide()
        self.p4.hide()
        self.angleRegion.setMovable(False)                      
        
        grid_resx = self.p.param('Data Processing', 'Grid Size Qx').value()
        grid_resy = self.p.param('Data Processing', 'Grid Size Qy').value()
        grid_resz = self.p.param('Data Processing', 'Grid Size Qz').value()

        self.busy = False

        tr = QTransform()  # prepare ImageItem transformation

        #Transform Detector View
        if self.ImgState == 1:
            tr.translate(self.xmin,self.zmin)
            tr.scale((self.xmax-self.xmin)/grid_resx,(self.zmax-self.zmin)/grid_resz)                    
        #HK
        elif projection == 0:
            tr.translate(-1*self.qrmax,-1*self.qrmax) 
            tr.scale((2*self.qrmax)/grid_resx,(2*self.qrmax)/grid_resy)                   
            #In the future we should transform for non-orthagonal coordinates here
        #HL    
        elif projection == 1:
            tr.translate(self.xmin,self.zmin) 
            tr.scale((self.xmax-self.xmin)/grid_resx,(self.zmax-self.zmin)/grid_resz)       
            
        #KL    
        elif projection == 2:
            tr.translate(self.ymin,self.zmin) 
            tr.scale((self.ymax-self.ymin)/grid_resy,(self.zmax-self.zmin)/grid_resz)      
            
        #3d
        elif projection == 3:    
            #This needs to be a Qt widget!!        
            importlib.reload(scene3d)
            scene3d.scene(self)

        self.imgLeft.setImage(self.showData,autoLevels=False) 
        self.imgLeft.setTransform(tr)
        self.qRegion.hide()
        self.statusLabel.setText('Vewing in Reciprocal lattice units, angluar integration is fixed, press detector view (CTRL+1) to go back')   
        
        self.progressBar.setValue(self.progressBar.maximum())
        #self.p3.autoRange()
        #self.update()   
        self.setupProjectionAxis()      
        self.hist.setImageItem(self.imgLeft)
        self.hideMasks() 
        self.showMasks()        
    
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

        self.pixMaps=dector_frame_to_hkl_frame(self,binning,second_det)  
        #these grids are coordinate grids that are used when selecting box profiles              
        #call the binning function
        
        #tempory for development
        importlib.reload(proj)
        self.xgrid,self.ygrid,self.zgrid,self.showData = proj.projection_binner(self,angles,image_data,progress_callback)       

    def make3dProjection(self):
        self.p.param('Data Processing', 'Select Projection').setValue(3)
        self.makeProjection()
        
    def detectorTransform(self):
        """Convert pixel coordinates to Reciprocal coordinates"""
        import time
        self.statusLabel.setText('Calculating reciprocal lattice coordinates (Please wait)..')
        self.update()
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
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
        res=int(self.p.param('Data Processing', 'View 2 Resolution').value())      

        #calculate the transformation[0]==h, [1]==k, [2]==l,[3]==sqrt(h^2 +k^2),[4]= intensity correciton
        self.pixMaps=dector_frame_to_hkl_frame(self,binning,second_det)
        #h_values = self.pixMaps[0].ravel()
        #k_values = self.pixMaps[1].ravel()
        l_values = self.pixMaps[2].ravel()
        hk_values = self.pixMaps[3].ravel()
        c_values = self.pixMaps[4]
    
        self.showData = self.image_stack.get_image(self.angleRegion.getRegion()[0], self.angleRegion.getRegion()[1])

        #create grid space for interpolation 
        self.grid_hk, self.grid_l=np.mgrid[np.min(hk_values):np.max(hk_values):res*1j,np.min(l_values):np.max(l_values):res*1j]        
    
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
        tr = QTransform()  # prepare ImageItem transformation:
        scaley = (np.max(l_values)-np.min(l_values))/res
        scalex = (np.max(hk_values)-np.min(hk_values))/res
        tr.translate(np.min(hk_values),np.min(l_values))
        tr.scale(scalex,scaley)       # scale horizontal and vertical axes

        
        self.imgLeft.setTransform(tr)    
        self.statusLabel.setText('Viewing in reciprocal lattice units, angluar integration is fixed, press Detector View to go back') 
        self.hideMasks() 
        self.showMasks()  

        #depending on the on the log intensity state, set ImgLeft to a logged showData
        if self.p.param('Data Processing', 'Log Intensity').value():
            self.logIntensity()  
        self.p3.autoRange()
        self.update()

    def theta2chiTransform(self):
        """Convert pixel coordinates to Reciprocal coordinates"""
        import time
        importlib.reload(chi2theta)
        self.statusLabel.setText('Calculating reciprocal lattice coordinates (Please wait)..')
        self.update()
        #IMG STATES: 0 = det view, 1 = trans det view, 2 - 2theta/chi, 3 = HK, 4=HL, 5=KL, 6=3d
        self.ImgState=2
        self.imgLeft.resetTransform()
        self.qRegion.hide()
        self.angleRegion.setMovable(False) 
        self.select_l_slice = True
        self.showData=None

        self.p3.getAxis('bottom').setGrid(255)
        self.p3.getAxis('left').setGrid(255)     
        self.p3.getAxis('bottom').setZValue(1)        
        self.p3.getAxis('left').setZValue(1)


        self.p3.getAxis('bottom').setLabel("2 THETA (Deg)")    
        self.p3.getAxis('left').setLabel("CHI (Deg)")        
                    
        binning=int(self.p.param('Data Processing', 'Binning (X by X)').value())
        second_det = self.p.param('Data Processing', 'Use 2nd detector').value()  
        res=int(self.p.param('Data Processing', 'View 2 Resolution').value())      

        #calculate the transformation[0]==h, [1]==k, [2]==l,[3]==sqrt(h^2 +k^2),[4]= intensity correciton
        self.pixMaps=chi2theta.dector_frame_to_theta2chi(self,binning,second_det)
        theta2_values = self.pixMaps[0].ravel()
        chi_values = self.pixMaps[1].ravel()
        c_values = self.pixMaps[2].ravel()
    
        self.showData = self.image_stack.get_image(self.angleRegion.getRegion()[0], self.angleRegion.getRegion()[1])

        #create grid space for interpolation 
        self.theta2grid, self.chigrid=np.mgrid[np.min(theta2_values):np.max(theta2_values):res*1j,np.min(chi_values):np.max(chi_values):res*1j]        
    
        #do 2d interpolation, and apply intensity corrections if needed
        if self.p.param('Data Processing', 'Apply Intensity Corrections').value():
            #we unravel each pixmap to get lists of coordinates (hk, l, correction)
            grid_i = interpolate.griddata((theta2_values,chi_values),np.ravel(c_values*np.rot90(self.showData,k=1)),(self.grid_theta2, self.grid_chi),method='linear',fill_value=0)
        else:
            grid_i = interpolate.griddata((theta2_values,chi_values),np.ravel(np.rot90(self.showData,k=1)),(self.grid_theta2, self.grid_chi),method='linear',fill_value=0)           
            

        #grid_i[grid_i<=0]=np.nan
        #replace showData with the 
        self.showData=grid_i
        self.gridl = self.grid_chi



        #set the image and scale the axis
        self.imgLeft.setImage(grid_i,autoLevels=False)
        #self.hist.setImageItem(self.imgLeft)
        tr = QTransform()  # prepare ImageItem transformation:
        tr.translate(np.min(theta2_values),np.min(chi_values)) 
        tr.scale((np.max(theta2_values)-np.min(theta2_values))/res,(np.max(chi_values)-np.min(chi_values))/res)          
        self.imgLeft.setTransform(tr)    

        self.statusLabel.setText('Viewing in reciprocal lattice units, angluar integration is fixed, press Detector View to go back') 
        self.hideMasks() 
        self.showMasks()  

        #depending on the on the log intensity state, set ImgLeft to a logged showData
        if self.p.param('Data Processing', 'Log Intensity').value():
            self.logIntensity()  
        self.p3.autoRange()
        self.p3.getViewBox().setAspectLocked(False)    
        self.update()
