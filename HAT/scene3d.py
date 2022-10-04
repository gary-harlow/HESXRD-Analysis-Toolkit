#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from PyQt6.QtWidgets import (QFileDialog)
#Export 3D data
def scene(hat):
    ''' This function outputs a 3d grid as an npz file'''
    data=hat.showData
    x = hat.xgrid
    y = hat.ygrid
    z = hat.zgrid  
    
    
    filename, _ = QFileDialog.getSaveFileName(hat,"Chose file name for 3d data", "","npz (*.npz);;All Files (*)")     
    if filename:  
       np.savez(filename, data=data,x=x,y=x,z=z)

