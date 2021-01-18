#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HESXRD Analysis Toolkit

This softwares allows one to extract meaningful data from
stacks of images such as those produced in high-energy-surface
x-ray diffraction experiemnts. 

Images can be addressed via image number, or angle if that is available.
Stacks can be converted to reciprocal space units, and rotations of stacks
(where there is angle information) can also be cut to produce an in-plane map. 

It is designed to be user-friendly and fairly quick so that it can
be easily used during an experiment. For there are some limitations in that the 
in-plane lattice vectors are expected to be 90 degrees.

For more general reciprocal space transformaitons abd binning
See the ESRF Binoculars software. 

@author: Sebastian Pfaff and Gary S. Harlow
"""

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

import interface

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = interface.Ram()
    sys.exit(app.exec_())
