#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys,os
import numpy as np
from scipy import interpolate

def add_to_divide(number, divisor):
    """ this function returns the int needed to be added to number to make it divisable by divisor"""
    tmp = number%divisor
    if tmp == 0:
        return 0
    return divisor - tmp               

def rebin(data,binlen):
    """rebins the data but pads the edge with zero if it isn't an exact divisor"""
    xlen = int(np.shape(data)[0])
    ylen = int(np.shape(data)[1])

    extra_x = add_to_divide(xlen,binlen)
    extra_y = add_to_divide(ylen,binlen)
    if extra_x != 0:         
        tmp = np.zeros((extra_x,ylen))
        data = np.concatenate((data,tmp),axis=0)
    if extra_y != 0: 
        tmp2 = np.zeros((xlen+extra_x,extra_y,))
        data = np.concatenate((data,tmp2),axis=1)    

    return data.reshape((xlen+extra_x)//binlen,binlen,(ylen+extra_y)//binlen,binlen).mean(axis=3).mean(axis=1)   



