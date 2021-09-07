# -*- coding: utf-8 -*-
#Plotting figures in matplotlib, eg. for publications 
#Note: this file is intended to be modified, and can be changed at runtime.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear
from mpl_toolkits.axisartist import Subplot
import matplotlib.ticker as ticker
from scipy import interpolate
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)


#inplane map
def plot_projection_hk(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    #the easiet way to do none cubic/rectangualr system to work in 
    #Q units and then covert when plotting. 
    def tr(h, k):
        h, k = np.asarray(h), np.asarray(k)
        return (np.sqrt(3)/2)*h*2.51,(k+0.5*h)*2.51

    def inv_tr(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return (x/(np.sqrt(3)/2))/2.51,(y-0.5*x)/2.51     

    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(20*cm,20*cm),dpi=600)

    grid_helper = GridHelperCurveLinear((tr, inv_tr),grid_locator1=FixedLocator([0.9,0.95,1,1.05,1.1]),grid_locator2=FixedLocator([-0.1,-0.05,0,0.05,0.1,1]))
    ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    grid_i = grid_i.T.astype(float)
    qmax = np.max(hat.grid_qr)  
    ax1.imshow(grid_i, extent=(-1*qmax,qmax,-1*qmax,qmax), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis',interpolation='bicubic')

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.8) 
    plt.xlabel('$H$ (RLU)')
    plt.ylabel('$K$ (RLU)')
    #plt.xlabel('$Q_{x}$ ($Å^{-1}$)')
    #plt.ylabel('$Q_{y}$ ($Å^{-1}$)')

    grid_i[grid_i <= 20] = np.nan  
    plt.xlim(2,2.5)
    plt.ylim(1.05,1.5)
    plt.savefig(outfile_name)
    plt.close()

def plot_projection_hl(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    print("Plotting h-l projection but nothing set")


def plot_projection_kl(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    print("Plotting k-l projection but nothing set")


def plot_projection_qrl(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    print("Plotting qr-l projection but nothing set")


def plot_transformed_detector(grid_hk,grid_l,grid_i,cmin,cmax,outfile_name):
    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    #see below for how to install fonts
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(20*cm,20*cm),dpi=600)
    ax1 = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax1)

    ax1.imshow(grid_i.T, extent=(grid_hk.min(),grid_hk.max(),grid_l.min(),grid_l.max()), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis')
    #ax1.grid(True,linewidth=0.05,color='w')
    ax1.set_aspect(1)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.xlabel('$\sqrt{H^{2}+{(rK)}^{2}}$ (RLU)')
    plt.ylabel('$L$ (RLU)')

    #plt.xlabel('$Q_{r}$ ($Å^{-1}$)')
    #plt.ylabel('$Q_{z}$ ($Å^{-1}$)')
    plt.savefig(outfile_name)
    plt.close()

