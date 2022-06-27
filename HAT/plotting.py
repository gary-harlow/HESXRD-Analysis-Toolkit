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

    #np.savez_compressed(outfile_name,data=grid_i)

    #return()
    b1 = hat.crystal.param('b₁').value()

    def tr(h, k):
        h, k = np.asarray(h), np.asarray(k)
        return (np.sqrt(3)/2)*h*b1,(k+0.5*h)*b1

    def inv_tr(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return (x/(np.sqrt(3)/2))/b1,(y-0.5*x)/b1   


    plt.rcParams.update({'font.size': 4})
    plt.rc('legend', fontsize=8, handlelength=2)
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(8*cm,8*cm),dpi=600)
    #grid_helper = GridHelperCurveLinear((tr, inv_tr),grid_locator1=FixedLocator([-1.04,-1,-0.096,-0.092]),grid_locator2=FixedLocator([-2,-1,0,0.96,1,1.04,1.08]))
    #grid_helper = GridHelperCurveLinear((tr, inv_tr),grid_locator1=FixedLocator([-3,-2,-1,0,1,2,3]),grid_locator2=FixedLocator([-3,-2,-1,0,1,2,3]))
    grid_helper = GridHelperCurveLinear((tr, inv_tr),grid_locator1=FixedLocator([2,2.05,1.95]),grid_locator2=FixedLocator([-1,-0.95,-1.05]))
    ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)
    #ax1 = Subplot(fig, 1, 1, 1,)
    fig.add_subplot(ax1)
    grid_i = grid_i.T.astype(float)
    grid_i[grid_i <= 1] = np.nan  

    tmp = ax1.imshow(grid_i, extent=(np.min(grid_h),np.max(grid_h),np.min(grid_k),np.max(grid_k)), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis',interpolation='gaussian')  
    #plt.colorbar(tmp)
    #xp = [0,1,-1,0,-1,1,1,2,2,3,2,-1,-2,3,2,1,0,-1,2,1,-2,-3,-2,-1,0,1,-2,-1]
    #yp = [1,0,0,-1,1,-1,1,1,0,-1,-1,-1,-1,-2,-2,-2,-2,-2,-3,-3,0,1,2,2,2,2,3,3]
    #p = tr(xp,yp)
    ax1.grid(True,linewidth=0.05,color='black')
    #ax1.plot(p[0],p[1], 'wo', mfc='none',markersize=5,linewidth=0.002)
    #ax1.set_facecolor('cornflowerblue')
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))

    #plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0.2) 

    #plt.xlabel('$Q_{x}$ ($Å^{-1}$)')
    #plt.ylabel('$Q_{y}$ ($Å^{-1}$)')
    #ax1.set_aspect(2.881/4.087)
    plt.xlabel('$H$ (RLU)')
    plt.ylabel('$K$ (RLU)')
    #plt.xlim(1.9,2.4)
    plt.ylim(-0.3,0.3)
    plt.xlim(4.2,4.65)

    plt.savefig(outfile_name)
    plt.close()
    print("Plotting complete")

def plot_projection_hl(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    print("Plotting h-l projection but nothing set")


def plot_projection_kl(hat, grid_h,grid_k,grid_i,cmin,cmax, outfile_name):  
    print("Plotting k-l projection but nothing set")

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
    ax1.grid(True,linewidth=0.05,color='w')
    ax1.set_aspect(0.8)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.xlabel('$|HK|$ (RLU)')
    plt.ylabel('$L$ (RLU)')

    #plt.xlabel('$Q_{r}$ ($Å^{-1}$)')
    #plt.ylabel('$Q_{z}$ ($Å^{-1}$)')
    plt.savefig(outfile_name)
    plt.close()


def plot_2theta_chi(grid_2theta,grid_chi,grid_i,cmin,cmax,outfile_name):
    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    #see below for how to install fonts
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(20*cm,20*cm),dpi=600)
    ax1 = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax1)
    ax1.imshow(grid_i.T, extent=(grid_2theta.min(),grid_2theta.max(),grid_chi.min(),grid_chi.max()), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.xlabel('$2 Theta (degrees)')
    plt.ylabel('$Chi (degrees)')
    plt.savefig(outfile_name)
    plt.close()

def plot_det_view(grid_i,cmin,cmax,outfile_name):
    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    #see below for how to install fonts
    plt.rcParams.update({'font.sans-serif': 'Arial'})
    cm = 1/2.54  # centimeters in inches
    fig = plt.figure(figsize=(20*cm,20*cm),dpi=600)
    ax1 = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax1)

    ax1.imshow(grid_i.T, origin='lower',vmin=cmin,vmax=cmax,cmap='viridis')
    ax1.grid(True,linewidth=0.05,color='w')
    ax1.set_aspect(0.8)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    plt.xlabel('X Pixel')
    plt.ylabel('Y Pixel')
    plt.savefig(outfile_name)
    plt.close()
