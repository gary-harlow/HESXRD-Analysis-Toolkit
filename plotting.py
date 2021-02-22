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


#inplane map
def plot_in_plane(grid_h,grid_k,grid_i,cmin,cmax, outfile_name):

    def tr(h, k):
        h, k = np.asarray(h), np.asarray(k)
        return (np.sqrt(3)/2)*h,(k+0.5*h)

    def inv_tr(x, y):
        x, y = np.asarray(x), np.asarray(y)
        return x/(np.sqrt(3)/2),(y-0.5*x)
        

    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    #see below for how to install fonts
    plt.rcParams.update({'font.sans-serif': 'Arial'})

    """Example plotting function for inplane map"""


    fig = plt.figure(figsize=(5,5),dpi=300)
    grid_helper = GridHelperCurveLinear((tr, inv_tr))
    ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)
    #ax1 = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax1)

    ax1.imshow(grid_i.T, extent=(grid_h.min(),grid_h.max(),grid_k.min(),grid_k.max()), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis')
    #ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
    #ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax1.grid(True,linewidth=0.1,color='black')
    ax1.set_aspect(1.)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    plt.xlabel('$H$ (RLU)')
    plt.ylabel('$K$ (RLU)')
    plt.xlim((-2,2))
    plt.ylim((-2.2,0))
    #plt.xlabel('$Q_{x}$ ($Å^{-1}$)')
    #plt.ylabel('$Q_{y}$ ($Å^{-1}$)')
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.2))

    plt.savefig(outfile_name)
    plt.close()


def plot_out_of_plane(grid_hk,grid_l,grid_i,cmin,cmax,outfile_name):
    plt.rcParams.update({'font.size': 8})
    plt.rc('legend', fontsize=8, handlelength=2)
    #see below for how to install fonts
    plt.rcParams.update({'font.sans-serif': 'Arial'})

    """Example plotting function for inplane map"""
    fig = plt.figure(figsize=(5,5),dpi=300)
    ax1 = Subplot(fig, 1, 1, 1)
    fig.add_subplot(ax1)

    ax1.imshow(grid_i.T, extent=(grid_hk.min(),grid_hk.max(),grid_l.min(),grid_l.max()), origin='lower',vmin=cmin,vmax=cmax,cmap='viridis')
    ax1.grid(True,linewidth=0.1,color='w')
    #ax1.set_aspect(1.)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.xlabel('$\sqrt{H^{2}+K^{2}}$ (RLU)')
    plt.ylabel('$L$ (RLU)')

    plt.xlabel('$Q_{r}$ ($Å^{-1}$)')
    plt.ylabel('$Q_{z}$ ($Å^{-1}$)')
    plt.savefig(outfile_name)
    plt.close()

#inplane map
def plot_in_plane2(grid_h,grid_k,grid_i,cmin,cmax, outfile_name):
    """Example plotting function for hexagonal coordinates"""

    #transformations to and from hexagonal coordinates
    trans= np.array([[-4/3,-2/3],[2/3,-2/3]])


    def tr2(h, k):
        bulk = [h,k]
        surf = np.linalg.solve(trans,bulk)
        return surf[0],surf[1]

    points_x, points_y=tr(points_h,points_k)
    #points_x,points_y = tr2(points_h2,points_k2)

    print("Binning Coordinates")
    H, xedges, yedges = np.histogram2d(points_x,points_y, weights=values, bins=2048)
    H = H.T #this is important!!
    X, Y = np.meshgrid(xedges[0:-1], yedges[0:-1])

    print("Plotting")
    H[H == 0] = np.nan #set values to measured to nan so they don't show in plot


    fig = plt.figure(figsize=(5,5),dpi=300)
    SMALL_SIZE = 4
    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)

    grid_helper = GridHelperCurveLinear((tr, inv_tr))
    ax1 = Subplot(fig, 1, 1, 1, grid_helper=grid_helper)
    fig.add_subplot(ax1)
    #grid_x, grid_y=tr(X,Y)
    #ax1.plot(xx, yy,'xr')
    ax1.pcolormesh(X,Y,
                H,vmin=50,vmax=25000,cmap='viridis')

    xx, yy = tr([-1,1,0,0,1,-1],[0,0,1,-1,-1,1])
    #ax1.plot(0, 0, 'ro',markersize=1,fillstyle='none')
    ax1.plot(xx, yy,'ro',markersize=15,fillstyle='none',markeredgewidth=0.5)
    #plt.text(0.08, -0.13, '(0,0)', fontsize=6,color='red')
    #plt.text(tr(1-0.09,-0.09), '(1,0)', fontsize=6,color='red')

    #xx2, yy2 = tr([1,-1,1,-1],[-2,-1,1,-1])
    xx2, yy2 = tr([1,-1,0,2,-2,0,-2,2,-1,-2,1,2,2],[1,-1,2,0,0,-2,1,-1,2,2,-2,-1,-2])
    #ax1.plot(xx2, yy2, 'bo',markersize=15,fillstyle='none',markeredgewidth=0.5)

    #ax1.plot(points_h,points_k)
    #plt.xlim((-1.7,1.7))
    #plt.ylim((-1.8,1.9))
    plt.grid(color='r',linestyle='-', linewidth=0.2)
    ax1.grid(True)
    plt.xlabel('$H_{surf}$ RLU')
    plt.ylabel('$K_{sur}$ RLU')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.6, hspace=0.6) 
    plt.savefig(outfile_name)

