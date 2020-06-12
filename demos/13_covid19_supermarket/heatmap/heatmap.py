#!/usr/bin/env python3
# https://towardsdatascience.com/simple-example-of-2d-density-plots-in-python-83b83b934f67

import sys
import time
import glob
import os
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from xml.dom.minidom import parse
import pandas as pd

from Utilities import read_subroom_walls
from Utilities import plot_geometry
from Utilities import geo_limits


def get_profile(data, xbins, ybins, Id,geometry_wall):
    """ make profile plots for density and velocity
    This method, discritises the geometry in regular grids.
    The value of every grid cell is the mean of <values> for points within each bin.
    Empty bins will be represented by 0.
    Note: <values> here means 1/A_i (density) or v_i (velocity) of Agent i

    :param data: pandas array containing IFD-Date (calculated by method I)
    :param xbins: Discretisation of the geometry in x-axis
    :param ybins: Discretisation of the geometry in y-axis
    :param Id: A number. Useful to produce movies
    :param geometry_wall: Geometry
    :returns: plots two figures: Density and Velocity profiles
    :rtype:

    """
    t = time.process_time()
    dx = xbins[1] - xbins[0]
    #print(dx)
    ret = stats.binned_statistic_2d(data['x'], data['y'], data['CD'], 'mean', bins=[xbins, ybins])
    prof = np.nan_to_num(ret.statistic.T)

    methods = ['none', 'nearest', 'bilinear', 'bicubic', 'spline16',
               'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
               'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

    for m in ['bicubic']:  # methods:
        fig, ax = plt.subplots(1, 1)
        im = ax.imshow(prof,
                         cmap=cm.jet,
                         interpolation=m,
                         origin='lower',
                         vmin=0, vmax=5,  # np.mean(data['d']) + np.mean(data['d']),
                         extent=[geominX, geomaxX, geominY, geomaxY])

        plot_geometry(ax, geometry_wall)
        ax.set_xlim(geominX,geomaxX)
        ax.set_ylim(geominY, geomaxY)
        # plot_peds(ax)
        #ax.set_title("fr: [{}-{}], dx={:.2f}, dy={:.2f}".format(begin_frame, end_frame, dx, dy))
        #figname = m + "_{:03d}_dx_{:.2f}_dy_{:.2f}_density.png".format(Id, dx, dy)
        figname="HeatMap.png"
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3.5%", pad=0.3)
        cb = plt.colorbar(im, cax=cax)
        cb.set_label('Distance Index', rotation=270, labelpad=15)
        fig.savefig(figname, dpi=300)
        elapsed_time = time.process_time() - t
        print("time: {:.3f} s".format(elapsed_time))
    return np.mean(prof), np.std(prof), np.max(prof), np.min(prof)

# Merge all the trajectories in to one file
def MergeTrajFile(path,name):
    filenames = os.listdir(path)
    result = name
    file = open(result, 'w+', encoding="utf-8")
    for filename in filenames:
        filepath = path + '/'
        filepath = filepath + filename
        for line in open(filepath, encoding="utf-8"):
            file.writelines(line)
        file.write('\n')
    file.close()

# Read trajectories as dataframe
def read_traj(traj_filename):
    df = pd.read_csv(traj_filename,
                     comment='#',sep='\t',
                     names=['i', 'f', 'x','y','z','A','B','ANGLE','COLOR','Inf','VC','VG','PI','CD'],
                     index_col=False)
    return df

#Get the min and max frame to identify the steady state, where the number of pedestrians in the room is constant
def GetSteadyState(traj_filename,MaxPedInRoom):
    begin_frame=0
    end_frame=0
    data=np.loadtxt(traj_filename)
    frames=np.unique(data[:,1])
    for frame in frames:
        NumIn=0
        dataf=data[data[:,1]==frame,:]
        for datafi in dataf:
            if datafi[2]>0:
                NumIn=NumIn+1
        if NumIn==MaxPedInRoom:
            end_frame=frame
            if begin_frame==0:
                begin_frame=frame
    return begin_frame,end_frame

if __name__ == '__main__':

    #Max number of pedestrians in the supermarket
    MaxPedInRoom=50

    # I need to read the document provided by mohcine
    dx = 0.2
    dy = 0.2
    print("dx: {:.2f}, fy: {:.2f}".format(dx, dy))

    # Read Geometry file
    geo_filename = 'supermarket_geo.xml'
    if not os.path.exists(geo_filename):
        sys.exit("{} does not exist".format(geo_filename))
    xml_date = open(geo_filename, "r")
    geo_xml = parse(xml_date)
    xml_date.close()
    geometry_wall = read_subroom_walls(geo_xml)
    geominX, geomaxX, geominY, geomaxY = geo_limits(geo_xml)
    geominX=0
    print("GeometrySize: X: ({:.2f},{:2f}), Y: ({:.2f},{:2f})".format(geominX, geomaxX, geominY, geomaxY))

    # Read Trajectories
    traj_dir = 'trajectories'
    traj_filename='traj_merge.txt'
    MergeTrajFile(traj_dir, traj_filename)
    if not os.path.exists(traj_filename):
        sys.exit("{} does not exist".format(traj_filename))
    begin_frame,end_frame=GetSteadyState(traj_filename,MaxPedInRoom)
    print("Begin steady: {}".format(begin_frame))
    print("End steady: {}".format(end_frame))

    #Data filter
    data = read_traj(traj_filename)
    if begin_frame is None or begin_frame < data['f'].iloc[0]:
        begin_frame = data['f'].iloc[0]
        print("Change begin frame to {}".format(begin_frame))
    if end_frame is None or end_frame > data['f'].iloc[-1]:
        print("Change end frame to {}".format(end_frame))
        end_frame = data['f'].iloc[-1]
    # ###################################################
    Id = 1
    xbins = np.arange(geominX, geomaxX + dx, dx)
    ybins = np.arange(geominY, geomaxY + dx, dy)
    data = data[(data['f'] >= begin_frame) & (data['f'] <= end_frame)]
    Mean, Std, Max, Min = get_profile(data, xbins, ybins, Id, geometry_wall)
    ###################################################
    print("{:.3f}  {:.3f}   {:.3f}   {:.3f}    {:.3f}".format(dx, Mean, Std, Max, Min))
    print("Max contact degree: {:.3f}".format(np.max(data["CD"])))
    print("Mean contact degree: {:.3f} +-{:.3f}".format(np.mean(data["CD"]), np.std(data["CD"])))
