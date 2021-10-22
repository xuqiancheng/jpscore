#!/usr/bin/env python3
# -*- coding: utf-8 -*-

' The influence of width of exit '

__author__ = 'Qiancheng'

import logging
import os
import subprocess
import sys
import glob
import shutil
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from os import path
from sys import argv, exit
from makeini import make_dir, SUCCESS, FAILURE

# ----------------------------------------------
if __name__ == "__main__":

    Widths=[]
    Width_datas={}
    files=glob.glob('*.txt')
    for file in files:
        if file.split('_')[0]=='Relation':
            Width=file.split('Width_')[1].split('.txt')[0]
            Widths.append(Width)
            data=np.loadtxt(file)
            Width_datas[Width]=data
        else:
            continue

    Widths=np.unique(Widths)
    Colors=['blue','green','red','cyan','magenta','black','brown']
    # plot figure
    figure_name='Relation_between_Tw_and_Number.png'
    plt.figure(dpi=600)
    color = 0
    for Width in Widths:
        result = Width_datas[Width]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 1], yerr=result[:, 2], fmt="o-",color=Colors[color], linewidth=1.2, capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 1], "o",color=Colors[color], label='Width=%s'%Width)
        color=color+1
    plt.legend(loc='best', numpoints=1)
    plt.xlabel("waiting time before removing agents (s)")
    plt.ylabel("Clogging times")
    plt.savefig("%s"%figure_name)

    figure_name = 'Relation_between_Tw_and_Evacuation_time.png'
    plt.figure(dpi=600)
    color = 0
    for Width in Widths:
        result = Width_datas[Width]
        result=np.array(result)
        plt.errorbar(result[:, 0], result[:, 3], yerr=result[:, 4], fmt="o-",color=Colors[color], linewidth=1.2, capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 3], "o",color=Colors[color], label='Width=%s' % Width)
        plt.legend(loc='best', numpoints=1)
        color = color + 1
    plt.xlabel("waiting time before removing agents (s)")
    plt.ylabel("Evacuation time (s)")
    plt.savefig("%s"%figure_name)

    exit(SUCCESS)
