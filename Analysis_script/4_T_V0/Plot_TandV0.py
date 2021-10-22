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

    V0s=[]
    V0_datas={}
    files=glob.glob('*.txt')
    for file in files:
        if file.split('_')[0]=='Relation':
            V0=file.split('V0_')[1].split('.txt')[0]
            V0s.append(V0)
            data=np.loadtxt(file)
            V0_datas[V0]=data
        else:
            continue

    V0s=np.unique(V0s)
    Colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 'brown']
    # plot figure
    figure_name = 'Relation_between_Ts_and_Number.png'
    plt.figure(dpi=600)
    color = 0
    for V0 in V0s:
        result = V0_datas[V0]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 1], yerr=result[:, 2], fmt="o-", color=Colors[color], linewidth=1.2,
                     capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 1], "o", color=Colors[color], label='V0=%s' % V0)
        color = color + 1
    plt.legend(loc='best', numpoints=1)
    plt.xlabel("Ts (s)")
    plt.ylabel("Clogging times")
    plt.savefig("%s" % figure_name)

    figure_name = 'Relation_between_Ts_and_Evacuation_time.png'
    plt.figure(dpi=600)
    color = 0
    for V0 in V0s:
        result = V0_datas[V0]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 3], yerr=result[:, 4], fmt="o-", color=Colors[color], linewidth=1.2,
                     capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 3], "o", color=Colors[color], label='V0=%s' % V0)
        plt.legend(loc='best', numpoints=1)
        color = color + 1
    plt.xlabel("Ts (s)")
    plt.ylabel("Evacuation time (s)")
    plt.savefig("%s" % figure_name)

    exit(SUCCESS)
