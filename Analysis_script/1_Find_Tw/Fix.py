#!/usr/bin/env python3
# -*- coding: utf-8 -*-

' Fix the clogging data '

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

# ==============
testname = "fix"
# ==============

# -----------------------------
logfile = "log_%s.txt" % testname
f = open(logfile, 'w')
f.close()
logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# ----------------------------------------------
HOME = path.expanduser("~")
DIR = path.dirname(path.realpath(argv[0]))
CWD = os.getcwd()


# ----------------------------------------------
if __name__ == "__main__":
    if CWD != DIR:
        logging.info("working dir is %s. Change to %s" % (CWD, DIR))
        os.chdir(DIR)

    os.chdir("../Jupedsim")
    TRUNK = os.getcwd()
    logging.info("change directory to %s" % TRUNK)
    executable = "%s/jpscore" % TRUNK
    if not path.exists(executable):
        executable="%s/jpscore.exe" % TRUNK #for windows
        if not path.exists(executable):
            logging.critical("executable <%s> does not exist yet." % executable)
            exit(FAILURE)

    # Generate trajectory and Clogging log
    os.chdir(DIR)
    logging.info("change directory to %s" % DIR)
    CLOGGINGLOG = "%s/cloggingLog" % DIR
    inifiles = glob.glob("inifiles/*.xml")

    # Move clogging log file to clogging log folder
    # Obtain seed_number factors clogging_times
    os.chdir('cloggingLog')
    logging.info("change directory to %s" % CLOGGINGLOG)
    Clogdata=[] #width,seed,case,clogging-times
    ClogLogs = glob.glob("*.txt")
    for ClogLog in ClogLogs:
        if ClogLog.find('CloggingLog') == -1:
            continue
        seed = ClogLog.split('seed_')[1].split('_')[0]
        width = ClogLog.split('geometry_')[1].split('_bottleneck')[0]
        Tw = ClogLog.split('Tw_')[1].split('.txt')[0]
        data = np.loadtxt(ClogLog)
        if data.size==5:
            clogging_times=1
        else:
            #clogging_times = np.max(data[:, 2])
            clogging_times = data[:,2].size
            # max = np.max(data[:,2])
            # for i in range(int(max)-1):
            #     Row=data[i+1,1]
            #     lastRow=data[i,1]
            #     if float(Row-lastRow)==float(Tw):
            #         max=max-1
            # clogging_times=max
        Clogdata.append([float(seed), float(width), float(Tw), float(clogging_times)])
    # obtain seed_number factors Evacuation_time
    Trajdata=[]
    Traj = ("%s/../trajectories" % CLOGGINGLOG)
    os.chdir(Traj)
    logging.info("change directory to %s" % Traj)
    TrajFiles = glob.glob("*.xml")
    for TrajFile in TrajFiles:
        if TrajFile.split('traj_') == -1:
            continue
        seed = TrajFile.split('seed_')[1].split('_')[0]
        width = TrajFile.split('geometry_')[1].split('_bottleneck')[0]
        Tw = TrajFile.split('Tw_')[1].split('.xml')[0]
        Trajdata.append([float(seed), float(width),float(Tw),0])
# Save data
    os.chdir(CWD)
    logging.info("change directory to %s" % CWD)
    Finaldata = []
    for i in Trajdata:
        found = 0
        for j in Clogdata:
            if (j[0] == i[0]) and (j[1] == i[1]) and (j[2] == i[2]):
                Finaldata.append([i[0], i[1], i[2], j[3], i[3]])  # seed,width,Tw,clogging_times,evacuation_time
                found = 1
        if found == 0:
            Finaldata.append([i[0], i[1], i[2], 0, i[3]])

    Finaldata = np.array(Finaldata)
    Widths = np.unique(Finaldata[:, 1])
    Tws = np.unique(Finaldata[:, 2])
    Width_datas = {}
    for Width in Widths:
        Width_data = []
        data = Finaldata[Finaldata[:, 1] == Width]
        for Tw in Tws:
            Tw_data = data[data[:, 2] == Tw]
            Width_data.append(
                [Tw, np.mean(Tw_data[:, 3]), np.std(Tw_data[:, 3]), np.mean(Tw_data[:, 4]), np.std(Tw_data[:, 4])])
        Width_datas[Width] = Width_data

        # Writing file
        result_file = 'Relation_between_Tw_and_Number_Width_%s.txt' % Width
        ff = open(result_file, 'w')
        ff.write(
            "# Tw\tnumber(mean)\tnumber(std)\tEva_time(mean)\tEva_time(std)\n")
        for data in Width_data:
            ff.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % (data[0], data[1], data[2], data[3], data[4]))
        ff.write("\n")
        ff.close()

    Colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'black', 'brown']
    # plot figure
    figure_name = 'Relation_between_Tw_and_Number.png'
    plt.figure(dpi=600)
    color = 0
    for Width in Widths:
        result = Width_datas[Width]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 1], yerr=result[:, 2], fmt="o-", color=Colors[color], linewidth=1.2,
                     capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 1], "o", color=Colors[color], label='Width=%s' % Width)
        color = color + 1
    plt.legend(loc='best', numpoints=1)
    plt.xlabel("waiting time before removing agents (s)")
    plt.ylabel("Clogging times")
    plt.savefig("%s" % figure_name)

    figure_name = 'Relation_between_Tw_and_Evacuation_time.png'
    plt.figure(dpi=600)
    color = 0
    for Width in Widths:
        result = Width_datas[Width]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 3], yerr=result[:, 4], fmt="o-", color=Colors[color], linewidth=1.2,
                     capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 3], "o", color=Colors[color], label='Width=%s' % Width)
        plt.legend(loc='best', numpoints=1)
        color = color + 1
    plt.xlabel("waiting time before removing agents (s)")
    plt.ylabel("Evacuation time (s)")
    plt.savefig("%s" % figure_name)

    exit(SUCCESS)
