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

# ==============
testname = "clogging_stepsize"
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

    # Generate inifiles for simulation
    logging.info("call makeini.py with -f %s/master_stepsize_ini.xml" % DIR)
    subprocess.call(["python", "makeini.py", "-f", "%s/master_stepsize_ini.xml" % DIR])
    # Specific the position of jpscore
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
    make_dir(CLOGGINGLOG)
    inifiles = glob.glob("inifiles/*.xml")
    for inifile in inifiles:
        if not path.exists(inifile):
            logging.critical("inifile <%s> does not exist" % inifile)
            exit(FAILURE)
        print(inifile)
        cmd = "%s --inifile=%s" % (executable, inifile)
        logging.info('----------------------------------\n start simulation with exe=<%s>' % cmd)
        subprocess.call([executable, "--inifile=%s" % inifile])
        logging.info('end simulation ...\n----------------------------------\n')

    # Move clogging log file to clogging log folder
    logfiles = glob.glob("inifiles/*.txt")
    for logfile in logfiles:
        if os.path.basename(logfile) != "log.txt.txt":
            shutil.move(logfile, CLOGGINGLOG)
            logging.info('Move %s to %s.' % (logfile, CLOGGINGLOG))


    # Obtain seed_number factors clogging_times
    os.chdir('cloggingLog')
    logging.info("change directory to %s" % CLOGGINGLOG)
    Clogdata=[] #width,seed,case,clogging-times
    ClogLogs = glob.glob("*.txt")
    for ClogLog in ClogLogs:
        if ClogLog.find('CloggingLog') == -1:
            continue
        seed = ClogLog.split('seed_')[1].split('_')[0]
        stepsize=ClogLog.split('stepsize_')[1].split('_update')[0]
        parallel=ClogLog.split('parallel_')[1].split('.txt')[0]
        data = np.loadtxt(ClogLog)
        if data.size==5:
            clogging_times=1
        else:
            clogging_times = np.max(data[:,2])
        Clogdata.append([float(seed), float(stepsize), float(parallel), float(clogging_times)])

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
        stepsize = TrajFile.split('stepsize_')[1].split('_update')[0]
        parallel = TrajFile.split('parallel_')[1].split('.xml')[0]
        tree = ET.parse(TrajFile)
        root = tree.getroot()
        header = root.find('header')
        frameRate = header.findtext('frameRate')
        for frame in root.iter('frame'):
            agent = frame.find('agent')
            if agent == None:
                continue
            Eframe = frame.attrib['ID']
        if float(Eframe) / float(frameRate)<2000:
            Trajdata.append([float(seed), float(stepsize),float(parallel),float(Eframe) / float(frameRate)])

    # Save data
    os.chdir(CWD)
    logging.info("change directory to %s" % CWD)
    Finaldata = []
    for i in Trajdata:
        found = 0
        for j in Clogdata:
            if (j[0] == i[0]) and (j[1] == i[1]) and (j[2]==i[2]):
                Finaldata.append([i[0], i[1], i[2], j[3], i[3]])#seed,width,case,clogging_times,evacuation_time
                found = 1
        if found == 0:
            Finaldata.append([i[0], i[1], i[2], 0, i[3]])

    Finaldata=np.array(Finaldata)
    Stepsizes=np.unique(Finaldata[:,1])
    Ups=np.unique(Finaldata[:,2])
    Up_datas={}
    for Up in Ups:
        Up_data=[]
        data=Finaldata[Finaldata[:,2]==Up]
        for Stepsize in Stepsizes:
            Stepsize_data=data[data[:,1]==Stepsize]
            Up_data.append([Stepsize,np.mean(Stepsize_data[:,3]),np.std(Stepsize_data[:,3]),np.mean(Stepsize_data[:,4]),np.std(Stepsize_data[:,4])])
        Up_datas[Up]=Up_data

        #Writing file
        result_file = 'Relation_between_stepsize_and_Number_parallel_%s.txt' % Up
        ff = open(result_file, 'w')
        ff.write(
            "# Stepsize(s)\tnumber(mean)\tnumber(std)\tEva_time(mean)\tEva_time(std)\n")
        for data in Up_data:
            ff.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % (data[0], data[1], data[2], data[3], data[4]))
        ff.write("\n")
        ff.close()


    Colors=['blue','green','red','cyan','magenta','black','brown']
    # plot figure
    figure_name='Relation_between_stepsize_and_Number.png'
    plt.figure(dpi=600)
    color = 0
    for Up in Ups:
        result = Up_datas[Up]
        result = np.array(result)
        plt.errorbar(result[:, 0], result[:, 1], yerr=result[:, 2], fmt="o-",color=Colors[color], linewidth=1.2, capthick=2, capsize=4)
        if Up==0:
            plt.plot(result[:, 0], result[:, 1], "o",color=Colors[color], label='Sequential update')
        else:
            plt.plot(result[:, 0], result[:, 1], "o", color=Colors[color], label='Parallel update')
        color=color+1
    plt.legend(loc='best', numpoints=1)
    plt.xlabel("Stepsize of update (s)")
    plt.ylabel("Clogging times")
    plt.savefig("%s"%figure_name)

    figure_name = 'Relation_between_stepsize_and_Evacuation_time.png'
    plt.figure(dpi=600)
    color = 0
    for Up in Ups:
        result = Up_datas[Up]
        result=np.array(result)
        plt.errorbar(result[:, 0], result[:, 3], yerr=result[:, 4], fmt="o-",color=Colors[color], linewidth=1.2, capthick=2, capsize=4)
        if Up==0:
            plt.plot(result[:, 0], result[:, 3], "o",color=Colors[color], label='Sequential update')
        else:
            plt.plot(result[:, 0], result[:, 3], "o", color=Colors[color], label='Parallel update')
        plt.legend(loc='best', numpoints=1)
        color = color + 1
    plt.xlabel("Stepsize of update (s)")
    plt.ylabel("Evacuation time (s)")
    plt.savefig("%s"%figure_name)

    exit(SUCCESS)
