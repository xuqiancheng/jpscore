#!/usr/bin/env python3
import glob
import logging
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import subprocess
import sys
import time
from os import path
from sys import argv, exit

utestdir='/home/qiancheng/Workspace/Simulations/jpscore/Utest'
sys.path.append(utestdir)
from utils import SUCCESS, FAILURE

#-------------------- DIRS ------------------------------
HOME = path.expanduser("~")
DIR= os.path.dirname(os.path.realpath(argv[0]))
CWD = os.getcwd()#CURRENT WORK DIR
print('HOME:',HOME)
print('DIR:',DIR)
print('CWD:',CWD)
#--------------------------------------------------------

if __name__=="__main__":
    if CWD!=DIR:
        os.chdir(DIR)

    os.chdir(utestdir)
    subprocess.call(["python","makeini.py","-f","%s/master_ini.xml"%DIR])
    os.chdir(DIR)

    #-------- get directory of the code TRUNK

    TRUNK='%s/../'%utestdir
    print(TRUNK)
    #get trajectory
    inifiles=glob.glob("inifiles/*.xml")
    executable="%s/bin/jpscore"%TRUNK
    if not path.exists(executable):
        exit(FAILURE)

    for inifile in inifiles:
        if not path.exists(inifile):
            exit(FAILURE)
        print(inifile)
        subprocess.call([executable,"--inifile=%s"%inifile])



