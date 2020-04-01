#!/usr/bin/env python3
"""
Test description
================
- Fundamental Diagram in 2D, test number 102
- Width = 1.8 m
- Length = 26.0 m
- Measurement area: X = [10, 16],  Y = [-0.9, 0.9]

Remarks
=======
TODO: Compare two "clouds" of points and return a number.

Source
======

"""

import os
import sys

_author__ = 'Qiancheng'

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


# ==============
testname = "validation_4"
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
    logging.info("call makeini.py with -f %s/master_ini.xml" % DIR)
    subprocess.call(["python", "makeini.py", "-f", "%s/master_ini.xml" % DIR])

    # Specific the position of jpscore
    os.chdir("/home/qiancheng/Workspace/Simulations/jpscore/bin")
    TRUNK = os.getcwd()
    logging.info("change directory to %s" % TRUNK)
    executable = "%s/jpscore" % TRUNK
    if not path.exists(executable):
        executable="%s/jpscore.exe" % TRUNK #for windows
        if not path.exists(executable):
            logging.critical("executable <%s> does not exist yet." % executable)
            exit(0)

    # Generate trajectory and Clogging log
    os.chdir(DIR)
    logging.info("change directory to %s" % DIR)
    inifiles = glob.glob("inifiles/*.xml")
    for inifile in inifiles:
        if not path.exists(inifile):
            logging.critical("inifile <%s> does not exist" % inifile)
            exit(0)
        print(inifile)
        cmd = "%s --inifile=%s" % (executable, inifile)
        logging.info('----------------------------------\n start simulation with exe=<%s>' % cmd)
        subprocess.call([executable, "--inifile=%s" % inifile])
        logging.info('end simulation ...\n----------------------------------\n')
    exit(1)

