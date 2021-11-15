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
testname = "run"
# ==============

# -----------------------------
logfile = "log_%s.txt" % testname
f = open(logfile, 'w')
f.close()
logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

if __name__ == "__main__":
    # Generate inifiles for simulation
    subprocess.call(["python", "makeini.py", "-f", "supermarket_ini.xml"])

    # Specific the position of jpscore
    TRUNK = "F:/workspace/jpscore_master/jpscore/build/bin/Release"
    executable = "%s/jpscore" % TRUNK
    if not path.exists(executable):
        executable="%s/jpscore.exe" % TRUNK #for windows
        if not path.exists(executable):
            logging.critical("executable <%s> does not exist yet." % executable)
            exit(0)
    print("exexutable:",executable)
    # Generate trajectory and Clogging log
    inifiles = glob.glob("inifiles/*.xml")
    for inifile in inifiles:
        if not path.exists(inifile):
            logging.critical("inifile <%s> does not exist" % inifile)
            exit(0)
        print(inifile)
        cmd = "%s --inifile=%s" % (executable, inifile)
        logging.info('----------------------------------\n start simulation with exe=<%s>' % cmd)
        subprocess.call([executable, "%s" % inifile])
        logging.info('end simulation ...\n----------------------------------\n')
    exit(1)

