#!/usr/bin/env python3
# -*- coding: utf-8 -*-

' Control the experiment '

import logging
from os import path
from sys import argv
import os


# ==============
testname = "control"
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
    EXP1 = 'scenario1'
    EXP2 = 'scenario2'
    EXP3 = 'scenario3'
    EXP4 = 'scenario4'


    # EXP1
    os.chdir("%s/%s" % (CWD, EXP1))
    logging.info("change directory to %s/%s" % (CWD, EXP1))
    logging.info("Start EXP1 simulation:")
    os.system("python run.py")
    logging.info("End EXP1")

    # EXP2
    os.chdir("%s/%s" % (CWD, EXP2))
    logging.info("change directory to %s/%s" % (CWD, EXP2))
    logging.info("Start EXP2 simulation:")
    os.system("python run.py")
    logging.info("End EXP2")

    # EXP3
    os.chdir("%s/%s" % (CWD, EXP3))
    logging.info("change directory to %s/%s" % (CWD, EXP3))
    logging.info("Start EXP3 simulation:")
    os.system("python run.py")
    logging.info("End EXP3")

    # EXP4
    os.chdir("%s/%s" % (CWD, EXP4))
    logging.info("change directory to %s/%s" % (CWD, EXP4))
    logging.info("Start EXP4 simulation:")
    os.system("python run.py")
    logging.info("End EXP4")


