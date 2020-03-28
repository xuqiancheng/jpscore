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
    SIM1 = "scenario_2_corridor_bidirection/script"
    SIM2 = "scenario_3_bottleneck/script"
    SIM3 = "scenario_4_four_exit_cross/script"

#Simulation1
    os.chdir("%s/%s" % (CWD,  SIM1))
    logging.info("change directory to %s/%s" % (CWD,  SIM1))
    logging.info("Start  SIM1 simulation:")
    os.system("python run_validation_x.py")
    logging.info("End SIM1")

#Simulation2
    os.chdir("%s/%s" % (CWD, SIM2))
    logging.info("change directory to %s/%s" % (CWD, SIM2))
    logging.info("Start SIM2 simulation:")
    os.system("python gettrajectory.py")
    logging.info("End SIM2")

#Simulation3
    os.chdir("%s/%s" % (CWD, SIM3))
    logging.info("change directory to %s/%s" % (CWD, SIM3))
    logging.info("Start SIM3 simulation:")
    os.system("python run_validation_4.py")
    logging.info("End SIM3")
