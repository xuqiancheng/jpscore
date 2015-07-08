#!/usr/bin/env python
"""
Test description
================
300 pedestrians are distributed in a room with two exits.
The pedestrians should prefer the nearest exit,
but some should (spontaneously) choose the second exit

Remarks
=======

Source
======
http://www.rimea.de/fileadmin/files/dok/richtlinien/r2.2.1.pdf
"""

import os
import sys
import glob
utestdir = os.path.abspath(os.path.dirname(os.path.dirname(sys.path[0])))
from sys import *
sys.path.append(utestdir)
from JPSRunTest import JPSRunTestDriver
from utils import *
import warnings

def run_rimea_test11(inifile, trajfile):
    files = glob.glob("trajectories/*_exit*")
    if len(files) == 0:
        logging.critical("%s exists with failure! Found no exit-files.", argv[0])
        exit(FAILURE)

    for f in files:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            d = np.loadtxt(f)

        if len(d) == 0:
            logging.critical("File %s is empty", f)
            logging.critical("%s exists with failure!", argv[0])
            exit(FAILURE)

        num_evacuated = max(d[:, 1]) # >0 ?
        logging.info("%d peds evacuated from exit <%s>",
                     num_evacuated, f.split(".dat")[0].split("_")[1])



if __name__ == "__main__":
    test = JPSRunTestDriver(11, argv0=argv[0], testdir=sys.path[0], utestdir=utestdir)
    test.run_test(testfunction=run_rimea_test11)
    logging.info("%s exits with SUCCESS" % (argv[0]))
    exit(SUCCESS)










