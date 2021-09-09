#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Author:Qiancheng

import os
import glob
import makeini
import subprocess

if __name__=='__main__':
    cwd=os.getcwd()
    Jupedsim="F:\workspace\jpscore_master\jpscore\\build\\bin\Release"
    makeini.main('%s\master_ini.xml'%cwd)
    inifiles = glob.glob("inifiles/*.xml")
    executable="%s/jpscore.exe" % Jupedsim
    for inifile in inifiles:
        subprocess.call([executable, "%s" % inifile])

