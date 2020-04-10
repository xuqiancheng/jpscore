#!/usr/bin/python
# import matplotlib
# matplotlib.use('macosx')
#You need to put traj.xml into folder /trajectories

from sys import exit
import numpy as np
import glob
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import shutil

cwd=os.getcwd()
Gap=1 # The time gap for avergae speed
vmin=0.1 # We see it as clogging if the average speed smaller than vmin

left=0
right=26
up=12
down=-6
vmax=3 # We need this to deal with the position changing of pedestrians becasue of the periodic case

def mkdir(path):
    if os.path.exists(path):
        shutil.rmtree(path)
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

if __name__=="__main__":
    trajxmls=glob.glob("%s/trajectories/*.xml" % cwd)
    FailFile = []
    Traj_data={}
    txt_folder='txt_average_speed'
    png_folder='png_average_speed'
    mkdir(txt_folder)
    mkdir(png_folder)
    for trajxml in trajxmls:
        tree=ET.parse(trajxml)
        trajxml=os.path.basename(trajxml)
        root=tree.getroot()
        frameRate=float(root.find('header').find('frameRate').text)
        agentsNumber=float(root.find('header').find('agents').text)
        frameGap = int(frameRate * Gap)
        frames=root.findall('frame')
        datas=[]
        for frame in frames:
            time=float(frame.attrib['ID'])
            agents=frame.findall('agent')
            if len(agents)!=agentsNumber:
                FailFile.append(trajxml)
                break
            for agent in agents:
                ID=float(agent.attrib['ID'])
                x=float(agent.attrib['x'])
                y=float(agent.attrib['y'])
                data=[time,ID,x,y]
                datas.append(data)
        datas=np.array(datas)

        ids=np.unique(datas[:,1])
        speeds=[]
        for id in ids:
            data=datas[datas[:,1]==id,:]
            data=data[::frameGap]
            i=0
            while i<len(data)-1:
                frame=(data[i][0]+data[i+1][0])/2
                x_gap=np.abs(data[i][2]-data[i+1][2])
                if x_gap>vmax*Gap:
                    x_gap=right-left-x_gap
                y_gap=np.abs(data[i][3]-data[i+1][3])
                if y_gap>vmax*Gap:
                    y_gap=up-down-y_gap
                speed=np.sqrt(x_gap*x_gap+y_gap*y_gap)/Gap
                speeds.append([id,frame,speed])
                i=i+1
        speeds=np.array(speeds)
        #speeds [id,frame,speed]

        ave_speed=[]
        frames=np.unique(speeds[:,1])
        for frame in frames:
            speed=speeds[speeds[:,1]==frame,:]
            mean=np.mean(speed[:,2])
            std=np.std(speed[:,2])
            ave_speed.append([frame/frameRate,mean,std])
        ave_speed=np.array(ave_speed)
        Traj_data[trajxml]=ave_speed
        #average_speed [time, mean, std]

        if ave_speed[-1][1]<vmin:
            FailFile.append(trajxml)
        # writing data
        result_file='average_speed_'+trajxml.split('.xml')[0]+'.txt'
        print(result_file)
        ff=open(result_file,'w')
        ff.write('# Frame\tSpeed_mean\tSpeed_std\n')
        for data in ave_speed:
            ff.write('%f\t%f\t%f\n'%(data[0],data[1],data[2]))
        ff.write("\n")
        ff.close()
        shutil.move(result_file,txt_folder)

        #Plot
        figure_name='average_speed_'+trajxml.split('.xml')[0]+'.png'
        print(figure_name)
        result=ave_speed
        #plt.errorbar(result[:, 0], result[:, 1], yerr=result[:, 2], fmt="o-",color='r', linewidth=1.2,capthick=2, capsize=4)
        plt.plot(result[:, 0], result[:, 1], "-", color='r', label='Average speed of all pedestrians')
        plt.legend(loc='best', numpoints=1)
        plt.xlabel("t (s)")
        plt.ylabel("Speed (m/s)")
        plt.savefig("%s" % figure_name,dpi=600)
        shutil.move(figure_name, png_folder)
        plt.clf()

    #writing failed simulation
    result_file = 'Failed_Simulation.txt'
    FailFile=np.unique(FailFile)
    ff = open(result_file, 'w')
    ff.write('# These simulations are failed\n')
    for file in FailFile:
        ff.write('%s\n'%file)
    ff.write("\n")
    ff.close()
    exit(0)

















