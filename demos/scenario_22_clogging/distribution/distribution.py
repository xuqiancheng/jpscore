import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import glob

# Loading Geometry file
Geo=glob.glob('*.xml')
tree=ET.parse(Geo[0])
root=tree.getroot()
plt.figure(1)
for polygon in root.iter('polygon'):
    GeoPoint = []
    for vertex in polygon.iter('vertex'):
        GeoPoint.append([float(vertex.attrib['px']),float(vertex.attrib['py'])])
    GeoPoint = np.array(GeoPoint)
    plt.plot(GeoPoint[:,0], GeoPoint[:,1],color='blue',linewidth=2)
for transition in root.iter('transition'):
    TraPoint = []
    for vertex in transition.iter('vertex'):
        TraPoint.append([float(vertex.attrib['px']),float(vertex.attrib['py'])])
    TraPoint = np.array(TraPoint)
    plt.plot(TraPoint[:, 0], TraPoint[:, 1], color='black',linestyle='--', linewidth=1)


#Loading position of clogging
dataname=glob.glob('*.txt')
data=np.loadtxt(dataname[0])
x=data[:,3]
y=data[:,4]
plt.scatter(data[:,3],data[:,4],edgecolors='black', linewidth='1',
            color='red', alpha=0.6, s=20, label='clogging' )
plt.title("Distribution of cloggings")
plt.legend(loc='upper left', numpoints=1)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.xlim(4,12)
plt.ylim(0,8)
plt.savefig("distribution_%s.png"% dataname[0].split('.txt')[0].split('_')[-1],dpi=300)


