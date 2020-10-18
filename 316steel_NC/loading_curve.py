# Yin Zhang , 17th Aug 2016
from odbAccess import *
import numpy as np
import math as mt
#import matplotlib as plt

#Reading output database

odbpath = 'J_8000.odb'

odb = openOdb(path=odbpath)

assembly = odb.rootAssembly

print 'Model data for ODB: ', odbpath

#Reading history output data

step = odb.steps['Step-1']

#print 'Node sets = ',odb.rootAssembly.nodeSets.keys()

#x_surface = assembly.nodeSets['_PickedSet18']

Node_Start = 1

Node_End = 9241

# U = displacement, F = loading force

U = np.zeros(1200)

F = np.zeros(1200)

for i in range(Node_Start,Node_End+21,21):
    node = 'Node PART-1-1.' + str(i)
    region = step.historyRegions[node]
    u1data = region.historyOutputs['U2'].data
    f1data = region.historyOutputs['RF2'].data
    j = 0
    for time,force in f1data:
        F[j] = F[j] + force
        j = j + 1

j = 0
for time,disp in u1data:
    U[j] = disp      
    j = j + 1

dispFile = open('fitting.dat','w')

for i in range(201):
    j = i - 1; 
    dispFile.write("%10.4E  %10.4E\n" % ((U[j]/20)*100,F[j]/1E6/400))
dispFile.close()
