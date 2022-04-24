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

U = np.zeros(201)

F = np.zeros(201)

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

for i in range(len(U)):
    j = i 
    dispFile.write("%10.4E  %10.4E\n" % ((U[j]/20)*100,F[j]/1E6/400*(1+U[j]/20)))
dispFile.close()
#######################################################################
#######################################################################
#######                                                         #######
#######            Lattice strain calculations!                 #######
#######                                                         #######
#######################################################################
#######################################################################
print('\nExtracting lattice strain along LD')

LD = []
LD.append(assembly.elementSets['FCC-LD200'])
LD.append(assembly.elementSets['FCC-LD220'])
LD.append(assembly.elementSets['FCC-LD111'])
LD.append(assembly.elementSets['FCC-LD311'])
LD.append(assembly.elementSets['FCC-LD331'])
LD.append(assembly.elementSets['BCC-LD200'])
LD.append(assembly.elementSets['BCC-LD110'])
LD.append(assembly.elementSets['BCC-LD211'])
LD.append(assembly.elementSets['BCC-LD321'])

fp = open('lattice_strain_LD.dat','w')
for i in odb.steps['Step-1'].frames:
#    print('\nExtracting from Frame:\t'+str(i.frameId))
    # calculate the average stress
    stress = 0.0
    x = i.fieldOutputs['S'].values
    for j in x:
        stress = stress + j.data[1]
    stress = stress/len(x)/1E6
    fp.write("%10.4E  " % stress)
    for j in LD:
        temp = 0.0
        field = i.fieldOutputs['SDV5'].getSubset(region=j).values
        for k in field:
            temp = temp + k.data
        temp = temp/len(field)
        fp.write("%10.6E  " % temp)
    fp.write("\n")

    tempF = 0.0
    numF = 0
    tempB = 0.0
    numB = 0
    for j in LD[0:4]:
        field = i.fieldOutputs['S'].getSubset(region=j).values
        for k in field:
            tempF = tempF + k.data[1]
            numF = numF + 1
    for j in LD[5:8]:
        field = i.fieldOutputs['S'].getSubset(region=j).values
        for k in field:
            tempB = tempB + k.data[1]
            numB = numB + 1
    print("%10.3f  %10.3f" % (tempF/numF/1E6,tempB/numB/1E6)) 

fp.close()


#print('\nExtracting lattice strain along TD')

TD = []
TD.append(assembly.elementSets['FCC-TD200'])
TD.append(assembly.elementSets['FCC-TD220'])
TD.append(assembly.elementSets['FCC-TD111'])
TD.append(assembly.elementSets['FCC-TD311'])
TD.append(assembly.elementSets['FCC-TD331'])
TD.append(assembly.elementSets['BCC-TD200'])
TD.append(assembly.elementSets['BCC-TD110'])
TD.append(assembly.elementSets['BCC-TD211'])
TD.append(assembly.elementSets['BCC-TD321'])

fp = open('lattice_strain_TD.dat','w')
for i in odb.steps['Step-1'].frames:
#    print('\nExtracting from Frame:\t'+str(i.frameId))
    # calculate the average stress
    stress = 0.0
    x = i.fieldOutputs['S'].values
    for j in x:
        stress = stress + j.data[1]
    stress = stress/len(x)/1E6
    fp.write("%10.4E  " % stress)
    for j in TD:
        temp = 0.0
        field = i.fieldOutputs['SDV4'].getSubset(region=j).values
        for k in field:
            temp = temp + k.data
        temp = temp/len(field)
        fp.write("%10.6E  " % temp)

    fp.write("\n")

fp.close()
