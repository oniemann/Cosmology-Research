import numpy as np
import os

path = 'testhalos/'

#input = raw_input("choose halo (no.npy)\n")
#halos = np.load(input+'.npy')




def swapindices(halos, haloname):
   halos.shape[1]
   #bigarray = np.arange(0,halos.shape[1]) + 1
   #halos[0] = bigarray
   print halos[0]
   partner = np.argmax(halos[2,1:50])+1
   
   #mass
   #print halos[2,0:10]
   #print halos[:,0]
   #print halos[:,4]
   
   halos[:,[0,partner]] = halos[:,[partner,0]]
   
   
   #print halos[:,0]
   #print halos[:,4]
   
   H_0 = 70
   x = halos[3,:]
   y = halos[4,:]
   z = halos[5,:]
   
   
   #Reconfigure the radii such that they reflect the new mass as target halo
   halos[1] = ((x-x[0])**2 + (y-y[0])**2 + (z-z[0])**2)**.5
   
   #print halos[0,0:10]
   #print halos[1,0:10]
   
   
   properorder = np.argsort(halos[1])
   haloshift = halos[:,properorder]
   print haloname
   np.save(path+haloname+'shift.npy', haloshift)




#MAIN: ########################################
files = []

for root, dirnames, filenames in os.walk(path):
   for file in filenames:
      files.append(file)

halofiles = [x for x in files if '.npy' in x]
trial = 0

for halofile in halofiles:
   if '.png' in halofile:
      continue
   trial += 1
   halo = path+halofile
   string = str(halofile)
   dot = string.find('.')
   haloname = string[0:dot]
   halos  = np.load(halo)
   swapindices(halos, haloname)
   print 'trial number:', trial

