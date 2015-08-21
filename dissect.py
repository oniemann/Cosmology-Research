from __future__ import division
from numpy import *
import numpy as np
h_cat = np.load("BolshoiSortedMass.npy")

#print "input format: dissect(lowx, highx, lowy, highy, lowz, highz)"
xcoord = h_cat[3,:]
ycoord = h_cat[4,:]
zcoord = h_cat[5,:]

def dissect(lowx, highx, lowy, highy, lowz, highz, lol):
	halolist = []
	haloradius = []
	mvir = []
	x = []
	y = []
	z = []
	vx = []
	vy = []
	vz = []
	mass = []
	issub = []
	for i in range(0, h_cat.shape[1]):
		if xcoord[i] > lowx and xcoord[i] < highx:
			if ycoord[i] > lowy and ycoord[i] < highy:
				if zcoord[i] > lowz and zcoord[i] < highz:
					halolist.append(i)
					mvir.append(h_cat[2,i])
					x.append(h_cat[3,i])  #""
					y.append(h_cat[4,i])  #""
					z.append(h_cat[5,i])  #""
					vx.append(h_cat[6,i])
					vy.append(h_cat[7,i])
					vz.append(h_cat[8,i])
					mass.append(h_cat[1,i])   #M_sun
					issub.append(h_cat[9,i])

	updatedlist = np.array([halolist, mass, mvir, x, y, z, vx, vy, vz, issub])
	np.save('quadrants/halos_0214_'+str(lol), updatedlist)

######################MAIN#############################
List = np.load("quadrants.npy")
for i in range(0, List.shape[0]):
	print "lowx:", List[i,0]
	print "highx:", List[i,3]
	dissect(List[i,0], List[i,3], List[i,1], List[i,4], List[i,2], List[i,5], i)
	print "percent: ", i / List.shape[0]
