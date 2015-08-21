import numpy as np

#normalized hubble parameter
h = .68
dest = "closehalos/"

def isolated_system(lbound, ubound, halos):
   isSub = halos[9,:]
   xc = halos[3,:]
   yc = halos[4,:]
   zc = halos[5,:]
   mass = halos[1,:]

	#chooses a halo which does not have a similar mass halo within 5Mpc
	#begins halfway through the list so to get a mean mass
   for i in range (lbound, ubound, omithalos):
      tooClose = False
      #print "before loop: ", tooClose
		#calculates the distance from the target halo and checks to see if it
		#is within a certain radius. If so, does not meet requirements and is
		#passed.
      #print "lbound: ", lbound
      #print "ubound: ", ubound

		#makes sure that there aren't any similar or greater masses
		#within 5Mpc/h (7.1Mpc) of target halo
      for j in range(lbound, h_cat.shape[1]):
         radius = ((xc[i]-xc[j])**2 + (yc[i]-yc[j])**2 + (zc[i]-zc[j])**2)**.5
         if (radius < 5 and i != j and  isSub[j] == -1):
            print radius
            tooClose = True
            break

      print "after loop: ", tooClose
		#checks if the target halo has already been analyzed
      beenDone = False
      for omit in range (0, len(omithalos)):
         if (i == omithalos[omit]):
            beenDone = True

      #print "beenDone: ", beenDone
      #print "tooClose: ", tooClose
      #print "Subhalo?: ", isSub[i]
      #print ""

		#if everything checks out, target halo is returned
      if (tooClose == False and isSub[i] == -1 and beenDone == False):
         return i

   return -1

def paired_system(lbound, ubound, halos):
   isSub = halos[9,:]
   xc = halos[3,:]
   yc = halos[4,:]
   zc = halos[5,:]
   mass = halos[1,:]

   for target in range(lbound, ubound):
      for partner in range(target+1, h_cat.shape[1]):
         radius = ((xc[target]-xc[partner])**2 + (yc[target]-yc[partner])**2 + (zc[target]-zc[partner])**2)**.5
         mdiff = np.absolute(mass[target]-mass[partner])

         if radius < 4:
            print "mdiff:", mdiff
            print "radius:", radius
            print "partner:", isSub[partner] == -1
            print "target: ", isSub[target] == -1
            print "samehalo?", target == partner
         
         if(radius < 2 and target != partner and isSub[partner] == -1 and isSub[target] == -1 and mdiff < mass[target]/10.0):
            print radius
            print "FOUND!!!"
            return target

   return -1


def find_halo(lbound, ubound, halos, omithalo, system):

   if (system[0] == 'i'):
      result = isolated_system(lbound,ubound, halos)

   else if (system[0] == 'p'):
      result = paired_system(lbound, ubound, halos)

   return result


def save_halo(radius, halos, halo, quadrant):
   #calculates the relative radius of all halos to the target
   mass = halos[1,:]
   x = halos[3,:]
   y = halos[4,:]
   z = halos[5,:]
   vx = halos[6,:]
   vy = halos[7,:]
   vz = halos[8,:]
   r_vir= halos[11,:]
   m_vir= halos[10,:]
   print vz
   
   rad = ((x-x[halo])**2 + (y-y[halo])**2 + (z-z[halo])**2)**.5

   #finds where the relative radius is less than inputted Mpc
   boundhalos = np.where(rad < radius)
   print "boundhalos:", boundhalos

   #constructs desired components
   hlist = np.arange(0, rad[boundhalos].shape[0])
   rf = rad[rad<radius]
   mf = mass[rad<radius]
   xf = x[rad<radius]
   yf = y[rad<radius]
   zf = z[rad<radius]
   vxf = vx[rad<radius]
   vyf = vy[rad<radius]
   vzf = vz[rad<radius]
   r_virf = r_vir[rad<radius]
   m_virf = m_vir[rad<radius]
   listy = [hlist, rf, mf, xf, yf, zf, vxf, vyf, vzf]
   #creates a single numpy array and saves it into "halos" directory 
   output = np.array(listy)
   print output.shape
   properorder = np.argsort(output[1])
   newhalos = output[:,properorder]
   newhalos[0] = np.arange(0, newhalos.shape[1])
   targg = str(newhalos[2,0])
   np.save("allhalos/mass"+targg[0]+targg[2:]+"quadrant"+str(quadrant)+".npy", newhalos)


def construct(radius, lmass, umass, halos, quadrant, system):
   lbound = ubound = 0
   halolist = []
   isSub = halos[9,:]

   #sets up lower and upper bounds of the list
   for i in range(0,halos.shape[1]):
      if (halos[1,i] < lmass):
         lbound = i
      elif (halos[1,i] < umass):
         ubound = i
      else:
         break

   #sweeps through all the indices  in order to find a halo which
   #matches necessary conditions and then saves them
   omithalo = []
   while (lbound < ubound):
      #print "lbound:", lbound
      #print "ubound:", ubound
      #print "to go", ubound-lbound
   
      halo = find_halo(lbound, ubound, halos, omithalo, system)

      if isSub[lbound] == -1:
         halo = lbound
      else:
         halo = -1

      lbound += 1   

      #signifies the whole list was searched to no avail
      if (halo == -1):
         continue
      else:
         omithalo.append(halo)

      save_halo(radius, halos, halo, quadrant)

      #continues to loop after conditional b/c there may be two higher
      #mass halos which share a fractional difference with each other      


####################MAIN#######################
system = raw_input("Looking for isolated systems or similar mass bound pairs (i or p)?\n")


#loops through every quadrant finding specifications included in "construct" arguments
for i in range(0, 20):
   print "BEGINNING OF TRIAL NUMBER", i+1, "###################################################"
   if (i < 10):
      halos = np.load("quadrants/halos_0214_0"+str(i)+".npy")
   else:
      halos = np.load("quadrants/halos_0214_"+str(i)+".npy")

   construct(10, 2*10**12, 10**13, halos, i, system)

