import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import os
import matplotlib.lines as lines
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from math import log


h = 0.68
H_0 = h * 100
path = 'ViableTargets/aclosehalos/'

#KINEMATIC FUNCTIONS
#returns the relative velocity, taking into account hubble expansion
def relative_velocity(xnew,ynew,znew,x,y,z,vx,vy,vz,system):
   vexpx = vx[system] + H_0 * np.abs((x[system]-xnew))
   vexpy = vy[system] + H_0 * np.abs((y[system]-ynew))
   vexpz = vz[system] + H_0 * np.abs((z[system]-znew))
   vrel = ((vexpx - vx[0])**2 + (vexpy - vy[0])**2 + (vexpz - vz[0])**2)**.5
    
   return vrel

#defines center of mass of system and returns components for new origin
def center_of_mass(rad,x,y,z,mass,system):
   totmass = np.sum(mass[0:system])
   compr = rad[0:system] * mass[0:system]
   compx = x[0:system] * mass[0:system]
   compy = y[0:system] * mass[0:system]
   compz = z[0:system] * mass[0:system]
   radnew = sum(compr) / totmass
   xnew = sum(compx) / totmass
   ynew = sum(compy) / totmass
   znew = sum(compz) / totmass

      
   return radnew, xnew, ynew, znew

def scale_mass(mass, system):
   maximum = np.log10(np.amax(mass))
   minimum = np.log10(np.amin(mass))
   ranges = maximum - minimum
   boundmass = 50**((np.log10(mass[:system]) - minimum) / ranges + 1)
   unboundmass = 50**((np.log10(mass[system:]) - minimum) / ranges + 1)

   #print "These are the bound masses:", boundmass
   return boundmass, unboundmass


def radial_velocity(x,y,z,vx,vy,vz,system):
   vxpol = vx + H_0 * (x - x[0])
   vypol = vy + H_0 * (y - y[0])
   vzpol = vz + H_0 * (z - z[0])
   numerator = vxpol*(x-x[0]) + vypol*(y-y[0]) + vzpol*(z-z[0])
   denominator = ((x-x[0])**2 + (y-y[0])**2 + (z-z[0])**2)**.5

   vradius = np.zeros(numerator.shape[0])
   vradius[1:] = numerator[1:]/denominator[1:]
#   try:
#      np.seterr(all = 'raise')
#      vradius[1:] = numerator[1:]/denominator[1:]
#
#   except FloatingPointError:
#      vradius[0] = 0
#
   #print "This is vrad:", vradius
   return vradius

########################BEGINNING OF ANALYSIS##############################

#TO DO: Figure out how to annotate the figures with their halo numbers!
def plot_velocities(halos, halofile, system):
   h = 0.68
   H_0 = h * 100
   listno = np.arange(0,system+200)
   rad = halos[1,:] / h
   mass = halos[2,:] / h
   x = halos[3,:] / h
   y = halos[4,:] / h
   z = halos[5,:] / h
   vx = halos[6,:]
   vy = halos[7,:]
   vz = halos[8,:]

   print "plot_velocities: "
   print rad[0:system]
 
   #print "What does system equal?:", system
   vrad = radial_velocity(x, y, z, vx, vy, vz, system)
  # print vrad[0:system]
   print system

   try:
      fig = plt.figure()
      ax1 = fig.add_subplot(111)
      xax = np.linspace(0, rad[system+200], 100)
      yax = H_0 * xax
      ax1.plot(xax, yax, color = 'c', linewidth = 2)
   
      bmass, ubmass = scale_mass(mass, system)
   
      target = ax1.scatter(0, vrad[0], s = bmass[0], c = 'black', alpha = 0.7)
      bound = ax1.scatter(rad[1:system], vrad[1:system], s = bmass[1:], c = 'green', alpha = 0.7)
      unbound = ax1.scatter(rad[system:system+200], vrad[system:system+200], s = ubmass[0:200], c = 'red', alpha = 0.7)
      for i in range(0,system):
         ax1.annotate('%d'%listno[i], xy = (rad[i], vrad[i]), xytext = (rad[i], vrad[i]), fontweight = 'bold')
   
      plt.xlim(xmin = 0)
      ax = plt.axes()
      ax.xaxis.grid()
      ax.yaxis.grid()
      ax1.set_xlabel('r [Mpc]')
      ax1.set_ylabel('velocity [km/s]')
      plt.legend((target, bound, unbound), ('Target Halo', 'Bound Halos', 'Unbound Halos'), scatterpoints = 1, loc = 'upper left', labelspacing = 1, fancybox=True, framealpha=0.5)
      massy = np.array_str(mass[0])
      plt.title('Target Halo Mass = ' + massy + " M_sun")
   ######################need to change to given directory!
      plt.savefig(path+halofile+'velocity.png')
      #plt.show()
      plt.close()
      #print "Plotting done!"
   except IndexError:
      print "Index out of bounds! Continuing..."


def plot_energy_ratio(Kinetic, Potential, halos, halofile, system):
   rad = halos[1,:] / h
   mass = halos[2,:] / h
   T = np.array(Kinetic)
   V = np.array(Potential)
   print "Kinetic", T
   print "Potential", V
   ratio = T / V
   #ratio[0] = 0
   print "Energy Ratio:"
   print rad[0:ratio.shape[0]]
   print ratio
   print ratio.shape[0]
   rat = np.array(np.abs(ratio))
   fig = plt.figure()
   ax1 = fig.add_subplot(111)
   ax2 = fig.add_subplot(111)
   bound, unbound = scale_mass(mass[2:], system+2) 
   #plt.clf()
   print bound
   size = ax1.scatter(rad[2:system+1], rat, s=bound, alpha=0.7)
   ratioo, = ax2.plot(rad[2:system+1], rat, 'darkmagenta')
   ax1.set_xlabel('r [Mpc]')
   ax1.set_ylabel('T/V')
   ax = plt.axes()
   ax.xaxis.grid()
   ax.yaxis.grid()

   plt.legend([ratioo, size],["T/V","Mass"], loc='upper left', scatterpoints=1)
   #plt.grid()
   listno = np.arange(2,rat.shape[0]+2)
   print rat
   print rad[2:system+1]
   print system
   for i in range(0,rat.shape[0]):
      ax1.annotate('%d'%listno[i], xy = (rad[i+2], rat[i]), xytext = (rad[i+2], rat[i]))
   #plt.xlim(0)
   #plt.ylim(0)
   massy = np.array_str(halos[2,0]/h)
   plt.title('Target Halo Mass = ' + massy + " M_sun")
   plt.savefig(path+halofile+'ratio.png')
   #plt.show()



def halo_scatter_plot(halos, halofile, system):
   x = halos[3,:] / h
   y = halos[4,:] / h
   z = halos[5,:] / h
   ub = system + 500

   plt.clf()

	#creates the plot/subplots
   fig = plt.figure(figsize = (12,12))
   ax = fig.add_subplot(111, projection = '3d')

	#extracts coordinates from bounded and unbounded numpy arrays
   ax.scatter(x[0], y[0], z[0], marker = '*', color = 'black', s=200, label = 'target halo')
   ax.scatter(x[1:system], y[1:system], z[1:system],marker = 'o', color='green', s=40, label = 'bounded halos')
   ax.scatter(x[system:ub], y[system:ub], z[system:ub],marker = 'x', color='red', s=40, label = 'unbounded halos')

	#sets labels
   ax.set_xlabel('x (Mpc)')
   ax.set_ylabel('y (Mpc)')
   ax.set_zlabel('z (Mpc)')

	#creates plot legend
   target_proxy = lines.Line2D([0],[0], linestyle = "none", c = 'black', marker = '*')
   bound_proxy = lines.Line2D([0],[0], linestyle = "none", c = 'green', marker = 'o')
   unbound_proxy = lines.Line2D([0],[0], linestyle = "none", c = 'red', marker = 'x')
   ax.legend([target_proxy, bound_proxy, unbound_proxy], ['Target Halo', 'Bound to System', 'Unbound'], numpoints = 1)

   mass = np.array_str(halos[2,0]/h)
   plt.title('Target Halo Mass = ' + mass + " M_sun")
   plt.savefig(path+halofile+'scatter.png')
   plt.close()


def halo_histogram(halos, halofile, system):
   rad = halos[1,:] / h
   try:
      bins = np.linspace(0,rad[5*system],25)
      plt.clf()
      plt.hist(rad[0:system], bins, alpha=0.5, label='bound')
      plt.hist(rad[0:5*system], bins, alpha=0.5, label='unbound')
      plt.legend(loc='upper right')
      mass = np.array_str(halos[2,0]/h)
      plt.xlabel('r (Mpc)')
      plt.ylabel('Halo Frequency')
      plt.title('Target Halo Mass = ' + mass + " M_sun")
      plt.savefig(path+halofile+'histogram.png')
      plt.close()

   except IndexError:
      print "Index out of bounds, continue..."


#############################MAIN FUNCTION##################################
#Produces energy output of every cluster step until the system is unbound
def calculate_energy(halos, halofile):
   h= 0.68
   virrad = halos[9,:]/h*10**-3
   mvir = halos[10,:]/ h
   rad = halos[1,:]/h     #Mpc
   mass = halos[2,:] / h    #M_sun
   x = halos[3,:] / h   #Mpc
   y = halos[4,:] / h   #Mpc
   z = halos[5,:] / h   #Mpc
   vx = halos[6,:]      #km/s
   vy = halos[7,:]      #km/s
   vz = halos[8,:]      #km/s
   G = 4.302 * 10**-9  #Mpc M_s^-1 (km/s)^2
   H_0 = h*100
   Potential = []
   Kinetic = []
   Energy = []
   E = 0
   system = 1
   hitvir = False
   clustradvir = 0

#######################################
   while (E <= 0):
      system += 1
      rnew, xnew, ynew, znew = center_of_mass(rad,x,y,z,mass,system)
      KE = 0
      PE = 0
       
      for i in range(0, system):
         targmass = mass[i]
         vrel = relative_velocity(xnew,ynew,znew,x,y,z,vx,vy,vz,i)
         KE += .5 * targmass * vrel**2
         "KE:", KE
         #print "targmass:", targmass
         #print "vrel:", vrel
   
         PE_comp = G * mass[i] * mass[0:i] / np.absolute((rad[i] - rad[0:i]))
         #print "Sums of PE", PE_comp
         PE -= np.sum(PE_comp)
         #print "PE:", PE
   
      E = KE + PE
      print "Potential: ", PE
      print "Kinetic: ", KE
      print "Total Energy: ", E
      print "Energy Ratio: ", np.abs(KE/PE)
      print "System: ", system
      print "Mass: ", mass[system]
      print "Radius: ", rad[system]
      print ""
   
      Kinetic.append(KE)
      Potential.append(PE)
      Energy.append(E)

      #Identifies radius at which the cluster is virialized
      if (hitvir == False and np.abs(KE/PE) <= 0.5):
         clustradvir = rad[system]
         hitvir = True
       
      if (system > 400):
         print "All halos bound"
         return -1, -1, -1, -1, -1, -1, -1
   #print "system: ", system
   number = system

   #Plots the velocity distribution and the |T/V| energy ratio as a function
   #of radius
   plot_velocities(halos, halofile, system)
   plot_energy_ratio(Kinetic, Potential, halos, halofile, system)
   #halo_scatter_plot(halos, halofile, system)
   #halo_histogram(halos,halofile,system)
   value = rad[number]
   totmass = np.sum(mass[0:number+1])
   print totmass
   halomass = mass[0]
   viradi = virrad[0]
   mvirr = mvir[0]

   return number, value, halomass, totmass, viradi, mvirr, clustradvir

#######################################################################


files = []

for root, dirnames, filenames in os.walk(path):
   for file in filenames:
      files.append(file)

print files

halofiles = [x for x in files if '.npy' in x]
trial = 0

#print halofiles
haloabundance = []
boundradius = []
halomass = []
totalmass = []
virrad = []
virmass = []
clusterrad = []

############all halos in directory
for halofile in halofiles:
   if '.png' in halofile:
      continue
   trial += 1
   print "##########################################################"
   print "file #:", trial
   print "##########################################################"
   halos = np.load(path+halofile)
   number, radius, hmass, totmass, rvir, mvirall, clustvirrad = calculate_energy(halos, halofile)
   if number == -1:
      continue
   halomass.append(hmass)
   haloabundance.append(number)
   boundradius.append(radius)
   totalmass.append(totmass)
   virrad.append(rvir)
   virmass.append(mvirall)
   clusterrad.append(clustvirrad)

#statistics = np.array([halomass,haloabundance,boundradius, totalmass, virrad, virmass, clusterrad])
#np.save("newallhalostatistics.npy", statistics)

#############single halo
#user = raw_input("choose halo\n")
#halos = np.load(path+user)
#calculate_energy(halos, user)
