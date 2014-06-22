#!/usr/bin/python

#  owner: rb

#  import statements
import sys
sys.path.append('/Users/rynebeeson/Desktop/GitHub/NTR14/src')
sys.path.append('/Users/rynebeeson/Desktop/GitHub/NTR14/src/aerospace_toolbox')


# GA picks integer parameters
#   Setup Object of Mission Class
#       Set SpaceCraft Object as Field of Mission
#       Perform Optimize Method of Mission
#       Perform ROI Method of Mission
#           Perform Revenue Method of Mission
#           Perform Cost Method of Mission
#
# ...... Continue GA until .....
#               max iterations  or
#               time expires    or
#               convergence



#  Setup GA Bounds and Parameters and Calls to Obj Functions
#  .........


#  GA Calls Mission
""" Can save results from each GA trial """
import missions
#  payload parameters
payload_mass   = 1000.0
payload_height = 2.0
payload_radius = 2.0
#  outer-loop optimizer call here
r = missions.main(1, 0, 0, payload_mass, payload_height, payload_radius)



#  plot the trajectory
from plotmodule import ArcPlot
trajplot = ArcPlot(1)
trajplot.flow2xyz(r)
trajplot.plot(color = 'b')



print ('')
print ('-------------------------------------')
print ('')



#  for reference, plot nso orbit
from spaceobjects import SpaceBody
nso = SpaceBody('LEO')
nso.set_sma(6900.0)
nso.set_i(28.5)
target = []
for i in range(90):
    nso.flow(60.0)
    target.append(nso.r)
trajplot.flow2xyz(target)
trajplot.plot(color = 'g')

print ('')
print ('-------------------------------------')
print ('')

#  for reference, plot the geo orbit
geo = SpaceBody('LEO')
geo.set_sma(42000.0)
geo.set_i(0.0)
target = []
for i in range(365):
    geo.flow(86400.0)
    target.append(geo.r)
trajplot.flow2xyz(target)
trajplot.plot(color = 'r')

trajplot.set_axis()
trajplot.show()




















#if __name__ == "__main__":
#	main()


"""
from spaceobjects import SpaceBody, SpaceCraft
from solarsystem_objects_parameters import moon
from math import pi


ntr = SpaceCraft('LEO')
ntr.add_tank()

print ntr.tank[0].name
print ntr.tank[0].mass
print ntr.mass

ntr.tank[0].add_fuel(50)

print ntr.tank[0].name
print ntr.tank[0].mass
print ntr.mass

print ntr.tanknames

print ntr.core.name

ntr.core.length = 1.0
print ntr.core.length

"""



"""
#  instantiate an object called earth and mars of the SpaceBody class
earth = SpaceBody('earth')
mars  = SpaceBody('mars')

#  simulate the earth and mars orbit for one year
r = []; r.append(earth.r)
m = []; m.append(mars.r)
for i in range(365):
    earth.flow(86400.0)
    mars.flow(86400.0)
    r.append(earth.r)
    m.append(mars.r)

#  plot the earth and mars for a one year period
from plotmodule import ArcPlot
oneplot = ArcPlot(1)
oneplot.flow2xyz(r)
oneplot.plot()
oneplot.flow2xyz(m)
oneplot.plot(color = 'r')



#  instantiate an object called spacecraft of the SpaceCraft class
#+ and give it earth's parameters
ntr = SpaceCraft()
ntr.set_parameters('earth')
ntr.name = 'ntr'
ntr.eom  = 'S2B'

ntr.transfer(mars.r, mars.v, 365 * 86400.0)

n = []; n.append(ntr.r)

ntr.v += ntr.dv1
for i in range(365):
    ntr.flow(86400.0)
    n.append(ntr.r)

oneplot.flow2xyz(n)
oneplot.plot(color = 'g')

n = []
ntr.v += ntr.dv2
for i in range(90):
    ntr.flow(86400.0)
    n.append(ntr.r)

from numpy.linalg import norm
print norm(ntr.dv1)
print norm(ntr.dv2)

oneplot.flow2xyz(n)
oneplot.plot(color = 'k')
oneplot.set_axis()
oneplot.show()
"""










