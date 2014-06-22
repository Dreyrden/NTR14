#
# Title:    solarsystem_objects_parameters
# Type:     Script
#
# Purpose:  Provides some basic parameters of solar system objects
#
# Inputs:   N/A
# Outputs:  N/A
#
# Original Date:    2014_04_31
# Original Author:  Ryne Beeson
#
# Versions / Comments
# Version: 0.01 / Original Code Created
#

# Units
units_absolute = {'Mass'  : 'kg',
                  'Length': 'km',
                  'Time'  : 's',
                  'Angle' : 'rad'}

units_derived = {'Force : kN',
                 'Speed : km/s'}

# Gravitational Constant
constant = {'G': 6.67384*1E-20} #(km^3/kg/s^2)

# List of Bodies Available
list = ['sun', 'earth', 'moon', 'mars']

# Sun
sun = {'Mass'  : 1.989*1E30,
       'Mu'    : 1.327*1E11,
       'Radius': 6.9599*1E5}

# Earth
earth = {'Mass'     : 5.97219*1E24,
         'g'        : 9.8066*1E-3,
         'Radius'   : 6.37812*1E3,
         'SemiMajor': 1.495978*1E8,
         'Eccentricity': 0.0167,
         'Inclination': 0.0}
earth['Mu'] = constant['G']*earth['Mass']

# Moon
moon = {'Mass'          : 7.3477*1E22,
        'Radius'        : 1.738*1E3,
        'SemiMajor'     : 3.844*1E5,
        'Eccentricity'  : 0.0549,
        'Inclination'   : 0.0898844564}
moon['Mu'] = constant['G']*moon['Mass']

# Earth-Moon Libration Points
earth_moon = {'L1': [ 3.217043538203148E+5,  0                   , 0],
              'L2': [ 4.442487738309655E+5,  0                   , 0],
              'L3': [-3.863465739338701E+5,  0                   , 0],
              'L4': [ 1.875281319000864E+5,  3.329001652147382E+5, 0],
              'L5': [ 1.875281319000864E+5, -3.329001652147382E+5, 0]}

# Mars
mars = {'Mass'          : 0.1074*earth['Mass'],
        'Radius'        : 0.532*earth['Radius'],
        'SemiMajor'     : 1.5237*earth['SemiMajor'],
        'Eccentricity'  : 0.0934,
        'Inclination'   : 0.032288591}
mars['Mu'] = constant['G']*mars['Mass']

















#