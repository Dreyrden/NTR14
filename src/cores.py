#  owner: rb

class Core:
    'Core class documentation string'

    #  initialization method
    def __init__(self, parent, thrust = 25000.0):
        
        """ DEBUG """
        self.DEBUG = False
    
        self.parent = parent
        self.thrust = thrust
        self.mass     = 0.0
        self.length   = 0.0
        self.diameter = 0.0
        self.set_thrust_derived_properties()
    
    #  a method for calculating the mass, height, and diameter of a
    #+ fiction core based on a given thrust
    def set_thrust_derived_properties(self, thrust = None):
        
        if self.DEBUG: print '<<< enter set_thrust... >>>'
        if thrust == None: thrust = self.thrust
        else: self.thrust = thrust
        #  retract current parent attributes
        if self.DEBUG: print self.parent.mass
        self.parent.mass   -= self.mass
        self.parent.length -= self.length
        if self.DEBUG: print self.parent.mass
        #  using Deason et. al. 13'
        #+ 25000lbf @ 2333 kg, 58cm length, 60cm diameter
        #+ fictional min: 1800kg, 48cm, 50cm
        mx0 = 1800.0
        mslope = (2333.0 - mx0) / 25000.0
        self.mass = mslope*thrust + mx0
        lx0 = 48.0
        lslope = (58.0 - lx0) / 25000.0
        self.length = lslope*thrust + lx0
        dx0 = 50.0
        dslope = (60.0 - dx0) / 25000.0
        self.diameter = dslope*thrust + dx0
        #  update parent attributes
        self.parent.mass   += self.mass
        self.parent.length += self.length
        if self.DEBUG: print self.parent.mass
        if self.DEBUG: print '<<< exit set_thrust... >>>'