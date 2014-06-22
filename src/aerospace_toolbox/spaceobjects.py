#  the SpaceBody class
#+ a generic space object that has the ability to be propagated using
#+ a kepler solver
class SpaceBody:
    'SpaceBody class documentation string'
    
    """
        
    FUTURE IMPROVEMENTS:
    2014_06_19: default ephemeris in set_parameters() are arbitrary
        
    """
    
    #  initialization method
    def __init__(self, body = 'earth'):
        #  call the set_parameters method
        #+ to initialize attributes to the specified body
        self.set_parameters(body)
        #  define the dynamical model used for propagating
        #+ the body
        self.eom = 'keplerian'
    
    def set_parameters(self, body, altitude = 400.0):
        import solarsystem_objects_parameters as ss
        #  setup parameters for the requested body if available
        #  first see if the user is requesting the 'Earth' or
        #+ if the body is not on the given list, then default
        #+ to the 'Earth'
        if body == 'earth' or body not in ss.list:
            self.name = 'earth'
            #  default classical orbital elements is
            #+ earth about sun
            self.sma    = ss.earth['SemiMajor']
            self.e      = ss.earth['Eccentricity']
            self.i      = ss.earth['Inclination']
            self.RAAN   = 0.0
            self.AOP    = 0.0
            self.f      = 0.0
            #  epoch (Julian Date), default January 3rd 2014
            #+ which is roughly when the Earth is at periapsis
            self.epoch = 2456660.500000 * 86400.0
            #  mass and mu of the body
            self.mass = ss.earth['Mass']
            self.mu   = ss.earth['Mu']
            #  radius of body
            self.radius = ss.earth['Radius']
            #  set the gravitational constant of 'self'
            #+ central body
            self.Mu = ss.sun['Mu']
        if body == 'sun':
            self.name = 'sun'
            #  default classical orbital elements is
            #+ earth about sun
            self.sma    = 0.0
            self.e      = 0.0
            self.i      = 0.0
            self.RAAN   = 0.0
            self.AOP    = 0.0
            self.f      = 0.0
            #  epoch (Julian Date), default January 3rd 2014
            #+ which is roughly when the Earth is at periapsis
            self.epoch = 2456660.500000 * 86400.0
            #  mass and mu of the body
            self.mass = ss.sun['Mass']
            self.mu   = ss.sun['Mu']
            #  radius of body
            self.radius = ss.sun['Radius']
            #  set the gravitational constant of 'self'
            #+ central body
            self.Mu = ss.sun['Mu']
        if body == 'moon':
            self.name = 'moon'
            #  default classical orbital elements is
            #+ earth about sun
            self.sma    = ss.moon['SemiMajor']
            self.e      = ss.moon['Eccentricity']
            self.i      = ss.moon['Inclination']
            self.RAAN   = 0.0
            self.AOP    = 0.0
            self.f      = 0.0
            #  epoch (Julian Date), default January 3rd 2014
            #+ which is roughly when the Earth is at periapsis
            self.epoch = 2456660.500000 * 86400.0
            #  mass and mu of the body
            self.mass = ss.moon['Mass']
            self.mu   = ss.moon['Mu']
            #  radius of body
            self.radius = ss.moon['Radius']
            #  set the gravitational constant of 'self'
            #+ central body
            self.Mu = ss.earth['Mu']
        if body == 'mars':
            self.name = 'mars'
            #  default classical orbital elements is
            #+ earth about sun
            self.sma    = ss.mars['SemiMajor']
            self.e      = ss.mars['Eccentricity']
            self.i      = ss.mars['Inclination']
            self.RAAN   = 0.0
            self.AOP    = 0.0
            self.f      = 0.0
            #  epoch (Julian Date), default January 3rd 2014
            #+ which is roughly when the Earth is at periapsis
            self.epoch = 2456660.500000 * 86400.0
            #  mass and mu of the body
            self.mass = ss.mars['Mass']
            self.mu   = ss.mars['Mu']
            #  radius of body
            self.radius = ss.mars['Radius']
            #  set the gravitational constant of 'self'
            #+ central body
            self.Mu = ss.sun['Mu']
        if body == 'LEO':
            self.name = 'leo_spacecraft'
            #  default classical orbital elements is
            #+ spacecraft about earth
            from math import pi
            self.sma    = ss.earth['Radius'] + altitude
            self.e      = 0.0
            self.i      = 28.5 * pi / 180.0
            self.RAAN   = 0.0
            self.AOP    = 0.0
            self.f      = 0.0
            #  epoch (Julian Date), default January 3rd 2014
            #+ which is roughly when the Earth is at periapsis
            self.epoch = 2456660.500000 * 86400.0
            #  mass and mu of the body
            self.mass = 1000.0
            self.mu   = 0.0
            #  set the gravitational constant of 'self'
            #+ central body
            self.Mu = ss.earth['Mu']
        #  place the classical orbital elements into an array
        self.coe = [self.sma, self.e, self.i, self.RAAN, \
                    self.AOP, self.f]
        #  calculate the r, v vectors for 'self'
        #+ based on the classical orbital elements 'coe'
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
        #  calculate additional useful orbital parameters
        #+ E: eccentricy anomaly
        #+ M: mean anomaly
        #+ T: orbital period
        #+ n: mean motion
        from auxiliary import f2E, E2M, orbitalperiod
        from math import pi
        self.E = f2E(self.e, self.f)
        self.M = E2M(self.e, self.E)
        self.T = orbitalperiod(self.Mu, self.sma)
        self.n = 2.0*pi / self.T
    
    #  method to reset sma
    def set_sma(self, sma):
        self.sma = sma
        #  reset sma in coe list
        self.coe[0] = sma
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
        #  reset the orbital period
        from auxiliary import orbitalperiod
        self.T = orbitalperiod(self.Mu, self.sma)
        #  reset the mean motion
        from math import pi
        self.n = 2.0*pi / self.T
    
    #  method to reset eccentricity
    def set_e(self, e):
        self.e = e
        #  reset e in coe list
        self.coe[1] = e
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
        #  reset
        #+ E: eccentricy anomaly
        #+ M: mean anomaly
        from auxiliary import f2E, E2M
        self.E = f2E(self.e, self.f)
        self.M = E2M(self.e, self.E)

    #  method to reset inclination
    def set_i(self, i):
        self.i = i
        #  reset i in coe list
        self.coe[2] = i
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
    
    #  method to reset RAAN
    def set_RAAN(self, RAAN):
        self.RAAN = RAAN
        #  reset RAAN in coe list
        self.coe[3] = RAAN
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
    
    #  method to reset AOP
    def set_AOP(self, AOP):
        self.AOP = AOP
        #  reset AOP in coe list
        self.coe[4] = AOP
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
    
    #  method to reset true anomaly
    def set_f(self, f):
        self.f = f
        #  reset f in coe list
        self.coe[5] = f
        #  reset the r and v vectors
        from auxiliary import coe2rv
        self.r, self.v = coe2rv(self.coe, self.Mu)
        #  reset
        #+ E: eccentricy anomaly
        #+ M: mean anomaly
        from auxiliary import f2E, E2M
        self.E = f2E(self.e, self.f)
        self.M = E2M(self.e, self.E)

    #  propagation method of body
    def flow(self, delta_time):
        if self.eom == 'keplerian':
            from kepler import kepler
            from auxiliary import E2f, E2M, coe2rv
            #  create a list of the necessary orbital parameters
            #+ for ephemeris() and order them in the default order
            ephemeris_list = [self.e, self.n, self.E]
            ephemeris_list = kepler(ephemeris_list, delta_time)
            #  update other relevant orbital parameters
            self.epoch += delta_time
            self.E = ephemeris_list[2]
            self.f = E2f(self.e, self.E)
            self.M = E2M(self.e, self.E)
            #  update coe list
            self.coe[5] = self.f
            #  update the r, v vectors
            self.r, self.v = coe2rv(self.coe, self.Mu)
        elif self.eom == 'P2B':
            from utilities import rv2ic
            from flow import P2B
            ic = rv2ic(self.r[0:2], self.v[0:2], self.Mu, STM = False)
            flow = P2B(ic, [0, delta_time])
            self.epoch += delta_time
            self.r[0:2] = flow['y'][-1][0:2]
            self.v[0:2] = flow['y'][-1][2:4]
        elif self.eom == 'P2B_varEqns':
            from utilities import rv2ic
            from flow import P2B
            ic = rv2ic(self.r[0:2], self.v[0:2], self.Mu, STM = True)
            flow = P2B(ic, [0, delta_time], eom = 'P2BP_varEqns')
            self.epoch += delta_time
            self.r[0:2] = flow['y'][-1][0:2]
            self.v[0:2] = flow['y'][-1][2:4]
            self.STM = flow['y'][-1][5:]
        elif self.eom == 'S2B':
            from utilities import rv2ic
            from flow import S2B
            ic = rv2ic(self.r, self.v, self.Mu, STM = False)
            flow = S2B(ic, [0, delta_time])
            self.epoch += delta_time
            self.r = flow['y'][-1][0:3]
            self.v = flow['y'][-1][3:6]
        elif self.eom == 'S2B_varEqns':
            from utilities import rv2ic
            from flow import S2B
            ic = rv2ic(self.r, self.v, self.Mu, STM = True)
            flow = S2B(ic, [0, delta_time], eom = 'S2BP_varEqns')
            self.epoch += delta_time
            self.r = flow['y'][-1][0:3]
            self.v = flow['y'][-1][3:6]
            self.STM = flow['y'][-1][6:]




#  the SpaceCraft class
#+ a generic spacecraft object that has the ability to be track the
#+ variational equations when propagating, so that the object can
#+ calculate maneuvers
#  the spacecraft class also adds thruster parameters to its fields
#+ (also known as attributes in Python-land)
class SpaceCraft(SpaceBody):
    'SpaceCraft class documentation string'

    #  initialization method
    def __init__(self, orbit = 'LEO'):
        #  call the set_parameters method
        #+ to initialize attributes to the specified body
        self.set_parameters(orbit)
        #  define the dynamical model used for propagating
        #+ the body
        self.eom = 'keplerian'
        #  initialize an empty array for storing tank objects
        #+ and tank names
        self.tank = []
        self.tanknames = []
        self.fuel = 0.0
        #  initialize a nuclear thermal core
        from cores import LEUCore as Core
        self.core = Core()
        #  initialize a nozzle
        from nozzles import Nozzle
        self.nozzle = Nozzle(self)
        #  initialize a structure
        #
        #  initialize a .....
    
    #  targeting method for transfers
    def transfer(self, target_position, target_velocity, flight_time, \
                 r_history = None, v_history = None, Minimize_Energy = False):
        #  giving lists to r_history and v_history arguments will
        #+ populate them with the r and v state histories
        #  Minimize_Energy = True will override 'flight_time' with the
        #  minimum energy time of flight given from lambert
        from utilities import rv2ic
        from shootingmodule import firstorder as shoot
        #  create the initial conditions vector
        #  select an appropriate eom model that includes the
        #+ variational equations
        print ('Setup Initial Conditions for Lambert ....')
        if len(target_position) == 2:
            ic = rv2ic(self.r[0:2], self.v[0:2], self.Mu, STM = True)
            temp_eom = 'P2BP_varEqns'
        elif len(target_position) == 3:
            ic = rv2ic(self.r, self.v, self.Mu, STM = True)
            temp_eom = 'S2BP_varEqns'
        #  Override 'flight_time' if Minimize_Energy = True
        if Minimize_Energy:
            from utilities import create_one_list
            from auxiliary import lambert
            ic0 = create_one_list([self.r, self.v], 0, 1)
            icf = create_one_list([target_position, target_velocity], 0, 1)
            #  find the minimum time for a lambert solution
            lambert_solution = lambert(ic0, icf, self.Mu, \
                                       FindMinEnergy = True, \
                                       NonDimUnits = False, \
                                       ScaleOutput = False)
            time_for_min_energy = lambert_solution['tm']
            #  use the minimum time as a flight time to compute
            #+ the lambert solution
            lambert_solution = lambert(ic0, icf, self.Mu, \
                                       TransferTime = time_for_min_energy, \
                                       NonDimUnits = False, \
                                       ScaleOutput = False)
            from numpy import array
            dv1_guess = list(lambert_solution['v1'] - array(ic0[2:4]))
            dv1_guess.append(0.0)
            flight_time = lambert_solution['time']
        #  apply the first order shooting method
        trajectory = shoot(ic, target_position, [0, flight_time], temp_eom, \
                           tol = 1E-1, iLimit = 30, \
                           damping = 0.3, guess = dv1_guess, \
                           print_status = True)
        #  store the delta-v impulsive maneuvers needed
        #+ for the given transfer
        if len(target_position) == 2:
            self.dv1 = trajectory['y'][0][2:4] - self.v
            self.dv2 = target_velocity - trajectory['y'][-1][2:4]
        elif len(target_position) == 3:
            self.dv1 = trajectory['y'][0][3:6] - self.v
            self.dv2 = target_velocity - trajectory['y'][-1][3:6]
        #  save trajectory states if desired
        from utilities import extract_elements
        if r_history != None:
            extract_elements(trajectory['y'], 0, 2, r_history)
        if v_history != None:
            extract_elements(trajectory['y'], 3, 5, v_history)
            
    #  set the tank model from the tank module
    def add_tank(self, tank = 'LH2_Tank'):
        if tank == 'LH2_Tank':
            from tanks import LH2Tank as Tank
            newtank = Tank(self)
            self.tank.append(newtank)
            self.tanknames.append(newtank.name)




















