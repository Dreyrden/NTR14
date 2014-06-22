#  owner: rb

def main(mission_type, launch_vehicle, payload_fairing,
         payload_mass, payload_height, payload_radius):
    #  unpack *parameters
    #  mission_type:    (integer)   [0, 4]
    #  launch_vehicle:  (integer)   [0, 15]
    #  payload_fairing  (integer)   [0, 3]
#    mission_type, launch_vehicle, payload_fairing = parameters
    #  convert mission type identifier to string
    if   mission_type == 0: mission_type = 'gto payload release, ' \
                                           'gto tank release'
    elif mission_type == 1: mission_type = 'gto payload release, ' \
                                           'nso tank release'
    elif mission_type == 2: mission_type = 'geo payload release, ' \
                                           'geo tank release'
    elif mission_type == 3: mission_type = 'geo payload release, ' \
                                           'gto tank release'
    elif mission_type == 4: mission_type = 'geo payload release, ' \
                                           'nso tank release'
    #  convert launch vehicle integer identifier to string
    if   launch_vehicle ==  0: launch_vehicle = 'AtlasV501'
    elif launch_vehicle ==  1: launch_vehicle = 'AtlasV511'
    elif launch_vehicle ==  2: launch_vehicle = 'AtlasV521'
    elif launch_vehicle ==  3: launch_vehicle = 'AtlasV531'
    elif launch_vehicle ==  4: launch_vehicle = 'AtlasV541'
    elif launch_vehicle ==  5: launch_vehicle = 'AtlasV701'
    elif launch_vehicle ==  6: launch_vehicle = 'AtlasV711'
    elif launch_vehicle ==  7: launch_vehicle = 'AtlasV721'
    elif launch_vehicle ==  8: launch_vehicle = 'AtlasV731'
    elif launch_vehicle ==  9: launch_vehicle = 'DeltaIVM42E'
    elif launch_vehicle == 10: launch_vehicle = 'DeltaIVM52E'
    elif launch_vehicle == 11: launch_vehicle = 'DeltaIVM54E'
    elif launch_vehicle == 12: launch_vehicle = 'Falcon9'
    elif launch_vehicle == 13: launch_vehicle = 'ArianeVECA'
    elif launch_vehicle == 14: launch_vehicle = 'ProtonM'
    elif launch_vehicle == 15: launch_vehicle = 'DeltaIVH'

    """ SETUP GA or BRUTE FORCE Parameters Here """
    #  initialize a mission of Mission class
    mission = Mission()
    #  set the spacecraft to have initial conditions of the
    #+ desired 'NSO' Nuclear Safe Orbit
    mission.nso_altitude = 1200.0
    mission.spacecraft.set_parameters('LEO', mission.nso_altitude)
    #  set the mission type
    mission.mission_type = mission_type
    #  set the launch vehicle
    mission.launch_vehicle = launch_vehicle
    #  set the payload fairing type
    #  some launch vehicles only have 1 choice for payload fairings,
    #+ others have up to 4
    mission.payload_fairing = payload_fairing
    #  set attributes for the payload
    mission.payload_mass   = payload_mass
    mission.payload_height = payload_height
    mission.payload_radius = payload_radius
    #  add a fuel tank to the spacecraft
    mission.spacecraft.add_tank()
    mission.spacecraft.tank[0].add_fuel(1000.0)
    #  reset the r, v, and m history vectors
    mission.r_history = []
    mission.v_history = []
    mission.m_history[0] = mission.spacecraft.mass
    #  set the equations of motion of the spacecraft to 'S2B'
    #+ for shooting methods
    mission.spacecraft.eom = 'S2B'

    """ TEMPORARY GEO TARGETING UNTIL TRANSFER METHOD IMPROVED 
        TO HIT COEs """
    #  initialize a mission target
    from spaceobjects import SpaceCraft
    mission.target = SpaceCraft()
    mission.target.set_sma(42000.0)
    mission.target.set_i(0.0)
    from math import pi
    mission.target.set_f(pi)

    """ Inner-Loop Optimization Call Here """
    #  import scipy.optimize module
    from scipy.optimize import minimize
    from scipy.optimize import fmin_l_bfgs_b
    from numpy import array
    #  parameters that scipy.optimize can modify
    #  nso_alititude:   (real)      [400, 2000] (km)    ic: 1200.0
    #
    #  perform inner-loop optimization
    """ For whatever reason, the minimize() interface is not
        working appropriately : 2014_06_22 
        , so using direct call to fmin_l_bfgs_b with 
        maxfun = 5. maxfun = 3 was sufficient for transfer 
        test problem """
#    profit = minimize(mission.optimize,
#                      ((nso_altitude), ),
#                      bounds = bnds, method = 'L-BFGS-B',
#                      options = {'disp': False, 'maxiter': 6})
    bnds = ((400.0, 2000.0), )
    initial_guess = ((mission.nso_altitude),)
    opt_result = fmin_l_bfgs_b(mission.optimize,
                               initial_guess,
                               bounds = bnds,
                               approx_grad = True,
                               epsilon = 100,
                               maxfun = 5)

    #  return inner-loop optimization result to
    #+ the global optimization routine
    print opt_result
    print opt_result[0]
    print opt_result[1]
#    return profit
    return mission.r_history







class Mission:
    'Mission class documentation string'

    #  initialization method
    def __init__(self):
        #  initialize a spacecraft for the mission
        from spaceobjects import SpaceCraft
        self.spacecraft = SpaceCraft()
        #  set the mission type
        self.mission_type = 'gto payload release, nso tank release'
        #  initialize r, v, and m history vectors to store
        #+ pertinent data about spacecraft history
        self.r_history = [self.spacecraft.r]
        self.v_history = [self.spacecraft.v]
        self.m_history = [self.spacecraft.mass]

    #  execute the specified mission
    def execute(self):
        #  import rocket equation for computing fuel usage
        from auxiliary    import tsiolkovsky
        from numpy.linalg import norm
        #  select and perform mission type
        if self.mission_type == 'gto payload release, gto tank release':
            pass
        elif self.mission_type == 'gto payload release, nso tank release':
            #  reset the r and v history
            self.r_history = []
            self.v_history = []
            #  calculate the first transfer phase from
            #+ NSO -> GTO
            self.spacecraft.transfer(self.target.r, self.target.v,
                                     0.1*86400.0, self.r_history,
                                     self.v_history,
                                     Minimize_Energy = True)
            #  based on the transfer calculation, compute the fuel burn
            fuel_used = self.spacecraft.mass - \
                        tsiolkovsky(norm(self.spacecraft.dv1),
                                    self.spacecraft.nozzle.Isp,
                                    self.spacecraft.mass)
            fuel_remaining = self.spacecraft.fuel - fuel_used
            print self.spacecraft.nozzle.Isp
            print ('Fuel Available: %8.5f ' %self.spacecraft.fuel)
            print ('Fuel Used: %8.5f' %fuel_used)
            print ('Fuel Remaining: %8.5f ' %fuel_remaining)
            #  flow NSO -> GTO trajectory and save state information
        elif self.mission_type == 'geo payload release, geo tank release':
            pass
        elif self.mission_type == 'geo payload release, gto tank release':
            pass
        elif self.mission_type == 'geo payload release, nso tank release':
            pass
                

    
            
    #  optimize the mission
    def optimize(self, nso_altitude):
        print ('\n ---- Starting a New Iteration ---- \n')
        #  set new mission.nso_altitude
        self.nso_altitude = nso_altitude[0]
        #  convert nso_altitude into nso_sma
        import solarsystem_objects_parameters as ss
        nso_sma = ss.earth['Radius'] + nso_altitude[0]
        print ('nso_altitude: %.12f ' %nso_altitude[0])
        print ('nso_sma:      %.12f ' %nso_sma)
        self.spacecraft.set_sma(nso_sma)
        #  execute the mission
        self.execute()
        #  calculate if feasible/infeasible
        pass
        #  call the revenue function, which calculates the
        #+ revenue generated by the mission
        from revenue import calculate_revenue
        calculate_revenue(self)
#        print self.revenue
        #  call the cost function, which calculates the cost
        #+ of the mission
        from cost import calculate_cost
        calculate_cost(self)
#        print self.cost
        #  add a 'profit' attribute to the mission
        self.profit = self.revenue - self.cost
        #  return the mission profit to the inner-loop solver
#        return self.profit
#        print ('\nprinting some very important information\n')
#        print self.spacecraft.dv1
#        print self.spacecraft.dv2
        from numpy import array
        from numpy.linalg import norm
        print ('Delta-V: %8.5f ' %(norm(array(self.spacecraft.dv1)))) #+ norm(array(self.spacecraft.dv2))
        return norm(array(self.spacecraft.dv1)) # + norm(array(self.spacecraft.dv2))









