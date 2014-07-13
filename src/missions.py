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
    #  set the ntr to have initial conditions of the
    #+ desired 'NSO' Nuclear Safe Orbit
    mission.nso_altitude = 1200.0
    mission.ntr.set_parameters('LEO', mission.nso_altitude)
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
    #  add a fuel tank to the ntr
    mission.ntr.add_tank()
    mission.ntr.tank[0].add_fuel(1000.0)
    #  reset the r, v, and m history vectors
#    mission.r_history = []
#    mission.v_history = []
    #  set the equations of motion of the ntr to 'S2B'
    #+ for shooting methods
    mission.ntr.eom = 'S2B'

    """ TEMPORARY GEO TARGETING UNTIL TRANSFER METHOD IMPROVED 
        TO HIT COEs """
    #  initialize a mission target
    from spaceobjects import SpaceCraft
    mission.target = SpaceCraft()
    mission.target.set_sma(42000.0)
    mission.target.set_i(0.0, 'Degrees')
    mission.target.set_f(180.0, 'Degrees')

    """ Inner-Loop Optimization Call Here """
    #  import scipy.optimize module
    from scipy.optimize import minimize
    from scipy.optimize import fmin_l_bfgs_b
    from numpy import array
    #  the decision_vector that scipy.optimize can modify
    #  nso_alititude:   (real)      [400, 2000]     (km)    ic: 1200
    #  core_thrust:     (real)      [5000, 25000]   (lbf)   ic: 12500
    #  perform inner-loop optimization
    """ For whatever reason, the minimize() interface is not
        working appropriately : 2014_06_22 
        , so using direct call to fmin_l_bfgs_b with 
        maxfun = 5. maxfun = 3 was sufficient for transfer 
        test problem """
#    profit = minimize(mission.optimize,
#                      ((nso_altitude),),
#                      bounds = bnds, method = 'L-BFGS-B',
#                      options = {'disp': False, 'maxiter': 6})
    bnds = ((400.0, 2000.0), )#(5000.0, 25000.0))
    initial_guess = array([mission.nso_altitude])#, \
#                           mission.core.thrust])
    opt_result = fmin_l_bfgs_b(mission.optimize,
                               initial_guess,
                               bounds = bnds,
                               approx_grad = True,
                               epsilon = 100,
                               maxfun = 5)
    #  store the decision vector and obj. func.
    mission.decision_vector = opt_result[0]
    mission.objfunc = opt_result[1]
    #  return inner-loop optimization result to
    #+ the global optimization routine
    print opt_result
    print opt_result[0]
    print opt_result[1]
#    return profit
    return mission







class Mission:
    'Mission class documentation string'

    #  initialization method
    def __init__(self):
        
        """ DEBUG """
        self.DEBUG = False
        
        if self.DEBUG: print '*** Mission Initialization ***'
        #  initialize a ntr for the mission
        from spaceobjects import NTR
        import solarsystem_objects_parameters as ss
        
        #  initializae an ntr spacecraft
        self.ntr = NTR()
        if self.DEBUG: print self.ntr.mass
        #  add a fuel tank to the ntr
        self.ntr.add_tank()
        self.ntr.tank[0].add_fuel(1000.0)
        if self.DEBUG: print self.ntr.mass
        #  set the equations of motion of the ntr to 'S2B'
        #+ for shooting methods
        self.ntr.eom = 'S2B'
        
        #  initialize global parameters
        #  set the generic global decision vector for the mission
        #+ and call the global optimizer setup
        #  the entries of global_decision_vector correspond with
        #  [0]: mission type
        #  [1]: launch vehicle
        #  [2]: payload fairing
        self.global_decision_vector = [1, 0, 0]
        self.set_global_decisions()
        
        #  initialize local parameters
        #  set the ntr to have initial conditions of a
        #+ generic 'NSO' Nuclear Safe Orbit
        self.nso_altitude = 1200.0
        self.ntr.set_sma(ss.earth['Radius'] + self.nso_altitude)
        #  set the thrust level of the ntr core
        self.ntr.core.thrust = 25000.0
        self.ntr.core.set_thrust_derived_properties()
        if self.DEBUG: print self.ntr.mass
        #  set a generic local decision vector for the mission
        self.local_decision_vector = [self.nso_altitude, self.ntr.core.thrust]
        
        #  set a generic target for the mission
        self.set_target()
        
        #  initialize a dictionary for storing relevant history information
        self.history = {}
        self.history['burnmass'] = []
    
    
    
    #  method for setting the payload
    def set_payload(self, name, mass, height, radius):
        #  set attributes for the payload
        self.payload_name   = name
        self.payload_mass   = mass
        self.payload_height = height
        self.payload_radius = radius



    def set_target(self):
        #  initialize a mission target
        from spaceobjects import SpaceCraft
        self.target = SpaceCraft()
        self.target.set_sma(42000.0)
        self.target.set_i(0.0, 'Degrees')
        self.target.set_f(180.0, 'Degrees')



    #  execute the specified mission
    def execute(self):
        #  import rocket equation for computing fuel usage
        from auxiliary    import tsiolkovsky
        from numpy.linalg import norm
        #  reset the r and v history
        self.ntr.r_history = []
        self.ntr.v_history = []
        #  select and perform mission type
        if self.mission_type == 'gto payload release, gto tank release':
            pass
        elif self.mission_type == 'gto payload release, nso tank release':
            #  calculate the first transfer phase from
            #+ NSO -> GTO
            self.ntr.transfer(self.target.r, self.target.v,
                              Minimize_Energy = True)
            #  based on the transfer calculation, compute the fuel burn
            fuel_used = self.ntr.mass - \
                        tsiolkovsky(norm(self.ntr.dv1),
                                    self.ntr.nozzle.Isp,
                                    self.ntr.mass)
            fuel_remaining = self.ntr.fuel - fuel_used
            if self.DEBUG: print ('ntr Isp: %f ' %self.ntr.nozzle.Isp)
            if self.DEBUG: print ('Fuel Available at Start: %8.5f ' %self.ntr.fuel)
            if self.DEBUG: print ('Fuel Used: %8.5f' %fuel_used)
            if self.DEBUG: print ('Fuel Remaining: %8.5f ' %fuel_remaining)
            #  store the burn history (i.e. fuel consumed at impulses)
            #+ and remove fuel from tanks
            self.history['burnmass'].append(fuel_used)
            self.ntr.tank[0].remove_fuel(fuel_used)
            if self.DEBUG: print ('CHECK -> Fuel Remaining: %f ' %self.ntr.fuel)
            
            #  flow NSO -> GTO trajectory and save state information
        elif self.mission_type == 'geo payload release, geo tank release':
            pass
        elif self.mission_type == 'geo payload release, gto tank release':
            pass
        elif self.mission_type == 'geo payload release, nso tank release':
            pass
       
       
       
    #  method for setting global optimization variables
    def set_global_decisions(self,
                             global_decision_vector = None,
                             mission_type = None,          # (integer)   [0, 4]
                             launch_vehicle = None,        # (integer)   [0, 15]
                             payload_fairing = None):      # (integer)   [0, 3]
        #  set the global decision vector if given as an argument
        if global_decision_vector != None:
            self.global_decision_vector = global_decision_vector
        #  set the global decision vector elements
        if mission_type    == None:
            mission_type    = self.global_decision_vector[0]
        if launch_vehicle  == None:
            launch_vehicle  = self.global_decision_vector[1]
        if payload_fairing == None:
            payload_fairing = self.global_decision_vector[2]
        #  convert mission type identifier to string
        if   mission_type == 0: self.mission_type = 'gto payload release, ' \
                                                    'gto tank release'
        elif mission_type == 1: self.mission_type = 'gto payload release, ' \
                                                    'nso tank release'
        elif mission_type == 2: self.mission_type = 'geo payload release, ' \
                                                    'geo tank release'
        elif mission_type == 3: self.mission_type = 'geo payload release, ' \
                                                    'gto tank release'
        elif mission_type == 4: self.mission_type = 'geo payload release, ' \
                                                    'nso tank release'
        #  convert launch vehicle integer identifier to string
        if   launch_vehicle ==  0: self.launch_vehicle = 'AtlasV501'
        elif launch_vehicle ==  1: self.launch_vehicle = 'AtlasV511'
        elif launch_vehicle ==  2: self.launch_vehicle = 'AtlasV521'
        elif launch_vehicle ==  3: self.launch_vehicle = 'AtlasV531'
        elif launch_vehicle ==  4: self.launch_vehicle = 'AtlasV541'
        elif launch_vehicle ==  5: self.launch_vehicle = 'AtlasV701'
        elif launch_vehicle ==  6: self.launch_vehicle = 'AtlasV711'
        elif launch_vehicle ==  7: self.launch_vehicle = 'AtlasV721'
        elif launch_vehicle ==  8: self.launch_vehicle = 'AtlasV731'
        elif launch_vehicle ==  9: self.launch_vehicle = 'DeltaIVM42E'
        elif launch_vehicle == 10: self.launch_vehicle = 'DeltaIVM52E'
        elif launch_vehicle == 11: self.launch_vehicle = 'DeltaIVM54E'
        elif launch_vehicle == 12: self.launch_vehicle = 'Falcon9'
        elif launch_vehicle == 13: self.launch_vehicle = 'ArianeVECA'
        elif launch_vehicle == 14: self.launch_vehicle = 'ProtonM'
        elif launch_vehicle == 15: self.launch_vehicle = 'DeltaIVH'
        #  set the payload fairing type
        #  some launch vehicles only have 1 choice for payload fairings,
        #+ others have up to 4
        self.payload_fairing = payload_fairing
        pass
    
    
    
    #  method for setting local optimization variables
    def set_local_decisions(self,
                            local_decision_vector = None,
                            nso_altitude = None,
                            core_thrust  = None):
        #  set the local decision vector if given as an argument
        if local_decision_vector != None:
            self.local_decision_vector = local_decision_vector
        #  import statement
        import solarsystem_objects_parameters as ss
        #  set the nso altitude
        if nso_altitude == None:
            self.nso_altitude = self.local_decision_vector[0]
        else:
            self.nso_altitude = nso_altitude
        #  set sma
        self.ntr.set_sma(ss.earth['Radius'] + self.nso_altitude)
        #  set the core thrust and derived parameters
        if core_thrust == None:
            self.ntr.core.thrust = self.local_decision_vector[1]
        else:
            self.ntr.core.thrust = core_thrust
        #  update relevant ntr attributes
        self.ntr.core.set_thrust_derived_properties()



    #  method call for global optimization
    def global_optimization(self, decision_vector = None):
        #  setup the global parameters if a decision vector is given
        if decision_vector != None:
            self.set_global_decisions(decision_vector)
        #  call the local optimizer
        self.local_optimization()



    #  method call for local optimization
    def local_optimization(self):
        """ Inner-Loop Optimization Call Here """
        #  import scipy.optimize module
        from scipy.optimize import minimize
        from scipy.optimize import fmin_l_bfgs_b
        from numpy import array
        #  the decision_vector that scipy.optimize can modify
        #  nso_alititude:   (real)      [400, 2000]     (km)    ic: 1200
        #  core_thrust:     (real)      [5000, 25000]   (lbf)   ic: 25000
        #  perform inner-loop optimization
        """ For whatever reason, the minimize() interface is not
            working appropriately : 2014_06_22
            , so using direct call to fmin_l_bfgs_b with
            maxfun = 5. maxfun = 3 was sufficient for transfer
            test problem """
        #    profit = minimize(mission.optimize,
        #                      ((nso_altitude),),
        #                      bounds = bnds, method = 'L-BFGS-B',
        #                      options = {'disp': False, 'maxiter': 6})
        bnds = ((400.0, 2000.0), (5000.0, 25000.0))
        initial_guess = array([self.nso_altitude, self.ntr.core.thrust])#, \
        #                           mission.core.thrust])
        opt_result = fmin_l_bfgs_b(self.optimize,
                                   initial_guess,
                                   bounds = bnds,
                                   approx_grad = True,
                                   epsilon = 100,
                                   maxfun = 5)
        #  store the decision vector and obj. func.
        self.local_decision_vector = opt_result[0]
        self.objfunc = opt_result[1]
        #  return inner-loop optimization result to
        #+ the global optimization routine
        if self.DEBUG: print opt_result
        if self.DEBUG: print opt_result[0]
        if self.DEBUG: print opt_result[1]
        #    return profit
        #    return mission


    
    #  optimize the mission
    def optimize(self, decision_vector):
        if self.DEBUG: print ('\n ---- Starting a New Iteration ---- \n')
        #  import statements
        from copy import copy, deepcopy
        import solarsystem_objects_parameters as ss
        from revenue import calculate_revenue
        from cost import calculate_cost
        from numpy import array
        from numpy.linalg import norm
        from time import sleep
    
        
        #  make a copy of the actual mission and its ntr to optimize on
        mission_copy     = deepcopy(self)
        mission_copy.set_local_decisions(decision_vector)
#        mission_copy.ntr = copy(self.ntr)
        if self.DEBUG: print ('--- Memory Locations ---')
        if self.DEBUG: print ('  mission level  ')
        if self.DEBUG: print mission_copy
        if self.DEBUG: print self
        if self.DEBUG: print ''
        if self.DEBUG: print ('  1 sub level  ')
        if self.DEBUG: print mission_copy.ntr
        if self.DEBUG: print self.ntr
        if self.DEBUG: print ''
        if self.DEBUG: print ('  1 obj sub level 1 attribute level  ')
        if self.DEBUG: print hex(id(mission_copy.ntr.mass))
        if self.DEBUG: print hex(id(self.ntr.mass))
        if self.DEBUG: print ''
        if self.DEBUG: print ('  2 obj sub levels  ')
        if self.DEBUG: print mission_copy.ntr.core
        if self.DEBUG: print self.ntr.core
        if self.DEBUG: print ''
        if self.DEBUG: print ('  2 obj ')
        if self.DEBUG: print ('-----------------------')
        if self.DEBUG: print mission_copy.ntr.mass
        #  set the decision_vector parameters
#        mission_copy.nso_altitude = decision_vector[0]
        if self.DEBUG: print ('self.nso_altitude is: %f ' %self.nso_altitude)
        if self.DEBUG: print ('copy.nso_altitude is: %f ' %mission_copy.nso_altitude)
        #  convert nso_altitude into nso_sma
#        nso_sma = ss.earth['Radius'] + decision_vector[0]
        if self.DEBUG: print ('nso_altitude: %.12f ' %decision_vector[0])
#        if self.DEBUG: print ('nso_sma:      %.12f ' %nso_sma)
#        mission_copy.ntr.set_sma(nso_sma)
        if self.DEBUG: print ('self.ntr.sma %f ' %self.ntr.sma)
        if self.DEBUG: print ('copy.ntr.sma %f ' %mission_copy.ntr.sma)
        #  set the core thrust and derived parameters
#        mission_copy.ntr.core.thrust = decision_vector[1]
#        mission_copy.ntr.core.set_thrust_derived_properties()
        if self.DEBUG: print ('mission_copy parameters')
        if self.DEBUG: print ('core thrust is: %f ' %mission_copy.ntr.core.thrust)
        if self.DEBUG: print ('core mass   is: %f ' %mission_copy.ntr.core.mass)
        if self.DEBUG: print ('ntr mass    is: %f ' %mission_copy.ntr.mass)
        if self.DEBUG: print ('self parameters')
        if self.DEBUG: print ('core thrust is: %f ' %self.ntr.core.thrust)
        if self.DEBUG: print ('core mass   is: %f ' %self.ntr.core.mass)
        if self.DEBUG: print ('ntr mass    is: %f ' %self.ntr.mass)
        
        if self.DEBUG: print 'sleep.....'
        if self.DEBUG: sleep(5.0)
        
        #  execute the mission
        mission_copy.execute()
        
        #  calculate if feasible/infeasible
        pass
        
        #  call the revenue function, which calculates the
        #+ revenue generated by the mission
        calculate_revenue(mission_copy)
        #  call the cost function, which calculates the cost
        #+ of the mission
        calculate_cost(mission_copy)
        #  add a 'profit' attribute to the mission
        mission_copy.profit = mission_copy.revenue - mission_copy.cost
        #  return the mission profit to the inner-loop solver
#        return mission_copy.profit
        if self.DEBUG: print ('Delta-V: %8.5f ' %(norm(array(mission_copy.ntr.dv1)))) #+ norm(array(self.ntr.dv2))
        
        return norm(array(mission_copy.ntr.dv1)) # + norm(array(mission_copy.ntr.dv2))



    #  method for evaluating a decision vector
    def evaluate(self,
                 global_decision_vector = None,
                 local_decision_vector  = None):
        #  import statements
        from revenue import calculate_revenue
        from cost import calculate_cost
        #  setting global decision vector
        if global_decision_vector == None: self.set_global_decisions()
        else: self.set_global_decisions(global_decision_vector)
        #  setting local decision vector
        if local_decision_vector == None: self.set_local_decisions()
        else: self.set_local_decisions(local_decision_vector)
        #  execute the mission
        self.execute()
        #  call the revenue function, which calculates the
        #+ revenue generated by the mission
        calculate_revenue(self)
        #  call the cost function, which calculates the cost
        #+ of the mission
        calculate_cost(self)
        #  add a 'profit' attribute to the mission
        self.profit = self.revenue - self.cost










