#----------------------------------------------------------#
#---------------------- flow Module -----------------------#
#----------------------------------------------------------#

def P2B(ic, tspan, \
        tstep = 0, \
        solver = 'rk4',    \
        eom = 'P2BP', \
        TU  = 1,            \
        trim_output = True, \
        print_status = False):
    '''
       2014_05_15
    '''
    # import statements
    from time import time
    from numpy import concatenate
    import integrators
    
    # calculate an appropriate tstep if not defined
    if tstep == 0: tstep = abs(tspan[1] - tspan[0])/1E2
    
    # select EOM model
    if eom == 'P2BP':
        from eom import P2BP as EOM
    elif eom == 'P2BP_varEqns':
        from eom import P2BP_varEqns as EOM
    else:
        print ('>> flow::P2B: EOM Model (%s) Not Available \n' %eom)
        return -1

    # start clock
    clockstart = time()

    # integrate states
    if solver == 'rk4':
        sol = integrators.rk4(EOM, tspan, tstep, ic)
    elif solver == 'rk4p':
        sol = integrators.rk4p(EOM, tspan, tstep, ic)
    else:
        print ('>> flow::P2B: Solver (%s) Not Available \n' %solver)
        return -1

    # print status if indicated
    if print_status:
        print ('flow::P2B: Number of %s steps taken: \t %-i' \
               %(solver, len(sol['x'])))
        print ('flow::P2B: Time step length was: \t \t %-5.1f (s)' \
               %(tstep*TU))
        print ('flow::P2B: %s simulation time was: \t %-5.1f (s)' \
               %(solver, time() - clockstart))
        print ('')

    # remove unwanted mu entries from sol['y']
    if eom == 'P2BP' and trim_output:
        for i in range(len(sol['x'])):
            sol['y'][i] = sol['y'][i][:-1]
    elif eom == 'P2BP_varEqns' and trim_output:
        for i in range(len(sol['x'])):
            sol['y'][i] = concatenate([sol['y'][i][:4], sol['y'][i][5:]])

    # return output as a dictionary with lists
    return sol


def S2B(ic, tspan, \
        tstep = 0, \
        solver = 'rk4',    \
        eom = 'S2BP', \
        TU  = 1,            \
        trim_output = True, \
        print_status = False):
    '''
        2014_06_09
        '''
    # import statements
    from time import time
    from numpy import concatenate
    import integrators
    
    # calculate an appropriate tstep if not defined
    if tstep == 0: tstep = abs(tspan[1] - tspan[0])/1E2
    
    # select EOM model
    if eom == 'S2BP':
        from eom import S2BP as EOM
    elif eom == 'S2BP_varEqns':
        from eom import S2BP_varEqns as EOM
    else:
        print ('>> Flow::P2B: EOM Model (%s) Not Available \n' %eom)
        return -1
    
    # start clock
    clockstart = time()
    
    # integrate states
    if solver == 'rk4':
        sol = integrators.rk4(EOM, tspan, tstep, ic)
    elif solver == 'rk4p':
        sol = integrators.rk4p(EOM, tspan, tstep, ic)
    else:
        print ('>> Flow::S2B: Solver (%s) Not Available \n' %solver)
        return -1
    
    # print status if indicated
    if print_status:
        print ('Flow::S2B: Number of %s steps taken: \t %-i' \
                %(solver, len(sol['x'])))
        print ('Flow::S2B: Time step length was: \t \t %-5.1f (s)' \
                %(tstep*TU))
        print ('Flow::S2B: %s simulation time was: \t %-5.1f (s)' \
                %(solver, time() - clockstart))
        print ('')
    
    # remove unwanted mu entries from sol['y']
    if eom == 'S2BP' and trim_output:
        for i in range(len(sol['x'])):
            sol['y'][i] = sol['y'][i][:-1]
    elif eom == 'S2BP_varEqns' and trim_output:
        for i in range(len(sol['x'])):
            sol['y'][i] = concatenate([sol['y'][i][:6], sol['y'][i][7:]])
    
    # return output as a dictionary with lists
    return sol

#----------------------------------------------------------#
#---------------------- flow Module -----------------------#
#----------------------------------------------------------#