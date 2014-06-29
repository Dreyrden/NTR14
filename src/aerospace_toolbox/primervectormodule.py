#-------------------------------------------------------------------#
#---------------------- Primer Vector Module -----------------------#
#-------------------------------------------------------------------#
"""
    module: primervectormodule.py
    date:   2014_05_13
    author: ryne beeson
"""
#  -------------------------
#      import statements
#  -------------------------


#  -----------------------
#       aux functions
#  -----------------------
#  extract the STM(tf, t0)
#  and block components of STM
def extractSTM(flow, i):
    #  import statements
    from numpy import array, matrix
    
    #  construct STM(t, t0)
    stm = array(flow['y'][i])

    #  strip the stm vector out of the flow
    #  this could be for a 2d or 3d vector
    #+ and either may have had 'mu' already removed
    if len(stm) == 21:
        stm = stm[5:]
    elif len(stm) == 20:
        stm = stm[4:]
    elif len(stm) == 43:
        stm = stm[7:]
    elif len(stm) == 42:
        stm = stm[6:]

    #  reshape the stm into a matrix
    if len(stm) == 16: stm = stm.reshape(4,4)
    if len(stm) == 36: stm = stm.reshape(6,6)
    stm = matrix(stm)
    
    # STM = [[M, N],
    #        [S, T]]
    if stm.shape[0] == 4:
        mblk = stm[0:2, 0:2]
        nblk = stm[0:2, 2:]
        sblk = stm[2:4, 0:2]
        tblk = stm[2:4, 2:]
    elif stm.shape[0] == 6:
        mblk = stm[0:3, 0:3]
        nblk = stm[0:3, 3:]
        sblk = stm[3:6, 0:3]
        tblk = stm[3:6, 3:]

    #  return the stm and the M, N Blocks
    return stm, mblk, nblk


def pslope(flow, p, start, end):
    # import statements
    from numpy.linalg import norm
    from numpy import array
    
    lx   = flow['x'][start] - flow['x'][end]
    pn1  = norm(array(p[start]))
    pn0  = norm(array(p[end]))
    pdot = (pn1 - pn0)/lx
    
    return pdot
    
    '''
    # time history
    lx = len(flow['x'])
    
    # generate primer vector magnitude history
    pnorm = [norm(array(p[0]))]
    for i in range(1, lx):
        pnorm.append(norm(array(p[i])))
    
    # find initial and final slope of primer vector magnitude
    p0dot = (pnorm[1]  - pnorm[0]) /(flow['x'][1]  - flow['x'][0])
    pfdot = (pnorm[-1] - pnorm[-2])/(flow['x'][-1] - flow['x'][-2])
    
    print p0dot
    print pfdot
    
    '''


# ------------------------
#      main functions
# ------------------------
def phistory(dv1, dv2, flow):
    # import statements
    from numpy.linalg import norm, solve
    from numpy import array, matrix
    
    # initial and final primer vector magnitude
    p0 = dv1/norm(dv1);
    pf = dv2/norm(dv2);
    
    # extract STM
    phi, mblk, nblk = extractSTM(flow, -1)
    
    # solve for p0 dot (i.e. pdot at time = t0)
    p0    = matrix(p0).T
    pf    = matrix(pf).T
    rhs   = pf - mblk*p0
    p0dot = solve(nblk, rhs)
    
    # solve for p(t) history
    p = [list(array(p0).reshape(len(p0),))]

    lx = len(flow['x'])
    for i in range(1,lx-1):
        phi, mblk, nblk = extractSTM(flow, i)
        pt = mblk*p0 + nblk*p0dot
        pt = list(array(pt).reshape(len(pt),))
        p.append(pt)
    
    # append pf
    p.append(list(array(pf).reshape(len(pf),)))

    return p


def pmaghistory(p):
    #  import statements
    from numpy.linalg import norm
    #  generate the primer vector magnitude history
    #+ when give a primer vector history
    pmag = []
    for item in p:
        pmag.append(norm(item))
    return pmag



# function for minimizing the sum of the squares of the
# primer vector slopes at the initial and final times.
# should be called with an appropriate function from
# scipy.optimize()
def opt_pterminal(x0, ic1, ic2, dt):
    '''
        1.) flow #1, #2 based on x0
        2.) run prussing_conway
        3.) run shooting method w/ variational eqns
        4.) generate delta-vs
        5.) generate primer history
        6.) extract slopes
        7.) compute cost function
        '''
    '''
        FUTURE WORK
        1.) add in intelligent step sizing for flow3 = shoot()
        '''
    # import statements
    from flow import P2B as flow
    from lambert import prussing_conway
    from shootingmodule import firstorder as shoot
    from numpy.linalg import norm
    from numpy import array, zeros, eye
    
    print ('Starting an iteration of opt_pterminal!')
    print ('x0 is: %e, %e ' %(x0[0], x0[1]))
    
    # function parameters
    ''' FW 1.) '''
    step  = 1E-2
    tol   = 1E-3
    
    # extract mu
    mu = ic1[4]
    
    # flow #1, #2 based on x0[0]
    tspan = [0, x0[0]]
    ministep = abs(x0[0]/2E2)
    print ('Midpoint-A of an iteration of opt_pterminal!')
    print ('tspan is: %e, %e ' %(tspan[0], tspan[1]))
    print ('ministep is: %e ' %ministep)
    flow1 = flow(ic1, tspan, ministep, trim_output = False, \
                 print_status = True, solver = 'rk4p')
    print ('Midpoint-B of an iteration of opt_pterminal!')
    flow2 = flow(ic2, tspan, ministep, trim_output = False)
    print ('Midpoint-C of an iteration of opt_pterminal!')
    # update ic1 and ic2
    ic1 = array(flow1['y'][-1])
    ic2 = array(flow2['y'][-1])
    # update of flow time based on x0
    tspan = [0, dt - x0[0] + x0[1]]
    # flow #2 based on x0[1]
    flow2 = flow(ic2, tspan, tstep = step, trim_output = False)
    # target position
    target = array(flow2['y'][-1][0:4])
    # solve for lambert arc
    print ('Midpoint-0 of an iteration of opt_pterminal!')
    arc = prussing_conway(ic1, target, mu, TransferTime = abs(tspan[1]))
    # setup initial conditions for flow of spacecraft with
    # variational equations
    ic = zeros(21)
    ic[0:2] = ic1[0:2]
    ic[2:4] = arc['v1']
    ic[4]   = mu
    ic[5:]  = eye(4).reshape(16,)
    # flow spacecraft using lambert arc
    print ('Midpoint-1 of an iteration of opt_pterminal!')
    flow3 = flow(ic, tspan, tstep = step, eom = 'P2BP_varEqns', \
                 trim_output = False)
    print ('Midpoint-2 of an iteration of opt_pterminal!')
    # position error
    poserror = flow3['y'][-1][0:2] - flow2['y'][-1][0:2]
    # apply shooting method to reduce error
    if norm(poserror) > tol:
        flow3 = shoot(ic, target[0:2], tspan, eom = 'P2BP_varEqns', \
                      tstep = step, damping = 1/3.0, \
                      tol = tol, trim_output = False)
        poserror = flow3['y'][-1][0:2] - flow2['y'][-1][0:2]
    print ('Midpoint-3 of an iteration of opt_pterminal!')
    # delta-vs
    dv1 =  array([flow3['y'][0][2],  flow3['y'][0][3]])  - \
           array([ic1[2], ic1[3]])
    dv2 = -array([flow3['y'][-1][2], flow3['y'][-1][3]]) + \
           array([flow2['y'][-1][2], flow2['y'][-1][3]])
    # generate primer vector history
    p = phistory(dv1, dv2, flow3)
    # extract the slopes of the primer vector
    dp0 = pslope(flow3, p,  1,  0)
    dp1 = pslope(flow3, p, -1, -2)
    # scalar cost function
    J = dp0**2 + dp1**2
    
    print ('Finished an iteration of opt_pterminal!')

    return J, flow1, flow2, flow3




# damped gradient based method for simple minimization
def pterminalcoast(param, step, damp, tol, res, iLimit, \
                    flow, p, start, end, \
                    ic1, ic2, target, TU):
    """
        currently hardcoded for t0 only..... 2014_05_13
    """
    
    DEBUG0 = False
    DEBUG1 = True
    DEBUG2 = False
    
    if DEBUG1: print ('')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('        TESTING.................        ')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('')
    
    # import statements
    from flow import P2B
    from auxiliary import lambert
    from numpy import array, zeros, eye, concatenate
    from shootingmodule import firstorder as shoot
    
    # initialize
    dt   = flow['x'][-1] - flow['x'][0]
    dx   = step
    i    = 1
    # set initial step direction based on slope
    if start == 1 and end == 0:
        if res < 0: dir = -1
        if res > 0: dir =  1
    
    solver = 'rk4p'
    tstep  = 1E-3
    mu     = ic1[4]
    
    if DEBUG0: print ('mu is: %5.2f \n' %mu)
    if DEBUG0: print ic1
    if DEBUG0: print ic2
    if DEBUG0: print target

    # use absolute value of residual
    res = abs(res)

    # while-loop unitl tolerance is met
    ''' modify while-loop statement to run first time through 
        (i.e. i = 0) like in shootingmodule::firstorder ? ? ?
        or maybe not necessary if p, flow is provided.....
        '''
    while abs(res) > tol:
        
        if DEBUG0: print ('\n Debug::Iteration is: %i \n' %i)
        
        # save old
        if i == 2: temp = res
        if i > 2:  res0 = temp; temp = res
        
        # flip the direction of step if necessary
        if i > 120:
            if abs(res) > abs(res0): dir = -1*dir
            # calculate damping needed
            elif abs(res) < abs(res0):
                slope = (res - res0)/step
                dx = -res/slope
                '''
                if res > 0: sgn = -1; else: sgn = 1
                dx = sgn*(res/slope)
                if abs(dx/step) < damp: step = step/2
                '''
                if abs(dx/step) < damp: step = dx/2
    
        # ---------------
        #      flows
        # ---------------
        # update target #1 and #2 start conditions
        # based on updated param
        '''
        if param > 0:
            tspan = [0, param]
        elif param < 0:
            tspan = [0, -param]
        '''
        tspan = [0, dir*step]
        ''' DEBUG '''
        if DEBUG2: print ('(dir*step) is: %f \n' %(dir*step))


        if DEBUG0: print ('\n Debug::flow1 \n')
        # target #1 flow and ic update
        ''' could get rid of these if(len()) lines by moving outside while()
            and using trim_output = False during while-loop
            '''
        if len(ic1) == 4: ic1 = concatenate((ic1.reshape(4,1), array([[mu]])))
        ''' choose a time step that is smaller than tspan!!!!'''
        if abs(dir*step) < 1E-3:
            ministep = abs(tspan[1] - tspan[0])/10
        else:
            ministep = 1E-4
        if DEBUG2: print ('delta-time is %e ' %tspan[1])
        if DEBUG2: print ('ministep is: %e ' %ministep)
        flow1  = P2B(ic1, tspan, tstep = ministep, trim_output = False)
        ''' DEBUG '''
        if DEBUG2: print ('flow(t)')
        if DEBUG2: print flow1['x'][-1] - flow1['x'][0]
        if DEBUG2: print ('flow1(0)....')
        if DEBUG2: print flow1['y'][0][0:4]
        if DEBUG2: print ('flow1(-1)....')
        if DEBUG2: print flow1['y'][-1][0:4]
    
        ic1 = array(flow1['y'][-1])
        
        ''' DEBUG '''
        if DEBUG2: print ('ic1....')
        if DEBUG2: print ic1[0:4]
        

        if DEBUG0: print ('\n Debug::flow2 \n')
        # target #2 flow and ic update
        ''' could get rid of these if(len()) lines by moving outside while()
            and using trim_output = False during while-loop
            '''
        if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
        flow2  = P2B(ic2, tspan, tstep = ministep, trim_output = False)
        ''' DEBUG '''
        if DEBUG2:  print ('flow2(-1)....')
        if DEBUG2: print flow2['y'][-1][0:4]
        
        ic2 = array(flow2['y'][-1])
        
        ''' DEBUG '''
        if DEBUG2: print ('ic2....')
        if DEBUG2: print ic2[0:4]
        
        
        # solve for new lambert arc
#        dt   = dt + param
#        lamb = prussing_conway(ic1, target, dt, mu, iLimit = 1500, tol = 1E-6);

        ''' need to update param here because target #1 and target #2
            may not have fully stepped the delta-time'''
        # update param
        param = param + (flow1['x'][-1] - flow1['x'][0])

        # skip lambert arc and just applying shooting method
        ic3      = zeros(21)
        ic3[0:2] = ic1[0:2]
        ic3[2:4] = flow['y'][0][2:4]
        ic3[4]   = mu
        # tack identity matrix on end for STM
        ic3[5:] = eye(4).reshape(16,)
        tspan = [0.0, dt - param]
        # apply shooting method
        if DEBUG0: print ('\n Debug::ShootingMethod \n')
        flow = shoot(ic3, target[0:2], tspan, eom = 'P2BP_varEqns', \
                     tstep = 1E-2, trim_output = False, \
                     print_status = False, \
                     damping = 1/3.0)
        if DEBUG2: print ('flow(t)')
        if DEBUG2: print flow['x'][-1] - flow['x'][0]



#        if DEBUG0: print ('\n Debug::flow3 \n')
#        # target #3 flow
#        flow  = P2B(solver, tspan, tstep, ic3, TU, \
#                    eom = 'P2BP_varEqns')

        if DEBUG0: print ('\n Debug::DeltaVs \n')
        # -----------------------
        #        Delta Vs
        # -----------------------
        dv1 =  array([flow['y'][0][2],  flow['y'][0][3]]) - ic1[2:4]
        dv2 = -array([flow['y'][-1][2], flow['y'][-1][3]]) + \
               target[2:4]
        
                    
        # initial primer vector condition
#        dv1 = lamb['v1'] - array([ic1[2], ic1[3]])
        # update velocity condition at target
        # and calculate final primer vector condition
        
#        if DEBUG0: print ('\n Debug::flow2 Update \n')

#        if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
#        flow2  = P2B(solver, tspan, tstep, ic2, TU, \
#                     eom = 'P2BP')
#        dv2 = -lamb['v2'] + array([flow2['y'][-1][2], flow2['y'][-1][3]])

        if DEBUG0: print ('\n Debug::Primer History \n')

        # primer vector history
        p = phistory(dv1, dv2, flow)
    
        if DEBUG0: print ('\n Debug::Primer Slope (Residual) \n')
    
        # calculate residual
        # (i.e. primer vector slope at terminal condition)
        res = pslope(flow, p, start, end)
#        res = abs(res)

        
        if DEBUG1:
            print('Iter: %4i \t Res: %9.5f \t Param: %9.6f \t' \
                  'Dir: %2i \t Step: %.6f \t dx: %.6f' \
                  %(i, res, param, dir, step, dx))
        
        
        # if limit reached -> break
        if i >= iLimit: break
        
        # update i
        i += 1

    if DEBUG2:
        print ('PVM::2: A BLOCK PRINT ------------ :)')
        print (' ~~~~ USING IC ~~~~~')
        print ('ic1....')
        print ic1[0:4]
        print ('ic2....')
        print ic2[0:4]
        print ('ic3....')
        print ic3[0:4]
        print (' ~~~~ USING flow ~~~~~')
        print ('flow1(-1)....')
        print flow1['y'][-1][0:4]
        print ('flow2(-1)....')
        print flow2['y'][-1][0:4]
        print ('flow3(0)....')
        print flow['y'][0][0:4]
    
    ''' update ic3 for output !!! '''
    ic3 = array(flow['y'][0][0:5])
        
    # initialze output dictionary
    ''' do the following if we request trim_output = True,
        otherwise should remove'''
    if len(ic1) == 4: ic1 = concatenate((ic1.reshape(4,1), array([[mu]])))
    if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
    if i == 1: ic3 = 0
    out = {'param':     param, \
           'tol':       tol, \
           'res':       res, \
           'iter':      i - 1, \
           'start':     start, \
           'end':       end, \
           'ic1':       ic1, \
           'ic2':       ic2, \
           'ic3':       ic3, \
           'dt':        dt, \
           'x':         flow['x'], \
           'y':         flow['y'], \
           'p':         p}

    if DEBUG1: print ('')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('        FINISHED   ------   FIN         ')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('')

    return out



'''
# terminal coast adjustment
#def pcoast(ic, fc, flow, p, mu):
def pcoast(flow, p):
    pass
'''
    

# damped gradient based method for simple minimization
def pterminalcoastU(ic, flow, target, \
                    step    = 1E-2,   \
                    damping = 10,     \
                    tol     = 1E-2,   \
                    iLimit  = 25,     \
                    start   = 0,      \
                    end     = 1,      \
                    param   = 0,      \
                    TU      = 1):
    """
        - editing to make it (U)niversal..... 2014_05_16
        - could hand function a list of burns with index numbers
        related to flow
        """
    
    DEBUG0 = False
    DEBUG1 = True
    DEBUG2 = False
    
    if DEBUG1: print ('')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('        TESTING.................        ')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('')
    
    # import statements
    from flow import P2B
    from auxiliary import lambert
    from numpy import array, zeros, eye, concatenate
    from shootingmodule import firstorder as shoot
    
    # initialize
    dt   = flow['x'][-1] - flow['x'][0]
    dx   = step
    i    = 1
    # set initial step direction based on slope
    if start == 1 and end == 0:
        if res < 0: dir = -1
        if res > 0: dir =  1
    
    solver = 'rk4p'
    tstep  = 1E-3
    mu     = ic1[4]
    
    if DEBUG0: print ('mu is: %5.2f \n' %mu)
    if DEBUG0: print ic1
    if DEBUG0: print ic2
    if DEBUG0: print target
    
    # use absolute value of residual
    res = abs(res)
    
    # while-loop unitl tolerance is met
    ''' modify while-loop statement to run first time through
        (i.e. i = 0) like in shootingmodule::firstorder ? ? ?
        or maybe not necessary if p, flow is provided.....
        '''
    while abs(res) > tol:
        
        if DEBUG0: print ('\n Debug::Iteration is: %i \n' %i)
        
        # save old
        if i == 2: temp = res
        if i > 2:  res0 = temp; temp = res
        
        # flip the direction of step if necessary
        if i > 120:
            if abs(res) > abs(res0): dir = -1*dir
            # calculate damping needed
            elif abs(res) < abs(res0):
                slope = (res - res0)/step
                dx = -res/slope
                '''
                    if res > 0: sgn = -1; else: sgn = 1
                    dx = sgn*(res/slope)
                    if abs(dx/step) < damping: step = step/2
                    '''
                if abs(dx/step) < damping: step = dx/2
        
        # ---------------
        #      flows
        # ---------------
        # update target #1 and #2 start conditions
        # based on updated param
        '''
            if param > 0:
            tspan = [0, param]
            elif param < 0:
            tspan = [0, -param]
            '''
        tspan = [0, dir*step]
        ''' DEBUG '''
        if DEBUG2: print ('(dir*step) is: %f \n' %(dir*step))
        
        
        if DEBUG0: print ('\n Debug::flow1 \n')
        # target #1 flow and ic update
        ''' could get rid of these if(len()) lines by moving outside while()
            and using trim_output = False during while-loop
            '''
        if len(ic1) == 4: ic1 = concatenate((ic1.reshape(4,1), array([[mu]])))
        ''' choose a time step that is smaller than tspan!!!!'''
        if abs(dir*step) < 1E-3:
            ministep = abs(tspan[1] - tspan[0])/10
        else:
            ministep = 1E-4
        if DEBUG2: print ('delta-time is %e ' %tspan[1])
        if DEBUG2: print ('ministep is: %e ' %ministep)
        flow1  = P2B(ic1, tspan, tstep = ministep, trim_output = False)
        ''' DEBUG '''
        if DEBUG2: print ('flow(t)')
        if DEBUG2: print flow1['x'][-1] - flow1['x'][0]
        if DEBUG2: print ('flow1(0)....')
        if DEBUG2: print flow1['y'][0][0:4]
        if DEBUG2: print ('flow1(-1)....')
        if DEBUG2: print flow1['y'][-1][0:4]
        
        ic1 = array(flow1['y'][-1])
        
        ''' DEBUG '''
        if DEBUG2: print ('ic1....')
        if DEBUG2: print ic1[0:4]
        
        
        if DEBUG0: print ('\n Debug::flow2 \n')
        # target #2 flow and ic update
        ''' could get rid of these if(len()) lines by moving outside while()
            and using trim_output = False during while-loop
            '''
        if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
        flow2  = P2B(ic2, tspan, tstep = ministep, trim_output = False)
        ''' DEBUG '''
        if DEBUG2:  print ('flow2(-1)....')
        if DEBUG2: print flow2['y'][-1][0:4]
        
        ic2 = array(flow2['y'][-1])
        
        ''' DEBUG '''
        if DEBUG2: print ('ic2....')
        if DEBUG2: print ic2[0:4]
        
        
        # solve for new lambert arc
        #        dt   = dt + param
        #        lamb = prussing_conway(ic1, target, dt, mu, iLimit = 1500, tol = 1E-6);
        
        ''' need to update param here because target #1 and target #2
            may not have fully stepped the delta-time'''
        # update param
        param = param + (flow1['x'][-1] - flow1['x'][0])
        
        # skip lambert arc and just applying shooting method
        ic3      = zeros(21)
        ic3[0:2] = ic1[0:2]
        ic3[2:4] = flow['y'][0][2:4]
        ic3[4]   = mu
        # tack identity matrix on end for STM
        ic3[5:] = eye(4).reshape(16,)
        tspan = [0.0, dt - param]
        # apply shooting method
        if DEBUG0: print ('\n Debug::ShootingMethod \n')
        flow = shoot(ic3, target[0:2], tspan, eom = 'P2BP_varEqns', \
                     tstep = 1E-2, trim_output = False, \
                     print_status = False, \
                     damping = 1/3.0)
        if DEBUG2: print ('flow(t)')
        if DEBUG2: print flow['x'][-1] - flow['x'][0]



        #        if DEBUG0: print ('\n Debug::flow3 \n')
        #        # target #3 flow
        #        flow  = P2B(solver, tspan, tstep, ic3, TU, \
        #                    eom = 'P2BP_varEqns')

        if DEBUG0: print ('\n Debug::DeltaVs \n')
        # -----------------------
        #        Delta Vs
        # -----------------------
        dv1 =  array([flow['y'][0][2],  flow['y'][0][3]]) - ic1[2:4]
        dv2 = -array([flow['y'][-1][2], flow['y'][-1][3]]) + \
         target[2:4]


        # initial primer vector condition
        #        dv1 = lamb['v1'] - array([ic1[2], ic1[3]])
        # update velocity condition at target
        # and calculate final primer vector condition

        #        if DEBUG0: print ('\n Debug::flow2 Update \n')

        #        if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
        #        flow2  = P2B(solver, tspan, tstep, ic2, TU, \
        #                     eom = 'P2BP')
        #        dv2 = -lamb['v2'] + array([flow2['y'][-1][2], flow2['y'][-1][3]])

        if DEBUG0: print ('\n Debug::Primer History \n')

        # primer vector history
        p = phistory(dv1, dv2, flow)

        if DEBUG0: print ('\n Debug::Primer Slope (Residual) \n')

        # calculate residual
        # (i.e. primer vector slope at terminal condition)
        res = pslope(flow, p, start, end)
        #        res = abs(res)


        if DEBUG1:
         print('Iter: %4i \t Res: %9.5f \t Param: %9.6f \t' \
               'Dir: %2i \t Step: %.6f \t dx: %.6f' \
               %(i, res, param, dir, step, dx))


        # if limit reached -> break
        if i >= iLimit: break

        # update i
        i += 1
    
    if DEBUG2:
        print ('PVM::2: A BLOCK PRINT ------------ :)')
        print (' ~~~~ USING IC ~~~~~')
        print ('ic1....')
        print ic1[0:4]
        print ('ic2....')
        print ic2[0:4]
        print ('ic3....')
        print ic3[0:4]
        print (' ~~~~ USING flow ~~~~~')
        print ('flow1(-1)....')
        print flow1['y'][-1][0:4]
        print ('flow2(-1)....')
        print flow2['y'][-1][0:4]
        print ('flow3(0)....')
        print flow['y'][0][0:4]
    
    ''' update ic3 for output !!!'''
    ic3 = array(flow['y'][0][0:5])
    
    # initialze output dictionary
    ''' do the following if we request trim_output = True,
        otherwise should remove'''
    if len(ic1) == 4: ic1 = concatenate((ic1.reshape(4,1), array([[mu]])))
    if len(ic2) == 4: ic2 = concatenate((ic2.reshape(4,1), array([[mu]])))
    if i == 1: ic3 = 0
    out = {'param':     param, \
        'tol':       tol, \
        'res':       res, \
        'iter':      i - 1, \
        'start':     start, \
        'end':       end, \
        'ic1':       ic1, \
        'ic2':       ic2, \
        'ic3':       ic3, \
        'dt':        dt, \
        'x':         flow['x'], \
        'y':         flow['y'], \
        'p':         p}
    
    if DEBUG1: print ('')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('        FINISHED   ------   FIN         ')
    if DEBUG1: print ('----------------------------------------')
    if DEBUG1: print ('')
    
    
    return out


















#-------------------------------------------------------------------#
#---------------------- Primer Vector Module -----------------------#
#-------------------------------------------------------------------#