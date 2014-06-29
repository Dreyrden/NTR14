#--------------------------------------------------------------#
#---------------------- Shooting Module -----------------------#
#--------------------------------------------------------------#

def firstorder(ic, target, tspan, eom, \
               solver  = 'rk4', \
               tstep   =     0, \
               damping =   0.1, \
               iLimit  =    25, \
               tol     =  1E-3, \
               TU      =     1, \
               guess   =  None, \
               print_status = False, \
               trim_output  = True):
    
    #  import statements
    from numpy.linalg import norm, solve
    from numpy import zeros, eye, concatenate
    from primervectormodule import extractSTM
    
    
    #  if eom argument is not an eom model with variational equations,
    #+ then turn on the appropriate model
    if eom == 'P2BP': eom = 'P2BP_varEqns'
    if eom == 'S2BP': eom = 'S2BP_varEqns'

    
    #  determine dimension of problem (i.e. planar or spatial)
    #+ and initialize the state vector and residual vector
    if len(target) == 2:
        from flow import P2B as flow
        # initialize a state vector
        if len(ic) == 21:
            s = ic
        if len(ic) == 5:
            s = zeros(21)
            s[0:5] = ic
            #  tack identity matrix on end for STM
            s[5:] = eye(4).reshape(16,1)
        #  residual vector
        res = s[0:2]
        #  make a copy of the velocity for computing/printing delta-v
        v0 = s[2:4].copy()
    elif len(target) == 3:
        from flow import S2B as flow
        #  initialize a state vector
        if len(ic) == 43:
            s = ic
        if len(ic) == 7:
            s = zeros(43)
            s[0:7] = ic
            #  tack identity matrix on end for STM
            s[7:] = eye(6).reshape(36,1)
        #  residual vector
        res = s[0:3]
        #  make a copy of the velocity for computing/printing delta-v
        v0 = s[3:6].copy()


    #  calculate an appropriate tstep if not defined
    if tstep == 0: tstep = abs(tspan[1] - tspan[0])/2E2

    #  if an initial guess for perturbation is given then use it,
    #+ otherwise perturb the shooting method in the direction of current
    #+ velocity, so as to generally seek a lower detla-v solution
    if len(target) == 2:
        if guess != None: s[2:4] += guess
        else: s[2:4] += 0.1*s[2:4]
    elif len(target) == 3:
        if guess != None: s[3:6] += guess
        else: s[3:6] += 0.1*s[3:6]
    

    #  set counter
    i   = 0
    #  iterate with a first-order shooting method
    while ((norm(res) > tol) and (i < iLimit + 1)) or (i == 0):
        
        #  solve for 'flow'
        solution = flow(s, tspan, tstep, solver, eom, TU, \
                        trim_output)
        #  extract the state transition matrix
        #  and relevant blocks for delta-v update
        stm, mblk, nblk = extractSTM(solution, -1)
        #  solve for delta-v update
        dv = solve(nblk, res)
        #+ update residual
        if   len(dv) == 2: res = (solution['y'][-1][0:2] - target)
        elif len(dv) == 3: res = (solution['y'][-1][0:3] - target)
        #  print iterative updates if desired
        if print_status and len(dv) == 2:
            print ('Iter: %4i   Res: %8e   ' \
                   'dv1x: %8.5f   dv1y: %8.5f ' \
                   %(i, norm(res), (s[2] - v0[0]), (s[3] - v0[1])))
        if print_status and len(dv) == 3:
            print ('Iter: %4i   Res: %8e   ' \
                  'dv1x: %8.5f   dv1y: %8.5f   dv1z: %8.5f ' \
                   %(i, norm(res), (s[3] - v0[0]), \
                    (s[4] - v0[1]), (s[5] - v0[2])))
        #  update initial velocity vector
        if   len(dv) == 2: s[2:4] = s[2:4] - damping*dv
        elif len(dv) == 3: s[3:6] = s[3:6] - damping*dv

        #  update counter
        i += 1

    return solution











#--------------------------------------------------------------#
#---------------------- Shooting Module -----------------------#
#--------------------------------------------------------------#