#---------------------------------------------------------------------------#
#---------------------- Astrodynamical Models Module -----------------------#
#---------------------------------------------------------------------------#

# ---------------------
#   Planar 2-Body EOM
# ---------------------
def P2BP(t, states):
    """
        
    """
    # import statements
    from numpy import zeros, ndarray
    
    mu = states[4]

    r = (states[0]**2 + states[1]**2)**(0.5)

    dstates = zeros(states.shape[0])

    # derivatives of states
    dstates[:2] = states[2:4]

    dstates[2] = -mu*states[0]/(r**3)
    dstates[3] = -mu*states[1]/(r**3)

    return dstates

# ---------------------
#   Planar 2-Body EOM
#           +
#  State Trans. Matrix
# ---------------------
def P2BP_varEqns(t, states):
    """
        input: 
                states: 21x1 or 1x21 array
        output:
                dstates: 21x1 array
    """
    # import statements
    from numpy import zeros, matrix, ndarray, array
    
    mu = states[4]
    
    r = (states[0]**2 + states[1]**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    # derivatives of states
    dstates[:2] = states[2:4]
    
    dstates[2] = -mu*states[0]/(r**3)
    dstates[3] = -mu*states[1]/(r**3)

    # terms for variational equations
    x = states[0]
    y = states[1]

    r3 = r**(-3.0)
    r5 = r**(-5.0)

    Uxx = mu*r3 - 3*mu*x*x*r5
    Uxy = -3*mu*x*y*r5
    Uyy = mu*r3 - 3*mu*y*y*r5

    # Df Matrix (i.e. grad(EOM))
    Df = matrix([[   0,    0, 1, 0], \
                 [   0,    0, 0, 1], \
                 [-Uxx, -Uxy, 0, 0], \
                 [-Uxy, -Uyy, 0, 0]])

    # reshape phi (STM) from array into matrix
    phi = states[5:]
    phi = phi.reshape(4,4)

    # matrix math to produce phi-dot
    phidot = Df*phi

    # reshape phi-dot into a 16x1 array and add to dstates
    phidot = phidot.reshape(16,)
    dstates[5:] = phidot

    # return dstates to caller
    return dstates

# ----------------------
#   Spatial 2-Body EOM
# ----------------------
def S2BP(t, states):
    """
        
        """
    # import statements
    from numpy import zeros, ndarray
    
    mu = states[6]
    
    r = (states[0]**2 + states[1]**2 + states[2]**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    # derivatives of states
    dstates[:3] = states[3:6]
    
    dstates[3] = -mu*states[0]/(r**3)
    dstates[4] = -mu*states[1]/(r**3)
    dstates[5] = -mu*states[2]/(r**3)
    
    return dstates

# ----------------------
#   Spatial 2-Body EOM
#           +
#  State Trans. Matrix
# ----------------------
def S2BP_varEqns(t, states):
    """
        input:
        states: 43x1 or 1x43 array
        output:
        dstates: 43x1 array
        """
    # import statements
    from numpy import zeros, matrix, ndarray, array
    
    mu = states[6]
    
    r = (states[0]**2 + states[1]**2 + states[2]**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    # derivatives of states
    dstates[:3] = states[3:6]
    
    dstates[3] = -mu*states[0]/(r**3)
    dstates[4] = -mu*states[1]/(r**3)
    dstates[5] = -mu*states[2]/(r**3)
    
    # terms for variational equations
    x = states[0]
    y = states[1]
    z = states[2]
    
    r3 = r**(-3.0)
    r5 = r**(-5.0)
    
    Uxx = mu*r3 - 3*mu*x*x*r5
    Uyy = mu*r3 - 3*mu*y*y*r5
    Uzz = mu*r3 - 3*mu*z*z*r5
    Uxy = -3*mu*x*y*r5
    Uxz = -3*mu*x*z*r5
    Uyz = -3*mu*y*z*r5
    
    # Df Matrix (i.e. grad(EOM))
    Df = matrix([[   0,    0,    0, 1, 0, 0], \
                 [   0,    0,    0, 0, 1, 0], \
                 [   0,    0,    0, 0, 0, 1], \
                 [-Uxx, -Uxy, -Uxz, 0, 0, 0], \
                 [-Uxy, -Uyy, -Uyz, 0, 0, 0], \
                 [-Uxz, -Uyz, -Uzz, 0, 0, 0]])
                 
    # reshape phi (STM) from array into matrix
    phi = states[7:]
    phi = phi.reshape(6,6)

    # matrix math to produce phi-dot
    phidot = Df*phi

    # reshape phi-dot into a 16x1 array and add to dstates
    phidot = phidot.reshape(36,)
    dstates[7:] = phidot

    # return dstates to caller
    return dstates

#---------------------------------------------------------------------------#
#---------------------- Astrodynamical Models Module -----------------------#
#---------------------------------------------------------------------------#