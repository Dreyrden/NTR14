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

    # reshape phi-dot into a 36x1 array and add to dstates
    phidot = phidot.reshape(36,)
    dstates[7:] = phidot

    # return dstates to caller
    return dstates

# -----------------------------------------
#   Planar Circular Restricted 3-Body EOM
# -----------------------------------------
def PCR3BP(t, states):
    """PCR3BP(t, states)
        
       INPUT
       ----------
       t        : 1x1 scalar
       states   : 8x1 list or array
                : [0]: x
                : [1]: y
                : [2]: vx
                : [3]: vy
                : [4]: mu1
                : [5]: mu2
                : [6]: x1 'x-coordinate of primary'
                : [7]: x2 'x-coordinate of secondary'
                
        ASSUMPTIONS: that primary and secondary are located on the x-axis
       
    """

    #  import statements
    from numpy import zeros, ndarray
    
    mu1 = states[4]
    mu2 = states[5]
    x1  = states[6]
    x2  = states[7]
    
    dx13 = states[0] - x1
    dx23 = states[0] - x2
    dy   = states[2]
    
    #  distance from primary
    r13 = (dx13**2 + dy**2)**(0.5)
    #  distance from secondary
    r23 = (dx23**2 + dy**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    #  derivatives of states
    dstates[:2] = states[2:4]
    
    dstates[2] =  2*states[3] + states[0] - mu1*dx13/(r13**3) - mu2*dx23/(r23**3)
    dstates[3] = -2*states[2] + states[1] - mu1*dy/(r13**3)   - mu2*dy/(r23**3)
    
    return dstates

# -----------------------------------------
#   Planar Circular Restricted 3-Body EOM
#                    +
#           State Trans. Matrix
# -----------------------------------------
def PCR3BP_varEqns(t, states):
    """PCR3BP_varEqns(t, states)
        
        INPUT
        ----------
        t        : 1x1 scalar
        states   : 24x1 list or array
                 : [0]: x
                 : [1]: y
                 : [2]: vx
                 : [3]: vy
                 : [4]: mu1
                 : [5]: mu2
                 : [6]: x1 'x-coordinate of primary'
                 : [7]: x2 'x-coordinate of secondary'
                 : [8-23]: Phi
                 
        ASSUMPTIONS: that primary and secondary are located on the x-axis
                 
        """
    
    #  import statements
    from numpy import zeros, ndarray
    
    mu1 = states[4]
    mu2 = states[5]
    x1  = states[6]
    x2  = states[7]
    
    dx13 = states[0] - x1
    dx23 = states[0] - x2
    dy   = states[2]
    
    #  distance from primary
    r13 = (dx13**2 + dy**2)**(0.5)
    #  distance from secondary
    r23 = (dx23**2 + dy**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    #  derivatives of states
    dstates[:2] = states[2:4]
    
    dstates[2] =  2*states[3] + states[0] - mu1*dx13/(r13**3) - mu2*dx23/(r23**3)
    dstates[3] = -2*states[2] + states[1] - mu1*dy/(r13**3)   - mu2*dy/(r23**3)
    
    # terms for variational equations
    x = states[0]
    y = states[1]
    
    zeta1 = (x - x1)
    zeta2 = (x - x2)
    
    tau1 = zeta1**2 + y**2
    tau2 = zeta2**2 + y**2
    
    psi1 = tau1**(-3.0/2.0)
    psi2 = tau2**(-3.0/2.0)
    
    phi1 = -3.0*tau1**(-5.0/2.0)
    phi2 = -3.0*tau2**(-5.0/2.0)
    
    Uxx = -1 + mu1*psi1 + mu1*phi1*zeta1**2 + \
               mu2*psi2 + mu2*phi2*zeta2**2
    
    Uxy = mu1*zeta1*y*phi1 + mu2*zeta2*y*phi2

    Uyy = -1 + mu1*psi1 + mu1*phi1*y**2 + \
               mu2*psi2 + mu2*phi2*y**2
    
    # Df Matrix (i.e. grad(EOM))
    Df = matrix([[   0,    0,  1, 0], \
                 [   0,    0,  0, 1], \
                 [-Uxx, -Uxy,  0, 2], \
                 [-Uxy, -Uyy, -2, 0]])
                 
    # reshape phi (STM) from array into matrix
    phi = states[8:]
    phi = phi.reshape(4,4)

    # matrix math to produce phi-dot
    phidot = Df*phi

    # reshape phi-dot into a 16x1 array and add to dstates
    phidot = phidot.reshape(16,)
    dstates[8:] = phidot

    # return dstates to caller
    return dstates

# ----------------------------------
#   Circular Restricted 3-Body EOM
# ----------------------------------
def CR3BP(t, states):
    """CR3BP(t, states)
        
        INPUT
        ----------
        t        : 1x1 scalar
        states   : 10x1 list or array
                 : [0]: x
                 : [1]: y
                 : [2]: z
                 : [3]: vx
                 : [4]: vy
                 : [5]: vz
                 : [6]: mu1
                 : [7]: mu2
                 : [8]: x1 'x-coordinate of primary'
                 : [9]: x2 'x-coordinate of secondary'
        
        ASSUMPTIONS: that primary and secondary are located on the x-axis
        
        """
    
    #  import statements
    from numpy import zeros, ndarray
    
    mu1 = states[6]
    mu2 = states[7]
    x1  = states[8]
    x2  = states[9]
    
    dx13 = states[0] - x1
    dx23 = states[0] - x2
    dy   = states[2]
    dz   = states[3]
    
    #  distance from primary
    r13 = (dx13**2 + dy**2 + dz**2)**(0.5)
    #  distance from secondary
    r23 = (dx23**2 + dy**2 + dz**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    #  derivatives of states
    dstates[:3] = states[3:6]
    
    dstates[3] =  2*states[4] + states[0] - mu1*dx13/(r13**3) - mu2*dx23/(r23**3)
    dstates[4] = -2*states[3] + states[1] - mu1*dy/(r13**3)   - mu2*dy/(r23**3)
    dstates[5] = -mu1*dz/(r13**3) - mu2*dz/(r23**3)
    
    return dstates

# ----------------------------------
#   Circular Restricted 3-Body EOM
#                 +
#        State Trans. Matrix
# ----------------------------------
def CR3BP_varEqns(t, states):
    """CR3BP_varEqns(t, states)
        
        INPUT
        ----------
        t        : 1x1 scalar
        states   : 46x1 list or array
                 : [0]: x
                 : [1]: y
                 : [2]: z
                 : [3]: vx
                 : [4]: vy
                 : [5]: vz
                 : [6]: mu1
                 : [7]: mu2
                 : [8]: x1 'x-coordinate of primary'
                 : [9]: x2 'x-coordinate of secondary'
                 : [10-45]: Phi
        
        ASSUMPTIONS: that primary and secondary are located on the x-axis
        
        """
    
    #  import statements
    from numpy import zeros, ndarray
    
    mu1 = states[6]
    mu2 = states[7]
    x1  = states[8]
    x2  = states[9]
    
    dx13 = states[0] - x1
    dx23 = states[0] - x2
    dy   = states[2]
    dz   = states[3]
    
    #  distance from primary
    r13 = (dx13**2 + dy**2 + dz**2)**(0.5)
    #  distance from secondary
    r23 = (dx23**2 + dy**2 + dz**2)**(0.5)
    
    dstates = zeros(states.shape[0])
    
    #  derivatives of states
    dstates[:3] = states[3:6]
    
    dstates[3] =  2*states[4] + states[0] - mu1*dx13/(r13**3) - mu2*dx23/(r23**3)
    dstates[4] = -2*states[3] + states[1] - mu1*dy/(r13**3)   - mu2*dy/(r23**3)
    dstates[5] = -mu1*dz/(r13**3) - mu2*dz/(r23**3)
    
    # terms for variational equations
    x = states[0]
    y = states[1]
    z = states[2]
    
    zeta1 = (x - x1)
    zeta2 = (x - x2)
    
    tau1 = zeta1**2 + y**2
    tau2 = zeta2**2 + y**2
    
    psi1 = tau1**(-3.0/2.0)
    psi2 = tau2**(-3.0/2.0)
    
    phi1 = -3.0*tau1**(-5.0/2.0)
    phi2 = -3.0*tau2**(-5.0/2.0)
    
    Uxx = -1 + mu1*psi1 + mu1*phi1*zeta1**2 + \
               mu2*psi2 + mu2*phi2*zeta2**2
    
    Uxy = mu1*zeta1*y*phi1 + mu2*zeta2*y*phi2

    Uxz = mu1*zeta1*z*phi1 + mu2*zeta2*z*phi2
    
    Uyy = -1 + mu1*psi1 + mu1*phi1*y**2 + \
               mu2*psi2 + mu2*phi2*y**2

    Uyz = mu1*y*z*phi1 + mu2*y*z*phi2

    Uzz = mu1*phi1*z**2 + mu2*phi2*z**2
    
    # Df Matrix (i.e. grad(EOM))
    Df = matrix([[   0,    0,     0,  1, 0, 0], \
                 [   0,    0,     0,  0, 1, 0], \
                 [   0,    0,     0,  0, 0, 1], \
                 [-Uxx, -Uxy,  -Uxz,  0, 2, 0], \
                 [-Uxy, -Uyy,  -Uyz, -2, 0, 0], \
                 [-Uxz, -Uyz,  -Uzz,  0, 0, 0]])
                 
    # reshape phi (STM) from array into matrix
    phi = states[10:]
    phi = phi.reshape(6,6)

    # matrix math to produce phi-dot
    phidot = Df*phi

    # reshape phi-dot into a 36x1 array and add to dstates
    phidot = phidot.reshape(36,)
    dstates[10:] = phidot

    # return dstates to caller
    return dstates

#---------------------------------------------------------------------------#
#---------------------- Astrodynamical Models Module -----------------------#
#---------------------------------------------------------------------------#