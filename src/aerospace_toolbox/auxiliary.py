#-------------------------------------------------------------------------#
#---------------------- Auxiliary Functions Module -----------------------#
#-------------------------------------------------------------------------#

#  Orbital Elements -> r, v Arrays
def coe2rv(oe, mu):
    #  import statements
    from math  import cos, sin, sqrt
    from numpy import matrix, array, zeros
    
    #  if 'a' is set to zero, then body is at rest in the
    #+ current frame of reference.
    #+ return a zero vector for r and v
    if oe[0] == 0.0: return zeros(3,1), zeros(3,1)
    
    a  = oe[0]
    e  = oe[1]
    i  = oe[2]
    Om = oe[3]
    om = oe[4]
    f  = oe[5]

    p  = a*(1 - e*e)
    r  = p/(1 + e*cos(f))
    rv = matrix([r*cos(f), r*sin(f),   0])
    vv = matrix([-sin(f),  e + cos(f), 0])
    vv = sqrt(mu/p)*vv

    c0 = cos(Om); s0 = sin(Om)
    co = cos(om); so = sin(om)
    ci = cos(i);  si = sin(i)

    R  = matrix([[c0*co - s0*so*ci, -c0*so - s0*co*ci,  s0*si],
                 [s0*co + c0*so*ci, -s0*so + c0*co*ci, -c0*si],
                 [so*si,             co*si,             ci]])

    ri = array(R*rv.T); ri = ri.reshape(3)
    vi = array(R*vv.T); vi = vi.reshape(3)
    
    return ri, vi

# r, v Arrays -> Orbital Elements
def rv2coe(rv, vv, mu):
    """ FUTURE WORK:
        1.) revisit trunctation occuring for arg > 1 : 2014_06_22
        """
    # import statements
    from math  import acos, pi, sqrt
    from numpy import ndarray, array, cross
    from numpy.linalg import norm
    
    # check if list or numpy.ndarray
    if type(rv) == list and len(rv) == 2: rv.append(0.0)
    if type(rv) == ndarray:
        rv = list(rv)
        if len(rv) == 2: rv.append(0.0)
    # check if list or numpy.ndarray
    if type(vv) == list and len(vv) == 2: vv.append(0.0)
    if type(vv) == ndarray:
        vv = list(vv)
        if len(vv) == 2: vv.append(0.0)

    rv = array(rv)
    vv = array(vv)

    k  = [0, 0, 1]
    hv = cross(rv, vv);
    nv = cross(array(k),  hv)
    n  = norm(nv)
    h2 = hv.dot(hv)
    v2 = vv.dot(vv)
    r  = norm(rv)
    ev = (1/mu)*((v2 - mu/r)*rv - (rv.dot(vv))*vv)
    p  = h2/mu
    e  = norm(ev)
    a  = p/(1 - e**2)
    i  = acos(hv[2]/sqrt(h2))
    Om = acos(nv[0]/n)

    # eps
    eps = 1E-6
    if nv[1] < - eps: Om = 2*pi - Om

    om = acos(nv.dot(ev)/n/e)
    if ev[2] < 0: om = 2*pi - om

    arg = ev.dot(rv)/e/r
    #  truncate acos() argument if it exceeds 1 and
    #  therefore violates acos() domain
    if norm(arg) > 1: arg = 1
    nu = acos(arg)
    if rv.dot(vv) < 0: nu = 2*pi - nu

    oe = [a, e, i, Om, om, nu]
    return oe

#  calculation of eccentric anomaly from true anomaly
def f2E(e, f):
    from math import sqrt, tan, atan, pi
    #  return the eccentric anomaly
    return 2.0*atan(sqrt((1 - e) / (1 + e)) * tan(f/2.0)) % (2.0*pi)

#  calculation of true anomaly from eccentric anomaly
def E2f(e, E):
    from math import sqrt, tan, atan, pi
    #  return the true anomaly
    return 2.0*atan(tan(E/2.0) * sqrt((1 + e) / (1 - e))) % (2.0*pi)

#  calculation of mean anomaly from eccentric anomaly
def E2M(e, E):
    from math import sin, pi
    #  return the mean anomaly
    return (E - e*sin(E)) % (2.0*pi)

#  orbital period
def orbitalperiod(mu, a):
    from math import pi, sqrt
    #  return the orbital period
    return 2.0*pi*sqrt((a**3)/mu)

#  calculate the mean angular rate
def meanmotion(mu, a):
    from math import pi
    return 2.0*pi / orbitalperiod(mu, a)

#  visviva
def visviva(mu, r, a):
    from math import sqrt
    #  return v
    return sqrt(mu*((2/r) - (1/a)))
             
#  hohmann transfer
def hohmann(mu, r1, r2):
    from math import sqrt
    dv1 = sqrt(mu/r1)*(sqrt(2*r2/(r1 + r2)) - 1)
    dv2 = sqrt(mu/r2)*(1 - sqrt(2*r1/(r1 + r2)))
    return dv1, dv2

#  tsiolkovsky rocket equation
def tsiolkovsky(dv, Isp, m0):
    #  dv:  (km/s)
    #  Isp: (s)
    from solarsystem_objects_parameters import earth
    from math import exp
    #  earth['g'] is in (km/s^2)
    return m0*exp(-dv/(Isp*earth['g']))

#  parabolic transfer time
def time_parabolic(s, c, mu, sgn):
    from math import sqrt
    
    tp = (sqrt(2.0)/3)*(s**(3.0/2) - sgn*((s - c)**(3.0/2)))/sqrt(mu)
    return tp

#  minimum energy transfer time
def time_minenergy(s, c, mu, theta):
    from math import asin, sin, sqrt, pi
    
    # (s - c)/s is an element of [0, 1]
    # therefore asin() is always real for the
    # arguement (s - c)/s.
    # specifically asin() is an element of [0, pi/2]
    beta = 2*asin(sqrt((s - c)/s))
    
    # flip beta if 180 deg <= theta < 360 deg
    if theta >= pi and theta < 2*pi: beta = -beta
    
    # calculate minimum energy transfer time
    tm = sqrt(s**(3.0)/8)*(pi - beta + sin(beta))/sqrt(mu)
    return tm

# lambert's Equation
def lamberteq(a, mu, dt, s, c, tm, theta):
    from math  import asin, sin, sqrt, pi
    from numpy.linalg import norm

    # alpha and beta are the rectilinear equivalent of
    # eccentric anomaly for transfer orbit
    alpha = 2*asin(sqrt(s/(2.0*a)))
    beta  = 2*asin(sqrt((s-c)/(2.0*a)))
    
    # perform necessary corrections
    if dt > tm: alpha = 2*pi - alpha
    if theta >= pi and theta < 2*pi: beta = -beta
    
    # residual of lambert's equation
    residual = sqrt(mu) * dt - \
               (a**(3.0/2)) * (alpha - beta - sin(alpha) + sin(beta))
    '''
    as of 2014_05_15, testing in Matlab showed that using norm()
    on residual and no abs() in while-loop for lambert() gives decent
    convergence, which can subsequently be converged to desired
    tolerance using a shooting method.
    >> should revisit while-loop iterative method to unconstrain this
    norm() function
    '''
    residual = norm(residual)
    
    return residual #, alpha, beta

# lambert Solver
def lambert(ic1, ic2, mu,                \
            TransferTime        = 0,     \
            FindTimeParabolic   = False, \
            FindMinEnergy       = False, \
            FindlambertArc      = False, \
            Retrograde          = False, \
            NonDimUnits         = True,  \
            ScaleOutput         = True,  \
            iLimit  = 500, \
            tol     = 1E-2):
    """
        lambert()
        
        Future Work
        1.) Need to eliminate v1, v2 input arguments / 2014_05_18
        2.) lambert should be able to use minimum time for 
            TransferTime without having to call lambert() twice
        3.) Should separate lambert functions into a separate module,
            which should eliminate r,v being set to 3d vectors from 2d 
            by rv2coe() call
        4.) scipy calls appear to be working fine, should eliminate 
            all the junk code floating around and cleanup
        5.) Shouldn't be rescaling if not using NonDimUnits
        6.) Need Spatial Capability!!!
    """
    
    ''' DEBUG_SCIPY '''
    DEBUG_SCIPY = True
    
    # import statements
    from math  import acos, asin, sin, pi, cos, sqrt
    from numpy import array
    from numpy.linalg import norm
    
    # initialze output dictionary
    out = {'a':     -1, \
           'res':   -1, \
           'iter':   0, \
           'alpha': -1, \
           'beta':  -1, \
           'tp':    -1, \
           'tm':    -1, \
           'v1':     0, \
           'v2':     0, \
           'DU':    -1, \
           'TU':    -1, \
           'MU':    -1, \
           'theta': -1}

    #  select lambert arc computation if neither
    #  parabolic time nor min. energy selected
    if not FindTimeParabolic and not FindMinEnergy:
        FindlambertArc = True

    #  initializations
    r1 = ic1[0:2]
    r2 = ic2[0:2]
    v1 = ic1[2:4]
    v2 = ic2[2:4]
    #  convert r, v to coe
    oe1 = rv2coe(r1, v1, mu)
    oe2 = rv2coe(r2, v2, mu)
    #  convert to array type
    r1 = array(r1[0:2])
    r2 = array(r2[0:2])
    v1 = array(v1[0:2])
    v2 = array(v2[0:2])

    #  nondim units
    G = 6.67384*1E-20 #(km^3/kg/s^2)
    DU = oe1[0]
    MU = mu/G
    TU = sqrt(1/(G*MU/DU**3))
    if NonDimUnits:
        mu = 1
        #  scaling
        oe1[0] = oe1[0]/DU
        oe2[0] = oe2[0]/DU
        r1 = r1/DU
        r2 = r2/DU
        v1 = v1*TU/DU
        v2 = v2*TU/DU
        TransferTime = TransferTime/TU

    # chord and semiperimeter
    chord = norm(r1 - r2)
    semip = (norm(r1) + norm(r2) + chord)/2

    # rotation of r1, 90 degrees
    temp = array([-r1[1], r1[0]])

    # find angle between r1 and r2
    theta  = acos(r1.dot(r2)/(norm(r1)*norm(r2)))
    theta2 = acos(temp.dot(r2)/(norm(r1)*norm(r2)))
    if theta2 >= pi/2: theta = 2*pi - theta

    # update out entry
    out['theta'] = theta

    # sign of sin(theta)
    if sin(theta) < 0:  sgn = -1
    if sin(theta) > 0:  sgn =  1
    if sin(theta) == 0: sgn =  0

    # parabolic transfer time
    if FindTimeParabolic or FindlambertArc:
        tp = time_parabolic(semip, chord, mu, sgn)
        out['tp'] = tp
    if FindTimeParabolic: return out

    # calculate the minimum energy transfer time (elliptic case)
    if FindMinEnergy or FindlambertArc:
        tm = time_minenergy(semip, chord, mu, theta)
        out['tm'] = tm
    if FindMinEnergy: return out

    #  let the user know if they selected an infeasbile
    #+ transfer time
    if TransferTime < tp:
        print ('lambert: The transfer time (%.2f) is too short ' \
               'for an Elliptic Orbit (< %.2f) --> Return' \
                %(TransferTime, tp))
        return out

    # calculate semi-major axis, alpha, and beta
    # scale tol if nondimunits == False
    if not NonDimUnits: tol = tol*DU

    ''' DEBUGGING for scipy.optimize 2014_05_18 '''
    if not DEBUG_SCIPY:
        # initialize
        res  = 10*tol
        a    = semip/2 #oe1[0]
        ainc = a*1E-3
        dir  = 1
        dx   = ainc
        damp = 10
        i    = 1

        # while-loop unitl tolerance is met
        while res > tol:
            # save old
            if i == 2: temp = res
            if i > 2:  res0 = temp; temp = res
        
            # flip the direction of ainc if necessary
            if i > 2 and res > res0: dir = -1*dir
            # calculate damping needed
            if i > 2 and res < res0:
                slope = (res - res0)/ainc
                if res > 0: dx = -res/slope
                if dx/ainc < damp: ainc = ainc/2
            # update a
            a = a + dir*ainc

            # calculate residual
            res = lamberteq(a, mu, TransferTime, \
                            semip, chord, tm, theta)

            # if limit reached -> break
            if i > iLimit: break
                
            # update i
            i += 1

    if DEBUG_SCIPY:
        ''' DEBUGGING for scipy.optimize 2014_05_18 '''
        from scipy.optimize import minimize
        a   = semip/2
        bnds = ((a, None),)
        res = minimize(lamberteq, array(a), \
                       (mu, TransferTime, semip, chord, tm, theta), \
                       bounds = bnds, method = 'L-BFGS-B')

        ''' DEBUGGING for scipy.optimize 2014_05_18 '''
#        print res
        a = res.x
        i = 0

    # alpha and beta are the rectilinear equivalent of
    # eccentric anomaly for transfer orbit
    alpha = 2*asin(sqrt(semip/(2.0*a)))
    beta  = 2*asin(sqrt((semip-chord)/(2.0*a)))

    # perform necessary corrections
    if TransferTime > tm: alpha = 2*pi - alpha
    if theta >= pi and theta < 2*pi: beta = -beta

    # switches for retrograde orbits.
    # will need post-shooting method to adjust
    if Retrograde == True: beta = -beta

    # compute delta-vs
    # use hohmann if beta ~0 or ~180 degrees
    # else
    # use default calculations for v1 and v2
    u1 = r1/norm(r1)
    u2 = r2/norm(r2)
    uc = (r2 - r1)/norm(r2 - r1)
    #if abs(pi - (alpha - beta)) < 1E-4 and oe1[1] < 1E-2 and oe2[1] < 1E-2:
    if abs(beta % pi) < 1E-4:
        dv1, dv2 = hohmann(mu, oe1[0], oe2[0])
        v1 = v1 + dv1*v1/norm(v1)
        v2 = v2 + dv2*v2/norm(v2)
    else:
        cota = cos(alpha/2)/sin(alpha/2)
        cotb = cos(beta/2)/sin(beta/2)
        A  = sqrt(mu/(4*a))*cota
        B  = sqrt(mu/(4*a))*cotb
        v1 = (B + A)*uc + (B - A)*u1
        v2 = (B + A)*uc - (B - A)*u2

    # compute semilatus rectum for output
    p = 4*a*(semip - norm(r1))*(semip - norm(r2))* \
        (sin((alpha + beta)/2)**2)/(chord**2)
    # compute orbital elements, output ecc and true anomaly
    oef = rv2coe(r1, v1, mu)
    f = oef[5]
    # if retrograde, mod() true anomaly (i.e. flip i pi)
    if abs(oef[2] - pi) < 1E-3: f = 2*pi - oef[5]
    psi = 2*pi - f
    ev = [cos(psi)*u1[0] - sin(psi)*u1[1], \
               sin(psi)*u1[0] + cos(psi)*u1[1]]
    # position of the vacant focus
    pf = [-2*a*oef[1]*ev[0], -2*a*oef[1]*ev[1]]

    # compile output
    out['a']     = a
    out['e']     = oef[1]
    out['ev']    = ev
    out['f']     = oef[5]
    out['T']     = orbitalperiod(mu, a)
    out['pf']    = pf
    ''' DEBUGGING scipy.optimize 2014_05_18 '''
    if DEBUG_SCIPY: out['res']   = res.fun
    if not DEBUG_SCIPY: out['res']   = res
    out['iter']  = i - 1
    out['alpha'] = alpha
    out['beta']  = beta
    out['time']  = TransferTime
    out['v1']    = v1
    out['v2']    = v2
    out['DU']    = DU
    out['TU']    = TU
    out['MU']    = MU
    # scale ouput if needed
    if ScaleOutput:
        out['a']  = out['a']*DU
        out['T']  = out['T']*TU
        out['tp'] = out['tp']*TU
        out['tm'] = out['tm']*TU
        out['v1'] = out['v1']*DU/TU
        out['v2'] = out['v2']*DU/TU

    return out


#-------------------------------------------------------------------------#
#---------------------- Auxiliary Functions Module -----------------------#
#-------------------------------------------------------------------------#