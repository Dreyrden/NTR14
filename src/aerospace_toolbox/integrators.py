#----------------------------------------------------------------#
#---------------------- ODE Solver Module -----------------------#
#----------------------------------------------------------------#


# ---------------------
#      RK4 - PRINT
# ---------------------
def rk4p(fhandle, xspan, h, states):
    """
        
     * COPIED FROM MATLAB rk4_print
     Really only working with arrays for states so no need to 
     check for transpose early on in function
       
     Title: [out] = rk4_print(fhandle, xspan, h, states)
     Type:  Function

     Purpose:   An Explicit 4th Order Runge-Kutta Method for Solution of
                Ordinary Differential Equations

     Inputs:
                fhandle : (e.g @myfunc)
                xspan   : 1x2 or 2x1 vector with start and end for indep.
                         variable
                h       : 1x1 scalar
                states  : nx1 column vector or 1xn row vector

     Outputs:
                out.x       : 1xm row vector
                out.states  : nxm matrix

     Original Date:   2014_03_09
     Original Author: Ryne Beeson
     
     Versions / Comments
     Version: 0.02 / Can Do Backwards Propagation
     Version: 0.01 / Original Code Created
     
    """
    
    # import statements
    from numpy import zeros
    
    ''' FUTURE WORK 
        1.) eliminate x/xf if-statements for printing
        '''
    
    # initializations
    print_inc  = 0.1
    print_frac = 0.1
    x  = xspan[0]
    xf = xspan[1]
    s  = zeros(states.shape[0])
    
    # convert all elements of the array states to floating point
    for i in range(states.shape[0]):
        s[i] = float(states[i])
    
    # dictionary output
    out = {'x': [x], 'y': [s]}

    # forward integration of states, else backwards
    if x < xf:
        #  increment x by h or (xf - x),
        #+ which ever is smaller
        if h <= (xf - x): x += h
        else: x += (xf - x)
        while x <= xf:

            k1 = fhandle(x, s)
            k2 = fhandle(x + h/2, s + h*k1/2)
            k3 = fhandle(x + h/2, s + h*k2/2)
            k4 = fhandle(x + h,   s + h*k3)

            s = s + h*(k1 + 2*k2 + 2*k3 + k4)/6
    
            out['x'].append(x)
            out['y'].append(s)

            #  if x is equal to the xf, then break from while-loop
            if x == xf: print ('Job Percent Completed: 100.0%'); break
            #  increment x by h or (xf - x),
            #+ which ever is smaller
            if h <= (xf - x): x += h
            else: x += (xf - x)
            #  backwards integration
            
            # print job percentage completed
            if x/xf > print_frac and xf != 0.0:
                print ('Job Percent Completed: %5.1f%%' %(100.0*x/xf))
                print_frac += print_inc
    elif x > xf:
        #  increment x by h or (x - xf),
        #+ which ever is smaller
        if h <= (x - xf): x -= h
        else: x -= (x - xf)
        while x >= xf:
    
            k1 = fhandle(x, s)
            k2 = fhandle(x - h/2, s - h*k1/2)
            k3 = fhandle(x - h/2, s - h*k2/2)
            k4 = fhandle(x - h,   s - h*k3)

            s = s - h*(k1 + 2*k2 + 2*k3 + k4)/6

            out['x'].append(x)
            out['y'].append(s)
            
            #  if x is equal to the xf, then break from while-loop
            if x == xf: print ('Job Percent Completed: 100.0%'); break
            #  increment x by h or (x - xf),
            #+ which ever is smaller
            if h <= (x - xf): x -= h
            else: x -= (x - xf)
            
            # print job percentage completed
            if abs(x/xf) > print_frac and xf != 0.0:
                print ('Job Percent Completed: %5.1f%%' %(100.0*abs(x/xf)))
                print_frac += print_inc

    # return output
    return out

# ---------------------
#          RK4
# ---------------------
def rk4(fhandle, xspan, h, states):
    """
        
        * COPIED FROM MATLAB rk4_print
        Really only working with arrays for states so no need to
        check for transpose early on in function
        
        Title: [out] = rk4_print(fhandle, xspan, h, states)
        Type:  Function
        
        Purpose:   An Explicit 4th Order Runge-Kutta Method for Solution of
        Ordinary Differential Equations
        
        Inputs:
        fhandle : (e.g @myfunc)
        xspan   : 1x2 or 2x1 vector with start and end for indep.
        variable
        h       : 1x1 scalar
        states  : nx1 column vector or 1xn row vector
        
        Outputs:
        out.x       : 1xm row vector
        out.states  : nxm matrix
        
        Original Date:   2014_03_09
        Original Author: Ryne Beeson
        
        Versions / Comments
        Version: 0.02 / Can Do Backwards Propagation
        Version: 0.01 / Original Code Created
        
        """

    # import statements
    from numpy import zeros
    
    # initializations
    x  = xspan[0]
    xf = xspan[1]
    s  = zeros(states.shape[0])
    
    # convert all elements of the array states to floating point
    for i in range(states.shape[0]):
        s[i] = float(states[i])
    
    # dictionary output
    out = {'x': [x], 'y': [s]}

    # forward integration of states, else backwards
    if x < xf:
        #  increment x by h or (xf - x),
        #+ which ever is smaller
        if h <= (xf - x): x += h
        else: x += (xf - x)
        while x <= xf:
            
            k1 = fhandle(x, s)
            k2 = fhandle(x + h/2, s + h*k1/2)
            k3 = fhandle(x + h/2, s + h*k2/2)
            k4 = fhandle(x + h,   s + h*k3)
            
            s = s + h*(k1 + 2*k2 + 2*k3 + k4)/6
            
            out['x'].append(x)
            out['y'].append(s)
            
            #  if x is equal to the xf, then break from while-loop
            if x == xf: break
            #  increment x by h or (xf - x),
            #+ which ever is smaller
            if h <= (xf - x): x += h
            else: x += (xf - x)
    #  backwards integration
    elif x > xf:
        #  increment x by h or (x - xf),
        #+ which ever is smaller
        if h <= (x - xf): x -= h
        else: x -= (x - xf)
        while x > xf:
            
            k1 = fhandle(x, s)
            k2 = fhandle(x - h/2, s - h*k1/2)
            k3 = fhandle(x - h/2, s - h*k2/2)
            k4 = fhandle(x - h,   s - h*k3)
            
            s = s - h*(k1 + 2*k2 + 2*k3 + k4)/6
            
            out['x'].append(x)
            out['y'].append(s)
            
            #  if x is equal to the xf, then break from while-loop
            if x == xf: break
            #  increment x by h or (x - xf),
            #+ which ever is smaller
            if h <= (x - xf): x -= h
            else: x -= (x - xf)
    # return output
    return out


#----------------------------------------------------------------#
#---------------------- ODE Solver Module -----------------------#
#----------------------------------------------------------------#