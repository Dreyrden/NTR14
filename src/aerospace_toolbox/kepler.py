#  laguerre_conway solver
def laguerre_conway(e, M):
    #  import statements
    from math  import sqrt, sin, cos
    from numpy import abs
    #  generate an initial guess for eccentric anomaly (E)
    En = (M*(1 - sin(M + e)) + (M + e)*sin(M)) / (1 + sin(M) - sin(M + e))
    n  = 4
    #  for-loop
    for i in range(1,20):
        f      =  M - En + e*sin(En)
        fdash  = -1 + e*cos(En)
        fddash = -e*sin(En)
        g = sqrt(((n - 1)**2) * (fdash**2) - n*(n - 1)*f*fddash)
        if fdash > 0:
            En1 = En - (n*f/(fdash + g))
        else:
            En1 = En - (n*f/(fdash - g))
        #  calculate error
        error = abs(En1 - En)
        En = En1
        if error <= 1E-4:
            break

    #  return the eccentric anomaly
    return En



#  solve kepler's equation using laguerre_conway
def kepler(list, delta_time, \
           index_dictionary = {'ecc': 0, \
                               'mean motion': 1, \
                               'eccentric anomaly': 2}):
    #  input delta_time in (secs)
    #  import statements
    from math import pi, sqrt, sin, tan, atan
    #  update each item in the list to the given epoch
    #  default index positions of elements in item come from
    #+ index_dictionary
    e_index = index_dictionary['ecc']
    n_index = index_dictionary['mean motion']
    E_index = index_dictionary['eccentric anomaly']
    #  check if list argument is a single list
    #+ or a list of list
    islist = False
    if type(list[0]) == type([]): islist = True
    if islist:
        for item in list:
            e   = item[e_index]
            rhs = item[E_index] - e*sin(item[E_index]) - \
                  item[n_index]*(-delta_time)
            E   = laguerre_conway(e, rhs)
            #  update entries of item
            item[E_index] = E % (2*pi)
    else:
        e   = list[e_index]
        rhs = list[E_index] - e*sin(list[E_index]) - \
              list[n_index]*(-delta_time)
        E   = laguerre_conway(e, rhs)
        #  update entries of item
        list[E_index] = E % (2*pi)
    
    return list



#  calculate ephemeris using laguerre_conway
def ephemeris(list, epoch, \
              index_dictionary = {'ecc': 0, \
                                  'mean motion': 1, \
                                  'eccentric anomaly': 2, \
                                  'current epoch': 3}):
    #  input epoch in (days)
    #  import statements
    from math import pi, sqrt, sin, tan, atan
    #  update each item in the list to the given epoch
    #  default index positions of elements in item come from
    #+ index_dictionary
    e_index = index_dictionary['ecc']
    n_index = index_dictionary['mean motion']
    E_index = index_dictionary['eccentric anomaly']
    c_index = index_dictionary['current epoch']
    #  check if list argument is a single list
    #+ or a list of list
    islist = False
    if type(list[0]) == type([]): islist = True
    if islist:
        for item in list:
            dt  = ( -epoch + item[c_index])
            e   = item[e_index]
            rhs = item[E_index] - e*sin(item[E_index]) - item[n_index]*dt
            E   = laguerre_conway(e, rhs)
            #  update entries of item
            item[c_index]  = epoch
            item[E_index] = E % (2*pi)
    else:
        dt  = ( -epoch + list[c_index])
        e   = list[e_index]
        rhs = list[E_index] - e*sin(list[E_index]) - list[n_index]*dt
        E   = laguerre_conway(e, rhs)
        #  update entries of item
        list[c_index] = epoch
        list[E_index] = E % (2*pi)

    return list

