#  converting length
def in2cm(inches):  return inches*2.54
def cm2in(cm):      return cm/2.54
def ft2m(ft):       return ft*0.3048
def m2ft(m):        return m/0.3048

#  converting mass
def lb2kg(lb): return lb*0.4535
def kg2lb(kg): return kg/0.4535

#  converting force
def lbf2N(lbf): return lbf*4.4482
def N2lbf(N):   return N/4.4482

#  converting temperature
def R2K(R): return R*5.0/9.0
def K2R(K): return K*9.0/5.0
def C2K(c): return c + 273.15
def K2C(K): return K - 273.15
def F2C(f): return (f - 32.0)*(5.0/9.0)
def C2F(c): return c*(9.0/5.0) + 32.0
def F2K(f): return (f + 459.67)*(5.0/9.0)

#  converting pressure
def atm2Pa(atm):  return atm*101325.0
def atm2psi(atm): return atm*14.7
def psi2atm(psi): return psi/14.7
def MPa2Pa(MPa):  return MPa*1E6
def Pa2MPA(Pa):   return Pa*1E-6
def psi2Pa(psi):  return psi*101325.0/14.7
def Pa2psi(Pa):   return Pa*14.7/101325.0
def psi2MPa(psi): return psi2Pa(psi)*1E-6
def MPa2psi(MPa): return Pa2psi(MPa*1E6)

#  main unit converter program
def converter(object, old_units, new_units):
    """ converter() should be passed a
        scalar or array to be converted, 
        the current unit identifier as a string, and
        a desired unit as a string.
        """
    
    degrees_list = ['deg', 'Deg', 'degrees', 'Degrees']
    radians_list = ['rad', 'Rad', 'radians', 'Radians']
    
    if old_units in degrees_list:
        if new_units in radians_list:
            from math import pi
            return object * pi / 180.0

    if old_units in radians_list:
        if new_units in degrees_list:
            from math import pi
            return object * 180.0 / pi