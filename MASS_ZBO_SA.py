'''
Created on Jul 2, 2013
@author: addison

Calculate the ZBO / Solar array 
power / mass requirements 
for an NTR ferry based 
based tank geometry.
'''
def calc_mass(A_tank,N,ZBO,Days,spow):
    '''
    Pass: Surface Area (m^2), N-layer MLI, ZBO (1,0), Days
    Return: Mass of Array (kg), Mass ZBO (kg), Prop. Mass Loss(kg)
    '''
    ##############################
    # Global Constants.
    ##############################
    Fs_m = 486.00             # MIN Solar Flux at Mars (W/m2) 
    
    ##############################
    # Heat transfer Constants
    ############################## 
    T_H = 302.00
    T_L = 20.00
        
    # empirical heat flux equation
    # as function of N-layer MLI, T_H, T_L
    Q_A = (1.8/N)*( 1.022*10**(-4)*( (T_H + T_L)/2 ) * (T_H - T_L )\
    + 1.67*10**(-11) * (T_H**(4.67) - T_L**(4.67)) )
    
    # Boil off Calculation
    dH_LH2= 446000.00                 # LH2 enthalpy of vap.  (J/kg)
    dH_LH2_day = (24*3600)/dH_LH2    # Boil rate  (kg/W-day)

    # Rate of boil-off (kg/day-m^2)
    rate_area =  dH_LH2_day *Q_A  # (kg/day-m2) 

    # heat loss / tank
    Q_tank =  Q_A *A_tank          # (W/tank) 
    
    # Q_loss total
    Q_tot =  Q_tank                # Q_loss total (W) 
    
    # Est. ZBO Power req.
    P_req = 100.00 * Q_tot            # Power req. (We) @100We/Wt
    
    # If No ZBO system, Power Req. = 0
    if ZBO == 0:
        P_req = 0
        
    ##############################
    # Calc. mass of ZBO system
    # using DRA estimates
    ##############################
    alpha_ZBO = 1100.0/7.37          # 1100 kg/ kWe ZBO
    M_ZBO = alpha_ZBO * P_req/1000.0
    
    P_SAreq = (P_req) + spow       # Power req. for Solar Array (W)
    eta= .2                        # Efficiency
    M_A = 455.0/25.0                   # Area density of solar array (kg/m^2)
    
    SA_area = P_SAreq/(Fs_m * eta) # Solar Array Area (m^2)
    M_array = SA_area*M_A          # Solar Array Mass (m^2)
    
    ##############################
    # Calc. Prop Mass Loss
    ##############################
    
    rate = A_tank * rate_area;     # Rate of boil-off kg/day
    LOSS = rate * Days;            # Total mass lost (kg)

    return (M_array,M_ZBO,LOSS)
"""
(Marray,Mzbo,Loss) = calc_mass(495.00,120.00,1,360.00)

print 'Solar Array Mass + %.3f kg' % (Marray)
print 'ZBO system Mass  + %.3f kg' % (Mzbo)
print 'Loss Prop. Mass  + %.3f kg' % (Loss)
"""