

def tank_power(mission, N, tank_area, ZBO, spow, days):
    '''
    mission, tank_area, ZBO,
    '''
    
    # Constants    
    solar_flux = 1340.00
    
    # Heat Transfer Constaints
    t_h = 302.00
    t_l = 20.00
    
    # Empirical Transfer Constants
    Q_A = (1.8/N)*( 1.022*10**(-4)*( (t_h + t_l)/2 ) * (t_h - t_l )\
    + 1.67*10**(-11) * (t_h**(4.67) - t_l**(4.67)) )
     
    print Q_A
    
    dh_LH2 = 446000.00
    dh_LH2_day = (24*3600)/dh_LH2
    
    rate_area = dh_LH2_day*Q_A
    
    Q_tank = Q_A*tank_area
    
    Q_total = Q_tank
    
    P_req = 100.00 * Q_total
    
    if ZBO == 0:
        P_req = 0
        
    alpha_ZBO = 1100.0/7.37
    M_ZBO = alpha_ZBO*P_req/1000.00
    
    P_SAreq = (P_req) + spow
    eta = 0.2
    M_A = 455.0/25.0
    
    SA_area = P_SAreq/(solar_flux * eta)
    M_array = SA_area*M_A
    
    rate = tank_area*rate_area
    loss = rate*days
    
    return (M_array, M_ZBO, loss)
    