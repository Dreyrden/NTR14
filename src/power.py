import math

def calculate_power(parent):
    '''
    mission, tank_area, ZBO,
   '''

    # System Trade-Study Mass/Power Requirements
    system_power = 0.0
    system_mass = 0.0
    parent.power += system_power

    orbital_parameter = 3.986E14
    earth_radius = 6378E3

    orbital_period = 2*math.pi()*math.pow(math.pow(parent.sma, 3)/orbital_parameter,1/2)
    eclipse_period = (orbital_period*math.asin(earth_radius/parent.sma))/180

    # Assuming Peak power tracking regulation scheme Xe = 0.6, Xd = 0.8
    Xe = 0.6
    Xd = 0.8
    solar_panel_power = (parent.power*eclipse_period/Xe + parent.power*(orbital_period-eclipse_period)/Xd)/Xd

    panel_efficiency = 0.28
    solar_flux = 1368.00
    inherent_degredation = 0.77
    sun_incidence_angle = 23.5
    power_beginning_of_life = panel_efficiency*solar_flux*inherent_degredation*math.cos(sun_incidence_angle)

    degredation_rate = 0.5
    lifetime = 15
    lifetime_degredation = math.pow((1-degredation_rate), lifetime)
    M_A = 455.0/25.0
    power_end_of_life = power_beginning_of_life*lifetime_degredation
    solar_panel_area = solar_panel_power/power_end_of_life
    solar_panel_mass = solar_panel_area*M_A

    parent.mass += array_mass + system_mass

