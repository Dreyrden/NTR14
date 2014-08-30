#  owner: rg
import math
import power

class Tank:
    'Tank class documentation string'

    #  initialization method
    def __init__(self):
        self.tanklength = 0.0
        self.radius = 0.0
        self.structural_mass = 0.0
        self.fuel = 0.0


class LH2Tank(Tank):
    'LH2Tank class documentation string'

    def __init__(self, parent):
        self.DEBUG = False
        self.name = 'LH2 Tank'
        self.tanklength = 0.0
        self.radius = 0.0
        self.structural_mass = 0.0
        self.fuel = 0.0
        self.power = 0.0
        self.zbo_mass = 0.0
        self.propellant_loss = 0.0
        #  add the parent object as an attribute
        #+ and update the parent's attributes
        self.parent = parent
        parent.mass += self.structural_mass + self.fuel
        parent.fuel += 0.0

    def add_fuel(self, fuel_to_add):
        self.fuel   += fuel_to_add
        self.parent.mass += fuel_to_add
        self.parent.fuel += fuel_to_add

    def remove_fuel(self, fuel_to_remove):
        self.fuel -= fuel_to_remove
        self.parent.mass -= fuel_to_remove
        self.parent.fuel -= fuel_to_remove

    def tank_size(self, radius, mli_layers, geometry = 'cylindrical'):
        '''Fuel must already have been added. Radius is the maximum constrained
        radius

        Usage:
        mission.ntr.add_tank()
        mission.ntr.tank[0].add_fuel(fuel_mass)
        mission.ntr.tank[0].tank_size(geometry, radius, mli_layers)
        '''

        self.radius = radius
        self.mli_layers = mli_layers

        # Define some constants
        hydrogen_density = 70.85
        tank_material_density = 2685
        stress_factor = 0.79
        welding_efficiency = 0.6
        max_working_stress = 234586466 # Lowest out of Ultimate Tensile Strength/1.65 and Yield Strength/1.33
        tank_pressure = 230000

        # Calculate total tank volume required
        total_tank_volume = (self.fuel)/hydrogen_density
        if self.DEBUG: print 'total_tank_volume', total_tank_volume

        if geometry == 'cylindrical':

            endcap_ellipsoid_ratio = 1.395

            # Endcaps -------------------------------------------------------------
            # Endcap Dimenions
            endcap_ellipsoid_height = self.radius/endcap_ellipsoid_ratio
            endcap_crown_radius = endcap_ellipsoid_ratio*self.radius
            eccentricity = math.pow((math.pow(self.radius,2)-math.pow(endcap_ellipsoid_height,2)),0.5)/self.radius
            endcap_volume = (2*math.pi*math.pow(self.radius,2)*endcap_ellipsoid_height)/3

            if self.DEBUG:
                print 'endcap_ellipsoid_height', endcap_ellipsoid_height
                print 'endcap_crown_radius', endcap_crown_radius
                print 'eccentricity', eccentricity
                print 'endcap_volume',  endcap_volume

            # Calculate Endcap thicknesses
            knuckle_thickness = (stress_factor*tank_pressure*self.radius)/(max_working_stress*welding_efficiency)
            crown_thickness = (tank_pressure*endcap_crown_radius)/(2*max_working_stress*welding_efficiency)
            equivalent_wall_thickness = (knuckle_thickness+crown_thickness)/2

            if self.DEBUG:
                print 'knuckle_thickness', knuckle_thickness
                print 'crown_thickness', crown_thickness
                print 'equivalent_wall_thickness', equivalent_wall_thickness

            # Endcap Surface Area
            endcap_surface_area = math.pow(self.radius,2) + ((math.pi*math.pow(endcap_ellipsoid_height,2)*math.log((1+eccentricity)/(1-eccentricity))))/(2*eccentricity)
            endcap_insulation_mass = endcap_surface_area*(0.78+0.015*self.mli_layers)
            design_factor = 2*endcap_ellipsoid_ratio + (1/math.pow((math.pow(endcap_ellipsoid_ratio,2)-1),0.5))*math.log((endcap_ellipsoid_ratio+math.pow(math.pow(endcap_ellipsoid_ratio,2)-1,0.5))/(endcap_ellipsoid_ratio-math.pow(math.pow(endcap_ellipsoid_ratio,2)-1,0.5)))
            endcap_mass = (math.pi*math.pow(self.radius,2)*equivalent_wall_thickness*design_factor*tank_material_density)/(2*endcap_ellipsoid_ratio) + endcap_insulation_mass

            if self.DEBUG:
                print 'endcap_surface_area', endcap_surface_area
                print 'endcap_insulation_mass', endcap_insulation_mass
                print 'design_factor', design_factor
                print 'endcap_mass', endcap_mass

            # Cylinder ------------------------------------------------------------
            cylinder_volume = total_tank_volume - 2*endcap_volume
            cylinder_length = cylinder_volume/(math.pi*math.pow(self.radius,2))
            cylinder_thickness = (tank_pressure*self.radius)/(max_working_stress*welding_efficiency)
            cylinder_surface_area = 2*math.pi*self.radius*cylinder_length
            cylinder_insulation_mass = cylinder_surface_area*(0.78+0.015*self.mli_layers)
            cylinder_mass = 2*math.pi*self.radius*cylinder_length*cylinder_thickness*tank_material_density + cylinder_insulation_mass

            if self.DEBUG:
                print 'cylinder_volume', cylinder_volume
                print 'cylinder_length', cylinder_length
                print 'cylinder_thickness', cylinder_thickness
                print 'cylinder_surface_area', cylinder_surface_area
                print 'cylinder_insulation_mass', cylinder_insulation_mass
                print 'cylinder_mass', cylinder_mass

            surface_area = cylinder_surface_area + 2*endcap_surface_area
            self.structural_mass += 2*endcap_mass + cylinder_mass #parent?
            self.tanklength += cylinder_length + 2*endcap_ellipsoid_height #parent?
            self.fuel_capacity = self.fuel #parent?
            #  parent mass update
            self.parent.mass += self.structural_mass


        elif geometry == 'spherical':

            self.radius = math.pow((total_tank_volume*3)/(4*math.pi),1/3)

            wall_thickness = (tank_pressure*self.radius)/(2*max_working_stress*welding_efficiency)
            sphere_mass = (4*math.pi*math.pow(self.radius,2)*wall_thickness*tank_material_density)

            surface_area = 4*math.pi*math.pow(self.radius,2)
            insulation_mass = surface_area*(0.78+0.015*self.mli_layers)

            self.structural_mass += sphere_mass + insulation_mass #parent?
            self.tanklength += 2*self.radius #parent?
            self.fuel_capacity = self.fuel #parent?
            #  parent mass update
            self.parent.mass += self.structural_mass

        else:
            print 'geometry selection error'

        zbo_power(self, surface_area)


    def zbo_power(self, surface_area):
        '''Function calculates

        removed days, spow (other power requirements), zbo as parameters. These
        are assigned/evaluated elsewhere'''
        #Constants

        # Heat Transfer Constaints
        t_h = 302.00
        t_l = 20.00

        # Empirical Transfer Constants
        Q_A = (1.8/self.mlilayers)*( 1.022*10**(-4)*( (t_h + t_l)/2 ) * (t_h - t_l )\
    + 1.67*10**(-11) * (t_h**(4.67) - t_l**(4.67)) )

        dh_LH2 = 446000.00
        dh_LH2_day = (24*3600)/dh_LH2

        rate_area = dh_LH2_day*Q_A

        Q_tank = Q_A*surface_area
        Q_total = Q_tank
        P_req = 100.00 * Q_total

        alpha_ZBO = 1100.0/7.37
        ZBO_mass = alpha_ZBO*P_req/1000.00

        rate = surface_area*rate_area

        self.zbo_mass += ZBO_mass #parent
        self.power += P_req #parent
        self.propellant_loss = rate #parent


#Testing Code
#tank = Tank()
#htank = LH2Tank(tank)
#htank.DEBUG = True
#htank.add_fuel(3000)
#htank.tank_size('cylindrical', 2, 120)
#print tank.fuel, tank.mass, tank.length, tank.radius
#print htank.fuel, htank.structural_mass, htank.tanklength

