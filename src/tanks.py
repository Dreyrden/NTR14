#  owner: rg

class Tank:
    'Tank class documentation string'

    #  initialization method
    def __init__(self):
        self.length = 0.0
        self.radius = 0.0


class LH2Tank(Tank):
    'LH2Tank class documentation string'

    def __init__(self, parent):
        self.name = 'LH2 Tank'
        self.tanklength = 0.0
        self.radius = 0.0
        self.structural_mass = 0.0
        self.fuel = 0.0
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

    def tank_size(self, radius, mli_layers):
        '''Fuel must already have been added. Radius is the maximum constrained
        radius

        Usage:
        mission.ntr.add_tank()
        mission.ntr.tank[0].add_fuel(fuel_mass)
        mission.ntr.tank[0].tank_size(max_radius)
        '''

        self.radius = radius
        self.mli_layers = mli_layers

        # Define some constants
        hydrogen_density = 70.85
        tank_material_density = 2685
        stress_factor = 0.79
        welding_efficiency = 0.6
        max_working_stress = 234586466 # Lowest out of Ultimate Tensile Strength/1.65 and Yield Strength/1.33
        endcap_ellipsoid_ratio = 1.395
        tank_pressure = 230000

        # Calculate total tank volume required
        total_tank_volume = (self.fuel)/hydrogen_density


        # Endcaps -------------------------------------------------------------
        # Endcap Dimenions
        endcap_ellipsoid_height = self.radius/ellipsoidal_ratio
        endcap_crown_radius = ellipsoidal_ratio*self.radius
        eccentricity = math.pow((math.pow(self.radius,2)-math.pow(endcap_ellipsoid_height,2)),0.5)/self.radius

        # Calculate Single Endcap Volume
        endcap_volume = (2*math.pi*math.pow(self.radius,2)*endcap_ellipsoid_height)/3

        # Calculate Endcap thicknesses
        knuckle_thickness = (stress_factor*tank_pressure*self.radius)/(max_working_stress*welding_efficiency)
        crown_thickness = (tank_pressure*endcap_crown_radius)/(2*max_working_stress*welding_efficiency)
        equivalent_wall_thickness = (knuckle_thickness+crown_thickness)/2

        # Endcap Surface Area
        endcap_surface_area = math.pow(self.radius,2) + ((math.pi*math.pow(endcap_ellipsoid_height,2)*math.log((1+eccentricity)/(1-eccentricity))))/(2*e1)
        endcap_insulation_mass = endcap_surface_area*(0.78+0.015*self.mli_layers)

        design_factor = 2*endcap_ellipsoid_ratio + (1/math.pow((math.pow(endcap_ellipsoid_ratio,2)-1),0.5))*math.log((endcap_ellipsoid_ratio+math.pow(math.pow(endcap_ellipsoid_ratio,2)-1,0.5))/(endcap_ellipsoid_ratio-math.pow(math.pow(endcap_ellipsoid_ratio,2)-1,0.5)))

        endcap_mass = (math.pi*math.pow(self.radius,2)*equivalent_wall_thickness*design_factor*tank_material_density)/(2*endcap_ellipsoid_ratio) + endcap_insulation_mass

        # Cylinder ------------------------------------------------------------
        cylinder_volume = total_tank_volume - 2*endcap_volume
        cylinder_length = cylinder_volumel/(math.pi*math.pow(self.radius,2))
        cylinder_thickness = (tank_pressure*self.radius)/(max_working_stress*welding_efficiency)
        cylinder_surface_area = 2*math.pi*self.radius*cylinder_length
        cylinder_insulation_mass = cylinder_surface_area*(0.78+0.015*self.mli_layers)
        cylinder_mass = 2*math.pi*self.radius*cylinder_length*cylinder_thickness*tank_material_density + cylinder_insulation_mass

        self.structural_mass += endcap_mass + cylinder_mass
        self.tanklength += cylinder_length + 2*endcap_ellipsoid_height
        self.parent.mass += self.structural_mass



