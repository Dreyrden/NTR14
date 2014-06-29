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
        self.length = 0.0
        self.radius = 0.0
        self.structural_mass = 150.0
        self.fuel = 1000.0
        #  add the parent object as an attribute
        #+ and update the parent's attributes
        self.parent = parent
        parent.mass += self.structural_mass + self.fuel
        parent.fuel += 1000.0

    def add_fuel(self, fuel_to_add):
        self.fuel   += fuel_to_add
        self.parent.mass += fuel_to_add
        self.parent.fuel += fuel_to_add

    def remove_fuel(self, fuel_to_remove):
        self.fuel -= fuel_to_remove
        self.parent.mass -= fuel_to_remove
        self.parent.fuel -= fuel_to_remove