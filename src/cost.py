#  owner: rb

def get_attributes(core):
    """
       [dict, dict] = get_attributes(Class)
       
       find what attributes the core() class currently has and set equal
       to the relevant dictionary key, otherwise give the dictionary key 
       a default value
    """

    #  dictionary for parameters and their default values
    parambook = {'enrichment': 93.0,
                'UO2': 0.0,
                'W': 0.0,
                'Graphite': 0.0,
                'ZrH1d8': 0.0,
                'Zirc4': 0.0,
                'ZrC': 0.0,
                'Gd203': 0.0,
                'ThO2': 0.0,
                'Re': 0.0,
                'Mo': 0.0,
                'BeO': 0.0,
                'Be': 0.0,
                'W_internal_shield': 0.0}
    #  dictionary for parameter default cost values
    rawcostbook = {'enrichment': 0.0,
                'UO2': 0.0,
                'W': 0.0,
                'Graphite': 2.0,
                'ZrH1d8': 132.0,
                'Zirc4': 132.0,
                'ZrC': 100.0,
                'Gd203': 146.0,
                'ThO2': 150.0,
                'Re': 1222.0,
                'Mo': 19.0,
                'BeO': 1195.0,
                'Be': 1195.0,
                'W_internal_shield': 7.0}
    
    #  search through the members of core() and if one matches a parambook{}
    #+ key, then assign the value from core.member to parambook{key}
    for member in core.__dict__:
        print member
        for key in parambook:
            print key
            if member == key: parambook[key] = core.__dict__[member]

    #  return the dictionary
    return(parambook, rawcostbook)


def reactor_tfu(core):

    #  pull out necessary attributes
    parambook, rawcostbook = get_attributes(core)
    
    #  dictionary to keep track of costs
    tfubook = {}

    #  constants or pseudo-constant (i.e. periodically adjusted) parameters
    #  UO2 to cover costs including
    #+  a.) fixed cost for facilities overhead (electricity, insurance, etc.)
    #+  b.) fixed cost for government oversight for the process
    #+  c.) cost for security at the final product processing
    #+  d.) handling / storage cost at the end-step processing facility
    #+  e.) overhead cost for the number of full-time engineering, test,
    #+      manufacturing that will be employed even when not working on steps
    
    """ 
        UO2_raw is a cost specific function parameterized on UO2 enrichment
        with a linear fit on the domain [19.75, 90.0]
    """
    #  Cost per Kilogram ($US2013 / kg) of RAW 19.75% LEU & 90% HEU
    leu_raw = 8478.0
    heu_raw = 40970.0
    slope = (heu_raw - leu_raw) / (90.0 - 19.75)
    rawcostbook['UO2'] = slope*(parambook['enrichment'] - 19.75) + leu_raw
    rawcostbook['W'] = rawcostbook['UO2'] / 3.0

    #  calculate the material cost for the reactor
    TFU = 0.0
    for key in parambook:
        tfubook[key] = parambook[key]*rawcostbook[key]
        TFU += tfubook[key]
    
    #  miscellaneous core cost from 'parametric model'
    misc = 0.0
    #+ (cost per kg, mass)
    #  W-Can Support
    misc += 7.0 * 5.0
    #  AL2024 PV(s)
    misc += 2.0 * 155.4
    #  Control Drums (#) B4C 120-deg segment
    num_control_drums = 10
    misc += 10 * 1000.0 * 149.2
    #  Control Sensors, Electrical Harness/Wiring
    misc += 10000.0 * 4.0
    #  Structure and Plenum and Duct Flanges
    misc += 5000.0 * 40.0
    tfubook['misc'] = misc
    TFU += misc

    #  margin: [0, +inf)
    margin = 0.3
    TFU = (1 + margin)*TFU

    #  discount tfu cost based on prior development efforts
    #+ (i.e. learned out)
    #  DDTE of Design: Low (0.4) "Minor" --> High (1.4) "New"
    DDTE_scaling = 1.0
    TFU *= DDTE_scaling
    #  Team Engineering Experience: Low (1.4) --> High (0.7)
    Experience_scaling = 1.0
    TFU *= Experience_scaling

    #  store the TFU and return the tfubook
    tfubook['TFU'] = TFU
    return tfubook




#  a function for calculating the mission cost
def calculate_cost(mission):
    
    #  the cost of the mission is added as
    #+ an attribute to the mission
    mission.cost = 0.5