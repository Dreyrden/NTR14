author: rb
date: 2014_06_22

the current code is a framework for the project. the majority of the structure
is there (i.e. outer-loop opt., inner-loop opt., mission and spacecraft classes, 
with modules in place for some subsystems of the spacecraft)

* the code can be run in a terminal using >> python ntrf14.py (the path to your
python may need to be set appropriately in ntrf.py). 
* you will also need to set the paths to your /src and /aerospace_toolbox directories
this can be upgraded in the near future using the python.sys library

the current code is only 
using the inner-loop opt and trivially finding the NSO altitude that reduces 
delta-v. as the cost and revenue functions are populated in the cost.py and 
revenue.py modules, then the inner-loop will optimmize for greatest mission 
profit. 

a run-down on the code structure,

1.) ntrf14.py setups up the outer-loop optimization problem and calls 
    mission.main()
2.) mission.main() instantiates a Mission class, which setups a spacecraft 
    and appropriate parameters for the inner-loop optimization. 
3.) the inner-loop optimizer calls the mission.optimize() method. it is within
    this method that the mission.execute(), cost(), and revenue() are all 
    called. 
    i.)  mission.execute(), executes the specific mission type
    ii.) cost() is owned by af
    iii.) revenue() is owned by rg

*.) other modules owned by rg and af will be attributes of the SpaceCraft class.
    an e.g. is the tanks.py modules, which has a parent Tank() class and child
    LH2Tank() class. These classes should have appropriate get, set, add, remove 
    methods that not only update themself, but also their parent SpaceCraft. 
    see SpaceCraft.add_tank() and tanks.Tank() for an example