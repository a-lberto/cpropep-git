# Cpropep is based on the theory presented by Gordon and McBride
# in the NASA report RP-1311. You can download a pdf version of
# this document at http://www.arocket.net/library/

# The thermodynamics data file thermo.dat coma also from McBride
# at the NASA Gleen Research center.

# Here is an example of an input file to be use by cpropep.
# Any line beginning by a '#' a space or a new_line is considered
# as a comment.

# This file should first contain a section named 'Propellant' which
# contain a list of all substance contain in the propellant. The 
# number refer to an element in the data file containing propellant
# information. In order to have a list of the substance, you could
# invoque the program like that:  'cpropep -p'

# There is two units that are support for ingredient quantity g (gram) or m (mole)

Propellant HTPB/KClO4/Al
+976      10 g

# You could then specify a list of problem to be solve. There is 4
# possible cases:

# TP for temperature-pressure fixed problem
# You have to specify the temperature and the pressure (of course)
# There is 4 pressure units (psi, kPa, atm and bar) and 3 temperature units (k, c and f)

TP
+chamber_pressure    20000 kPa  
+chamber_temperature 673 k

# HP for enthalpy-pressure fixed problem. It use the enthalpy of
# the propellant describe at the beginning.

# Only the chamber pressure shoud be specified. The temperature of
# the product will be the adiabatic flame temperature.

#HP
#+chamber_pressure 20.4 atm   # 136 atm

# FR is used to compute frozen performance.
# You have to specify the chamber pressure and an exit condition.
# This condition could be one of the following three:

# exit_pressure:         pressure at the exit.
# supersonic_area_ratio: exit to throat area for an area after the nozzle
# subsonic_area_ratio:   exit to throat area for an area before any nozzle

#FR
#+chamber_pressure      20.4 atm #136 atm
#+exit_pressure         1   atm
#+supersonic_area_ratio 10
#+subsonic_area_ratio   5

# EQ is used to compute shifting equilibrium performance.
# The options are the same as for frozen.

#EQ
#+chamber_pressure      20.4 atm #136 atm
#+exit_pressure         1   atm
#+supersonic_area_ratio 10 
#+subsonic_area_ratio   5
