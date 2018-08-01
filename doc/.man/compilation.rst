*************
Changelog
*************

New functionality present in version 3.4.0 of ProtoMS.

===========
GCAP simulations
===========

Added functionality to perform grand canonical alchemical perturbation (GCAP) simulations, where ligand perturbations can be performed with GCMC water sampling. These can be set up using protoms.py -s gcap_single or protoms.py -s gcap_dual. These can either be performed at a range of B values, known as surface-GCAP (two-dimensional simulation in B and lambda) or at the equilibrium B value, which is single-GCAP. These new simulation methods come with a new analysis tool ``calc_gcap_surface.py`` which is able to calculate the free energy from the two-dimensional GCAP surface.

===========
Grand canonical integration
===========

More automation has been introduced to ``calc_gci.py`` to calculate the appropriate number of steps for the titration fit. Water occupancy at the equilibrium B value is returned. 

===========
Free energy calculation
===========

Updates have been made to free energy calculations. ``calc_ti_decomposed.py`` allows the different contributions to the free energy changes to be considered. This is useful to understand where changes in relative free energies are arising from. ``calc_dg_cycle.py`` allows the cycle closure and free energy changes of various legs to be determined for a set of free energy perturbations. 

===========
Clustering
===========

``calc_clusters.py`` has been changed to prevent two waters observed in the same frame being placed in the same frame. This previously allowed for cluster occupancies of greater than 100%.
