*************
Change Log
*************

This section describess the new and updated functionality present in version 3.4.0 of ProtoMS.

===========================
GCAP Simulations
===========================

ProtoMS is now able to perform grand canonical alchemical perturbation (GCAP) simulations, where ligand perturbations can be performed with GCMC water sampling. These can be set up using ``protoms.py -s gcap_single`` or ``protoms.py -s gcap_dual``, for single and dual topology calculations, respectively. These can either be performed at a range of :math:`B` values, known as surface-GCAP (two-dimensional simulation in :math:`B` and :math:`\lambda`) or at the equilibrium :math:`B` value, which is single-GCAP. These new simulation methods come with a new analysis tool ``calc_gcap_surface.py`` which is able to calculate the free energy from the two-dimensional GCAP surface.

===========================
Grand Canonical Integration
===========================

More automation has been introduced to ``calc_gci.py`` to calculate the appropriate number of steps for the titration fit. The water occupancy at the equilibrium :math:`B` value is also returned. 

===========================
Free Energy Calculations
===========================

Updates have been made to the analysis of free energy calculations. ``calc_ti_decomposed.py`` allows the different contributions to the free energy changes to be extracted. This is useful to understand where changes in relative free energies are arising from. ``calc_dg_cycle.py`` allows the cycle closure and free energy changes of various legs to be determined for a set of free energy perturbations. 

Additionally, a 'partial' implementation of dual topology calculations is now available. If the change between two ligands in a relative free energy calculation is too large for single topology, but dual topology would not be appropriate, it is now possible to apply softcore potentials to only the changing groups, whilst leaving the common structure unchanged, reducing noise in the calculated free energy.

===========================
protomslib
===========================
The underlying class structures and many of the fundamental functions of the Python layer of ProtoMS have now been restructured into an importable module, ``protomslib``. This also allows for support of Python :math:`\geq` 3.5 (retaining support for Python 2.7) with the setup and analysis scripts.

===========================
Restarting Simulations
===========================
The writing of restart files in ProtoMS has now been significantly improved, such that manual editing of the ``.cmd`` file is no longer required to restart a simulation which was terminated early or did not complete. Now, the simulation can be restarted where it left off by simply re-running the original command (something like ``$PROTOMSHOME/protoms3 run.cmd``). ProtoMS will read the restart file written out and will know where to continue the simulation from.

===========================
Water Clustering
===========================

The clustering of GCMC waters with ``calc_clusters.py`` has been changed to prevent two waters observed in the same frame from being placed in the same cluster. This previously allowed for cluster occupancies of greater than 100%, and has now been corrected. The clustering should also now be less sensitive to the distance cutoff used.
