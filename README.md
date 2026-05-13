# Microscale-Serotonin-Dynamics
This repository contains the code corresponding to the paper "Compartmental-reaction diffusion framework for microscale dynamics of extracellular serotonin in brain tissue", by M. Pelz, S. Janusonis, and G. Handy published in SIAM Life Sciences (2026).

The numerical simulations were completed in Julia and can be found in the subfolder 'run-numerical-simulatinos'. The file run_simuation_2D.jl contains the main numerical algorithm

The subfolder 'reproduce-paper-figures' used Matlab and compressed data files to reproduce all main figures of the paper. The comments at the top of these figures also contain instructions as to which Julia files (and in which order) to run to complete the corresponding numerical simulations from scratch.
