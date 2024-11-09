# EOS_for_rotating_neutron_star_analysis
The repository contains programs that implement the algorithms developed by Read et al. and Greif et al. for a piecewise polytropic EOS and a speed of sound parametrization EOS. Also, it contains some codes I used to analyze the output of the program rns from Stergioulas.

-eos_generatorsCS.py implements the algorithm for a speed of sound paramettrization EOS. It checks that the speed of sound does not pass the speed of light. One can freely modify the values from a1 to a5 to get different kinds of EOS from soft to stiff

-eos_generatorPP.py implements the algorithm for a piecewise polytropic EOS. It checks that the speed of sound does not pass the speed of light. One can freely modify the values of the adiabatic index and of the of the energy density intervals to get different kinds of EOS from soft to stiff

-rns_input_reader.py reads the output from rns to produce arrays with the quantities needed for the analysis of universal relations for a rotating neutron stars.


