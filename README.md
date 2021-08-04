# steadyState_inflowModels
Software implementation of static (steady state) inflow model for a rotor in hover (BEMT) and forward flight (Glauert approximation). 

Based on Leishman's book and Blade Element Theory

Leishman, Gordon J. Principles of helicopter aerodynamics with CD extra. Cambridge university press, 2006.


## Main MT and BET files
bet_main_single.m
- Loads (thrust, torque, power, etc) for a given rotor ata given condition (eomga, rho, lambda_c, etc)

bet_main_omega.m
- Thrust and torque as a function of RPM

bet_main_advratio.m
- Inflow, induced inflow and aoa as a function of advance ratio
- Thrust and torque as a function of advance ratio

mt_main_advratio.m
- Inflow, induced inflow and aoa as a function of advance ratio
- Thrust and torque as a function of advance ratio

bet_main_coaxial.m
- Upper and lower results (thrust, rpm, torque, power, etc) as a funciton of eta_T 
- Takes about 2 hours to complete

mt_coax_power.m 
- Induced velocity ratio as a function of thrust ratio
- Single and coaxial power as a function of thrust ratio

bet_main_NASALangleyWT.m
- Comparison between WT, CFD and BET results for NASA Lengley WT tests

