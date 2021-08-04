function bet_main_single()

    clear all
    close all
    format compact
    % format short
    format long
    % format bank
    clc
    
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 
    
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi; 
    
    omega       = 3100 * rpm2rads;         

    % rotortype   = 'KDE6213XF185_KDECF305DP';
    rotortype   = 'KDE6213XF185_KDECF245DP';
    modeltype   = 'interp1';
    % modeltype   = 'quadratic';
    [kde_T, kde_Q, kde_P] = kde_rotor_TQP(omega, rotortype, modeltype);
    
    rotortype   = 'KDECF245DP';
    modeltype   = 'quadratic';
    [ane_T, ane_Q, ane_P] = anecha_rotor_TQP(omega, rotortype, modeltype);

    % rotortype   = 'KDECF305DP';
    rotortype   = 'KDECF245DP';               
    [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
    omega       = omega;
    CT_target   = NaN;              % CT_target is only needed for mu != 0
    blade_st    = blade_model(...
        rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);

    % Calculate single rotor loads
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    bet_st = bet_forces(blade_st);
    bet_st = bet_forces_add_total(bet_st, false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    str1        = 'bet_main_single';
    str2        = rotortype;
    legend_cell = {rotortype};
    bet_plot_single(bet_st, str1, str2, legend_cell);
                       
    close all;   
end

