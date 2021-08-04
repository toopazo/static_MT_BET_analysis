function mt_bet_main()
    clear all
    clc
    close all
    format compact
    format short
    
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 

    % Hover
    bet_main_single
    
    bet_main_omega    

    % Hover coaxial
    
    mt_main_coaxial
    
    bet_main_coaxial
    
    return
    
    % Forward flight
    mt_main_advratio
    
    bet_main_advratio
    
    bet_main_NASALangleyWT
     
    
end

