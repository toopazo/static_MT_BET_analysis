clear all
close all
format compact
format short
clc

% Add functions to the path
addpath('mt_bet');  
addpath('mt_bet_data');       
addpath('mt_bet_coaxial_rotor');   
addpath('mt_bet_single_rotor');   
addpath('mt_bet_plot');       

% Radius of lower rotor wake at upper rotor location
coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10];
% Radius of upper rotor wake at lower rotor location
coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0];
% Upper and lower rotor type 
rotortype   = 'KDECF245DP';

for coaxu_r0 = coaxu_r0_arr
    for coaxl_r0 = coaxl_r0_arr
        coaxu_r0
        coaxl_r0
        % Calculate results
        bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0);
    end
end

% Radius of lower rotor wake at upper rotor location
coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10];
% Radius of upper rotor wake at lower rotor location
coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0];
% Upper and lower rotor type 
rotortype   = 'KDECF245TP';

for coaxu_r0 = coaxu_r0_arr
    for coaxl_r0 = coaxl_r0_arr
        % Calculate results
        bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0);
    end
end
