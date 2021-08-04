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

% Radius of lower rotor wake at upper rotor location (plane)
coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10];
% Radius of upper rotor wake at lower rotor location (plane)
coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0];
% Inflow velocity of lower rotor wake at upper rotor location (plane)
coaxu_l0_arr = [0.00];
% Inflow velocity of upper rotor wake at lower rotor location (plane)
coaxl_l0_arr = [0.04, 0.05, 0.06, 0.07, 0.08];
% Upper and lower rotor type 
rotortype   = 'KDECF245DP';

for coaxu_r0 = coaxu_r0_arr
    for coaxl_r0 = coaxl_r0_arr
        for coaxu_r0 = coaxu_l0_arr
            for coaxl_r0 = coaxl_l0_arr
                coaxu_r0
                coaxl_r0
                coaxu_l0
                coaxl_l0
                % Calculate results
                bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0, coaxu_l0, coaxl_l0);
            end
        end
    end
end

% Radius of lower rotor wake at upper rotor location (plane)
coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10];
% Radius of upper rotor wake at lower rotor location (plane)
coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0];
% Inflow velocity of lower rotor wake at upper rotor location (plane)
coaxu_l0_arr = [0.00];
% Inflow velocity of upper rotor wake at lower rotor location (plane)
coaxl_l0_arr = [0.04, 0.05, 0.06, 0.07, 0.08];
% Upper and lower rotor type 
rotortype   = 'KDECF245TP';

for coaxu_r0 = coaxu_r0_arr
    for coaxl_r0 = coaxl_r0_arr
        for coaxu_r0 = coaxu_l0_arr
            for coaxl_r0 = coaxl_l0_arr
                coaxu_r0
                coaxl_r0
                coaxu_l0
                coaxl_l0
                % Calculate results
                bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0, coaxu_l0, coaxl_l0);
            end
        end
    end
end
