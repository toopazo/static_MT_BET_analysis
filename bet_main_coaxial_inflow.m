function bet_main_coaxial_inflow()
    
    % % format short
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot');       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear all
    close all

    rotortype   = 'KDECF245DP';               
    % rotortype   = 'KDECF245TP';     
    % rotortype   = 'KDECF245HP';  
    % rotortype   = 'KDECF305DP';
    
    coax_thrust = 45;
    % coax_thrust = 40;

    coaxu_r0 = 10;     % Radius of lower rotor wake at upper rotor location
    coaxl_r0 = 0.8;     % Radius of upper rotor wake at lower rotor location

    str1 = 'bet_main_coaxial_inflow';
    str2 = rotortype;

    filename = [
        'img/'            ...
        str1              ...
        '_' rotortype     ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '.mat' ...
        ];   

    read_from_file = true;
    % read_from_file = false;
    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_thrust_lc_surface(...
            coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end

    bet_plot_thrust_lc_surface(str1, str2, db);
    
end
