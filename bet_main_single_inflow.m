function bet_main_single_inflow()
    
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data'); 
    addpath('mt_bet_single_rotor');         
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear all
    close all

    rotortype   = 'KDECF245DP';  
    % rotortype   = 'KDECF245TP';     
    % rotortype   = 'KDECF245HP';  
    % rotortype   = 'KDECF305DP';

    tot_thrust  = 35;
    % coax_thrust = 40;

    str1 = 'bet_main_single_inflow';
    str2 = rotortype;
    filename = ['img/' str1 '_' str2];

    % read_from_file = true;
    read_from_file = false;
    if read_from_file
        st = load(filename);
        db = st.db;
    else
        db = bet_single_thrust_lc_surface(tot_thrust, rotortype);
        save(filename, 'db');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bet_plot_thrust_lc_surface(str1, str2, db);

end
