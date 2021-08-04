function bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0)

    % clear all
    % close all
    % format compact
    % format short
    % clc
    
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;       

    read_from_file = false;
    % read_from_file = true;
    
    tic
    % The program section to time.     

    % rotortype   = 'KDECF245DP';               
    % rotortype   = 'KDECF245TP';     
    % rotortype   = 'KDECF245HP';  
    % rotortype   = 'KDECF305DP';
    
    % coaxu_r0 = 5.0;     % Radius of lower rotor wake at upper rotor location
    % coaxl_r0 = 0.5;     % Radius of upper rotor wake at lower rotor location

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 35;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 40;
    end

    filename = [
        'img/bet_coax_eta_thrust' ...
        '_' rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '.mat' ...
        ];   

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end
    
    dcollpitch = 0;
    str1 = 'bet_coax_eta_thrust';
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 
    close all;          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 45;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 50;
    end

    filename = [
        'img/bet_coax_eta_thrust' ...
        '_' rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '.mat' ...
        ];   

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end

    dcollpitch = 0;
    str1 = 'bet_coax_eta_thrust';
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 
    close all;                           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 55;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 60;
    end

    filename = [
        'img/bet_coax_eta_thrust' ...
        '_' rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '.mat' ...
        ];   

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end

    dcollpitch = 0;
    str1 = 'bet_coax_eta_thrust';
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 
    close all;                           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    toc
end
