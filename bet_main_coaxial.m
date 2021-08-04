function bet_main_coaxial(...
        rotortype   , ...
        coaxu_r0    , ...
        coaxl_r0    , ...
        coaxu_l0    , ...
        coaxl_l0      ...
    )

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
    
    str1 = 'bet_coax_eta_thrust';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 35;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 40;
    end

    str2 = [
        rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '_' num2str(coaxu_l0) '_' num2str(coaxl_l0) ...
    ];
    filename = ['img/' str1 '_' str2];

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end
    
    dcollpitch = 0;
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 
    close all;          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 45;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 50;
    end

    str2 = [
        rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '_' num2str(coaxu_l0) '_' num2str(coaxl_l0) ...
    ];
    filename = ['img/' str1 '_' str2];

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end

    dcollpitch = 0;
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 
    close all;                           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    if strcmp(rotortype, 'KDECF245DP')
        coax_thrust = 55;
    end
    if strcmp(rotortype, 'KDECF245TP')
        coax_thrust = 60;
    end

    str2 = [
        rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '_' num2str(coaxu_l0) '_' num2str(coaxl_l0) ...
    ];
    filename = ['img/' str1 '_' str2];

    if read_from_file
        st = load(filename);
        db = st.db;          
    else
        db = bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0);
        save(filename, 'db');
    end

    dcollpitch = 0;
    bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1, str2); 
    close all;                           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    toc
end
