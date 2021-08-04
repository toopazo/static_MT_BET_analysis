function bet_main_omega()

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
    
    rotortype   = 'KDECF245DP';
    % rotortype   = 'KDECF245TP';
    % rotortype   = 'KDECF245HP';

    % omegarpm_arr = [2000, 2500, 3000, 3500, 3890, 4000, 4038, 4500];
    % omegarpm_arr = [2000, 2500, 3000, 3500, 4000, 4500];
    omegarpm_arr = linspace(2000, 4500, 10);
    npoints = length(omegarpm_arr);

    % T and Q vs RPM for anecha and BET
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:npoints
        omega   = omegarpm_arr(i) * rpm2rads;
        fprintf('[bet_main_omega] omega %d RPM \n', omegarpm_arr(i));
    
        % rotortype   = 'KDE6213XF185_KDECF305DP';
        % rotortype   = 'KDE6213XF185_KDECF245DP';
        % rotortype   = 'KDE6213XF185_KDECF245TP';
        % modeltype   = 'interp1';
        % modeltype   = 'quadratic';
        % [kde_T, kde_Q, kde_P] = kde_rotor_TQP(omega, rotortype, modeltype);
        % kde_T_arr(i) = kde_T;
        % kde_Q_arr(i) = kde_Q;
        
        % Anecha experimental results
        modeltype   = 'quadratic';
        [ane_T, ane_Q, ane_P] = anecha_rotor_TQP(omega, rotortype, modeltype);
        ane_T_arr(i) = ane_T;
        ane_Q_arr(i) = ane_Q;

        % BET results
        [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype); 
        omega       = omega;
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_st    = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);           
        bet_st = bet_forces(blade_st);
        bet_st = bet_forces_add_total(bet_st, false);
        bet_T_arr(i) = bet_st.total.T;           
        bet_Q_arr(i) = bet_st.total.Q;       
        
        % Other results
        rhoAvtip2_arr(i) = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;
        rhoAvtip3_arr(i) = blade_st.rho * blade_st.rotArea * blade_st.Vtip^3;
    end    
    
    % RMS error
    % err_T_arr = bet_T_arr - kde_T_arr; 
    % err_Q_arr = bet_Q_arr - kde_Q_arr; 
    err_T_arr = bet_T_arr - ane_T_arr;
    err_Q_arr = bet_Q_arr - ane_Q_arr;
    
    mean_err_T = mean(err_T_arr);
    mean_err_Q = mean(err_Q_arr);
    std_err_T = std(err_T_arr); % sqrt(mean( err_T_arr.*err_T_arr ) )
    std_err_Q = std(err_Q_arr); % sqrt(mean( err_Q_arr.*err_Q_arr ) )
    
    fprintf('stdev err_T %.4f \n', std_err_T);
    fprintf('stdev err_Q %.4f \n', std_err_Q);
    fprintf('mean err_T %.4f \n', mean_err_T);
    fprintf('mean err_Q %.4f \n', mean_err_Q);
    
    perr_T_arr = (err_T_arr) ./ ane_T_arr .* 100;
    perr_Q_arr = (err_Q_arr) ./ ane_Q_arr .* 100;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot
    st.omegarpm_arr = omegarpm_arr;
    st.ane_T_arr = ane_T_arr;
    st.ane_Q_arr = ane_Q_arr;
    st.bet_T_arr = bet_T_arr;
    st.bet_Q_arr = bet_Q_arr;
    st.err_T_arr = err_T_arr;
    st.err_Q_arr = err_Q_arr;
    st.mean_err_T = mean_err_T;
    st.mean_err_Q = mean_err_Q;
    st.std_err_T = std_err_T;
    st.std_err_Q = std_err_Q;
    st.perr_T_arr = perr_T_arr;
    st.perr_Q_arr = perr_Q_arr;
    str1 = 'bet_main_omega';
    str2 = rotortype;
    bet_main_omega_plot(st, str1, str2);

    % Save data to file
    table_thrust = bet_T_arr;
    table_omega = omegarpm_arr .* rpm2rads;
    table_power = bet_Q_arr .* table_omega; 
    table_torque = bet_Q_arr;   
    table_rhoAvtip2 = rhoAvtip2_arr;
    table_rhoAvtip3 = rhoAvtip3_arr;
    T = table(...
        table_thrust(:)     , ...
        table_omega(:)      , ...
        table_power(:)      , ...
        table_torque(:)     , ...
        table_rhoAvtip2(:)  , ...
        table_rhoAvtip3(:)    ...
    );
    writetable(T, ['img/bet_main_omega_' rotortype '.txt']);
    % writetable(struct2table(table_coax), 'img/table_txt')
end

