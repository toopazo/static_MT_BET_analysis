function db = bet_single_thrust_lc_surface(tot_thrust, rotortype)
   
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;    

    % thrust_arr = [35, 40, 45, 50, 55, 60];
    thrust_arr = [tot_thrust];
    % lc_arr = 0:0.025:0.1;
    lc_arr = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
    for i = 1:length(thrust_arr)
        thrust = thrust_arr(i);
        for j = 1:length(lc_arr)           
            [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
            lambda_c = lc_arr(j);

            fprintf('[bet_single_thrust_lc_surface] lambda_c %d T %d \n', ...
            lambda_c, thrust);

            % Get rotor and initial guess of omega
            omega       = 5000 * rpm2rads;     ;
            CT_target   = NaN;      % CT_target is only needed for mu != 0
            blade_st    = blade_model(...
                rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);

            % Find omega that matches desired thrust
            blade_st = bet_single_match_thrust(blade_st, thrust);

            % Calculate BET
            bet_st = bet_forces(blade_st);
            bet_st = bet_forces_add_total(bet_st, false);

            % Average result for lower blade
            bet_omega_arr(i, j) = bet_st.blade_st.omega;
            bet_thrust_arr(i, j) = bet_st.total.T;
            bet_torque_arr(i, j) = bet_st.total.Q;
            bet_power_arr(i, j) = bet_st.total.P; 
            bet_lambdai_arr(i, j) = bet_st.total.lambda_i;
            bet_lambdac_arr(i, j) = bet_st.total.lambda_c;
            bet_lambda_arr(i, j) = bet_st.total.lambda;
            % Distribution along r
            bet_lcr_arr(i, j, :) = bet_st.total.lcr;
            bet_lir_arr(i, j, :) = bet_st.total.lir;
            bet_lr_arr(i, j, :) = bet_st.total.lr;
            bet_aoa_arr(i, j, :) = bet_st.total.aoar;

            % Calculate MT inflow for single rotor
            % mt_st = mt_single_axial_climb(blade_st, thrust, lambda_c);
            mt_st = mt_single_axial_climb(...
                bet_st.blade_st, ...
                bet_thrust_arr(i, j), ...   % Use blade thrust
                bet_lambdac_arr(i, j));     % Use blade lamda_c
            mt_lambda_arr(i, j) = mt_st.mt_lambda;
            mt_lambdai_arr(i, j) = mt_st.mt_lambda_i;
            mt_v_arr(i, j) = mt_st.mt_v;
            mt_power_arr(i, j) = mt_st.mt_power_i;
        end

        mt_lambda_h(i) = mt_lambdai_arr(i, find(lc_arr==0));
        mt_power_h(i) = mt_power_arr(i, find(lc_arr==0));
    end

    % db.r_arr = linspace(0, 1, length(bet_lr_arr(j, :)));
    db.thrust_arr = thrust_arr;
    db.lc_arr = lc_arr;
    
    db.bet_omega_arr = bet_omega_arr;
    db.bet_thrust_arr = bet_thrust_arr;
    db.bet_torque_arr = bet_torque_arr;
    db.bet_power_arr = bet_power_arr;
    db.bet_lambdac_arr = bet_lambdac_arr;
    db.bet_lambdai_arr = bet_lambdai_arr;
    db.bet_lambda_arr = bet_lambda_arr;

    db.bet_lcr_arr = bet_lcr_arr;
    db.bet_lir_arr = bet_lir_arr;
    db.bet_lr_arr = bet_lr_arr;
    db.bet_aoa_arr = bet_aoa_arr;
    
    db.mt_lambda_arr = mt_lambda_arr;
    db.mt_lambdai_arr = mt_lambdai_arr;
    db.mt_v_arr = mt_v_arr;
    db.mt_power_arr = mt_power_arr;
    
    db.mt_lambda_h = mt_lambda_h;
    db.mt_power_h = mt_power_h;
end


