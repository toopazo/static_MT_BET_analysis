function ...
    [...
        db            ...
    ] = bet_coax_eta_thrust(...
        coax_thrust , ... 
        rotortype   , ...  
        coaxu_r0    , ... 
        coaxl_r0      ...
    )

    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Given:')
    disp('  constant total thrust')
    disp('  theta0_u = collpitch_u')
    disp('  theta0_l = collpitch_l')
    disp('  T_u      = W*(eta_T)/(1+eta_T)')
    disp('  T_l      = W - T_u')
    % disp('  eta_T    = eta_T_arr(i)')
    disp('Find:')
    disp('  omega_u and omega_l')     
    
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;    
                       
    [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
    omega       = 3000 * rpm2rads;
    CT_target   = NaN;      % CT_target is only needed for mu != 0

    collpitch0  = collpitch;
    
    % filename = ['img/bet_coax_eta_T_arr_' rotortype '.mat'];
    filename = [
        'img/bet_coax_eta_thrust' ...
        '_' rotortype  ...
        '_' num2str(coax_thrust) ...
        '_' num2str(coaxu_r0) '_' num2str(coaxl_r0) ...
        '.mat' ...
        ];
    
    % Hover case for KDECF305DP
    % omega       = 3890 * rpm2rads;
    % db.blade_h_st   = blade_model(...
    %     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
    
    % Hover case for KDECF245DP
    db.blade_h_st   = blade_model(...
        rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
    T_target        = coax_thrust/2;
    db.blade_h_st   = bet_find_omega(db.blade_h_st, T_target);
        
    db.bet_h_st     = bet_forces(db.blade_h_st);        
    db.bet_h_st     = bet_forces_add_total(db.bet_h_st, false);
    db.T_h          = db.bet_h_st.total.T; 
    db.Q_h          = db.bet_h_st.total.Q; 
    db.P_h          = sqrt(2)*2*db.bet_h_st.total.P;
    db.W            = 2*db.T_h;
    db.P_0          = db.W/(2*db.blade_h_st.rho*db.blade_h_st.rotArea);           

    db.coaxu_r0     = coaxu_r0;
    db.coaxl_r0     = coaxl_r0;

    % db.dcollpitch_arr  = -1.0:0.5:+2.0;
    db.dcollpitch_arr  = [0];
    % db.eta_T_arr       = 0.5:0.1:4;
    db.eta_T_arr       = 0.5:0.5:4;
    
    % db.blade_h_st
    % db

    for i = 1:length(db.dcollpitch_arr)
        dcollpitch  = db.dcollpitch_arr(i)
        collpitch_u = deg2rad( collpitch0 - dcollpitch );
        collpitch_l = deg2rad( collpitch0 + dcollpitch );
    
        db.collpitch_u_arr(i) = collpitch_u;
        db.collpitch_l_arr(i) = collpitch_l;
    
        % Get blade model
        blade_u_st = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch_u, omega, CT_target);
        blade_l_st = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch_l, omega, CT_target);

        % Set blades for coax and sbs cases
        bladecoax_u_st  = blade_u_st;
        bladecoax_l_st  = blade_l_st; 
        bladesbs_u_st   = blade_u_st;
        bladesbs_l_st   = blade_l_st;

        % Add coaxial r0 parameters to both coax blades
        bladecoax_u_st.coaxu_r0 = coaxu_r0;
        bladecoax_u_st.coaxl_r0 = coaxl_r0;
        bladecoax_l_st.coaxu_r0 = coaxu_r0;
        bladecoax_l_st.coaxl_r0 = coaxl_r0;

        for j = 1:length(db.eta_T_arr)
            
            % Calculate desired upper and lower thrust
            eta_T       = db.eta_T_arr(j)
            T_u_target  = db.W*(eta_T)/(1+eta_T);
            T_l_target  = db.W - T_u_target;  
                                
            % Find RPM that matches desired upper and lower thrust
            [bladecoax_u_st, bladecoax_l_st] = bet_coax_match_thrust(...
                bladecoax_u_st, bladecoax_l_st, T_u_target, T_l_target);
            omegacoax_u = bladecoax_u_st.omega;
            omegacoax_l = bladecoax_l_st.omega;                                  
            
            % Calculate coaxial loads
            [bet_u_st, bet_l_st] = bet_coax_forces(...
                bladecoax_u_st, bladecoax_l_st);

            % Save coax results to db
            db.coax.bet_u_st(i, j) = bet_u_st;
            db.coax.bet_l_st(i, j) = bet_l_st;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find RPM that matches desired upper and lower thrust
            [bladesbs_u_st, bladesbs_l_st] = bet_sbs_match_thrust(...
                bladesbs_u_st, bladesbs_l_st, T_u_target, T_l_target);
            omegasbs_u = bladesbs_u_st.omega;
            omegasbs_l = bladesbs_l_st.omega;            
            
            % Calculate side-by-side loads
            [bet_u_st, bet_l_st] = bet_sbs_forces(...
                bladesbs_u_st, bladesbs_l_st);
            
            % Save sbs results to db
            db.sbs.bet_u_st(i, j) = bet_u_st;
            db.sbs.bet_l_st(i, j) = bet_l_st;
        end
    end    
end
