function [blade_u_st, blade_l_st] = bet_coax_match_thrust(...
    blade_u_st, blade_l_st, T_u_target, T_l_target)

    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;  

    disp('[bet_coax_match_thrust] Running fsolve ...')

    options = optimoptions(@fsolve);
    options = optimoptions(options, ...
        'Display','off', 'TolFun', 10^-5, 'OptimalityTolerance', 10^-10');
    % options = optimoptions(options, ...
    %     'Display','off', 'TolFun', 10^-2, 'OptimalityTolerance', 10^-4');

    % Special case, zero upper thrust => zero upper rpm
    if T_u_target == 0
        xguess = [blade_l_st.omega];
        x0 = fsolve(@funzero_lower, xguess, options);
        x0 = [0; x0];
    end

    % Special case, zero lower thrust => zero lower rpm
    if T_l_target == 0
        xguess = [blade_u_st.omega];
        x0 = fsolve(@funzero_upper, xguess, options);
        x0 = [x0; 0];
    end

    % General case
    if T_u_target ~= 0 && T_l_target ~= 0
        xguess = [blade_u_st.omega; blade_l_st.omega];
        x0 = fsolve(@funzero_coax, xguess, options);
        x0 = [x0(1); x0(2)];
    end

    % Calculate error
    err = funzero_coax(x0);
    err_u = err(1);
    err_l = err(2);
    
    % Check error magnitude
    if norm(err) > 10^(-4)
        err
        blade_u_st
        blade_l_st
        error('norm(err) > 10^-4')
    end
    
    % Apply changes
    blade_u_st.omega      = x0(1);
    blade_u_st.Vtip       = blade_u_st.omega*blade_u_st.R;
    blade_l_st.omega      = x0(2);
    blade_l_st.Vtip       = blade_l_st.omega*blade_l_st.R;  

    % Calculate BET loads
    [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st);                
    T_u = bet_u_st.total.T;
    Q_u = bet_u_st.total.Q;
    lambda_u = bet_u_st.total.lambda;
    T_l = bet_l_st.total.T;     
    Q_l = bet_l_st.total.Q;
    lambda_l = bet_l_st.total.lambda;

    fprintf('[bet_coax_match_thrust] omega_u %.4f, omega_l %.4f \n', ...
        blade_u_st.omega * rads2rpm, blade_l_st.omega * rads2rpm);
    fprintf('[bet_coax_match_thrust] lambda_u %.4f, lambda_l %.4f \n', ...
        lambda_u, lambda_l);  
    fprintf('[bet_coax_match_thrust] T_u %.4f, T_l %.4f \n', ...
        T_u, T_l);
    fprintf('[bet_coax_match_thrust] err_u %.4f, err_l %.4f \n', ...
        err_u, err_l);
                
    return

    function err = funzero_coax(x)
        omega_u = x(1);
        omega_l = x(2);

        [err_u, err_l] = bet_coax_match_thrust_fzero(...
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target);
        err = [err_u, err_l];
    end
    function err = funzero_upper(x)
        omega_u = x(1);
        omega_l = 0; % x(2);

        [err_u, err_l] = bet_coax_match_thrust_fzero(...
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target);
        err = [err_u];
    end
    function err = funzero_lower(x)
        omega_u = 0; % x(1);
        omega_l = x(1); % x(2);

        [err_u, err_l] = bet_coax_match_thrust_fzero(...
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target);
        err = [err_l];
    end
end
