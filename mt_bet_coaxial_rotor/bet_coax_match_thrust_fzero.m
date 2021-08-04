function [err_u, err_l] = bet_coax_match_thrust_fzero(...
        omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)

    % if omega_u < 0
    %     error('omega_u < 0')
    % end
    % if omega_l < 0
    %     error('omega_l < 0')
    % end
    omega_u = sqrt(omega_u * omega_u);
    omega_l = sqrt(omega_l * omega_l);

    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;      

    blade_u_st.omega      = omega_u;
    blade_u_st.Vtip       = blade_u_st.omega * blade_u_st.R;
    blade_l_st.omega      = omega_l;
    blade_l_st.Vtip       = blade_l_st.omega * blade_l_st.R;         
    
    % Calculate BET loads
    [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st);                
    T_u = bet_u_st.total.T;
    Q_u = bet_u_st.total.Q;
    lambda_u = bet_u_st.total.lambda;
    T_l = bet_l_st.total.T;     
    Q_l = bet_l_st.total.Q;
    lambda_l = bet_l_st.total.lambda;
            
    % extra_err_u = 0;
    % extra_err_l = 0;
    % if x(1) < 0
    %     extra_err_u = x(1) * x(1) * 100000;
    % end
    % if x(2) < 0
    %     extra_err_l = x(2) * x(2) * 100000;
    % end
    % err_u = T_u - T_u_target + extra_err_u;
    % err_l = T_l - T_l_target + extra_err_l;
    err_u = T_u - T_u_target;
    err_l = T_l - T_l_target;

    err = [err_u; err_l];

    % fprintf('[bet_coax_match_thrust] ------------------------------------ \n');
    % fprintf('[bet_coax_match_thrust] omega_u %.4f, omega_l %.4f \n', ...
    %     blade_u_st.omega * rads2rpm, blade_l_st.omega * rads2rpm);
    % fprintf('[bet_coax_match_thrust] lambda_u %.4f, lambda_l %.4f \n', ...
    %     lambda_u, lambda_l);  
    % fprintf('[bet_coax_match_thrust] T_u %.4f, T_l %.4f \n', ...
    %     T_u, T_l);
    % fprintf('[bet_coax_match_thrust] err_u %.4f, err_l %.4f \n', ...
    %     err_u, err_l);

    % lambda_c = bet_u_st.total.lambda_c
    % lambda_i = bet_u_st.total.lambda_i
    % lambda = bet_u_st.total.lambda
    % % rpm_u = x(1) * 30 / pi
    % rpm_l = x(2) * 30 / pi
    % % T_u
    % T_l
end