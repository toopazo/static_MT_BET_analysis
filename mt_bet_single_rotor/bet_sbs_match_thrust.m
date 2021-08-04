function [blade_u_st, blade_l_st] = bet_sbs_match_thrust(...
    blade_u_st, blade_l_st, T_u_target, T_l_target)

    disp('[bet_sbs_match_thrust] Running fsolve ...')

    options = optimset('Display','off');
    xguess = [blade_u_st.omega; blade_l_st.omega];
    x0 = fsolve(@funzero, xguess, options);
    
    err = funzero(x0);
    if norm(err) > 10^(-4)
        blade_u_st
        blade_l_st
        error('norm(err) > 10^-4')
    end      
    fprintf('[bet_sbs_match_thrust] err_u %.4f, err_l %.4f \n', err(1), err(2));
    
    % Apply changes and return
    blade_u_st.omega      = x0(1);
    blade_u_st.Vtip       = blade_u_st.omega*blade_u_st.R;
    blade_l_st.omega      = x0(2);
    blade_l_st.Vtip       = blade_l_st.omega*blade_l_st.R;     

    function err = funzero(x)
        blade_u_st.omega      = x(1);
        blade_u_st.Vtip       = blade_u_st.omega*blade_u_st.R;
        blade_l_st.omega      = x(2);
        blade_l_st.Vtip       = blade_l_st.omega*blade_l_st.R;        
        
        [bet_u_st, bet_l_st] = bet_sbs_forces(blade_u_st, blade_l_st);
        T_u = bet_u_st.total.T;
        Q_u = bet_u_st.total.Q;
        T_l = bet_l_st.total.T;        
        Q_l = bet_l_st.total.Q;
                
        err_u = T_u - T_u_target;
        err_l = T_l - T_l_target;
        err = [err_u; err_l];
    end
end
