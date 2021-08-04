function [blade_st] = bet_single_match_thrust(blade_st, thrust)

    disp('[bet_single_match_thrust] Running fsolve ...')

    options = optimset('Display','off');
    xguess = blade_st.omega;
    x0 = fsolve(@funzero, xguess, options);
    
    err = funzero(x0);
    if norm(err) > 10^(-4)
        blade_st
        error('norm(err) > 10^-4')
    end      
    fprintf('[bet_single_match_thrust] err %.4f \n', err);
    
    % Apply changes and return
    blade_st.omega    = x0(1);
    blade_st.Vtip     = blade_st.omega*blade_st.R;

    function err = funzero(x)
        blade_st.omega    = x; % sqrt(x*x);
        blade_st.Vtip     = blade_st.omega*blade_st.R;
        
        bet_st    = bet_forces(blade_st);
        bet_st    = bet_forces_add_total(bet_st, false);
        
        bet_thrust = bet_st.total.T;                
        err = bet_thrust - thrust;
        err = err;
    end
end
