function [bet_u_st, bet_l_st] = bet_sbs_forces(blade_u_st, blade_l_st)
    
    % Non-iterative claculation of coax rotor    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    blade_u_st.position     = 'upper';
%    blade_u_st.lambda_infl  = x;                       % Lower on upper rotor inflow  
    bet_u_st                = bet_forces(blade_u_st);
    bet_u_st                = bet_forces_add_total(bet_u_st, false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    blade_l_st.position     = 'lower';
%    blade_l_st.lambda_infl  = bet_u_st.total.lambda;   % Upper on lower rotor inflow
    bet_l_st                = bet_forces(blade_l_st);
    bet_l_st                = bet_forces_add_total(bet_l_st, false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
end

