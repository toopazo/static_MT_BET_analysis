function Cd0 = bet_calculate_Cd0_c81table(blade_st)

    % Blade properties are assummed to be constant
    alpha       = blade_st.TPP_alpha;
    CT_target   = blade_st.CT_target;
    mu          = blade_st.mu;
    lambda_c    = blade_st.lambda_c;
    omega       = blade_st.omega;
    R           = blade_st.R;
    rho         = blade_st.rho;
    A           = blade_st.rotArea;
    Vtip        = blade_st.Vtip;
    Nb          = blade_st.Nb;
    Cd0         = blade_st.Cd0;
    sigma       = mean(blade_st.sigma_arr);

    % Calculate load according to BET
    bet_st = bet_forces(blade_st);
    
    y_arr = bet_st.blade_st.y_arr;
    dy    = y_arr(2) - y_arr(1);
    y_arr = y_arr + dy/2;
    y_arr = y_arr(1:end-1);
    L_arr = bet_st.L_arr;
    D_arr = bet_st.D_arr;
    
    % Net produced loads along azimuth
    Lr_arr      = bet_forces_mean_along_dpsi(L_arr);
    Dr_arr      = bet_forces_mean_along_dpsi(D_arr);    
    
    % Po = (1/8) * sigma * Cd0 * (rho * A * Vtip3);
    integral_Dy = sum(Dr_arr .* y_arr);
    Po = omega * Nb * integral_Dy;
    Cd0 = Po / ( (1/8) * sigma * (rho * A * Vtip^3) );
end

