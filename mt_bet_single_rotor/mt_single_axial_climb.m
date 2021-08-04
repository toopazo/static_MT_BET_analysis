function mt_st = mt_single_axial_climb(blade_st, thrust, lambda_c)

    rhoAvtip2 = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;
    % rhoAvtip3 = blade_st.rho * blade_st.rotArea * blade_st.Vtip^3;

    % MT inflow for single rotor
    Vtip = blade_st.omega * blade_st.R;
    CT = thrust / rhoAvtip2;
    mu = 0;
    TPP_alpha = 0;
    % lambda_c = 0;
    
    mt_st.mt_lambda = mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c);
    mt_st.mt_lambda_i = mt_st.mt_lambda - lambda_c;
    mt_st.mt_v = mt_st.mt_lambda * Vtip;
    mt_st.mt_power_i = 1.15 * thrust * mt_st.mt_v;  

    % if lambda_c == 0
    %     omega = blade_st.omega
    %     thrust
    %     mt_st.mt_lambda
    %     mt_st.mt_lambda_i
    %     mt_st.mt_v
    %     mt_st.mt_power_i
    % end
end


