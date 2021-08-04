function bet_st = bet_forces_add_total(bet_st, verbose)
    
    % Unpack bet_forces results
    blade_st        = bet_st.blade_st;
    psi_arr         = bet_st.psi_arr;
    r_arr           = bet_st.r_arr;
    T_arr           = bet_st.T_arr;
    Q_arr           = bet_st.Q_arr;
    P_arr           = bet_st.P_arr;

    lambda_c_arr    = bet_st.lambda_c_arr;
    lambda_i_arr    = bet_st.lambda_i_arr;
    lambda_arr      = bet_st.lambda_arr;
    alpha_arr       = bet_st.alpha_arr;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Net produced thrust, torque and power
    % Net produced loads along azimuth
    Tr_arr      = bet_forces_mean_along_dpsi(T_arr);
    Qr_arr      = bet_forces_mean_along_dpsi(Q_arr);
    Pr_arr      = bet_forces_mean_along_dpsi(P_arr);
    
    % A blade produces lift at all r, but only during one psi at a time      
    % Therefore, the following is wrong
    %   T = sum(T_arr, 'all'); 
    %   Q = sum(Q_arr, 'all');
    % And the this is right
    T = sum(Tr_arr);    % Sum of avgerage of T(psi, r) along psi
    Q = sum(Qr_arr);    % Sum of avgerage of Q(psi, r) along psi    
    P = Q * blade_st.omega;

    % percT = T / blade_st.Thover * 100;
    % percQ = Q / ( blade_st.Pavail / blade_st.omega ) * 100;
    % percP = P / blade_st.Pavail * 100;
    
    CT = T / ( blade_st.rho * blade_st.rotArea * blade_st.Vtip^2 );
    CQ = Q / ( blade_st.rho * blade_st.rotArea * blade_st.Vtip^2 * blade_st.R );
    CP = CQ;    % CQ == CP    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inflow    
    lcr_arr     = bet_forces_mean_along_dpsi(lambda_c_arr);
    lir_arr     = bet_forces_mean_along_dpsi(lambda_i_arr);
    lr_arr      = bet_forces_mean_along_dpsi(lambda_arr);
    aoa_arr     = bet_forces_mean_along_dpsi(alpha_arr);
    
    % lambda_c    = bet_st.blade_st.lambda_c;
    lambda_c    = mean(lcr_arr);    % Avgerage of { avgerage of li(psi, r) along psi } along r
    lambda_i    = mean(lir_arr);    % Avgerage of { avgerage of li(psi, r) along psi } along r
    lambda      = mean(lr_arr);     % Avgerage of { avgerage of li(psi, r) along psi } along r
    % lambda      = lambda_c + lambda_i;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Net dissymmetry of lift moment
    [dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr] = bet_forces_dolM(bet_st);
    
    mean_dolMx = mean(dolMxpsi_arr);
    mean_dolMy = mean(dolMypsi_arr);
    mean_dolMz = mean(dolMzpsi_arr);
             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
    % Add to total bet_st
    bet_st.total.T = T; 
    bet_st.total.Q = Q;
    bet_st.total.P = P;
    
    % bet_st.total.percT = percT; 
    % bet_st.total.percQ = percQ;
    % bet_st.total.percP = percP;  
      
    bet_st.total.CT = CT;
    bet_st.total.CQ = CQ;      
    bet_st.total.CP = CP;   
    
    % Average along r
    bet_st.total.lambda_c   = lambda_c;
    bet_st.total.lambda_i   = lambda_i;
    bet_st.total.lambda     = lambda;

    % distribution along r
    bet_st.total.lcr        = lcr_arr;
    bet_st.total.lir        = lir_arr;
    bet_st.total.lr         = lr_arr;
    % bet_st.total.lcr        = lambda_c .* lir_arr ./ lir_arr;
    % bet_st.total.lcr(isnan(bet_st.total.lcr)) = 0;
    % bet_st.total.lr         = bet_st.total.lcr + bet_st.total.lir;
    bet_st.total.aoar       = aoa_arr;

    % dissymmetry of lift
    bet_st.total.mean_dolMx = mean_dolMx;
    bet_st.total.mean_dolMy = mean_dolMy;
    bet_st.total.mean_dolMz = mean_dolMz;   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % if blade_st.omega < 0
    %     blade_st.omega
    %     bet_st.total
    %     error('blade_st.omega < 0')
    % end
    % if bet_st.total.lambda_i < 0
    %     blade_st.omega
    %     bet_st.total
    %     error('bet_st.total.lambda_i < 0')
    % end
    % if bet_st.total.T < 0
    %     blade_st.omega
    %     bet_st.total
    %     error('bet_st.total.T < 0')
    % end

    % 4) Print
    if verbose    
        % Print blade type
        fprintf('  %s \n  ', blade_st.rotortype);
        
        % Angular speed
        fprintf('  Omega = %.2f rad/s, Vtip = %.2f m/s \n', ...
            blade_st.omega, blade_st.Vtip );
      
        % Print net produced Thrust and Torque 
        fprintf('  T %.2f N, Q %.2f Nm, P = %.2f \n', T, Q, P);
        % fprintf('  percThover = %.2f, percPhover = %.2f \n', percT, percP);    
        fprintf('  CT %.4f, CQ %.4f \n', CT, CQ);   
        
        % Print net dissymmetry of lift moment
        fprintf('  Mean dolM about psi [0, 2 pi] is [%.2f; %.2f; %.2f] \n', ...
            mean_dolMx, mean_dolMy, mean_dolMz);  
    end
end
