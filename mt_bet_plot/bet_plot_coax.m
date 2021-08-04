function bet_plot_coax(bet_u_st, bet_l_st, str1, str2, legend_cell)
    
    bet_st = bet_u_st;
    % Unpack bet_forces results
    blade_st        = bet_st.blade_st;
    psi_arr         = bet_st.psi_arr;
    r_arr           = bet_st.r_arr;
    T_arr           = bet_st.T_arr;
    Q_arr           = bet_st.Q_arr;
    P_arr           = bet_st.P_arr;
    lambda_i_arr    = bet_st.lambda_i_arr;    
    alpha_arr       = bet_st.alpha_arr;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi    = bet_st.dpsi;
    dr      = bet_st.dr;
    % T_arr(psi, r) is a partition of the continuous thrust
    % and represents the contribution of the ENTIRE section dr 
    % whose center is localted in (psi, r)
    % A true distribution is therefore T_arr(psi, r)/dr
    % because it is "nsections" independent
    % The same is true for Q_arr and P_arr, but NOT for lambda_i_arr    
    psi_arr         = psi_arr;
    T_arr           = T_arr ./ dr;
    Q_arr           = Q_arr ./ dr;
    P_arr           = P_arr ./ dr;
    lambda_i_arr    = lambda_i_arr;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1) Calculate net produced Thrust and Torque
    Tr_arr      = bet_forces_mean_along_dpsi(T_arr);
    Qr_arr      = bet_forces_mean_along_dpsi(Q_arr);
    lir_arr     = bet_forces_mean_along_dpsi(lambda_i_arr);    
    aoar_arr    = bet_forces_mean_along_dpsi(alpha_arr);
        
    fig = figure(1);
    subplot(2, 1, 1)
    hold on;
    grid on;
    plot(r_arr, Tr_arr, '-*');
    arg = [ 'Mean T(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Thrust N');
    ylim([-50, 300]);   

    fig = figure(2);
    subplot(2, 1, 1)
    hold on;
    grid on;
    plot(r_arr, Qr_arr, '-*');
    arg = [ 'Mean Q(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Torque Nm');
    ylim([-5, 15]);   
    
    fig = figure(3);
    subplot(2, 1, 1)
    hold on;
    grid on;
    plot(r_arr, lir_arr, '-*');
    arg = [ 'Mean \lambda_i(r) distribution about \psi [0, 2 pi]' ]; 
    title(arg);
    xlabel('Radius');
    ylabel('Induced inflow');   
    ylim([-0.05, 0.10]);   
    
    fig = figure(4);
    subplot(2, 1, 1)
    hold on;
    grid on;
    plot(r_arr, rad2deg(aoar_arr), '-*');
    arg = [ 'Mean aoa(r) distribution about \psi [0, 2 pi]' ]; 
    title(arg);
    xlabel('Radius');
    ylabel('Angle of attack deg');       
    ylim([-60, 10]);
    
    %
    % Now, the same, but for the lower blade
    %
             
    bet_st = bet_l_st;
    % Unpack bet_forces results
    blade_st        = bet_st.blade_st;
    psi_arr         = bet_st.psi_arr;
    r_arr           = bet_st.r_arr;
    T_arr           = bet_st.T_arr;
    Q_arr           = bet_st.Q_arr;
    P_arr           = bet_st.P_arr;
    lambda_i_arr    = bet_st.lambda_i_arr;    
    alpha_arr       = bet_st.alpha_arr;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dpsi    = bet_st.dpsi;
    dr      = bet_st.dr;
    % T_arr(psi, r) is a partition of the continuous thrust
    % and represents the contribution of the ENTIRE section dr 
    % whose center is localted in (psi, r)
    % A true distribution is therefore T_arr(psi, r)/dr
    % because it is "nsections" independent
    % The same is true for Q_arr and P_arr, but NOT for lambda_i_arr    
    psi_arr         = psi_arr;
    T_arr           = T_arr ./ dr;
    Q_arr           = Q_arr ./ dr;
    P_arr           = P_arr ./ dr;
    lambda_i_arr    = lambda_i_arr;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1) Calculate net produced Thrust and Torque
    Tr_arr      = bet_forces_mean_along_dpsi(T_arr);
    Qr_arr      = bet_forces_mean_along_dpsi(Q_arr);
    lir_arr     = bet_forces_mean_along_dpsi(lambda_i_arr);    
    aoar_arr    = bet_forces_mean_along_dpsi(alpha_arr);
        
    fig = figure(1);
    subplot(2, 1, 2)
    hold on;
    grid on;
    plot(r_arr, Tr_arr, '-*');
    arg = [ 'Mean T(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Thrust N');
    ylim([-50, 300]);   

    fig = figure(2);
    subplot(2, 1, 2)
    hold on;
    grid on;
    plot(r_arr, Qr_arr, '-*');
    arg = [ 'Mean Q(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Torque Nm');
    ylim([-5, 15]);           
    
    fig = figure(3);
    subplot(2, 1, 2)
    hold on;
    grid on;
    plot(r_arr, lir_arr, '-*');
    arg = [ 'Mean \lambda_i(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Induced inflow');      
    ylim([-0.05, 0.10]);  
    
    fig = figure(4);
    subplot(2, 1, 2)
    hold on;
    grid on;
    plot(r_arr, rad2deg(aoar_arr), '-*');
    arg = [ 'Mean aoa(r) distribution about \psi [0, 2 pi]' ];
    title(arg);
    xlabel('Radius');
    ylabel('Angle of attack deg'); 
    ylim([-60, 10]);
    
    %
    % Now, save it all
    %                    
    
    nfig = 1;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 2;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 3;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 4;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northwest');    
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

end
