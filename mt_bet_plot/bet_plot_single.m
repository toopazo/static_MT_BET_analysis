function bet_plot_single(bet_st, str1, str2, legend_cell)

%    persistent cnt
%    if isempty(cnt)
%        cnt = 0;
%    end
%    if reset_cnt
%        cnt = 0;
%    end
%    
%    cnt = cnt + 1;
%    switch cnt
%        case 1
%            linspec = '-b*';
%        case 2
%            linspec = '-r*';
%        case 3
%            linspec = '-y*';
%        otherwise
%            disp('cnt is out of bound')
%            cnt
%            linspec = NaN;
%    end
    linspec = '-b*';

    % Unpack bet_forces results
    blade_st        = bet_st.blade_st;
    psi_arr         = bet_st.psi_arr;
    r_arr           = bet_st.r_arr;
    T_arr           = bet_st.T_arr;
    Q_arr           = bet_st.Q_arr;
    P_arr           = bet_st.P_arr;
    lambda_i_arr    = bet_st.lambda_i_arr;    
    alpha_arr       = bet_st.alpha_arr;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % 1) Calculate net produced Thrust and Torque
    Tr_arr      = bet_forces_mean_along_dpsi(T_arr);
    Qr_arr      = bet_forces_mean_along_dpsi(Q_arr);
    lir_arr     = bet_forces_mean_along_dpsi(lambda_i_arr);    
    aoar_arr    = bet_forces_mean_along_dpsi(alpha_arr);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 1;
    fig = figure(nfig);
    hold on;
    grid on;
    plot(r_arr, Tr_arr, '-*');
    arg = [ 'Mean T(r) distribution about \psi [0, 2 pi]' ]; %, total T = ' num2str(T) ];
    title(arg);
    xlabel('Radius');
    ylabel('Thrust N');
    
    omegarpm = blade_st.omega*30/pi;
    Ttot = round(bet_st.total.T, 2);
    text(0.05, Ttot*1.2, ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N'])           
    
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 2;
    fig = figure(nfig);
    hold on;
    grid on;
    plot(r_arr, Qr_arr, '-*');
    arg = [ 'Mean Q(r) distribution about \psi [0, 2 pi]' ]; %, total Q = ' num2str(Q) ];
    title(arg);
    xlabel('Radius');
    ylabel('Torque Nm');
    
    omegarpm = blade_st.omega*30/pi;
    Qtot = round(bet_st.total.Q, 2);
    text(0.05, Qtot*1.2, ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
    
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);          
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 3;
    fig = figure(nfig);
    hold on;
    grid on;
    plot(r_arr, lir_arr, '-*');
    arg = [ 'Mean \lambda_i(r) distribution about \psi [0, 2 pi]' ]; %, total Q = ' num2str(Q) ];
    title(arg);
    xlabel('Radius');
    ylabel('Induced inflow');

    omegarpm = blade_st.omega*30/pi;
    lambda_i = round(bet_st.total.lambda_i, 4);
    text(0.3, lambda_i, ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);

    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(4);
    hold on;    
    [X, Y] = meshgrid(r_arr, rad2deg(psi_arr));    
    surf(X, Y, lambda_i_arr);        
    surf(X, Y, lambda_i_arr.*0); %), 'FaceAlpha',0.5,'EdgeColor','none')  % plot zero surface
    arg = [ '\lambda_i(\psi, r) distribution' ];
    title(arg);
    xlabel('Radius');
    ylabel('\psi deg');
    zlabel('induced inflow');
    view(35, 35);
    grid on;     

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 5;
    fig = figure(nfig);
    hold on;    
    grid on;       
    [X, Y] = plot3D_polar(nfig, linspec, psi_arr, r_arr, lambda_i_arr);       
    arg = [ '\lambda_i(\psi, r) distribution' ];
    title(arg);
    xlabel('x axis');
    ylabel('y axis');
    zlabel('induced inflow');  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 6;
    fig = figure(nfig);
    hold on;  
    grid on;  
    [X, Y] = plot3D_polar(nfig, linspec, psi_arr, r_arr, T_arr);    
    arg = [ 'Thrust(\psi, r) distribution' ];
    title(arg);
    xlabel('x axis');
    ylabel('y axis');
    zlabel('Thrust');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 7;
    [dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr] = bet_forces_dolM(bet_st);    
    bet_plot_dolM(nfig, linspec, psi_arr, dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig1 = 8;
    nfig2 = 9;
    blade_plot(nfig1, nfig2, blade_st);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(10);
    hold on;
    grid on;
    plot(r_arr, rad2deg(aoar_arr), '-*');
    arg = [ 'Mean aoa(r) distribution about \psi [0, 2 pi]' ]; %, total T = ' num2str(T) ];
    title(arg);
    xlabel('Radius');
    ylabel('Angle of attack deg');       

    % Save figures
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nfig = 4;
    fig = figure(nfig);
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 5;
    fig = figure(nfig);
    % legend(legend_cell, 'Location','northwest');    
    plot3D_rotorAxes(nfig);       
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 6;
    fig = figure(nfig);
    % legend(legend_cell, 'Location','northwest');    
    plot3D_rotorAxes(nfig);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename); 
     
    nfig = 7;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northwest');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
    nfig = 8;
    fig = figure(nfig);
    legend(legend_cell, 'Location','northeast');
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);

    nfig = 9;
    fig = figure(nfig);
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);    
    
    nfig = 10;
    fig = figure(nfig);
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
%    nfig = 11;
%    fig = figure(nfig);
%    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
%    saveas(fig, filename);     
                 
end
