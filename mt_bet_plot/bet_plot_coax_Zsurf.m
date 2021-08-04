function bet_plot_coax_Zsurf(rotortype, db)

    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;  
    
    % [red green blue] matrix
    C = transpose(db.Pcoax_arr);
    for i = 1:size(C, 1)
        for j = 1:size(C, 2)
           Cr(i, j, :) = [0.8 0.0 0.0];
           Cg(i, j, :) = [0.0 0.8 0.0];
           Cb(i, j, :) = [0.0 0.0 0.8];
           Ck(i, j, :) = [0.8 0.8 0.8];
           Cy(i, j, :) = [0.8 0.8 0.0];
        end    
    end

    fig = figure(15);
    hold on;
    x = db.dcollpitch_arr;
    y = db.eta_Tcoax_arr(1, :);
    [X, Y] = meshgrid(x, y);
    Z = transpose(db.Pcoax_arr) ./ db.P_h;
    surf(X, Y, Z, Cr);
    text(x(end) + 0.1, y(end)+0.1, Z(end, end)+0.1, 'coax');
    Z = transpose(db.Psbs_arr) ./ db.P_h;
    surf(X, Y, Z, Cb);
    text(x(end) + 0.1, y(end)+0.1, Z(end, end)+0.1, 'sbs');
    Z = transpose(db.kint_arr);
    surf(X, Y, Z, Cg);
    text(x(end) + 0.1, y(end)+0.1, Z(end, end)+0.1, 'k_{int}');
    % colorbar
    title('coax and sbs power');
    xlabel(['Differential pitch \theta ' char(176)]);
    ylabel('Thrust ratio \eta_{T}');
    zlabel('Normalized power');
    zlim([0, 2.2]);
    view(35, 35);
    grid on;               
    
    fig = figure(16);
    hold on;
    x = db.dcollpitch_arr;
    y = db.eta_Tcoax_arr(1, :);
    [X, Y] = meshgrid(x, y);
    Z = transpose(db.Qcoax_u_arr.*rads2rpm);
    surf(X, Y, Z, Cr);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'coax\_u');
    Z = transpose(db.Qcoax_l_arr.*rads2rpm);
    surf(X, Y, Z, Cr);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'coax\_l');
    Z = transpose(db.Qsbs_u_arr.*rads2rpm);
    surf(X, Y, Z, Cb);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'sbs\_u');
    Z = transpose(db.Qsbs_l_arr.*rads2rpm);
    surf(X, Y, Z, Cb);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'sbs\_l');
    % colorbar
    title('coax and sbs rotor torque');
    xlabel(['Differential pitch \theta ' char(176)]);
    ylabel('Thrust ratio \eta_{T}');
    % zlabel('RPM ratio \eta_{\Omega}');
    zlabel('Torque Nm');
    view(35, 35);
    grid on;     
    
    fig = figure(17);
    hold on;
    x = db.dcollpitch_arr;
    y = db.eta_Tcoax_arr(1, :);
    [X, Y] = meshgrid(x, y);
    Z = transpose(db.eta_Qcoax_arr);
    surf(X, Y, Z, Cr);
    text(x(end) + 0.1, y(end)+0.1, Z(end, end)+0.1, 'coax');
    Z = transpose(db.eta_Qsbs_arr);
    surf(X, Y, Z, Cb);
    text(x(end) + 0.1, y(end)+0.1, Z(end, end)+0.1, 'sbs');
    Z = transpose(db.eta_Qsbs_arr ./ db.eta_Qsbs_arr);  % Plane at Z == 1
    surf(X, Y, Z, Ck);
    % colorbar
    title('coax and sbs torque ratio \eta_{\Omega}');
    xlabel(['Differential pitch \theta ' char(176)]);
    ylabel('Thrust ratio \eta_{T}');
    zlabel('Torque ratio \eta_{\Omega}');
    view(35, 35);
    grid on;               

    fig = figure(18);
    hold on;
    x = db.dcollpitch_arr;
    y = db.eta_Tcoax_arr(1, :);
    [X, Y] = meshgrid(x, y);
    Z = transpose(db.omegacoax_u_arr.*rads2rpm);
    surf(X, Y, Z, Cr);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'coax\_u');
    Z = transpose(db.omegacoax_l_arr.*rads2rpm);
    surf(X, Y, Z, Cr);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'coax\_l');
    Z = transpose(db.omegasbs_u_arr.*rads2rpm);
    surf(X, Y, Z, Cb);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'sbs\_u');
    Z = transpose(db.omegasbs_l_arr.*rads2rpm);
    surf(X, Y, Z, Cb);
    text(x(end) + 0.01, y(end)+0.01, Z(end, end)+0.1, 'sbs\_l');
    % surf(X, Y, transpose(omegacoax_u_arr .* 0));  % Plane at Z == 1
    % colorbar
    title('coax and sbs rotor speed');
    xlabel(['Differential pitch \theta ' char(176)]);
    ylabel('Thrust ratio \eta_{T}');
    % zlabel('RPM ratio \eta_{\Omega}');
    zlabel('Rotor speed \Omega RPM');
    view(35, 35);
    grid on;  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str1 = 'bet_coax_eta_T_arr';
    str2 = rotortype;  
    nfig_arr = [15, 16, 17, 18];
    for nfig = nfig_arr
        fig = figure(nfig);
        plot_save_nfig_arr(str1, str2, nfig, true);
        
        view(80, 10);
        filename = ['img/' str1 '_' str2 '_' num2str(nfig) 'b.jpg'];
        saveas(fig, filename);
    end
end
    
