function bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1, str2)       

    dcollpitch_arr = db.dcollpitch_arr;
    collpitch_u_arr = db.collpitch_u_arr;
    collpitch_l_arr = db.collpitch_l_arr;
    P_h = db.P_h;
    W = db.W;
    coaxu_r0 = db.coaxu_r0;
    coaxl_r0 = db.coaxl_r0;
    rho = db.blade_h_st.rho; 
    rotArea = db.blade_h_st.rotArea;
    R = db.blade_h_st.R;

    % i = dcollpitch 
    % j = eta_thrust
    
    i = find(db.dcollpitch_arr == dcollpitch);
    for j = 1:length(db.eta_T_arr)
        bet_u_st = db.coax.bet_u_st(i, j);
        bet_l_st = db.coax.bet_l_st(i, j);

        Tcoax_u_arr(j) = bet_u_st.total.T;
        Tcoax_l_arr(j) = bet_l_st.total.T;
        eta_Tcoax_arr(j) = Tcoax_u_arr(j) ./ Tcoax_l_arr(j);

        lambdacoax_u_arr(j) = bet_u_st.total.lambda;
        lambdacoax_l_arr(j) = bet_l_st.total.lambda;
        eta_lambdacoax_arr(j) = lambdacoax_u_arr(j) ./ lambdacoax_l_arr(j);
        vu_arr(j) = lambdacoax_u_arr(j) * bet_u_st.blade_st.Vtip;
        vl_arr(j) = lambdacoax_l_arr(j) * bet_l_st.blade_st.Vtip;
        % distribution along r
        lrcoax_u_arr(j, :) = bet_u_st.total.lcr + bet_u_st.total.lir;
        lrcoax_l_arr(j, :) = bet_l_st.total.lcr + bet_l_st.total.lir;
        
        omegacoax_u_arr(j) = bet_u_st.blade_st.omega;
        omegacoax_l_arr(j) = bet_l_st.blade_st.omega;
        eta_omegacoax_arr(j) = omegacoax_u_arr(j) ./ omegacoax_l_arr(j);
        
        Qcoax_u_arr(j) = bet_u_st.total.Q;
        Qcoax_l_arr(j) = bet_l_st.total.Q;
        eta_Qcoax_arr(j) = Qcoax_u_arr(j) ./ Qcoax_l_arr(j);
        
        Pcoax_u_arr(j) = bet_u_st.total.P;
        Pcoax_l_arr(j) = bet_l_st.total.P;
        Pcoax_arr(j) = Pcoax_u_arr(j) + Pcoax_l_arr(j);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bet_u_st = db.sbs.bet_u_st(i, j);
        bet_l_st = db.sbs.bet_l_st(i, j);

        Tsbs_u_arr(j) = bet_u_st.total.T;
        Tsbs_l_arr(j) = bet_l_st.total.T;
        eta_Tsbs_arr(j) = Tsbs_u_arr(j) ./ Tsbs_l_arr(j);
        
        lambdasbs_u_arr(j) = bet_u_st.total.lambda;
        lambdasbs_l_arr(j) = bet_l_st.total.lambda;
        eta_lambdasbs_arr(j) = lambdasbs_u_arr(j) ./ lambdasbs_l_arr(j);

        omegasbs_u_arr(j) = bet_u_st.blade_st.omega;
        omegasbs_l_arr(j) = bet_l_st.blade_st.omega;
        eta_omegasbs_arr(j) = omegasbs_u_arr(j) ./ omegasbs_l_arr(j);
        
        Qsbs_u_arr(j) = bet_u_st.total.Q;
        Qsbs_l_arr(j) = bet_l_st.total.Q;
        eta_Qsbs_arr(j) = Qsbs_u_arr(j) ./ Qsbs_l_arr(j);

        Psbs_u_arr(j) = bet_u_st.total.P;
        Psbs_l_arr(j) = bet_l_st.total.P;
        Psbs_arr(j) = Psbs_u_arr(j) + Psbs_l_arr(j);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kint_arr(j)   = Pcoax_arr(j) ./ Psbs_arr(j);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mt_vind_arr(j) = sqrt(W/(2*rho*pi*R^2));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assuming only upper on lower interference and perfect contraction
        % vi_u = np.sqrt(m1_thrust / (2 * density * area))
        % vc_l = vi_u
        % vh_l = np.sqrt(m2_thrust / (2 * density * area))
        % vi_l = vh_l * (
        %         - (0.5 * vc_l / vh_l) + np.sqrt((0.5 * vc_l / vh_l) ** 2 + 1)
        % )
        
        % MT inflow for upper rotor
        Vtip = omegacoax_u_arr(j) * R;
        CT = Tcoax_u_arr(j) / (rho * rotArea * Vtip^2);
        mu = 0;
        TPP_alpha = 0;
        lambda_c = 0;
        mt_lambdacoax_u_arr(j) = mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c);
        mt_vu_arr(j) = mt_lambdacoax_u_arr(j) * Vtip;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MT inflow for lower rotor
        Vtip = omegacoax_l_arr(j) * R;
        CT = Tcoax_l_arr(j) / (rho*rotArea*Vtip^2);
        mu = 0;
        TPP_alpha = 0;
        lambda_c = mt_lambdacoax_u_arr(j);
        mt_lambdacoax_l_arr(j) = mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c);
        mt_vl_arr(j) = mt_lambdacoax_l_arr(j) * Vtip;

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % MT inflow for lower rotor using extra RPM
        % omega_l = omegacoax_u_arr(j)*1.2;
        % csvp = db;
        % m1_thrust = Tcoax_u_arr(j);
        % omega_u = omegacoax_u_arr(j)
        % m2_thrust = Tcoax_l_arr(j)
        % [thrust, a_thrust] = mt_coaxial_lower_rotor_model_jack(...
        %     omega_l, csvp, m1_thrust, omega_u, m2_thrust)
    end

    fprintf('[bet_plot_coax_dcollpitch] dcollpitch %.4f, i %d, P_h %.4f \n', ...
        dcollpitch, i, P_h);
    fprintf('[bet_plot_coax_dcollpitch] collpitch_u %.4f, collpitch_l %.4f \n', ...
        rad2deg(collpitch_u_arr(i)), rad2deg(collpitch_l_arr(i)));

    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;  
    
    % Save data to file
    table_m1_thrust = Tcoax_u_arr;
    table_m2_thrust = Tcoax_l_arr;    
    table_eta_trust = eta_Tcoax_arr;
    table_m1_omega = omegacoax_u_arr;
    table_m2_omega = omegacoax_l_arr;
    table_eta_omega = eta_omegacoax_arr;
    table_m1_power = Qcoax_u_arr .* table_m1_omega;
    table_m2_power = Qcoax_l_arr .* table_m2_omega;   
    table_coax_power = Pcoax_arr;
    T = table(...
        table_m1_thrust(:)     , ...
        table_m2_thrust(:)     , ...
        table_eta_trust(:)     , ...
        table_m1_omega(:)      , ...
        table_m2_omega(:)      , ...
        table_eta_omega(:)     , ...
        table_m1_power(:)      , ...
        table_m2_power(:)      , ...
        table_coax_power(:)      ...
    );
    
    filename = ['img/' str1 '_' str2 '.txt'];
    writetable(T, filename);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 1;
    figure(nfig);
    subplot(2, 1, 1)        
    hold on;
    grid on;
    plot(eta_Tcoax_arr, Tcoax_u_arr, '-r*');
    plot(eta_Tcoax_arr, Tcoax_l_arr, '-b*');
    W_arr = Tcoax_u_arr + Tcoax_l_arr;
    plot(eta_Tcoax_arr, W_arr, '-k*');
    legend('coax u', 'coax l', 'coax W', 'Location','NorthEast');          
    xlabel('\eta_T');
    ylabel('Thrust N'); 
    if strcmp(rotortype, 'KDECF305DP')
        ylim([0, 300]);
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        xlim([0, 5]);    
        ylim([0, 70]); yticks([0, 10, 20, 30, 40, 50, 60, 70]);
    end
              
    subplot(2, 1, 2)
    hold on;
    grid on;
    plot(eta_Tsbs_arr, Tsbs_u_arr, 'r-*');
    plot(eta_Tsbs_arr, Tsbs_l_arr, 'b-*');
    W_arr = Tsbs_u_arr + Tsbs_l_arr;
    plot(eta_Tsbs_arr, W_arr, '-k*');        
    legend('sbs u', 'sbs l', 'sbs W', 'Location','NorthEast');  
    xlabel('\eta_T');
    ylabel('Thrust N'); 
    if strcmp(rotortype, 'KDECF305DP')
        ylim([0, 300]);
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        xlim([0, 5]);    
        ylim([0, 70]);  yticks([0, 10, 20, 30, 40, 50, 60, 70]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 2;
    figure(nfig);
    subplot(2, 1, 1)
    hold on;
    grid on;
    yyaxis left
    plot(eta_Tcoax_arr, omegacoax_u_arr.*rads2rpm, 'r-*');        
    plot(eta_Tcoax_arr, omegacoax_l_arr.*rads2rpm, 'b-*');       
    xlabel('\eta_T');
    ylabel('\Omega RPM');        
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tcoax_arr(1) - 0.3, omegacoax_u_arr(1).*rads2rpm, ...
            [num2str(dcollpitch) char(176)]);    
        text(eta_Tcoax_arr(1) - 0.3, omegacoax_l_arr(1).*rads2rpm, ...
            [num2str(dcollpitch) char(176)]); 
        ylim([2000, 7000]);
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 4000]);
        xlim([0, 5]);
    end
    yyaxis right
    plot(eta_Tcoax_arr, eta_omegacoax_arr, 'k-*');       
    legend('coax u', 'coax l', '\eta_{\Omega}', 'Location','SouthEast');     
    ylabel('\eta_{\Omega}');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k'; 
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tcoax_arr(end) + 0.1, eta_omegacoax_arr(end), ...
            [num2str(dcollpitch) char(176)]); 
        ylim([0, 3]);
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 2]);
        xlim([0, 5]);
    end
    
    subplot(2, 1, 2)
    hold on;
    grid on;    
    yyaxis left    
    plot(eta_Tsbs_arr, omegasbs_u_arr.*rads2rpm, 'r-*');           
    plot(eta_Tsbs_arr, omegasbs_l_arr.*rads2rpm, 'b-*');
    legend('sbs u', 'sbs l', 'Location','southeast');     
    xlabel('\eta_T');
    ylabel('\Omega RPM');        
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tsbs_arr(1) - 0.3, omegasbs_u_arr(1).*rads2rpm, ...
            [num2str(dcollpitch) char(176)]);     
        text(eta_Tsbs_arr(1) - 0.3, omegasbs_l_arr(1).*rads2rpm, ...
            [num2str(dcollpitch) char(176)]); 
        ylim([2000, 7000]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 4000]);
        xlim([0, 5]);
    end    
    yyaxis right
    plot(eta_Tsbs_arr, eta_omegasbs_arr, 'k-*');
    legend('sbs u', 'sbs l', '\eta_{\Omega}', 'Location','southeast');     
    ylabel('\eta_{\Omega}');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k'; 
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tsbs_arr(end) + 0.1, eta_omegasbs_arr(end), ...
            [num2str(dcollpitch) char(176)]);
        ylim([0, 3]);
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 2]);
        xlim([0, 5]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 3;
    figure(nfig);
    subplot(2, 1, 1)
    hold on;
    grid on;    
    plot(omegacoax_u_arr.*rads2rpm, Tcoax_u_arr, 'r-*');        
    plot(omegacoax_l_arr.*rads2rpm, Tcoax_l_arr, 'b-*');       
    legend('coax u', 'coax l', 'Location','northwest');    
    xlabel('\Omega RPM');
    ylabel('Thrust N');  
    if strcmp(rotortype, 'KDECF305DP')
        text(omegacoax_u_arr(end).*rads2rpm + 10, Tcoax_u_arr(end) + 20, ...
            [num2str(dcollpitch) char(176)]);      
        text(omegacoax_l_arr(end).*rads2rpm - 10, Tcoax_l_arr(end) - 20, ...
            [num2str(dcollpitch) char(176)]);      
        xlim([2000, 7000]);
        ylim([0, 300]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        xlim([0, 4000]);
        ylim([0, 70]); yticks([0, 10, 20, 30, 40, 50, 60, 70]);
    end
    
    subplot(2, 1, 2)
    hold on;
    grid on;    
    plot(omegasbs_u_arr.*rads2rpm, Tsbs_u_arr, 'r-*');       
    plot(omegasbs_l_arr.*rads2rpm, Tsbs_l_arr, 'b-*');            
    legend('sbs u', 'sbs l', 'Location','northwest'); 
    xlabel('\Omega RPM');
    ylabel('Thrust N');     
    if strcmp(rotortype, 'KDECF305DP')
        text(omegasbs_u_arr(end).*rads2rpm + 10, Tsbs_u_arr(end) + 20, ...
            [num2str(dcollpitch) char(176)]);        
        text(omegasbs_l_arr(end).*rads2rpm - 10, Tsbs_l_arr(end) - 20, ...
            [num2str(dcollpitch) char(176)]);      
        xlim([2000, 7000]);
        ylim([0, 300]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        xlim([0, 4000]);
        ylim([0, 70]); yticks([0, 10, 20, 30, 40, 50, 60, 70]);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 4;
    figure(nfig);
    hold on;
    grid on;
    yyaxis left
    % plot(eta_Tcoax_arr, Pcoax_arr ./ P_h, 'r-*');
    % plot(eta_Tsbs_arr, Psbs_arr ./ P_h, 'b-*');
    plot(eta_Tcoax_arr, Pcoax_arr, 'r-*');
    plot(eta_Tsbs_arr, Psbs_arr, 'b-*');
    xlabel('\eta_T');
    ylabel('Power W');    
    if strcmp(rotortype, 'KDECF305DP')
        ylim([0, 2000]); % ylim([0, 2.0]);    
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 750]); yticks([0, 250, 500, 750]);    
        xlim([0, 5]);
    end    
    yyaxis right
    plot(eta_Tcoax_arr, kint_arr, 'k-*');
    legend('coax', 'sbs', 'k_{int}', 'Location','northeast');
    ylabel('k_{int}');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';  
    if strcmp(rotortype, 'KDECF305DP')
        text(0.1, kint_arr(1), [num2str(dcollpitch) char(176)]);    
        ylim([0, 3.0]);    
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 3.0]);    
        xlim([0, 5]);
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 5;
    figure(nfig);
    subplot(2, 1, 1)
    hold on;
    grid on;
    yyaxis left
    plot(eta_Tcoax_arr, Qcoax_u_arr, 'r-*');
    plot(eta_Tcoax_arr, Qcoax_l_arr, 'b-*');                  
    xlabel('\eta_T');
    ylabel('Torque Nm');
    if strcmp(rotortype, 'KDECF305DP')
%        text(eta_Tcoax_arr(end) + 0.1, Qcoax_l_arr(end) + 0.1, ...
%            [num2str(dcollpitch) char(176)]);    
        ylim([0, 10]);   
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 2.0]);    
        xlim([0, 5]);
    end        
    yyaxis right
    plot(eta_Tcoax_arr, eta_Qcoax_arr, 'k-*');         
    legend('coax u', 'coax l', '\eta_Q', 'Location','north');     
    ylabel('\eta_Q');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tcoax_arr(end) + 0.1, eta_Qcoax_arr(end) + 0.1, ...
            [num2str(dcollpitch) char(176)]);    
        ylim([0, 5]);   
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 5.0]);    
        xlim([0, 5]);
    end  
    
    subplot(2, 1, 2)
    hold on;
    grid on;
    yyaxis left  
    plot(eta_Tsbs_arr, Qsbs_u_arr, 'r-*');
    plot(eta_Tsbs_arr, Qsbs_l_arr, 'b-*');
    xlabel('\eta_T');
    ylabel('Torque Nm');
    if strcmp(rotortype, 'KDECF305DP')
%        text(eta_Tsbs_arr(end) + 0.1, Qsbs_l_arr(end) + 0.1, ...
%            [num2str(dcollpitch) char(176)]);     
        ylim([0, 10]);   
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 2.0]);    
        xlim([0, 5]);
    end        
    yyaxis right
    plot(eta_Tsbs_arr, eta_Qsbs_arr, 'k-*');   
    legend('sbs u', 'sbs l', '\eta_Q', 'Location','north');     
    ylabel('\eta_Q');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    if strcmp(rotortype, 'KDECF305DP')
        text(eta_Tsbs_arr(end) + 0.1, eta_Qsbs_arr(end) + 0.1, ...
            [num2str(dcollpitch) char(176)]);       
        ylim([0, 5]);   
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        ylim([0, 5]);    
        xlim([0, 5]);
    end            
    
    pause(30*0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 6;
    figure(nfig);
    subplot(2, 1, 1)
    hold on;
    grid on;
    % yyaxis left
    % plot(eta_Tcoax_arr, lambdacoax_u_arr, 'r-*');
    % plot(eta_Tcoax_arr, lambdacoax_l_arr, 'b-*');  
    % plot(eta_Tcoax_arr, mt_lambdacoax_u_arr, 'r-.*');
    % plot(eta_Tcoax_arr, mt_lambdacoax_l_arr, 'b-.*');
    % plot(eta_Tcoax_arr, 0.0708 * ones(size(mt_lambdacoax_l_arr)), 'k--*');
    plot(eta_Tcoax_arr, vu_arr, 'r-*');
    plot(eta_Tcoax_arr, vl_arr, 'b-*');  
    plot(eta_Tcoax_arr, mt_vu_arr, 'r-.*');
    plot(eta_Tcoax_arr, mt_vl_arr, 'b-.*');    
    plot(eta_Tcoax_arr, mt_vind_arr, 'k--*');
    xlabel('\eta_T');
    % ylabel('\lambda = \lambda_c + \lambda_i');
    ylabel('v = v_c + v_i');
    if strcmp(rotortype, 'KDECF305DP')
%        text(eta_Tcoax_arr(end) + 0.1, lambdacoax_l_arr(end) + 0.1, ...
%            [num2str(dcollpitch) char(176)]);    
        ylim([0, 10]);   
        xlim([0, 4.5]);
    end
    if strcmp(rotortype, 'KDECF245DP')
        % ylim([0, 0.15]);
        ylim([0, 15]);
        xlim([0, 5]);
    end        
    legend('coax u', 'coax l', 'mt coax u', 'mt coax l', 'mt single', ...
        'Location','northeast');  

    % lambda_mtu = mt_lambdacoax_u_arr(1)
    % lambda_mtl = mt_lambdacoax_l_arr(1)
    % lambda_mtu =
    %     0.0659
    % lambda_mtl =
    %     0.1016
    
    subplot(2, 1, 2)
    hold on;
    grid on;
    for j = 1:length(db.eta_T_arr)
        plot(linspace(0, 1, length(lrcoax_u_arr(j, :))), lrcoax_u_arr(j, :), '-r');
        plot(linspace(0, 1, length(lrcoax_l_arr(j, :))), lrcoax_l_arr(j, :), '-b');
    end
    xlabel('Normalized radius');
    ylabel('BET \lambda = \lambda_c + \lambda_i');
    % ylim([0, 2]);    
    % xlim([0, 5]);
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig_arr = [1, 2, 3, 4, 5, 6];
    savefig = false; 
    plot_save_nfig_arr(str1, str2, nfig_arr, savefig);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (length(dcollpitch_arr) > 1) && (length(eta_T_arr) > 1)
        bet_plot_coax_Zsurf(rotortype, db);
    end
end
    
