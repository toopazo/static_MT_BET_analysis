function [...
        Pcoaxn_i_arr , ...
        eta_T_arr   , ...
        kint_arr      ...
    ] = mt_main_coaxial()
    clear all
    clc
    close all
    format compact
    % format short
    
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 

    W = 35;
    rotortype = 'KDECF245DP';
    [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
    
    density     = rho;
    radius      = 0.6223 / 2;           % 24.5in = 0.6223m
    area        = pi * radius * radius; % m2
    chord       = 4 / 100;              % blade chord  
    nblades     = 2;
    solidity    = (nblades * chord) / (pi * radius)  ;
    Cd0         = 0.013640; 
    
    P0          = ( W^(3/2) )/( sqrt(2*rho*area) )
    
    motornum    = 12;
    thrust_lvl  = W;
    casestr     = ['motornum' num2str(motornum) '_nblades' num2str(nblades) '_tlv' num2str(thrust_lvl)];
    
    eta_T_arr = 0:0.001:4;
    for i = 1:length(eta_T_arr)
        eta_T = eta_T_arr(i);
        T_u  = W*(eta_T)/(1+eta_T);
        T_l  = W - T_u;       
        
        % Get omega_u and omega_l
        motornum = motornum;
        nblades = nblades;
        thrust_lvl = thrust_lvl;
        [omega_u, omega_l] = ...
            mt_find_omega_a_thrust(T_u, T_l, motornum, nblades, thrust_lvl);
        
        modeltype = 'leishman2006aerodynamic';
        mt_st = mt_coax_power_plus_a_thrust(...
            T_u         , ...
            T_l         , ...
            omega_u     , ...
            omega_l     , ...        
            radius      , ...
            density     , ...
            area        , ...
            solidity    , ...
            Cd0         , ...
            modeltype     ...
        );
        
        % coax        
        Pcoaxu_i  = mt_st.Pcoaxu_i;
        Pcoaxl_i  = mt_st.Pcoaxl_i;
        Pcoax_i   = mt_st.Pcoax_i;
        Pcoaxu_o  = mt_st.Pcoaxu_o;
        Pcoaxl_o  = mt_st.Pcoaxl_o;
        Pcoax_o   = mt_st.Pcoax_o;
        Pcoax     = mt_st.Pcoax;
        % sbs
        Psbsu_i   = mt_st.Psbsu_i;
        Psbsl_i   = mt_st.Psbsl_i;
        Psbs_i    = mt_st.Psbs_i;
        Psbsu_o   = mt_st.Psbsu_o;
        Psbsl_o   = mt_st.Psbsl_o;
        Psbs_o    = mt_st.Psbs_o;
        Psbs      = mt_st.Psbs;
        % upper and lower ratios
        kint      = mt_st.kint;
        eta_T     = mt_st.eta_T;
        eta_v     = mt_st.eta_v;  
        if (eta_T - eta_T_arr(i)) > 10^-4
            error('Error in eta_T')
        end                
        
        Pcoax_i_arr(i)  = Pcoax_i;
        Psbs_i_arr(i)   = Psbs_i;
        kint_arr(i)     = kint;
        eta_v_arr(i)    = eta_v;
        
        k               = 1.15;
        Pcoaxu_arr(i)   = Pcoaxu_i + Pcoaxu_o;
        Pcoaxl_arr(i)   = Pcoaxl_i + Pcoaxl_o;
        Pcoax_arr(i)    = Pcoax;
        
        Psbsu_arr(i)    = Psbsu_i + Psbsu_o;
        Psbsl_arr(i)    = Psbsl_i + Psbsl_o;
        Psbs_arr(i)     = Psbs;        
        
%        if eta_T == 1
%            eta_T
%            eta_v
%            Pcoax_i / P0
%            Psbs_i / P0
%        end 
    end
    Pcoaxn_i_arr = Pcoax_i_arr ./ P0;

    [val, ind] = min(Pcoax_i_arr);
    min_Pcoax = val
    min_eta_T = eta_T_arr(ind)
    min_kint = kint_arr(ind)

    [val, ind] = min(Pcoax_arr);
    min_Pcoax = val
    min_eta_T = eta_T_arr(ind)
    
    % min_eta_omega = eta_omega(ind)
    % min_eta_torque = eta_torque(ind)
    % return

    fig = figure(1);
    hold on;
    grid on;
    plot(eta_T_arr, eta_v_arr, '-', 'LineWidth', 2);
    % plot(eta_T_arr, 1./(eta_T_arr-1), '-');
    % plot(eta_T_arr, sqrt(eta_T_arr), '-');
    title('Induced velocity ratio', 'FontSize', 18)
    ylabel('\eta_v = v_u / v_l', 'FontSize', 18);
    xlabel('\eta_T = T_u / T_l', 'FontSize', 18); 
    xlim([0, 5]);  
    ylim([0, 5]); 
    % axis square
    set(fig,'units','centimeters','position', [0, 0, 14, 10]); 
    % pause(60)
    saveas(fig, ['img/mt_main_coaxial_' casestr '_1.jpg']);

    fig = figure(2);
    hold on;
    grid on;
    title('Induced coaxial power', 'FontSize', 18)
    plot(eta_T_arr, Psbs_i_arr / P0, '-', 'LineWidth', 2);
    plot(eta_T_arr, Pcoax_i_arr / P0, '-', 'LineWidth', 2);
    ylabel('Power W', 'FontSize', 18);
    xlabel('\eta_T = T_u / T_l', 'FontSize', 18); 
    xlim([0, 5]);  
    ylim([0.5, 1.5]); 
    % ylim([0, 750]); 
    % yticks([0, 250, 500, 750])
    yyaxis right;
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k'; 
    plot(eta_T_arr, kint_arr, '-', 'LineWidth', 2);
    ylabel('k_{int}', 'FontSize', 18);
    legend('MT side-by-side', 'MT coax', 'k_{int}', 'Location', 'SouthEast')
    % axis square
    ylim([0.5, 1.5]); 
    set(fig,'units','centimeters','position', [0, 0, 14, 10]); 
    % pause(60)
    % return
    saveas(fig, ['img/mt_main_coaxial_' casestr '_2.jpg']);
    
    fig = figure(3);
    hold on;
    grid on;
    plot(eta_T_arr, Psbs_arr, '-', 'Color', 'black');
    plot(eta_T_arr, Pcoax_arr, '-', 'Color', 'blue');
    plot(eta_T_arr, Pcoaxu_arr, '-', 'Color', 'red');
    plot(eta_T_arr, Pcoaxl_arr, '-', 'Color', 'green');
    title('Coaxial power', 'FontSize', 18)
    ylabel('Power W', 'FontSize', 18);
    xlabel('\eta_T = T_u / T_l', 'FontSize', 18); 
    xlim([0, 5]);  
    ylim([0, 750]);  
    yticks([0, 250, 500, 750])
    legend('MT side-by-side', 'MT coax total', 'MT coax upper', 'MT coax lower', 'FontName', 'TimesNewRoman')
    % axis square
    set(fig,'units','centimeters','position', [0, 0, 14, 10]); 
    pause(60*0)
    saveas(fig, ['img/mt_main_coaxial_' casestr '_3.jpg']);  

    close all

    % Tu = W * eta_T / (1 + eta_T)
    % Tl = W - Tu
end