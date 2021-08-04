function bet_main_advratio()

    clear all
    close all
    format compact
    clc
    
    % Add functions to the path
    % addpath([ bpath '/' 'mt_bet']);    
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 

    nfig1 = 1;
    nfig2 = 2;
    figure(nfig1);
    figure(nfig2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Find omega that generates 15*9.8/4 = 36.75 N thurst
    omega_h = 1700 * pi / 30;     
    modeltype = 'interp1';
    [T, Q, P] = kde_rotor_TQP(omega_h, 'KDE8218XF120_KDECF305DP', modeltype);  
    Thover = T;    

    % Calculate thrust, torque and inflow at hover    
    rotortype   = 'KDECF305DP';
    % rotortype   = 'KDECF245DP';               
    [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
    omega       = omega_h;
    CT_target   = NaN;              % CT_target is only needed for mu != 0
    blade_st    = blade_model(...
        rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
    % Add CT_target for mu != 0 inflow calculations
    rhoAVtip2   = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;
    blade_st.CT_target  = Thover / rhoAVtip2;
    
    bet_st      = bet_forces(blade_st);
    bet_st      = bet_forces_add_total(bet_st, false);    
    lambda_c    = bet_st.total.lambda_c;
    lambda_i    = bet_st.total.lambda_i;
    lambda      = bet_st.total.lambda;
    lambda_h    = lambda;

    lambda_c_arr    = linspace(0, 0.1, 4);
    mu_arr          = linspace(0, 0.2, 10);
    
    for j = 1:length(lambda_c_arr)
        fprintf('lambda_c iteration %d/%d \n', j, length(lambda_c_arr))
        for i = 1:length(mu_arr)                   
        
            blade_st = blade_model(...
                rotortype, rho, lambda_c_arr(j), mu_arr(i), collpitch, ...
                omega, Thover / rhoAVtip2);
            
            % Calculate thrust, torque and inflow
            bet_st = bet_forces(blade_st);
            bet_st = bet_forces_add_total(bet_st, false);
            
            % Save variables
            TPP_alpha_arr(i)    = blade_st.TPP_alpha;
            T_arr(i)            = bet_st.total.T;
            Q_arr(i)            = bet_st.total.Q;
            lambda_arr(i)       = bet_st.total.lambda;
            lambda_i_arr(i)     = bet_st.total.lambda_i;        
        end     
        
         % Save label for this case
        switch j
            case 1
                lbl_1 = ['\lambda_c = ' num2str(round(blade_st.lambda_c, 2))];
            case 2
                lbl_2 = ['\lambda_c = ' num2str(round(blade_st.lambda_c, 2))];
            case 3
                lbl_3 = ['\lambda_c = ' num2str(round(blade_st.lambda_c, 2))];
            case 4
                lbl_4 = ['\lambda_c = ' num2str(round(blade_st.lambda_c, 2))];
            case 5       
                lbl_5 = ['\lambda_c = ' num2str(round(blade_st.lambda_c, 2))];
        end
       
        figure(nfig1)
        subplot(3, 1, 1)
        hold on;
        grid on;
        plot(mu_arr, lambda_arr./lambda_h, '-*');
        % plot(mu_arr, lambda_arr./lambda_h, '-*');   
        ylabel('\lambda / \lambda_h');
        xlabel('Advance ratio \mu');
        ylim([0, 3]);

        subplot(3, 1, 2)
        hold on;
        grid on;
        plot(mu_arr, lambda_i_arr./lambda_h, '-*');
        ylabel('\lambda_i / \lambda_h');
        xlabel('Advance ratio \mu');  
        ylim([0, 2]);
        
        subplot(3, 1, 3)
        hold on;
        grid on;
        plot(mu_arr, rad2deg(TPP_alpha_arr), '-*');
        % ylabel('Angle of attack \alpha');
        ylabel('\alpha_{TPP} deg');
        xlabel('Advance ratio \mu');
        ylim([0, 90]);
        
        figure(nfig2)
        subplot(2, 1, 1)
        hold on;
        grid on;
        plot(mu_arr, T_arr, '-*');
        ylabel('Thrust N');
        xlabel('Advance ratio \mu');
        ylim([0, 50]);
        
        subplot(2, 1, 2)
        hold on;
        grid on;
        plot(mu_arr, Q_arr, '-*');
        ylabel('Torque Nm');
        xlabel('Advance ratio \mu'); 
        ylim([0, 3]); 
       
    end
    
    % Add legend and save
    str1 = 'bet_main_advratio';
    str2 = rotortype;

    fig = figure(nfig1);
    sgtitle('BET at constant climb inflow');
    subplot(3, 1, 2);
    legend(lbl_1, lbl_2, lbl_3, lbl_4, 'Location', 'best'); % , lbl_5);
    set(fig,'units','centimeters','position',[0, 0, 20, 10]);
    filename = ['img/' str1 '_' str2 '_' num2str(nfig1) '.jpg'];
    saveas(fig, filename);

    fig = figure(nfig2);
    sgtitle('BET at constant climb inflow');
    subplot(2, 1, 1);
    legend(lbl_1, lbl_2, lbl_3, lbl_4, 'Location', 'north'); % , lbl_5);
    set(fig,'units','centimeters','position',[0, 0, 20, 10]);
    filename = ['img/' str1 '_' str2 '_' num2str(nfig2) '.jpg'];    
    saveas(fig, filename);

    close all;
    return
end
