function mt_main_advratio()
    clear
    clc
    close all

    % Add functions to the path
    addpath('mt_bet');

    rpm2rads    = (pi/30);  % 60 RPM = 1 RPsecond = 1 Hz = 2*pi rad/s

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_c    = 0.0;
    mu          = 0.0;
    
    % Find omega that generates 15*9.8/4 = 36.75 N thurst
    omega = 2100 * pi / 30;     
    modeltype = 'interp1';
    % [T, Q, P] = kde_rotor_TQP(omega, 'KDE8218XF120_KDECF305DP', modeltype);
    [T, Q, P] = kde_rotor_TQP(omega, 'KDE6213XF185_KDECF245DP', modeltype);
    Thover = T;
       
    % rotortype   = 'KDECF305DP';
    rotortype   = 'KDECF245DP';               
    [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
    omega       = omega;
    CT_target   = NaN;              % CT_target is only needed for mu != 0
    blade_st    = blade_model(...
        rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);     
    
    CT_hover    = Thover / (blade_st.rho * blade_st.rotArea * blade_st.Vtip^2);
    blade_st.CT_target  = CT_hover;   
    Cd0         = bet_calculate_Cd0_c81table(blade_st);

    [T, Q, lambda_i] = mt_forces(blade_st);

    lambda      = blade_st.lambda_c + lambda_i;
    lambda_h    = lambda;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda_c_arr = linspace(0, 0.1, 4);
    mu_arr = linspace(0, 0.2, 10);
    nifg1 = 1;
    nifg2 = 2;
    fig1 = figure(nifg1);
    fig2 = figure(nifg2);

    for j = 1:length(lambda_c_arr)
        for i = 1:length(mu_arr)           
            
            blade_st    = blade_model(...
                rotortype, rho, lambda_c_arr(j), mu_arr(i), collpitch, omega, CT_target);     
            blade_st.CT_target  = CT_hover;   
            blade_st.Cd0        = Cd0;    
            
            [T, Q, lambda_i] = mt_forces(blade_st);    
            lambda = blade_st.lambda_c + lambda_i;
            
            T_arr(i) = T;
            Q_arr(i) = Q;
            lambda_i_arr(i) = lambda_i;

            lambda_arr(i) = lambda;
            TPP_alpha_arr(i) = blade_st.TPP_alpha;        
            
            % Save label for this case
            switch j
                case 1
                    lbl_1 = ['\lambda_c = ' num2str(blade_st.lambda_c)];
                case 2
                    lbl_2 = ['\lambda_c = ' num2str(blade_st.lambda_c)];
                case 3
                    lbl_3 = ['\lambda_c = ' num2str(blade_st.lambda_c)];
                case 4
                    lbl_4 = ['\lambda_c = ' num2str(blade_st.lambda_c)];
                case 5       
                    lbl_5 = ['\lambda_c = ' num2str(blade_st.lambda_c)];
            end                        
        end

        figure(nifg1)
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
        ylabel('\alpha_{TPP} deg');
        xlabel('Advance ratio \mu');      
        ylim([0, 90]);
        
        figure(nifg2)
        subplot(2, 1, 1)
        hold on;
        grid on;
        plot(mu_arr, T_arr, '-*');
        ylabel('Thrust N');
        xlabel('Advance ratio \mu');
        % ylim([0, 10^5]);
        ylim([0, 50]);
        
        subplot(2, 1, 2)
        hold on;
        grid on;
        plot(mu_arr, Q_arr, '-*');
        ylabel('Torque Nm');
        xlabel('Advance ratio \mu');
        arg = '$Torque = \frac{Pi \ + \ Pc \ + \ Po \ + \ Pp}{\Omega}$';
        text(0.25, 2.8*10^5, arg, 'Interpreter', 'latex')
        ylim([0, 3]);   
        % ylim([0, 5*10^5]);
              
    end
    
    % Add legend and save
    str1 = 'mt_main_advratio';
    str2 = rotortype;

    fig1 = figure(nifg1);
    sgtitle('MT at constant climb inflow')
    subplot(3, 1, 2)
    legend(lbl_1, lbl_2, lbl_3, lbl_4, 'Location', 'best') % , lbl_5);
    set(fig1,'units','centimeters','position',[0, 0, 20, 10])
    % pause(10)
    filename = ['img/' str1 '_' str2 '_' num2str(nifg1) '.jpg'];
    saveas(fig1, filename)

    fig2 = figure(nifg2);
    sgtitle('MT at constant climb inflow')
    subplot(2, 1, 1)
    legend(lbl_1, lbl_2, lbl_3, lbl_4, 'Location', 'best') % , lbl_5);
    set(fig2,'units','centimeters','position',[0, 0, 20, 10])
    % pause(20)
    filename = ['img/' str1 '_' str2 '_' num2str(nifg2) '.jpg'];
    saveas(fig2, filename);

    close all
end




