function bet_main_omega_plot(st, str1, str2)

    % Loadcel precision
    % According to https://www.ati-ia.com/products/ft/ft_models.aspx?id=Delta
    % and model US-150-600
    % https://www.convertunits.com/from/lbf/to/N
    % the max error in thrust is 1/16 = 0.0625 lbf = 0.2780 N
    % http://convert-to.com/conversion/torque/convert-lbf-in-to-n-m.html
    % the max error in torque is 1/16 = 0.0625 lbf-in =  0.01 Nm

    % Experimentally determined error are about
    % Error in thrust is 1.00Nm
    error_thrust_arr = 1.0 * ones(size(st.ane_T_arr));
    % Error in torque is 0.03Nm
    error_torque_arr = 0.03 * ones(size(st.ane_Q_arr));
    % Error in RPM is 30 RPM
    error_rpm_arr = 30 * ones(size(st.omegarpm_arr));
    % Error in power is 0.03Nm * 30rpm = 0.03 * 3.14 rad/s = 0.1 W
    error_power_arr = 0.1;

    plot_error = true;

    nfig = 1;   
    fig = figure(nfig);
    subplot(2, 1, 1)
    hold on;
    grid on;
    set(gca,'FontSize', 14)
    % plot(omegarpm_arr, kde_T_arr, '-*');
    % plot(st.omegarpm_arr, st.ane_T_arr, '-*');
    plot(st.omegarpm_arr, st.bet_T_arr, '-b');
    errorbar(st.omegarpm_arr, st.ane_T_arr, error_thrust_arr, '-r')
    errorbar(st.omegarpm_arr, st.ane_T_arr, error_rpm_arr, '-r', 'horizontal')
    xlabel('Rotor speed RPM');
    ylabel('Thrust N');  
    xlim([2000, 4500]); 
    ylim([0, 80]); 

    if plot_error
        yyaxis right
        plot(st.omegarpm_arr, st.perr_T_arr, '-k');
        % legend('KDE', 'AneCha', 'BET', 'error', 'Location', 'northwest');     
        % legend('AneCha', 'BET', 'error', 'Location', 'northwest');
        ylabel('Error %');
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k'; 
        ylim([-10, 10]); 

    %    ind = find(omegarpm_arr==3890);
    %    omegaind = omegarpm_arr(ind);
    %    val = round(kde_T_arr(ind), 1);
    %    text(2600, 17, ['KDE: ' num2str(omegaind) ' RPM, ' num2str(val) ' N'])
    %    val = round(bet_T_arr(ind), 1);
    %    text(2600, 12, ['BET: ' num2str(omegaind) ' RPM, ' num2str(val) ' N'])
        % 3890 rpm 126.5090 N 4.7044 Nm
        text(2800, -5, ['std dev error: ' num2str(st.std_err_T) ' N'])
        text(2800, -8, ['   mean error: ' num2str(st.mean_err_T) ' N'])
    end

    subplot(2, 1, 2)
    hold on;
    grid on;
    set(gca,'FontSize', 14)
    % plot(omegarpm_arr, kde_Q_arr, '-*');
    % plot(st.omegarpm_arr, st.ane_Q_arr, '-*');
    plot(st.omegarpm_arr, st.bet_Q_arr, '-b');
    errorbar(st.omegarpm_arr, st.ane_Q_arr, error_torque_arr, '-r')
    errorbar(st.omegarpm_arr, st.ane_Q_arr, error_rpm_arr, '-r', 'horizontal')
    xlabel('Rotor speed RPM');
    ylabel('Torque Nm');
    xlim([2000, 4500]); 
    ylim([0, 3]); 
    
    if plot_error
        yyaxis right
        plot(st.omegarpm_arr, st.perr_Q_arr, '-k');
        % legend('KDE', 'AneCha', 'BET', 'error', 'Location', 'northwest');     
        % legend('AneCha', 'BET', 'error', 'Location', 'northwest');
        legend('BET', 'Exp. uncertainty', 'Exp. uncertainty', 'Perc. error', 'Location', 'northwest');
        ylabel('Error %');
        ax = gca;
        ax.YAxis(1).Color = 'k';
        ax.YAxis(2).Color = 'k'; 
        ylim([-10, 10]); 

    %    ind = find(omegarpm_arr==3890);
    %    omegaind = omegarpm_arr(ind);
    %    val = round(kde_Q_arr(ind), 1);
    %    text(2600, 17, ['KDE: ' num2str(omegaind) ' RPM, ' num2str(val) ' Nm'])
    %    val = round(bet_Q_arr(ind), 1);
    %    text(2600, 12, ['BET: ' num2str(omegaind) ' RPM, ' num2str(val) ' Nm'])
        % text(2800, 8, ['std dev error: ' num2str(st.std_err_Q) ' Nm'])
        % text(2800, 5, ['   mean error: ' num2str(st.mean_err_Q) ' Nm'])
        text(3300, -2, ['std dev error: ' num2str(st.std_err_Q) ' Nm'])
        text(3300, -5, ['  mean error: ' num2str(st.mean_err_Q) ' Nm'])
    else
        legend('BET', 'Exp. uncertainty', 'Exp. uncertainty', 'Location', 'northwest');
    end

    fig = figure(nfig);
    filename = ['img/' str1 '_' str2 '.jpg'];
    saveas(fig, filename);  
    
    close all;   
end

