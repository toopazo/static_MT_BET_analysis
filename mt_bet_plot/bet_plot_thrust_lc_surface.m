function bet_plot_thrust_lc_surface(str1, str2, db)

    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;    

    % Unpack structure
    thrust_arr = db.thrust_arr;
    lc_arr = db.lc_arr;
    bet_omega_arr = db.bet_omega_arr;
    bet_thrust_arr = db.bet_thrust_arr;
    bet_torque_arr = db.bet_torque_arr;
    bet_power_arr = db.bet_power_arr;
    bet_lambdac_arr = db.bet_lambdac_arr;
    bet_lambdai_arr = db.bet_lambdai_arr;
    bet_lambda_arr = db.bet_lambda_arr;

    bet_lcr_arr = db.bet_lcr_arr;
    bet_lir_arr = db.bet_lir_arr;
    bet_lr_arr = db.bet_lr_arr;
    bet_aoa_arr = db.bet_aoa_arr;
    
    mt_lambda_arr = db.mt_lambda_arr;
    mt_lambdai_arr = db.mt_lambdai_arr;
    mt_v_arr = db.mt_v_arr;
    mt_power_arr = db.mt_power_arr;


    tdiff = diff(bet_thrust_arr);
    if norm(tdiff .* tdiff) > 1
        error('norm(tdiff .* tdiff) > 1')
    end
    mt_thrust = bet_thrust_arr(1)
    mt_lambda_h = db.mt_lambda_h
    mt_power_h = db.mt_power_h
    % mt_lambda_h =
    %     0.0708
    % mt_power_h =
    %     417.0932

    % tot_thrust = thrust_arr(1);
    tot_thrust = bet_thrust_arr(1);
    % tot_thrust = mt_thrust(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 1;
    fig = figure(nfig);
    hold on
    grid on

    nfig = 2;
    fig = figure(nfig);
    hold on
    grid on

    % Define colors 
    n = 6; 
    % colors = jet(n);
    % colors = hot(n);
    % colors = cool(n);
    colors = parula(n);
    % colors = [...
    %     [1 1 0]; ... % yellow
    %     [1 0 1]; ... % magenta	
    %     [0 1 1]; ... % cyan
    %     [1 0 0]; ... % red
    %     [0 1 0]; ... % green
    %     [0 0 1]; ... % blue
    %     [1 1 1]; ... % white
    %     [0 0 0]  ... % black
    % ];

    pcnt = 0;
    for i = 1:length(thrust_arr)
        for j = 1:length(lc_arr)   
            pcnt = pcnt + 1;

            r_arr = linspace(0, 1, length(bet_lr_arr(i, j, :)));
            r_arr = squeeze(r_arr);

            nfig = 1;
            fig = figure(nfig);
            inflow_lims = [-0.02, 0.18];

            subplot(2, 1, 1)
            hold on;
            grid on;
            y_arr = bet_lcr_arr(i, j, :);
            plot(r_arr, y_arr(:), '-*', 'Color', colors(pcnt,:))
            ylabel('\lambda_c')
            xlabel('Normalized radius')
            ylim(inflow_lims)

            subplot(2, 1, 2)
            hold on;
            grid on;
            y_arr = bet_lir_arr(i, j, :);
            plot(r_arr, y_arr(:), '-*', 'Color', colors(pcnt,:))
            ylabel('\lambda_i')
            xlabel('Normalized radius')
            ylim(inflow_lims)

            nfig = 2;
            fig = figure(nfig);

            subplot(2, 1, 1)
            hold on;
            grid on;
            y_arr = bet_lr_arr(i, j, :);
            plot(r_arr, y_arr(:), '-*', 'Color', colors(pcnt,:))
            ylabel('\lambda')
            xlabel('Normalized radius')
            ylim(inflow_lims)

            subplot(2, 1, 2)
            hold on;
            grid on;
            y_arr = bet_aoa_arr(i, j, :);
            plot(r_arr, y_arr(:) * 180 / pi, '-*', 'Color', colors(pcnt,:))
            ylabel('\alpha deg');
            xlabel('Normalized radius')
            ylim([-20, +10])
        end
    end

    nfig = 1;
    fig = figure(nfig);
    sgtitle(['T ' num2str(round(tot_thrust, 2)) ' N'])
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, filename);      

    nfig = 2;
    fig = figure(nfig);
    sgtitle(['T ' num2str(round(tot_thrust, 2)) ' N'])
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, filename);    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dcollpitch = 0;
    % str1 = 'bet_coax_eta_thrust';
    % bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1); 

    nfig = 3;
    fig = figure(nfig);
    hold on
    grid on
    % plot(lc_arr ./ mt_lambda_h, bet_lambdai_arr ./ mt_lambda_h, '-*r')
    % plot(lc_arr ./ mt_lambda_h, mt_lambdai_arr ./ mt_lambda_h, '-*g')
    plot(bet_lambdac_arr ./ mt_lambda_h, bet_lambdai_arr ./ mt_lambda_h, '-*r')
    plot(bet_lambdac_arr ./ mt_lambda_h, mt_lambdai_arr ./ mt_lambda_h, '-*g')
    xlabel('\lambda_c / \lambda_h^{MT}')
    ylabel('Mean \lambda_i / \lambda_h^{MT}')
    ylim([0, 2])
    yyaxis right
    hold on;
    grid on;
    % plot(lc_arr ./ mt_lambda_h, bet_power_arr ./ mt_power_h, '-*b')
    % plot(lc_arr ./ mt_lambda_h, mt_power_arr ./ mt_power_h, '-*k')
    plot(bet_lambdac_arr ./ mt_lambda_h, bet_power_arr ./ mt_power_h, '-*b')
    plot(bet_lambdac_arr ./ mt_lambda_h, mt_power_arr ./ mt_power_h, '-*k')
    xlabel('\lambda_c / \lambda_h^{MT} ')
    ylabel('P / P_h^{MT}')
    ylim([0, 4])
    % xlim([0, 1.5])
    legend('BET ind inflow', 'MT ind inflow', ...
        'BET tot power', 'MT ind power', ...
        'Location', 'NorthWest');
    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k'; 

    title([...
        'T ' ...
        num2str(round(tot_thrust, 2)) ...
        ' N, \lambda_h^{MT} '...
        num2str(round(mt_lambda_h, 4)) ...
        ' , P_h^{MT} '...
        num2str(round(mt_power_h, 2)) ...
        ' W' ...
    ]);

    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, filename);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 4;
    fig = figure(nfig);
    hold on
    grid on
    plot(lc_arr ./ mt_lambda_h, bet_omega_arr * rads2rpm, '-*')
    xlabel('\lambda_c / \lambda_h^{MT} ')
    ylabel('\Omega^{BET} RPM')
    ylim([0, 5500])

    title(['Thrust ' num2str(round(tot_thrust, 2)) ' N'])
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, filename);  

    close all;          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    
    return

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nfig = 4;
    fig = figure(nfig);
    hold on;
    grid on;

    if norm(diff(db.mt_lambda_h)) < 10^-8
        disp('norm(diff(db.mt_lambda_h)) < 10^-8')
        mt_lambda_h = db.mt_lambda_h(1)
    else
        error('norm(diff(db.mt_lambda_h)) >= 10^-8')
        % for i = 1:length(db.thrust_arr)
        %     Z(i, :) = db.bet_lambdai_arr(i, :) ./ db.mt_lambda_h(i);
        % end
    end

    x = db.thrust_arr;
    y = db.lc_arr ./ mt_lambda_h;
    [X, Y] = meshgrid(x, y);
    % Z = transpose(db.mt_v_arr); %  ./ db.P_h;        
    Z = transpose(db.bet_lambdai_arr ./ mt_lambda_h);
    surf(X, Y, Z, Cr);
    Z = transpose(db.mt_lambdai_arr ./ mt_lambda_h);
    surf(X, Y, Z, Cg);

    % colorbar
    % title('coax and sbs power');
    xlabel('Thrust N');
    ylabel('\lambda_c / \lambda_h^{MT}');
    % zlabel('Induced velocity m/s');
    zlabel('\lambda / \lambda_h^{MT}');
    zlim([0, 2]);
    legend('BET', 'MT')
    view(35, 35);
    grid on;    
    
    filename = ['img/' str1 '_' str2 '_' num2str(nfig)];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, [filename '.jpg']);  
    savefig(fig, [filename '.fig']);  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfig = 5;
    fig = figure(nfig);
    hold on;
    grid on;

    if norm(diff(db.mt_lambda_h)) < 10^-8
        disp('norm(diff(db.mt_lambda_h)) < 10^-8')
        mt_lambda_h = db.mt_lambda_h(1)
    else
        error('norm(diff(db.mt_lambda_h)) >= 10^-8')
    end
    % if norm(diff(db.mt_power_h)) < 10^-8
    %     disp('norm(diff(db.mt_power_h)) < 10^-8')
    %     mt_power_h = db.mt_lambda_h(1)
    % else
    %     error('norm(diff(db.mt_power_h)) >= 10^-8')
    % end

    x = db.thrust_arr;
    y = db.lc_arr ./ mt_lambda_h;
    [X, Y] = meshgrid(x, y);
    % Z = transpose(db.mt_v_arr); %  ./ db.P_h;        
    Z = transpose(db.bet_power_arr);
    surf(X, Y, Z, Cr);
    Z = transpose(db.mt_power_arr);
    surf(X, Y, Z, Cg);

    % colorbar
    % title('coax and sbs power');
    xlabel('Thrust N');
    ylabel('\lambda_c / \lambda_h^{MT}');
    zlabel('Power W');
    legend('BET', 'MT')
    % zlim([0, 4]);
    view(35, 35);
    grid on;    
    
    filename = ['img/' str1 '_' str2 '_' num2str(nfig)];
    set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 8], ...
        'PaperUnits', 'centimeters', 'PaperSize', [12, 8]);
    saveas(fig, [filename '.jpg']);  
    savefig(fig, [filename '.fig']);  

    close all

end
