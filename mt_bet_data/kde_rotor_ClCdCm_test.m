function kde_rotor_ClCdCm_test()

    clear all
    close all
    format compact
    clc
    
    mach_arr = [0.0000 0.2849 0.2850 0.3649 0.3650 0.6000];
    aoa_arr = -180:1:+180;
    for i=1:length(aoa_arr)
        aoa = aoa_arr(i);
        mach = mach_arr(end);
        [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, deg2rad(aoa), 'NACA5407_Station6');
        
        Cl_arr(i) = Cl;
        Cd_arr(i) = Cd;
        Cm_arr(i) = Cm;
    end    
    
    rotortype = 'NACA5407';
    [y_arr, chord_arr, theta_arr] = kde_rotor_geometry([rotortype 'shape']);
    NACA5407x_arr = chord_arr - 1/4;
    NACA5407y_arr = theta_arr;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(1);
    title([rotortype ': CP assumed to be at 1/4 chord' ]);
    hold on;
    grid on;
    axis square
    axis equal
    plot(NACA5407x_arr, NACA5407y_arr, 'k-*');
    daoa = 22.5;
    alpha_arr = -180:daoa:(+180-daoa);
    for i=1:length(alpha_arr)
        aoa = alpha_arr(i);
        mach = mach_arr(end);
        [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, deg2rad(aoa), 'NACA5407_Station6');
        P = 1;       
        Clx = -P*cos(deg2rad(aoa));     Cly = -P*sin(deg2rad(aoa));
        plot([0 Clx], [0 Cly], 'r-*');
        % lbl = ['aoa ' num2str(aoa) ' Cl ' num2str(round(Cl, 2))];
        lbl = ['Cl ' num2str(round(Cl, 2))];
        text(Clx+0.05, Cly-0.05, lbl)
        % plot(Cdy, Cdy, 'g-*');
        % plot(Cmy, Cmy, 'b-*');
        
    end    
    xlim([-1 1])
    ylim([-1 1])
    xlabel('x axis m')
    ylabel('y axis m')    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(2);
    title([rotortype ': CP assumed to be at 1/4 chord' ]);
    hold on;
    grid on;
    axis square
    axis equal
    plot(NACA5407x_arr, NACA5407y_arr, 'k-*');
    daoa = 22.5;
    alpha_arr = -180:daoa:(+180-daoa);
    for i=1:length(alpha_arr)
        aoa = alpha_arr(i);
        mach = mach_arr(end);
        [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, deg2rad(aoa), 'NACA5407_Station6');
        Cl = Cd; P = 1;       
        Clx = -P*cos(deg2rad(aoa));     Cly = -P*sin(deg2rad(aoa));
        plot([0 Clx], [0 Cly], 'r-*');
        % lbl = ['aoa ' num2str(aoa) ' Cl ' num2str(round(Cl, 2))];
        lbl = ['Cd ' num2str(round(Cl, 2))];
        text(Clx+0.05, Cly-0.05, lbl)
        % plot(Cdy, Cdy, 'g-*');
        % plot(Cmy, Cmy, 'b-*');
        
    end    
    xlim([-1 1])
    ylim([-1 1])
    xlabel('x axis m')
    ylabel('y axis m')    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(3);
    title([rotortype ': CP assumed to be at 1/4 chord' ]);
    hold on;
    grid on;
    axis square
    axis equal
    plot(NACA5407x_arr, NACA5407y_arr, 'k-*');
    daoa = 22.5;
    alpha_arr = -180:daoa:(+180-daoa);
    for i=1:length(alpha_arr)
        aoa = alpha_arr(i);
        mach = mach_arr(end);
        [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, deg2rad(aoa), 'NACA5407_Station6');
        Cl = Cm; P = 1;       
        Clx = -P*cos(deg2rad(aoa));     Cly = -P*sin(deg2rad(aoa));
        plot([0 Clx], [0 Cly], 'r-*');
        % lbl = ['aoa ' num2str(aoa) ' Cl ' num2str(round(Cl, 2))];
        lbl = ['Cm ' num2str(round(Cl, 2))];
        text(Clx+0.05, Cly-0.05, lbl)
        % plot(Cdy, Cdy, 'g-*');
        % plot(Cmy, Cmy, 'b-*');
        
    end    
    xlim([-1 1])
    ylim([-1 1])
    xlabel('x axis m')
    ylabel('y axis m')            

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(4);
    % sgtitle(rotortype);
    subplot(3, 1, 1)
    hold on;
    grid on;
    set(gca,'FontSize', 14)
    plot(aoa_arr, Cl_arr, 'r');
    xlabel('Angle of attack deg')
    ylabel('Cl')
    % yyaxis right
    % plot(aoa_arr, Cla_arr, 'b-*');
    % xlabel('Angle of attack deg')
    % ylabel('Cla')
    
    subplot(3, 1, 2)
    hold on;
    grid on;
    set(gca,'FontSize', 14)
    plot(aoa_arr, Cd_arr, 'r');
    xlabel('Angle of attack deg')
    ylabel('Cd')
    
    subplot(3, 1, 3)
    hold on;
    grid on;
    set(gca,'FontSize', 14)
    plot(aoa_arr, Cm_arr, 'r');   
    xlabel('Angle of attack deg')
    ylabel('Cm')
    % legend('interp1', 'quadratic', 'Location', 'northwest')    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure(5);
    hold on;
    grid on;
    aoa_arr = [0, 5, 10];
    mach_arr = 0:0.1:0.6;
    % mach_arr = [0.0000 0.2849 0.2850 0.3649 0.3650 0.6000 ];
    Cl_arr = [];
    for i=1:length(aoa_arr)
        for j=1:length(mach_arr)
            mach = mach_arr(j);
            aoa = aoa_arr(i);
            [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, deg2rad(aoa), 'NACA5407_Station6');
            Cl_arr(i, j) = Cl;
        end  
    end      
    plot(mach_arr, Cl_arr(1, :), 'r-*');
    plot(mach_arr, Cl_arr(2, :), 'g-*');
    plot(mach_arr, Cl_arr(3, :), 'b-*');
    xlabel('Mach number')
    ylabel('Cl')
    legend('0 deg', '5 deg', '10 deg')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str1 = 'kde_rotor_ClCdCm_test';
    str2 = rotortype;
    
    nfig = 1;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
    nfig = 2;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
    nfig = 3;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
    nfig = 4;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);      
    
    nfig = 5;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);        
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);     
    
    close all
end

