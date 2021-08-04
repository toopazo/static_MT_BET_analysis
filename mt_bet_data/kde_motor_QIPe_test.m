function kde_motor_QIPe_test()

    clear all
    close all
    format compact
    format short
    clc   
    
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;         
    
    npoints         = 5;
    motortype       = 'KDE8218XF';
    % motortype       = 'KDE6213XF';
    throttle_arr    = linspace(0, 1, npoints);   
    voltage_arr     = linspace(10, 60, npoints);
    omega_arr       = 0:5:1000;
    for i=1:npoints
        throttle = throttle_arr(i);
        voltage = voltage_arr(i);
        for j=1:length(omega_arr)
            omega = omega_arr(j);
            % [Q, I, V, P, e] = kde_motor_QIPe(omega, voltage, motortype, 'voltage');
            [Q, I, V, P, e] = kde_motor_QIPe(omega, throttle, motortype, 'throttle');
            Qarr(j) = Q;
            Iarr(j) = I;
            Varr(j) = V;
            Parr(j) = P;
            earr(j) = e;
        end

        fig = figure(1);
        sgtitle(motortype);
        subplot(3, 1, 1)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Qarr, '-*');
        xlabel('angvel RPM')
        ylabel('Torque Nm')

        subplot(3, 1, 2)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Iarr, '-*');
        xlabel('angvel RPM')
        ylabel('Current A')

        subplot(3, 1, 3)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Varr, '-*');
        xlabel('angvel RPM')
        ylabel('Voltage V')  
        for iN = 1:npoints
            % legendCell{iN} = num2str(voltage_arr(iN));
            legendCell{iN} = num2str(throttle_arr(iN));
        end
        legend(legendCell, 'Location', 'northwest');          
        
        fig = figure(2);
        sgtitle(motortype);
        subplot(2, 1, 1)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Parr, '-*');
        xlabel('angvel RPM')
        ylabel('Shaft Power W')

        subplot(2, 1, 2)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, earr, '-*');
        xlabel('angvel RPM')
        ylabel('Motor efficiency')
        for iN = 1:npoints
            % legendCell{iN} = num2str(voltage_arr(iN));
            legendCell{iN} = num2str(throttle_arr(iN));
        end
        legend(legendCell, 'Location', 'northwest');        
    
    end
    
    str1 = 'kde_motor_QIPe_test';
    str2 = motortype;
    
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
    
    close all
end
