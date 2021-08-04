function kde_rotor_TQP_test()
    clear all
    clc
    format long
    close all

    % Test the replacement of interp1 by a quadratic fit   
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
    uomega_arr = [00.00 1920  2500  3240  3980  4520  5260  5880 ] * (pi / 30);
    urotortype = 'KDE6213XF185_KDECF185DP'
    umodeltype = 'interp1'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'r-*')
    umodeltype = 'quadratic'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'b-*');
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') 
    uomega_arr = [00.00 1680  2280  2820  3480  3900  4380  4800 ] * (pi / 30);
    urotortype = 'KDE6213XF185_KDECF245DP'
    umodeltype = 'interp1'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'r-*')
    umodeltype = 'quadratic'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'b-*');
    
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')     
    % uomega_arr = [00.00 1140  1620  2040  2520  2940  3240  3600 ] * (pi / 30);
    uomega_arr = [00.00 1660  2220  2880  3420  3900  4380  4920 ] * (pi / 30);
    urotortype = 'KDE8218XF120_KDECF305DP'
    umodeltype = 'interp1'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'r-*')
    umodeltype = 'quadratic'
    [Tp, Qp, Pp] = plot_data(uomega_arr, urotortype, umodeltype, 'b-*');    
    
    str1 = 'kde_rotor_TQP_test';
    str2 = '185DP_245DP_305DP';
    
    nfig = 1;
    fig = figure(nfig);
    % set(fig,'units', 'centimeters', 'position', [0, 0, 20, 10]);       
    legend('interp1', 'quadratic', 'Location', 'northwest') 
    filename = ['img/' str1 '_' str2 '_' num2str(nfig) '.jpg'];
    saveas(fig, filename);  
    
    close all    
    
    function [Tp, Qp, Pp] = plot_data(omega_arr, rotortype, modeltype, linespce)
        
        rpm2rads = pi / 30;
        rads2rpm = 30 / pi; 
        
        for i=1:length(omega_arr)
            omega = omega_arr(i);        
            [T, Q, P] = kde_rotor_TQP(omega, rotortype, modeltype);
            Trotor_arr(i) = T;
            Qrotor_arr(i) = Q;
            Protor_arr(i) = P;
        end

        x = omega_arr;
        y = Trotor_arr;
        n = 2;
        [Tp, S]         = polyfit(x,y,n);  
        Terr = S.normr
        % [y_fit, delta]  = polyval(Tp,x,S);
        
        x = omega_arr;
        y = Qrotor_arr;
        n = 2;
        [Qp, S]         = polyfit(x,y,n); 
        Qerr = S.normr
        % [y_fit, delta]  = polyval(Qp,x,S);
        
        x = omega_arr;
        y = Protor_arr;
        n = 2;
        [Pp, S]         = polyfit(x,y,n);
        Perr = S.normr
        % [y_fit, delta]  = polyval(Pp,x,S);
        
        fig = figure(1);
        % sgtitle(rotortype);
        subplot(3, 1, 1)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Trotor_arr, linespce);
        xlabel('angvel RPM')
        ylabel('Thrust N')
        
        subplot(3, 1, 2)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Qrotor_arr, linespce);
        xlabel('angvel RPM')
        ylabel('Torque Nm')
        
        subplot(3, 1, 3)
        hold on;
        grid on;
        plot(omega_arr.*rads2rpm, Protor_arr, linespce);
        xlabel('angvel RPM')
        ylabel('Power W')    

%        % Plot polyfit
%        for i=1:length(omega_arr)
%            omega = omega_arr(i);  
%        
%            omega = abs(omega);  
%            omega2 = omega.*omega;       

%            abc = Tp;
%            a = abc(1);
%            b = abc(2);
%            c = abc(3);
%            Tpolyfit = a*omega2 + b*omega + c;
%            
%            abc = Qp;
%            a = abc(1);
%            b = abc(2);
%            c = abc(3);
%            Qpolyfit = a*omega2 + b*omega + c;
%            
%            abc = Pp;
%            a = abc(1);
%            b = abc(2);
%            c = abc(3);
%            Ppolyfit = a*omega2 + b*omega + c; 

%            pfit_Trotor_arr(i) = Tpolyfit;
%            pfit_Qrotor_arr(i) = Qpolyfit;
%            pfit_Protor_arr(i) = Ppolyfit;
%        end    
    end    
end

