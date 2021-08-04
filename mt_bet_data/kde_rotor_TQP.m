function [T, Q, P] = kde_rotor_TQP(...
    omega       , ... % 
    rotortype   , ... %
    modeltype     ... %
    )
    
    persistent data  
    if isempty(data)   
        rpm2rads   = pi/30;
        rads2rpm   = 30/pi;
        in2cm      = 2.54;     % 1in = 2.54cm
        ft2cm      = 30.48;    % 1ft = 30.48cm
        cm2m       = 1/100;
        
        effcy      = 0.80;  % Assumed efficiency for KDE ESC and Motors (check kde_motor_QIPe_test__2.jpg)
        
        % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde6213xf-185
        % https://www.kdedirect.com/collections/multi-rotor-propeller-blades/products/kde-cf185-dp                
        data.KDE6213XF_CF185.name       = "KDE6213XF-185 motor and KDE-CF185-DP @ 8S voltage";
        data.KDE6213XF_CF185.voltage    = 30.8;     % 8S = 30.8V to 34.8 V
        data.KDE6213XF_CF185.throttle   = [00.00 25.00 37.50 50.00 62.50 75.00 87.50 100.0]./100;
        data.KDE6213XF_CF185.current    = [00.00 00.90 01.80 03.10 05.10 07.90 11.70 15.50];            % A
        data.KDE6213XF_CF185.power      = [00.00 31.00 62.00 107.0 177.0 274.0 407.0 539.0].*effcy;     % W
        data.KDE6213XF_CF185.thrust     = [00.00 04.12 07.26 11.57 16.48 22.36 30.11 38.34];            % N
        data.KDE6213XF_CF185.angvel     = [00.00 1920  2500  3240  3980  4520  5260  5880 ].*rpm2rads;  % rad/s
        data.KDE6213XF_CF185.torque     = data.KDE6213XF_CF185.power ./ data.KDE6213XF_CF185.angvel;
        data.KDE6213XF_CF185.torque(1)  = 00.00;    % Solve 0/0 = NaN problem

        % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde6213xf-185
        % https://www.kdedirect.com/collections/multi-rotor-propeller-blades/products/kde-cf245-dp
        % https://cdn.shopify.com/s/files/1/0496/8205/files/
        % KDE_Direct_XF_CF_Brushless_Performance_Testing_KDE6213XF-185.pdf?665
        data.KDE6213XF_CF245.name       = "KDE6213XF-185 motor and KDE-CF245-DP @ 8S voltage";
        data.KDE6213XF_CF245.voltage    = 30.8;     % 8S = 30.8V to 34.8 V
        data.KDE6213XF_CF245.throttle   = [00.00 25.00 37.50 50.00 62.50 75.00 87.50 100.0]./100;
        data.KDE6213XF_CF245.current    = [00.00 01.80 03.90 07.70 013.7 019.7 029.2 038.8];            % A
        data.KDE6213XF_CF245.power      = [00.00 62.00 135.0 267.0 476.0 685.0 1016  1350 ].*effcy;     % W
        data.KDE6213XF_CF245.thrust     = [00.00 09.02 16.77 26.67 39.32 50.50 64.63 79.73];            % N
        data.KDE6213XF_CF245.angvel     = [00.00 1680  2280  2820  3480  3900  4380  4800 ].*rpm2rads;  % rad/s
        data.KDE6213XF_CF245.torque     = data.KDE6213XF_CF245.power ./ data.KDE6213XF_CF245.angvel;
        data.KDE6213XF_CF245.torque(1)  = 00.00;    % Solve 0/0 = NaN problem

%        % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde8218xf-120
%        % https://www.kdedirect.com/collections/multi-rotor-propeller-blades/products/kde-cf305-dp
%        data.KDE8218XF_CF305.name       = "KDE8218XF-120 motor and KDE-CF305-DP @ 8S voltage";
%        data.KDE8218XF_CF305.voltage    = 30.8;     % 8S = 30.8V to 34.8 V
%        data.KDE8218XF_CF305.throttle   = [00.00 25.00 37.50 50.00 62.50 75.00 87.50 100.0]./100;
%        data.KDE8218XF_CF305.current    = [00.00 01.70 04.10 08.70 14.60 22.80 32.60 44.20];            % A
%        data.KDE8218XF_CF305.power      = [00.00 59.00 142.0 302.0 508.0 793.0 1134  1538 ].*effcy;     % W
%        data.KDE8218XF_CF305.thrust     = [00.00 10.59 20.79 34.91 50.60 68.35 86.69 109.7];            % N
%        data.KDE8218XF_CF305.angvel     = [00.00 1140  1620  2040  2520  2940  3240  3600 ].*rpm2rads;  % rad/s
%        data.KDE8218XF_CF305.torque     = data.KDE8218XF_CF305.power ./ data.KDE8218XF_CF305.angvel;
%        data.KDE8218XF_CF305.torque(1)  = 00.00;    % Solve 0/0 = NaN problem
        
        % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde8218xf-120
        % https://www.kdedirect.com/collections/multi-rotor-propeller-blades/products/kde-cf305-dp
        data.KDE8218XF_CF305.name       = "KDE8218XF-120 motor and KDE-CF305-DP @ 12S voltage";
        data.KDE8218XF_CF305.voltage    = 46.2;     % 12S = 46.2V to 52.2 V

        data.KDE8218XF_CF305.throttle   = [00.00 25.00 37.50 50.00 62.50 75.00 87.50 100.0]./100;
        data.KDE8218XF_CF305.current    = [00.00 03.10 07.30 15.00 26.30 39.00 58.40 77.10];            % A
        data.KDE8218XF_CF305.power      = [00.00 161.0 381.0 783.0 1372  2035  3048  4024 ].*effcy;     % W
        data.KDE8218XF_CF305.thrust     = [00.00 21.77 40.21 66.39 96.20 123.8 161.1 199.8];            % N
        data.KDE8218XF_CF305.angvel     = [00.00 1660  2220  2880  3420  3900  4380  4920 ].*rpm2rads;  % rad/s
        data.KDE8218XF_CF305.torque     = data.KDE8218XF_CF305.power ./ data.KDE8218XF_CF305.angvel;
        data.KDE8218XF_CF305.torque(1)  = 00.00;    % Solve 0/0 = NaN problem        

    end
    
    omega = abs(omega); 
    omega2 = omega.*omega;   

    T = NaN;
    Q = NaN; 
    P = NaN;

    if strcmp(modeltype, 'interp1')
        switch rotortype
            case 'KDE6213XF185_KDECF185DP'
                angvel_arr  = data.KDE6213XF_CF185.angvel;
                thrust_arr  = data.KDE6213XF_CF185.thrust;
                torque_arr  = data.KDE6213XF_CF185.torque;
                power_arr   = data.KDE6213XF_CF185.power;
                P = interp1(angvel_arr, power_arr, omega, 'linear');
            case 'KDE6213XF185_KDECF245DP'
                angvel_arr  = data.KDE6213XF_CF245.angvel;
                thrust_arr  = data.KDE6213XF_CF245.thrust;
                torque_arr  = data.KDE6213XF_CF245.torque;
                power_arr   = data.KDE6213XF_CF245.power;
            case 'KDE8218XF120_KDECF305DP'
                angvel_arr  = data.KDE8218XF_CF305.angvel;
                thrust_arr  = data.KDE8218XF_CF305.thrust;
                torque_arr  = data.KDE8218XF_CF305.torque;
                power_arr   = data.KDE8218XF_CF305.power;                           
            otherwise
                rotortype
                error('Unrecognized rotortype')
        end 
        T = interp1(angvel_arr, thrust_arr, omega, 'linear');
        Q = interp1(angvel_arr, torque_arr, omega, 'linear');
        P = interp1(angvel_arr, power_arr, omega, 'linear');          
    end 
    
    if strcmp(modeltype, 'quadratic')
        % Check kde_rotor_TQP_test.m
        % for more info
        
        % Due to polyfit, some values of P = ax2 + bx + c were going negative
        % at values of x = 1000 RPM, that is why abs(P) was added
        switch rotortype
            case 'KDE6213XF185_KDECF185DP'
                Tp  = [0.000102258581588  -0.001580898616871   0.171712520586930];
                Qp  = [0.000001448687037   0.000317623787623   0.003479420754307];
                Pp  = [0.001663184872615  -0.302823815853915   6.715407982612643];
            case 'KDE6213XF185_KDECF245DP'
                Tp  = [0.000327734408560  -0.008052086345186   0.110101785370729];
                Qp  = [0.000009331395472  -0.000220228085597   0.012167581471604];
                Pp  = [0.006859613056730  -1.290662285012739  19.174913498116599];
            case 'KDE8218XF120_KDECF305DP'
                % @ 8S voltage
                % Tp  = [0.000805810278059  -0.019115751572062   0.578189922634889];
                % Qp  = [0.000023705473638   0.000102352177260   0.023295938649086];
                % Pp  = [0.013679189229579  -1.895399103022898  27.212836053358288];
                % @ 12S voltage
                Tp  = [0.000777342382212  -0.010897398956374   0.061642566031112];
                Qp  = [0.000025691619747  -0.000046150086881  -0.003557371563189];
                Pp  = [0.019694910562288  -3.684000544184555  52.426122573969295];
            otherwise
                rotortype
                error('Unrecognized rotortype')
        end      
        T   = abs( Tp(1)*omega2 + Tp(2)*omega + Tp(3) );
        Q   = abs( Qp(1)*omega2 + Qp(2)*omega + Qp(3) );
        P   = abs( Pp(1)*omega2 + Pp(2)*omega + Pp(3) );
    end
      
end

