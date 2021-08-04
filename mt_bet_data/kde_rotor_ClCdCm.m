function [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, aoa, rotortype)

    persistent NACA5407_Station6
    if isempty(NACA5407_Station6)
        % filename = [bpath '/' 'NACA5407_Station6'];
        filename = 'mt_bet_data/NACA5407_Station6';
        st = load(filename);
        NACA5407_Station6 = st.NACA5407_Station6;
    end        

    % mach_arr = [0.0000 0.2849 0.2850 0.3649 0.3650 0.6000];
    % mach = 0.3649;
    if mach > 0.6000
        mach = 0.6000;
    end
    aoa = rad2deg(aoa);
%    if aoa > 180
%        aoa = 180;
%    end
%    if aoa < -180
%        aoa = -180;
%    end    

    if strcmp(rotortype, 'NACA5407')
        rotortype = 'NACA5407_Station6';
    end
    if strcmp(rotortype, 'KDECF245DP')
        rotortype = 'NACA5407_Station6';
    end
    if strcmp(rotortype, 'KDECF245TP')
        rotortype = 'NACA5407_Station6';
    end
    if strcmp(rotortype, 'KDECF245HP')
        rotortype = 'NACA5407_Station6';
    end
    if strcmp(rotortype, 'KDECF305DP')
        rotortype = 'NACA5407_Station6';
    end
    if strcmp(rotortype, 'NACA5407_Station6')
        rpm2rads   = pi/30;
        rads2rpm   = 30/pi;
        in2cm      = 2.54;     % 1in = 2.54cm
        ft2cm      = 30.48;    % 1ft = 30.48cm
        cm2m       = 1/100;
        
        use_fixed_mach = true;
        if use_fixed_mach
            % [Cl, Cd, Cm] @ mach = 0.28
            col = 3;
            
            Xq = aoa;
            X = NACA5407_Station6.ANGLE_Cl(:, col);
            V = NACA5407_Station6.Cl(:, col);
            Cl = interp1(X, V, Xq, 'spline');   
            
            Xq = aoa;
            X = NACA5407_Station6.ANGLE_Cd(:, col);
            V = NACA5407_Station6.Cd(:, col);
            Cd = interp1(X, V, Xq, 'spline'); 
            
            Xq = aoa;
            X = NACA5407_Station6.ANGLE_Cm(:, col);
            V = NACA5407_Station6.Cm(:, col);
            Cm = interp1(X, V, Xq, 'spline');                
        else
            % Vq = interp2(X,Y,V,Xq,Yq);        
            Xq = mach;
            Yq = aoa;
            X = NACA5407_Station6.MACH_Cl;
            Y = NACA5407_Station6.ANGLE_Cl;
            V = NACA5407_Station6.Cl;
            Cl = interp2(X, Y, V, Xq, Yq, 'linear');
                    
            Xq = mach;
            Yq = aoa;
            X = NACA5407_Station6.MACH_Cd;
            Y = NACA5407_Station6.ANGLE_Cd;
            V = NACA5407_Station6.Cd;
            Cd = interp2(X, Y, V, Xq, Yq, 'linear');     
            
            Xq = mach;
            Yq = aoa;
            X = NACA5407_Station6.MACH_Cm;
            Y = NACA5407_Station6.ANGLE_Cm;
            V = NACA5407_Station6.Cm;
            Cm = interp2(X, Y, V, Xq, Yq, 'linear');   
        end
        
        % Date: Sept30th 
        % bet_main_omega.m produces a very nice fit as shown in 
        % static_RotorInflowModel/img/bet_main_omega_KDECF245DP_1.jpg
        % Thrust predicted by BET shows an error below 0.5% when compared
        % with experimental data (quadratic fit of AneCha data)
        % However the torque calculated by BET is under-predicting by -10%
        
        % Important parameters are
        %   rotortype   = 'KDECF245DP';
        %   rho         = 1.1388;
        %   lambda_c    = 0.0;
        %   mu          = 0.0;
        %   collpitch   = -deg2rad(0.0);

        % This are the sames parameters as in 
        %   rotortype = 'KDECF245DP'
        %   [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype)
        
        % By manually tunning Cd I was able to bring the torque error 
        % down to -2.5% % and keep the thrust error at 0.2%      
        Cd = Cd*2.0;
        % Check results in bet_main_omega_KDECF245DP_1.jpg
    end        
   
    if (isnan(Cl) == true) || (isfinite(Cl) == false) || (isreal(Cl) == false)
        Cl 
        mach 
        aoa 
        rotortype
        error('Cl is NaN, Inf or Complex')
    end
    if (isnan(Cd) == true) || (isfinite(Cd) == false) || (isreal(Cd) == false)
        Cd 
        mach 
        aoa 
        rotortype
        error('Cd is NaN, Inf or Complex')
    end
    if (isnan(Cm) == true) || (isfinite(Cm) == false) || (isreal(Cm) == false)
        Cm 
        mach 
        aoa 
        rotortype
        error('Cm is NaN, Inf or Complex')
    end      
end

