function blade_st = blade_model(rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
    
    rpm2rads    = (pi/30);  % 60 RPM = 1 RPsecond = 1 Hz = 2*pi rad/s
    in2m        = (0.0254); 
    ft2m        = (0.3048);
    lb2kg       = (0.453592);
    
    Vsound      = 340.2941;

    % In the case of hover => lambda_c = mu = 0
    %   lambda_c = mu*tan(TPP_alpha)
    %   TPP_alpha = atan(lambda_c/mu)
    %   atan(0/0) = NaN
    %   atan(1/0) = pi/2 = 90 deg
    %   tan(atan(1/0)) = 1.633123935319537e+16
    %   0*tan(atan(1/0)) = 0
    %   blade_st.TPP_alpha  = atan((blade_st.lambda_c + 10^-6)/(blade_st.mu));

    % Blade geometry (pitch, chord, solidity)
    blade_st.rotortype = rotortype;
        
    known_rotortype = false;
    if strfind(rotortype,'constant_twist_blade') == 1
        known_rotortype = true;
        
        % UH-60A                    
        [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry();    
    
        
        % 1) General rotor properties 
        blade_st.Nb         = Nb;
        blade_st.R          = R;   
        blade_st.rotArea    = pi*blade_st.R^2;    
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;
        
        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;    
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan( (blade_st.lambda_c + 10^-6) / (blade_st.mu) ); 
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);
        
        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        blade_st.Cla        = Cla;
        blade_st.Cl0        = 0.0;
        blade_st.Cd0        = Cd0;
        blade_st.d1         = 0;
        blade_st.d2         = 0;   
        % blade_st.c81table.get_ClCdCm = NaN;
        
        % 5) Pitch, chord and solidity distribution
        % nsections = 16 => 
        %                   len(psi_arr) = 17
        %                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000        
        nsections           = 16;
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections;
        blade_st.y_arr      = linspace(0, R, nsections + 1); 
        theta0str = extractBetween(rotortype,'constant_twist_blade_','rad');
        if isempty(theta0str)
            theta0 = deg2rad(22.0335);    % deg2rad(22.0);
        else
            theta0 = str2double(theta0str);
        end      
        blade_st.theta_arr  = theta0 .* blade_st.ones_arr + collpitch;                     
        blade_st.chord_arr  = chord .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);      
        
    end    
    if strcmp(rotortype, 'ideally_twisted_blade')    
        known_rotortype = true;
        
        % UH-60A                    
        [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry();             
        
        % 1) General rotor properties       
        blade_st.Nb         = Nb;
        blade_st.R          = R;   
        blade_st.rotArea    = pi*blade_st.R^2;    
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;
        
        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;    
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan( (blade_st.lambda_c + 10^-6) / (blade_st.mu) ); 
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);
        
        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        blade_st.Cla        = Cla;
        blade_st.Cl0        = 0.0;
        blade_st.Cd0        = Cd0;
        blade_st.d1         = 0;
        blade_st.d2         = 0;   
        % blade_st.c81table.get_ClCdCm = NaN;
        
        % 5) Pitch, chord and solidity distribution
        % nsections = 16 => 
        %                   len(psi_arr) = 17
        %                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000        
        nsections           = 16;   
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections; 
        blade_st.y_arr      = linspace(0, R, nsections + 1); 
        itb_theta0          = deg2rad(15.020);     % deg2rad(15.2);
        blade_st.theta_arr  = itb_theta0 .* blade_st.ones_arr ./ linspace(0, R, nsections+1);
        blade_st.theta_arr  = blade_st.theta_arr + collpitch;
        theta_max           = deg2rad(70);
        ind_arr             = find(blade_st.theta_arr >= theta_max);
        for ind = ind_arr
            blade_st.theta_arr(ind) = theta_max;
        end                     
        blade_st.chord_arr  = chord .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);              
    end
    if strcmp(rotortype, 'linearly_twisted_blade')  
        known_rotortype = true;    

        % UH-60A                    
        [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry();    

        ltb_theta0          = deg2rad(37.560);     % deg2rad(35.0);
        ltb_theta1          = deg2rad(17.000);     % deg2rad(18.0);             
        
        % 1) General rotor properties        
        blade_st.Nb         = Nb;
        blade_st.R          = R;   
        blade_st.rotArea    = pi*blade_st.R^2;    
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;
        
        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;    
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan( (blade_st.lambda_c + 10^-6) / (blade_st.mu) ); 
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);
        
        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        blade_st.Cla        = Cla;
        blade_st.Cl0        = 0.0;
        blade_st.Cd0        = Cd0;
        blade_st.d1         = 0;
        blade_st.d2         = 0;   
        % blade_st.c81table.get_ClCdCm = NaN;
        
        % 5) Pitch, chord and solidity distribution
        % nsections = 16 => 
        %                   len(psi_arr) = 17
        %                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000        
        nsections           = 16;   
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections; 
        blade_st.y_arr      = linspace(0, R, nsections + 1); 
        blade_st.theta_arr  = linspace(ltb_theta0, ltb_theta1, nsections + 1) + collpitch;
        blade_st.chord_arr  = chord .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);
    end
    if strcmp(rotortype, 'KDECF245DP') || strcmp(rotortype, 'KDECF245TP') || strcmp(rotortype, 'KDECF245HP')  
        known_rotortype     = true;        

        [y_arr, chord_arr, theta_arr] = kde_rotor_geometry(rotortype);

        if strcmp(rotortype, 'KDECF245DP') 
            nblades = 2;
        end
        if strcmp(rotortype, 'KDECF245TP') 
            nblades = 3;
        end
        if strcmp(rotortype, 'KDECF245HP')  
            nblades = 6;
        end      

        % 1) General rotor properties       
        blade_st.Nb         = nblades;
        blade_st.R          = y_arr(end);
        blade_st.rotArea    = pi*blade_st.R^2;
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;

        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;       
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan( (blade_st.lambda_c + 10^-6) / (blade_st.mu) ); 
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);

        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        use_c81table = true;
        if use_c81table
            blade_st.Cla        = NaN; % 2*pi*0.9488;  % clbar = 0.5987281447 => theta approx 5.5deg   
            blade_st.Cl0        = NaN; % 0.6158;       % Cl at zero angle of attack
            blade_st.Cd0        = NaN; % blade_st.Cla * deg2rad(5) / 29;;
            blade_st.d1         = NaN; % 0;
            blade_st.d2         = NaN; % 0;   
            blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm;
        else
            blade_st.Cla        = 2*pi*0.9488;  % clbar = 0.5987281447 => theta approx 5.5deg   
            blade_st.Cl0        = 0.6158;       % Cl at zero angle of attack
            blade_st.Cd0        = blade_st.Cla * deg2rad(5) / 29;;
            blade_st.d1         = 0;
            blade_st.d2         = 0;
            % blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm;
        end
        % 5) Pitch, chord and solidity distribution
        nsections           = length(y_arr) - 1;
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections;      
        blade_st.y_arr      = y_arr;
        blade_st.r_cutoff   = 0.15;
        blade_st.theta_arr  = theta_arr + collpitch;
        blade_st.chord_arr  = chord_arr .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);
    end    
    if strcmp(rotortype, 'KDECF305DP')  
        known_rotortype     = true;        
                       
        [y_arr, chord_arr, theta_arr] = kde_rotor_geometry(rotortype);
           
        % 1) General rotor properties      
        blade_st.Nb         = 2;
        blade_st.R          = y_arr(end);   
        blade_st.rotArea    = pi*blade_st.R^2;    
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;
        
        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan( (blade_st.lambda_c + 10^-6) / (blade_st.mu) ); 
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);
        
        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        blade_st.Cla        = NaN; % 2*pi*0.9488;  % clbar = 0.5987281447 => theta approx 5.5deg   
        blade_st.Cl0        = NaN; % 0.6158;       % Cl at zero angle of attack
        blade_st.Cd0        = NaN; % blade_st.Cla * deg2rad(5) / 29;;
        blade_st.d1         = NaN; % 0;
        blade_st.d2         = NaN; % 0;   
        blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm;
        
        % 5) Pitch, chord and solidity distribution
        nsections           = length(y_arr) - 1;
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections;      
        blade_st.y_arr      = y_arr;
        blade_st.theta_arr  = theta_arr + collpitch;
        blade_st.chord_arr  = chord_arr .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);
    end    
    if strcmp(rotortype, 'NACA5407')  
        known_rotortype     = true;	
      
        [y_arr, chord_arr, theta_arr] = kde_rotor_geometry(rotortype);
           
        % 1) General rotor properties    
        blade_st.Nb         = 2;
        blade_st.R          = y_arr(end);   
        blade_st.rotArea    = pi*blade_st.R^2;    
        blade_st.omega      = omega;
        blade_st.Vtip       = blade_st.omega*blade_st.R;
        blade_st.Vsound     = Vsound;
        blade_st.mach       = blade_st.Vtip / blade_st.Vsound;
        
        % 2) Airflow properties
        blade_st.rho        = rho;            % density in kg/m3        
        blade_st.normf      = blade_st.rho * blade_st.rotArea * blade_st.Vtip^2;    
        blade_st.lambda_c   = lambda_c;
        blade_st.mu         = mu;
        blade_st.CT_target  = CT_target;
        blade_st.TPP_alpha  = atan((blade_st.lambda_c + 10^-6)/(blade_st.mu));
        blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha);
        
        % 3) Blade motion  
        blade_st.b0         = 0;
        blade_st.b1c        = 0;
        blade_st.b1s        = 0; 
        
        % 4) Blade lift and drag
        blade_st.Cla        = NaN; % 2*pi*0.9488; % clbar = 0.5987281447    
        blade_st.Cl0        = NaN; % 0.6158;       % Cl at zero angle of attack
        blade_st.Cd0        = NaN; % blade_st.Cla * deg2rad(5) / 29;
        blade_st.d1         = NaN; % 0;
        blade_st.d2         = NaN; % 0;   
        blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm;
        
        % 5) Pitch, chord and solidity distribution
        nsections           = length(y_arr) - 1;
        blade_st.ones_arr   = ones(1, nsections + 1);
        blade_st.nsections  = nsections;      
        blade_st.y_arr      = y_arr;
        blade_st.theta_arr  = theta_arr + collpitch;
        blade_st.chord_arr  = chord_arr .* blade_st.ones_arr;
        blade_st.sigma_arr  = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R);
    end       
    if known_rotortype == false
        rotortype
        error('Unknown rotortype')
    end
    
    % Case name
    blade_st.casename   = ['mu' num2str(blade_st.mu) '_lc' num2str(blade_st.lambda_c)];
    blade_st.casename   = strrep(blade_st.casename,'.','p');
    
end
