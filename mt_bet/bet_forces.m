function bet_st = bet_forces(blade_st)
    
%    @book{leishman2006principles,
%      title={Principles of helicopter aerodynamics with CD extra},
%      author={Leishman, Gordon J},
%      year={2006},
%      publisher={Cambridge university press}
%    }

    % blade_st parameters
    Nb          = blade_st.Nb;          % Number of blades
    R           = blade_st.R;           % Blade length
    omega       = blade_st.omega;       % Blade angular speed        
    y_arr       = blade_st.y_arr;       % Blade sectional location distribution
    theta_arr   = blade_st.theta_arr;   % Blade pitch distribution
    sigma_arr   = blade_st.sigma_arr;   % Blade solidity distribution
    chord_arr   = blade_st.chord_arr;   % Blade chord distribution
    rho         = blade_st.rho;         % Fluid density
    Vsound      = blade_st.Vsound;      % Speed of sound
    Cla         = blade_st.Cla;         % Coeff of drag term associated with alpha^1
    Cl0         = blade_st.Cl0;         % Coeff of drag term associated with alpha^0
    Cd0         = blade_st.Cd0;         % Coeff of drag term associated with alpha^0
    d1          = blade_st.d1;          % Coeff of drag term associated with alpha^1
    d2          = blade_st.d2;          % Coeff of drag term associated with alpha^2   
    b0          = blade_st.b0;          % Blade conning angle
    b1c         = blade_st.b1c;         % Blade flapping associated with cos(psi)
    b1s         = blade_st.b1s;         % Blade flapping associated with sin(psi)
    lambda_c    = blade_st.lambda_c;    % Normalized axial flow velocity
    mu          = blade_st.mu;          % Normalized TPP flow velocity
    CT_target   = blade_st.CT_target;
    nsections   = blade_st.nsections;
    r_cutoff    = blade_st.r_cutoff;

    dr      = 1.0 / nsections;
    num_dr  = nsections;
    % blade_st.r_arr    = 0     1     2     3
    % r_arr             =   0.5   1.5   2.5 
    r_arr = linspace(dr, 1, num_dr) - dr/2;
    
    dpsi        =  deg2rad(360) / ( nsections );
    num_dpsi    = nsections;
    psi_arr     = linspace(0, 2*pi-dpsi, num_dpsi);
    % psi_arr = [ 0, 22.5, 45.0, .. 315.0, 337.5 ]
    % psi_arr = deg2rad( linspace(0, 360-22.5, 16) );
    % num_dpsi = size(psi_arr, 2);
    
    L_arr = zeros(num_dpsi, num_dr);
    D_arr = zeros(num_dpsi, num_dr);
    T_arr = zeros(num_dpsi, num_dr);
    Q_arr = zeros(num_dpsi, num_dr);
    P_arr = zeros(num_dpsi, num_dr);
    
    lambda_c_arr = zeros(num_dpsi, num_dr);
    lambda_i_arr = zeros(num_dpsi, num_dr);
    lambda_arr = zeros(num_dpsi, num_dr);
        
    for i = 1:num_dpsi
        psi = psi_arr(i);
        % fprintf('psi %.2f/%.2f \n', rad2deg(psi), rad2deg(psi_arr(end)));
        
        if (mu == 0) && (i >= 2)
            % disp('Copy BET along r around psi to avoid redundant calculations')')
            for j = 1:num_dr
                L_arr(i, j) = L_arr(1, j);
                D_arr(i, j) = D_arr(1, j);
                T_arr(i, j) = T_arr(1, j);
                Q_arr(i, j) = Q_arr(1, j);
                P_arr(i, j) = P_arr(1, j);
                
                lambda_c_arr(i, j) = lambda_c_arr(1, j);
                lambda_i_arr(i, j) = lambda_i_arr(1, j);
                lambda_arr(i, j) = lambda_arr(1, j);
                alpha_arr(i, j) = alpha_arr(1, j); 
            end
            continue
        end
        
        % disp('Proceed to calculate BET along r')
        for j = 1:num_dr
        
            % Normalized sectional location
            r = r_arr(j);
            
            % Sectional location
            y   = r * R;            
            dy  = dr * R;            
            
            % 1) Calculate sectional blade geometry
            chord   = interp1(y_arr, chord_arr, y, 'linear');            
            sigma   = interp1(y_arr, sigma_arr, y, 'linear');
            theta   = interp1(y_arr, theta_arr, y, 'linear');
            
            % 2) Calculate blade motion
            beta    =  b0 + b1c*cos(psi) + b1s*sin(psi);
            betadot = omega * (-b1c*sin(psi) + b1s*cos(psi) );     
            
            % Modification of lambda_c for coax rotors 
            lambda_c = blade_st.lambda_c;
            if isfield(blade_st,'coax_pos')
                lambda_c = bet_forces_coax_pos(blade_st, r, r_arr);
            end

            % 3) Calculate inflow            
            if (mu <= 0.01)     % inflow according to BEMT = (BET + MT + mu=0 )
                if isfield(blade_st,'c81table')
                    % It is necessary to simultaneously solve for [lambda, UT, UP, UR, Cl]
                    %   [UT, UP, UR] = bet_UT_UP_UR(...
                    %       r, lambda, mu, beta, betadot, psi, omega, R);
                    %   U       = sqrt(UT^2 + UP^2 + UR^2);            
                    %   mach    = U / Vsound;
                    %   phi     = atan2(UP, UT); % UP / UT ;                    
                    %   alpha   = theta - phi;                                     
                    %   [Cl, Cd, Cm] = blade_st.c81table.get_ClCdCm(...
                    %       mach, alpha, blade_st.rotortype);     
                    % Such that BEMT assumption holds true and
                    %   err = ( 0.5*Cl*r ) - ( 4*F*lambda*(lambda - lambda_c) ) 
                    %   is zero or very close to zero                    
                    [lambda, lambda_i] = bet_forces_inflow_bemt_c81table(...
                        r, mu, beta, betadot, psi, omega, R, ...
                        Vsound, theta, sigma, lambda_c, blade_st);    
                else
                    % It is assumed that Cl = Cl0 + Cla*alpha 
                    [lambda, lambda_i] = bet_forces_inflow_bemt(...
                        sigma, Cla, Cl0, theta, r, lambda_c); 
                end
            end
            if (mu > 0.01)      % inflow according to BET and Glauert
                kx = 1.2;
                % tan(TPP_alpha) = lambda_c / mu
                TPP_alpha = atan(lambda_c / (mu+10^-6));
                % lambda_i from MT
                Vtip        = NaN;  % Only used when mu == 0
                lambda_MT   = mt_inflow(...
                    CT_target, mu, TPP_alpha, Vtip, lambda_c);
                lambda_0    = lambda_MT - lambda_c;
                % Glauert approximation
                lambda_i    = lambda_0 * ( 1 + kx*r*cos(psi) );
                % total inflow
                lambda      = lambda_c + lambda_i;
            end
            
            % 4) Calculate velocities UT UP UR
            % UT        = vel along the TPP
            % UP        = vel perpendicular to the TPP
            % UR        = vel radial to blade section                   
            [UT, UP, UR] = bet_UT_UP_UR(...
                r, lambda, mu, beta, betadot, psi, omega, R);
                
            % 5) Sectional angle of attack = theta - phi
            % phi = atan(UP / UT) 
            % At small angles => phi =  UP / UT
            %   UT = omega y = omega R r 
            %   UP = Vc + vi
            %   phi = ( Vc + vi ) / ( omega R r ) = lambda / r
            phi     = atan2(UP, UT); % lambda / r;  
            alpha   = theta - phi;                          
                    
            % 6) Calculate sectional [Cl, Cd, Cm]
            if isfield(blade_st,'c81table')
                U       = sqrt(UT^2 + UP^2 + UR^2);            
                mach    = U / Vsound;
                            
                [Cl, Cd, Cm] = blade_st.c81table.get_ClCdCm(...
                    mach, alpha, blade_st.rotortype);
            else
                Cl = Cl0 + Cla*alpha;
                Cd = Cd0 + d1*alpha + d2*alpha^2;
                Cm = 0;
            end                     
                        
            % 7) Calculate sectional dT dQ dP
            % [dT, dQ, dP] = bet_dT_dQ_dP_using_Cl(...
            %     UT, UP, UR, chord, Cl, Cd, Nb, rho, omega, y, dy);
            % 7.1) Calculate sectional dL dD
            dL = 0.5 * rho * chord * Cl * ( UT^2 + UP^2 ) * dy ;
            dD = 0.5 * rho * chord * Cd * ( UT^2 + UP^2 ) * dy ;
            % 7.2) Calculate sectional dT dQ dP
            dT = Nb * dL;
            dQ = Nb * (phi*dL + dD) * y;
            dP = dQ * omega;
            
            % 8) Save local section results
            % NOTE: T_arr(psi, r) is a partition of the continuous thrust
            % and represents the contribution of the ENTIRE section dr 
            % whose center is localted at (psi, r)
            % A true distribution is therefore T_arr(psi, r)/dr
            % because it is "dr" -or "nsections"- independent
            % The same is true for Q_arr and P_arr, but NOT for 
            % lambda_c_arr, lambda_i_arr or lambda_arr
            if (r <= r_cutoff) || (omega == 0)
                dL          = 0;
                dD          = 0;
                dT          = 0;
                dQ          = 0;
                dP          = 0;  
                
                lambda_c    = 0;  % it is an extrinsic variable
                lambda_i    = 0;
                lambda      = 0;
                alpha       = 0;
            end
            L_arr(i, j) = dL;
            D_arr(i, j) = dD;
            T_arr(i, j) = dT;
            Q_arr(i, j) = dQ;
            P_arr(i, j) = dP; 
            
            lambda_c_arr(i, j) = lambda_c; 
            lambda_i_arr(i, j) = lambda_i; 
            lambda_arr(i, j) = lambda;                 
            alpha_arr(i, j) = alpha;     
        end        
    end       
    
    % Put results into bet_st
    bet_st.blade_st     = blade_st;
    bet_st.psi_arr      = psi_arr;
    bet_st.r_arr        = r_arr;
    bet_st.dpsi         = psi_arr;    
    bet_st.dr           = dr;
    
    bet_st.L_arr        = L_arr;
    bet_st.D_arr        = D_arr;
    bet_st.T_arr        = T_arr;
    bet_st.Q_arr        = Q_arr;
    bet_st.P_arr        = P_arr;
    
    bet_st.lambda_c_arr = lambda_c_arr;
    bet_st.lambda_i_arr = lambda_i_arr;
    bet_st.lambda_arr   = lambda_arr;
    bet_st.alpha_arr    = alpha_arr;
    
    % bet_st = bet_forces_add_total(bet_st, false);    
end
