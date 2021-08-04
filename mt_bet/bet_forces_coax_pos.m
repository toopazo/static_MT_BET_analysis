function lambda_c = bet_forces_coax_pos(blade_st, r, r_arr)

    r_cutoff    = blade_st.r_cutoff;
    lambda_c    = blade_st.lambda_c;    % Normalized axial flow velocity
    mu          = blade_st.mu;          % Normalized TPP flow velocity

    % r0 is dependant on h/R, but the exact form seems to 
    % depend on rotor design and operating point
    %
    % NASA Langley WT test was at 
    %   h/R = 0.65 
    %   eta_T = 1.35 (climb), 1.41 (hover) and 1.21 (ff)
    % 
    % coaxu_r0 = 1.70; % KDECF305DP, NASA Langeley WT and Jason's CFD
    % coaxl_r0 = 0.70; % KDECF305DP, NASA Langeley WT and Jason's CFD
    if (blade_st.coax_pos == 'upper')   % lower on upper influence
        % Radius of lower rotor downwash at upper rotor location
        r0 = blade_st.coaxu_r0;
        if (r <= 1)
            % inflow 1
            lambda_l = mean(blade_st.coax_lambda(r_arr>r_cutoff));
            % inflow 2
            % lambda_l = interp1([0; r_arr(:)],  ...
            %     [0; blade_st.coax_lambda(:)], (r/r0)*r_arr(end), 'linear');
            % inflow 3
            % lambda_l = interp1(r_arr, ...
            %     blade_st.coax_lambda, r, 'linear');
            
            lambda_infl = (1/r0)*(1/r0)*lambda_l;
            if isnan(lambda_infl)
                error('lambda_u is NaN')
            end
            lambda_c = blade_st.lambda_c + lambda_infl;
        end
    end
    if (blade_st.coax_pos == 'lower')   % upper on lower influence
        % Radius of upper rotor downwash at lower rotor location
        r0 = blade_st.coaxl_r0;
        if (r <= r0)
            % inflow 1
            lambda_u = mean(blade_st.coax_lambda(r_arr>r_cutoff));
            % inflow 2
            % lambda_u = interp1(r_arr,  ...
            %     blade_st.coax_lambda, (r/r0)*r_arr(end), 'linear');
            % inflow 3
            % lambda_u = interp1(r_arr, ...
            %     blade_st.coax_lambda, r, 'linear');
            
            lambda_infl = (1/r0)*(1/r0)*lambda_u;
            if isnan(lambda_infl)
                error('lambda_l is NaN')
            end
            lambda_c = blade_st.lambda_c + lambda_infl;
        end
    end
end
