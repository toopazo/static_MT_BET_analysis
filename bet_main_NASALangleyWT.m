function bet_main_NASALangleyWT()

    clear all
    close all
    format compact
    % format short
    format long
    % format bank
    clc
    
    % Add functions to the path
    addpath('mt_bet');  
    addpath('mt_bet_data');       
    addpath('mt_bet_coaxial_rotor');   
    addpath('mt_bet_plot'); 
    
    rpm2rads = pi / 30;
    rads2rpm = 30 / pi;   
    
%    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');    
%    for val = [0.00, 1.00, 1.6, 1.8, 2.00, 2.20, 2.40]
%        collpitch = - deg2rad(val);
%        fprintf('collpitch %.4f deg \n', rad2deg(collpitch));
%        coax_wt_cfd_bet_climb(collpitch, true);
%        coax_wt_cfd_bet_hover(collpitch, true);
%        coax_wt_cfd_bet_ff(collpitch, true);
%    end

    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    % Compare BET results with: Dragonfly_PhaseB_CFD_Update_14x22_Comparisons.pptx
    coax_wt_cfd_bet_climb(true);
    coax_wt_cfd_bet_hover(true);
    coax_wt_cfd_bet_ff(true);
    
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    kde_bet_at_3890rpm();       

    function coax_wt_cfd_bet_climb(collpitch, verbose)       
        rotortype   = 'KDECF305DP';
        [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
        
        % blade_u_st
        omega       = 3889 * rpm2rads; 
        lambda_c    = 2.29 / omega * 0.38735;       
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_u_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % blade_u_st                    
        
        % blade_l_st
        omega       = 4092 * rpm2rads;
        lambda_c    = 2.29 / omega * 0.38735;        
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_l_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % blade_l_st        
        
        % Simulate coax rotor
        [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st);
        T_u = bet_u_st.total.T;
        Q_u = bet_u_st.total.Q;
        T_l = bet_l_st.total.T;
        Q_l = bet_l_st.total.Q;        

        % Add error with WT data
        Terr = abs(T_u - 126.1) + abs(T_l - 89.5);
        Qerr = abs(Q_u - 4.7) + abs(Q_l - 4.4);
            
        if verbose                     
            fprintf('T_u %.4f, T_l %.4f, Terr %.4f \n', ...
                bet_u_st.total.T, bet_l_st.total.T, Terr);
            fprintf('Q_u %.4f, Q_l %.4f, Qerr %.4f \n', ...
                bet_u_st.total.Q, bet_l_st.total.Q, Qerr);
            fprintf('lambda_i_u %.4f, lambda_i_l %.4f \n', ...
                bet_u_st.total.lambda_i, bet_l_st.total.lambda_i);
                
            % Add WT, CFD and BET results
            fig = figure(1);
            subplot(2, 1, 1)
            text(0.05, 270, ' WT: 3889 RPM, 126.1 N');
            text(0.05, 230, 'CFD: 3889 RPM, 125.3 N');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Ttot = round(bet_u_st.total.T, 1);
            text(0.05, 190, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            subplot(2, 1, 2)
            text(0.05, 180, ' WT: 4092 RPM, 89.5 N');
            text(0.05, 140, 'CFD: 4092 RPM, 100.1 N');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Ttot = round(bet_l_st.total.T, 1);
            text(0.05, 100, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            
            fig = figure(2);
            subplot(2, 1, 1)
            text(0.05, 12, ' WT: 3889 RPM, 4.7 Nm');
            text(0.05, 10, 'CFD: 3889 RPM, 4.7 Nm');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Qtot = round(bet_u_st.total.Q, 1);
            text(0.05, 8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            subplot(2, 1, 2)
            text(0.05, 9, ' WT: 4092 RPM, 4.4 Nm');
            text(0.05, 7, 'CFD: 4092 RPM, 4.6 Nm');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Qtot = round(bet_l_st.total.Q, 1);
            text(0.05, 5, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            
            fig = figure(3);
            subplot(2, 1, 1)
            omegarpm = blade_u_st.omega*rads2rpm;
            lambda_i = round(bet_u_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);  
            subplot(2, 1, 2)
            omegarpm = blade_l_st.omega*rads2rpm;
            lambda_i = round(bet_l_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);
            
            legend_cell = {rotortype};
            bet_plot_coax(bet_u_st, bet_l_st, ...
                'bet_main_NASALangleyWT', 'coax_climb', legend_cell);
            
            close all;  
        end
    end
    
    function coax_wt_cfd_bet_hover(collpitch, verbose)
        rotortype   = 'KDECF305DP'; 
        [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
        
        % blade_u_st
        omega       = 3890 * rpm2rads;        
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_u_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % blade_u_st                    
        
        % blade_l_st
        omega       = 4038 * rpm2rads;        
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_l_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % blade_l_st
        
        % Simulate coax rotor
        [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st);
        T_u = bet_u_st.total.T;
        Q_u = bet_u_st.total.Q;
        T_l = bet_l_st.total.T;
        Q_l = bet_l_st.total.Q; 
        
        Terr = abs(T_u - 128.6) + abs(T_l - 94.9);
        Qerr = abs(Q_u - 4.5) + abs(Q_l - 4.4);
            
        if verbose                     
            fprintf('T_u %.4f, T_l %.4f, Terr %.4f \n', ...
                bet_u_st.total.T, bet_l_st.total.T, Terr);
            fprintf('Q_u %.4f, Q_l %.4f, Qerr %.4f \n', ...
                bet_u_st.total.Q, bet_l_st.total.Q, Qerr);
            fprintf('lambda_i_u %.4f, lambda_i_l %.4f \n', ...
                bet_u_st.total.lambda_i, bet_l_st.total.lambda_i);
        
            % Add WT, CFD and BET results to plots
            fig = figure(1);
            subplot(2, 1, 1)
            text(0.05, 270, ' WT: 3890 RPM, 128.6 N');
            text(0.05, 230, 'CFD: 3890 RPM, 140.6 N');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Ttot = round(bet_u_st.total.T, 1);
            text(0.05, 190, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            subplot(2, 1, 2)
            text(0.05, 180, ' WT: 4038 RPM, 94.9 N');
            text(0.05, 140, 'CFD: 4038 RPM, 110.0 N');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Ttot = round(bet_l_st.total.T, 1);
            text(0.05, 100, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            
            fig = figure(2);
            subplot(2, 1, 1)
            text(0.05, 12,  ' WT: 3890 RPM, 4.5 Nm');
            text(0.05, 10, 'CFD: 3890 RPM, 4.9 Nm');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Qtot = round(bet_u_st.total.Q, 1);
            text(0.05, 8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            subplot(2, 1, 2)
            text(0.05, 9, ' WT: 4038 RPM, 4.4 Nm');
            text(0.05, 7, 'CFD: 4038 RPM, 4.8 Nm');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Qtot = round(bet_l_st.total.Q, 1);
            text(0.05, 5, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            
            fig = figure(3);
            subplot(2, 1, 1)
            omegarpm = blade_u_st.omega*rads2rpm;
            lambda_i = round(bet_u_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);  
            subplot(2, 1, 2)
            omegarpm = blade_l_st.omega*rads2rpm;
            lambda_i = round(bet_l_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);
            
            legend_cell = {rotortype};
            bet_plot_coax(bet_u_st, bet_l_st, ...
                'bet_main_NASALangleyWT', 'coax_hover', legend_cell);
            
            close all;  
        end
    end

    function coax_wt_cfd_bet_ff(collpitch, verbose)                       
        rotortype   = 'KDECF305DP';
        [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
        
        lambda_c    = 0.0;
        mu          = 0.038;   
        
        % blade_u_st
        omega       = 3926 * rpm2rads;        
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_u_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % Add CT_target for mu != 0 inflow calculations
        Tu          = 120.5;
        rhoAVtip2   = blade_u_st.rho * blade_u_st.rotArea * blade_u_st.Vtip^2;
        blade_u_st.CT_target = Tu / rhoAVtip2;
        % blade_u_st        
                
        
        % blade_l_st
        omega       = 4385 * rpm2rads;        
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_l_st  = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        % Add CT_target for mu != 0 inflow calculations
        Tl          = 120; %99.2;
        rhoAVtip2   = blade_l_st.rho * blade_l_st.rotArea * blade_l_st.Vtip^2;
        blade_l_st.CT_target = Tl / rhoAVtip2;
        % blade_l_st
        
        % Calculate coax rotor
        [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st);
        T_u = bet_u_st.total.T;
        Q_u = bet_u_st.total.Q;
        T_l = bet_l_st.total.T;
        Q_l = bet_l_st.total.Q; 
        
        Terr = abs(T_u - 120.5) + abs(T_l - 99.2);
        Qerr = abs(Q_u - 4.3) + abs(Q_l - 4.4);
            
        if verbose                     
            fprintf('T_u %.4f, T_l %.4f, Terr %.4f \n', ...
                bet_u_st.total.T, bet_l_st.total.T, Terr);
            fprintf('Q_u %.4f, Q_l %.4f, Qerr %.4f \n', ...
                bet_u_st.total.Q, bet_l_st.total.Q, Qerr);
            fprintf('lambda_i_u %.4f, lambda_i_l %.4f \n', ...
                bet_u_st.total.lambda_i, bet_l_st.total.lambda_i);
        
            % Add WT, CFD and BET results
            fig = figure(1);
            subplot(2, 1, 1)
            text(0.05, 270, ' WT: 3926 RPM, 120.5 N');
            text(0.05, 230, 'CFD: 3926 RPM, 135.7 N');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Ttot = round(bet_u_st.total.T, 1);
            text(0.05, 190, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            subplot(2, 1, 2)
            text(0.05, 180, ' WT: 4385 RPM, 99.2 N');
            text(0.05, 140, 'CFD: 4385 RPM, 131.7 N');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Ttot = round(bet_l_st.total.T, 1);
            text(0.05, 100, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Ttot) ' N']);
            
            fig = figure(2);
            subplot(2, 1, 1)
            text(0.05, 12, ' WT: 3926 RPM, 4.3 Nm');
            text(0.05, 10, 'CFD: 3926 RPM, 4.6 Nm');
            omegarpm = blade_u_st.omega .* rads2rpm;
            Qtot = round(bet_u_st.total.Q, 1);
            text(0.05, 8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            subplot(2, 1, 2)
            text(0.05, 9, ' WT: 4385 RPM, 4.6 Nm');
            text(0.05, 7, 'CFD: 4385 RPM, 5.3 Nm');
            omegarpm = blade_l_st.omega .* rads2rpm;
            Qtot = round(bet_l_st.total.Q, 1);
            text(0.05, 5, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(Qtot) ' Nm']);
            
            fig = figure(3);
            subplot(2, 1, 1)
            omegarpm = blade_u_st.omega*rads2rpm;
            lambda_i = round(bet_u_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);  
            subplot(2, 1, 2)
            omegarpm = blade_l_st.omega*rads2rpm;
            lambda_i = round(bet_l_st.total.lambda_i, 4);
            text(0.3, lambda_i*0.8, ...
                ['BET: ' num2str(omegarpm) ' RPM, ' num2str(lambda_i)]);
            
            legend_cell = {rotortype};
            bet_plot_coax(bet_u_st, bet_l_st, ...
                'bet_main_NASALangleyWT', 'coax_ff', legend_cell);
            
            close all;
        end  
    end
        
    function kde_bet_at_3890rpm(collpitch)
        omega       = 3890 * rpm2rads;
        
        rotortype   = 'KDE8218XF120_KDECF305DP'
        % modeltype   = 'interp1'
        % [T, Q, P]   = kde_rotor_TQP(omega, rotortype, modeltype)
        modeltype   = 'quadratic'
        [T, Q, P]   = kde_rotor_TQP(omega, rotortype, modeltype)   
        
        rotortype   = 'KDECF305DP';
        [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype);
        omega       = omega;
        CT_target   = NaN;              % CT_target is only needed for mu != 0
        blade_st    = blade_model(...
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
        
        % Simulate single rotor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        bet_st = bet_forces(blade_st);
        bet_st = bet_forces_add_total(bet_st, true);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        str1        = 'bet_main_NASALangleyWT';
        str2        = rotortype;
        legend_cell = {rotortype};
        bet_plot_single(bet_st, str1, str2, legend_cell);
        
        close all;
    end
end

