function [bet_u_st, bet_l_st] = bet_coax_forces(blade_u_st, blade_l_st)
    
    ones_arr    = ones(length(blade_l_st.y_arr)-1, 1);

    % Non-iterative claculation of coax rotor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    blade_u_st.coax_pos     = 'upper';
    % Lower on upper rotor inflow
    % Lower on upper interference 
    blade_u_st.coax_lambda  = bet_l_st.coaxu_l0 .* ones_arr;
    bet_u_st                = bet_forces(blade_u_st);
    bet_u_st                = bet_forces_add_total(bet_u_st, false);
    lambda_u                = bet_u_st.total.lr;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    blade_l_st.coax_pos     = 'lower';
    % Upper on lower rotor inflow
    % Upper on lower interference
    blade_l_st.coax_lambda  = bet_l_st.coaxl_l0 .* ones_arr;
    bet_l_st                = bet_forces(blade_l_st);
    bet_l_st                = bet_forces_add_total(bet_l_st, false);
    lambda_l                = bet_l_st.total.lr;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    return

    % Iterative claculation of coax rotor
    lambda_l    = 0.0 * ones_arr;
    T_u_prev    = 0;
    T_l_prev    = 0;
        
    niter       = 100;
    for i = 1:niter
        
        % Simulate coax rotor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        blade_u_st.coax_pos     = 'upper';
        % Lower on upper rotor inflow
        % Lower on upper interference 
        blade_u_st.coax_lambda  = lambda_l;
        bet_u_st                = bet_forces(blade_u_st);
        bet_u_st                = bet_forces_add_total(bet_u_st, false);
        lambda_u                = bet_u_st.total.lr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        blade_l_st.coax_pos     = 'lower';
        % Upper on lower rotor inflow
        % Upper on lower interference
        blade_l_st.coax_lambda  = lambda_u;
        bet_l_st                = bet_forces(blade_l_st);
        bet_l_st                = bet_forces_add_total(bet_l_st, false);
        lambda_l                = bet_l_st.total.lr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
        
        % Check for terminating condition
        delta_T_u   = abs(bet_u_st.total.T - T_u_prev);
        delta_T_l   = abs(bet_l_st.total.T - T_l_prev);
        delta_T     = delta_T_u + delta_T_l;
        T_u_prev    = bet_u_st.total.T;
        T_l_prev    = bet_l_st.total.T;       

    %    fprintf('[bet_coax_forces] i %d, delta_T %.4f N, lambda_u %.4f, lambda_l %.4f \n', ...
    %        i, delta_T, bet_u_st.total.lambda, bet_l_st.total.lambda);
        if delta_T < 10^(-2)
            % arg = ['[bet_coax_forces] i ' num2str(i) ' delta ' num2str(delta)];
            % disp(arg);
            % fprintf('[bet_coax_forces] i %d, delta %.4f N, lambda_u %.4f, lambda_l %.4f \n', ...
            %     i, delta, bet_u_st.total.lambda, bet_l_st.total.lambda);
            break;
        end
        if i >= (niter-1)
            error('Too many iterations without convergence')
        end
    end            
    % if i < (2)
    %     error('Too few iterations for convergence')
    % end
end
