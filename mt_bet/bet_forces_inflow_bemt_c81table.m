function [lambda, lambda_i] = bet_forces_inflow_bemt_c81table(...
    r, mu, beta, betadot, psi, omega, R, Vsound, theta, sigma, lambda_c, blade_st)

    if r <= blade_st.r_cutoff
        lambda_i    = 0;
        lambda      = lambda_c + lambda_i;
        return
    end

    options = optimoptions(@fsolve, 'Display', 'off', 'MaxFunctionEvaluations', 1500);
    xguess = lambda_c + 0.1;
    [x0, fval, exitflag, output] = fsolve(@funzero, xguess, options);
    
    lambda      = x0;
    lambda_i    = lambda - lambda_c; 
        
    % err = funzero(lambda);
    % if norm(err) > 10^(-4)
    if exitflag < 1 % 
        options = optimoptions(@fsolve, 'Display', 'iter', 'MaxFunctionEvaluations', 1500);
        xguess = 0.1;
        [x0, fval, exitflag, output] = fsolve(@funzero, xguess, options);

        format long
        r
        x0
        fval
        exitflag
        output
        lambda
        error('exitflag >= 1')
    end   

    function err = funzero(x)
        lambda = x;
        % lambda = sqrt(x*x);
    
        [UT, UP, UR] = bet_UT_UP_UR(r, lambda, mu, beta, betadot, psi, omega, R);
        U       = sqrt(UT^2 + UP^2 + UR^2);            
        mach    = U / Vsound;
        phi     = atan2(UP, UT); % atan(UP / UT) ;                    
        alpha   = theta - phi;                 
        
        [Cl, Cd, Cm] = blade_st.c81table.get_ClCdCm(...
            mach, alpha, blade_st.rotortype);
                   
        % leishman2006principles says that Prandtl tip loss factor is a good 
        % approximation for helicopter blades, but not so much
        % for propellers
        % So, I assume it means that is not good for "smaller" blades
        % Or most likely, not good for blades with high pitch angles
        if phi > 0
            f = (blade_st.Nb/2) * (1-r)/(r*phi);  
            F = (2/pi) * acos(exp(-f));        
        else
            F = 1;
            if phi ~= 0
                phi_deg = phi * 180 / pi
                disp('phi is negative')
            end
        end
        % F = 1;

        % BEMT 0.5*Cl*r^2 dr = 4*F*lambda*(lambda - lambda_c) r dr 
        err = ( 0.5*sigma*Cl*r )  - ( 4*F*lambda*(lambda - lambda_c) );

        if (isnan(err) == true) || (isfinite(err) == false) || (isreal(err) == false)
            omega
            lambda 
            [UT, UP, UR] 
            phi 
            alpha
            [Cl, Cd, Cm]
            F
            err
            error('err is NaN, Inf or Complex')
        end

        if phi < 0
            omega
            lambda 
            [UT, UP, UR] 
            phi 
            alpha
            [Cl, Cd, Cm]
            F
            err

            phi_deg = phi * 180 / pi
            error('err is negative')
        end

    end
end
