function lambda = mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c)
%    @book{leishman2006principles,
%      title={Principles of helicopter aerodynamics with CD extra},
%      author={Leishman, Gordon J},
%      year={2006},
%      publisher={Cambridge university press}
%    }

    % lambda = lambda_c + lambda_i
    % lambda_c    = Vc / (omega*blade_st.R);   % vertical vel  
    % lambda_h    = sqrt( CT/2 )  
    % mu          = Vf / (omega*blade_st.R);   % horizontal vel
    
    if isnan(CT)
        error('CT is NaN')
    end
    if isnan(TPP_alpha)
        error('TPP_alpha is NaN')
    end
    
    if mu == 0.0
        Vc          = Vtip * lambda_c;
        
        % CT = T / (rho*A*Vtip^2)
        % vh = sqrt(T/(2*rho*A))  
        % vh = sqrt(0.5 * CT * Vtip^2)                  
        vh          = Vtip * sqrt(CT/2);
        vi_vh       = - (0.5*Vc/vh) + sqrt( (0.5*Vc/vh)^2 + 1 );
        vi          = vi_vh * vh;
        
        lambda_i    = vi / Vtip;
        
        lambda      = lambda_c + lambda_i;
        return 
    else
        options = optimset('Display','off');
        xguess = 10;
        x0 = fsolve(@funzero, xguess, options);
        err = funzero(x0);
        
        lambda = x0;
        return
    end
    
    function err = funzero(x)
        err = x - mu*tan(TPP_alpha) - (CT)/(2*sqrt(mu^2 + x^2));
    end
end
