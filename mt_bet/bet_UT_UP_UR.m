function [UT UP UR] = bet_UT_UP_UR(r, lambda, mu, beta, betadot, psi, omega, R)
%    @book{leishman2006principles,
%      title={Principles of helicopter aerodynamics with CD extra},
%      author={Leishman, Gordon J},
%      year={2006},
%      publisher={Cambridge university press}
%    }

    % UT        = vel along the TPP
    % UP        = vel perpendicular to the TPP
    % UR        = vel radial to blade section
    UT = (omega*R) * ( r + mu*sin(psi) ) ;
    UP = (omega*R) * ( lambda + r*betadot/omega + mu*beta*cos(psi) ) ;
    UR = (omega*R) * ( mu*cos(psi) ) ;

    if omega == 0
        UP = (omega*R) * ( lambda + 0 + mu*beta*cos(psi) ) ;
    end

    if isnan(UT)
        error('UT is NaN')
    end
    if isnan(UP)
        error('UP is NaN')
    end
    if isnan(UR)
        error('UR is NaN')
    end
    
    if isfinite(UT) == false
        error('UT is not finite')
    end
    if isfinite(UP) == false
        error('UP is not finite')
    end
    if isfinite(UR) == false
        error('UR is not finite')
    end

    if isreal(UT) == false
        error('UT is not real')
    end
    if isreal(UP) == false
        error('UP is not real')
    end
    if isreal(UR) == false
        error('UR is not real')
    end

end

