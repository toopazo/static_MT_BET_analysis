function [dT dQ dP] = bet_dT_dQ_dP_using_Cl(UT, UP, UR, c, Cl, Cd, Nb, rho, omega, y, dy)
%    @book{leishman2006principles,
%      title={Principles of helicopter aerodynamics with CD extra},
%      author={Leishman, Gordon J},
%      year={2006},
%      publisher={Cambridge university press}
%    }

    % Simplifying assumptions
    % UR is always zero in hover
    % UR impact on lift is ignored in forward flight
    % UR impact on drag is ignored in forward flight
    % phi = UP / UT (small angle approx)
    % dD * phi = 0 (or negligible)    
    phi = UP / UT;
    
    % Derived equations
    % lambda = phi * r is valid for small angle phi
    % alpha = theta - phi
    % dr = dy / R
    
    % rho   = air density
    % c     = blade chord     
    % Cl    = coeff of lift
    % Cd    = coeff of drag
    % c*dy  = reference area    
    
    dL = 0.5 * rho * c * Cl * ( UT^2 + UP^2 ) * dy ;
    dD = 0.5 * rho * c * Cd * ( UT^2 + UP^2 ) * dy ;

    dT = Nb * dL;
    dQ = Nb * (phi*dL + dD) * y;
    dP = dQ * omega;    

end

