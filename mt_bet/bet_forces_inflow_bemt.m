function [lambda, lambda_i] = bet_forces_inflow_bemt(sigma, Cla, Cl0, theta, r, lambda_c)
%    @book{leishman2006principles,
%      title={Principles of helicopter aerodynamics with CD extra},
%      author={Leishman, Gordon J},
%      year={2006},
%      publisher={Cambridge university press}
%    }

    % lambda = lambda_c + lambda_i
    val = ( sigma*Cla/16 - lambda_c/2 )^2 + ( sigma*Cla*theta*r/8 + sigma*Cl0*r/8);
    lambda = sqrt( val ) - ( sigma*Cla/16 - lambda_c/2 );

    lambda_i = lambda - lambda_c;     
    
    if isnan(lambda)
        lambda
        error('lambda is NaN')
    end
    if isreal(lambda) == false
        lambda
        error('lambda is complex')
    end
end
