function [T, Qt, Pt, Pi, Pc, Po, Pp] = mt_T_Q_P(rho, A, sigma, omega, R, Vc, vi, Cd0)

   % (lambda, blade_st)

%    mu          = blade_st.mu;
%    lambda_c    = blade_st.lambda_c;
%    omega       = blade_st.omega;
%    R           = blade_st.R;
%    rho         = blade_st.rho;
%    sigma       = mean(blade_st.sigma_arr);
%    A           = blade_st.rotArea;
%    Cd0         = blade_st.Cd0;

%    lambda_i = lambda - lambda_c;    
%    
    Vtip    = (omega*R);
    Vtip2   = (omega*R)^2;
    Vtip3   = (omega*R)^3;
%    Vc      = lambda_c * Vtip;
%    vi      = lambda_i * Vtip;
%    
%    (rho, A, sigma, omega, R, Vc, vi, Cd0)
    
    % lambda_i = CT / 2*sqrt(lambda^2 + mu^2)
    % 2*sqrt(lambda^2 + mu^2)*lambda_i = T / rho*A*Vtip2
    % T = 2*rho*A*Vtip2*sqrt(lambda^2 + mu^2)*lambda_i
    % T = 2 * rho * A * Vtip2 * sqrt(lambda^2 + mu^2) * lambda_i;
    T = rho*A*(Vc+vi)*(2*vi);
    
    % P = T(Vc+vi) = T * lambda * omega * R = (T * lambda * R) * omega
    % Q = T * lambda * R;
    % P = Q * omega;

    Pi = T * vi;
    Pc = T * Vc;
    % Po = (1/8) * rho * Nb * c * Cd0 * omega^3 * R^4
    Po = (1/8) * sigma * Cd0 * (rho * A * Vtip3);
    Pp = 0;
    Pt = Pi + Pc + Po + Pp;
    Qt = Pt / omega;
end
