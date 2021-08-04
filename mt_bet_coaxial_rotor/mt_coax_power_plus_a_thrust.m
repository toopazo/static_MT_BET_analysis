function mt_st = mt_coax_power_plus_a_thrust(...
        T_u         , ...
        T_l         , ...
        omega_u     , ...
        omega_l     , ...        
        radius      , ...
        density     , ...
        area        , ...
        solidity    , ...
        Cd0         , ...
        modeltype     ...
    )
    
    
    % Thrust ratio
    eta_T = T_u / T_l;
    
    % omega ratio
    eta_omega = omega_u / omega_l;
    
    % Vtip
    Vtip_u = omega_u * radius;
    Vtip_l = omega_l * radius;
    
    k   = 1.1;
    W   = T_u + T_l;

    v_l_div_v_u = NaN;
    kint = NaN;
    if strcmp(modeltype, 'leishman2006aerodynamic')
        % Induced velocity
        v_l_div_v_u = (1/2) * (...
            sqrt( (eta_T^3 + 4*eta_T^2 + 8*eta_T + 4)/(eta_T) ) - eta_T - 2 ...
            );     
        eta_v = 1 / v_l_div_v_u;
        
        % eta_T
        % eta_v
        % eta_T = 1 => eta_v = 1.780776406404415 => 1/eta_v = 0.5616

        % Induced power
        P0      = W^(3/2) / sqrt(2*density*area);
        Pcoax_i = P0 * sqrt( (1+eta_T) / (eta_T)) * ( eta_v / (1+eta_v) );    
        Psbs_i  = P0 * ( 1 + eta_T^(3/2) ) / ( (1+eta_T)^(3/2) );   
        kint    = Pcoax_i / Psbs_i;    
    end
    
    Pcoax_i   = Pcoax_i;
    Psbs_i    = Psbs_i;
    kint      = kint;
    
    % Upper and lower iduced power for coax
    v_u         = sqrt(T_u / (2 * density * area) );
    w_u         = 2 * v_u;
    Pcoaxu_i    = T_u * v_u;
    v_l         = v_l_div_v_u * v_u;
    Pcoaxl_i    = T_l * (v_u + v_l);    
    if (Pcoax_i - Pcoaxu_i - Pcoaxl_i) > 10^-4
        error('Error in upper and lower iduced coax power')
    end  
    
    % Upper and lower induced power for sbs
    Psbsu_i = (1 / sqrt(2 * density * area)) * T_u^(3/2);
    Psbsl_i = (1 / sqrt(2 * density * area)) * T_l^(3/2);
    if (Psbs_i - Psbsu_i - Psbsl_i) > 10^-4
        error('Error in upper and lower iduced sbs power')
    end 
    
    % Profile power
    % Po        = (solidity * Cd0 / 8) * density * area * (Vtip_u^3 + Vtip_l^3);
    Pu_o        = (solidity * Cd0 / 8) * density * area * Vtip_u^3;
    Pl_o        = (solidity * Cd0 / 8) * density * area * Vtip_l^3;
    Pcoaxu_o    = Pu_o;
    Pcoaxl_o    = Pl_o;
    Pcoax_o     = Pcoaxu_o + Pcoaxl_o; % Po;
    Psbsu_o     = Pu_o;
    Psbsl_o     = Pl_o;
    Psbs_o      = Psbsu_o + Psbsl_o; % Po;
    
    % Total power
    Pcoax = k * Pcoax_i + Pcoax_o;
    Psbs  = k * Psbs_i + Psbs_o;
    
    % coax        
    mt_st.Pcoaxu_i  = Pcoaxu_i;
    mt_st.Pcoaxl_i  = Pcoaxl_i;
    mt_st.Pcoax_i   = Pcoax_i;
    mt_st.Pcoaxu_o  = Pcoaxu_o;
    mt_st.Pcoaxl_o  = Pcoaxl_o;
    mt_st.Pcoax_o   = Pcoax_o;
    mt_st.Pcoax     = Pcoax;
    % sbs
    mt_st.Psbsu_i   = Psbsu_i;
    mt_st.Psbsl_i   = Psbsl_i;
    mt_st.Psbs_i    = Psbs_i;
    mt_st.Psbsu_o   = Psbsu_o;
    mt_st.Psbsl_o   = Psbsl_o;
    mt_st.Psbs_o    = Psbs_o;
    mt_st.Psbs      = Psbs;
    % upper and lower ratios
    mt_st.kint      = kint;
    mt_st.eta_T     = eta_T;
    mt_st.eta_v     = eta_v;
        
end

