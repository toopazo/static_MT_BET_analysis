function [rho, lambda_c, mu, collpitch] = kde_rotor_defaults(rotortype)

%    rho     = 1.225;            % density at NASA Langley WT
%    rho     = 1.2115900277;     % density at KDE
%    rho     = 1.1849;           % density at AneCha
%    rho     = 1.1388;           % density at AneCha @(400 msl, 25 deg)

% collpitch   = -deg2rad(0.0);
% collpitch   = -deg2rad(1.2);
% collpitch   = -deg2rad(2.5);

    if strcmp(rotortype, 'KDECF305DP')
        rotortype   = 'KDECF305DP';
        rho         = 1.225;              % density at NASA Langley WT
        lambda_c    = 0.0;
        mu          = 0.0;           
        collpitch   = -deg2rad(2.5);
        % omega       = omega;
        % CT_target   = NaN;              % CT_target is only needed for mu != 0
        
        % blade_st    = blade_model(...
        %     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);
    end
    if strcmp(rotortype, 'KDECF245DP')
        % rotortype   = 'KDECF305DP';
        rotortype   = 'KDECF245DP';
        % rho         = 1.2115900277;     % density at KDE
        % rho         = 1.1849;           % density at AneCha
        rho         = 1.1388;           % density at AneCha @(400 msl, 25 deg)
        lambda_c    = 0.0;
        mu          = 0.0;   
        % collpitch   = -deg2rad(1.2);
        collpitch   = -deg2rad(0.0);
        % omega       = omega;
        % CT_target   = NaN;              % CT_target is only needed for mu != 0
        
        % blade_st    = blade_model(...
        %     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);      
    end
    if strcmp(rotortype, 'KDECF245TP')
        % rotortype   = 'KDECF305DP';
        rotortype   = 'KDECF245DP';
        % rho         = 1.2115900277;     % density at KDE
        % rho         = 1.1849;           % density at AneCha
        rho         = 1.1388;           % density at AneCha @(400 msl, 25 deg)
        lambda_c    = 0.0;
        mu          = 0.0;   
        % collpitch   = -deg2rad(1.2);
        collpitch   = -deg2rad(0.0);
        % omega       = omega;
        % CT_target   = NaN;              % CT_target is only needed for mu != 0
        
        % blade_st    = blade_model(...
        %     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);      
    end
    if strcmp(rotortype, 'KDECF245HP')
        % rotortype   = 'KDECF305DP';
        rotortype   = 'KDECF245DP';
        % rho         = 1.2115900277;     % density at KDE
        % rho         = 1.1849;           % density at AneCha
        rho         = 1.1388;           % density at AneCha @(400 msl, 25 deg)
        lambda_c    = 0.0;
        mu          = 0.0;   
        % collpitch   = -deg2rad(1.2);
        collpitch   = -deg2rad(0.0);
        % omega       = omega;
        % CT_target   = NaN;              % CT_target is only needed for mu != 0
        
        % blade_st    = blade_model(...
        %     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target);      
    end        

end

