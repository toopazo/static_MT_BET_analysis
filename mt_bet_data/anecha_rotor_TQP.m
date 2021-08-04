function [T, Q, P] = anecha_rotor_TQP(...
    omega       , ... % 
    rotortype   , ... %
    modeltype     ... %
    )
      
    if strcmp(modeltype, 'quadratic')
        % [thrust_vs_rpm] m1, motornum 111, nblades 2, afit 0.0003353, lsqe 18.552769
        % [thrust_vs_rpm] m1, motornum 111, nblades 3, afit 0.0004221, lsqe 19.409315
        % [thrust_vs_rpm] m1, motornum 111, nblades 6, afit 0.0005811, lsqe 37.446083
        a_thrust_nblades2 = 0.000335;
        a_thrust_nblades3 = 0.000422;
        a_thrust_nblades6 = 0.000581;

        % [thrust_vs_rpm] m1, motornum 111, nblades 2, solidity 0.08184088347968264
        % [thrust_vs_rpm] m1, motornum 111, nblades 3, solidity 0.12276132521952397
        % [thrust_vs_rpm] m1, motornum 111, nblades 6, solidity 0.24552265043904795
        solidity_nblades2 = 0.08184088347968264;
        solidity_nblades3 = 0.12276132521952397;
        solidity_nblades6 = 0.24552265043904795;

        % [thrust_vs_rpm] m1, motornum 111, nblades 2, density 1.1849
        % [thrust_vs_rpm] m1, motornum 111, nblades 3, density 1.1849
        % [thrust_vs_rpm] m1, motornum 111, nblades 6, density 1.1849
        % density_nblades2 = 1.1849;
        % density_nblades3 = 1.1849;
        % density_nblades6 = 1.1849;
        density_nblades2 = 1.1388;
        density_nblades3 = 1.1388;
        density_nblades6 = 1.1388;

        % [thrust_vs_rpm] m1, motornum 111, nblades 2, area 0.304151164328273
        % [thrust_vs_rpm] m1, motornum 111, nblades 3, area 0.304151164328273
        % [thrust_vs_rpm] m1, motornum 111, nblades 6, area 0.304151164328273
        area_nblades2 = 0.304151164328273;
        area_nblades3 = 0.304151164328273;
        area_nblades6 = 0.304151164328273;

        rho = density_nblades2;
        A = area_nblades2;
        R = 0.6223 / 2;     % 24.5in = 0.6223m
        k = 1.15;
        Cd0 = 0.013640;
        Vtip3 = (omega*R)^3;

        switch rotortype
            case 'KDECF245DP'
                sigma = solidity_nblades2;
                aT = a_thrust_nblades2;

                T = aT * omega^2;
                P = k*(T^(3/2))/sqrt(2*rho*A) + (sigma*Cd0/8)*(rho*A*Vtip3);
                Q = P / omega;
                return
            case 'KDECF245TP'
                sigma = solidity_nblades3;
                aT = a_thrust_nblades3;
                
                T = aT * omega^2;
                P = k*(T^(3/2))/sqrt(2*rho*A) + (sigma*Cd0/8)*(rho*A*Vtip3);
                Q = P / omega;
                return
            case 'KDECF245HP'
                sigma = solidity_nblades6;
                aT = a_thrust_nblades6;
                
                T = aT * omega^2;
                P = k*(T^(3/2))/sqrt(2*rho*A) + (sigma*Cd0/8)*(rho*A*Vtip3);
                Q = P / omega;
                return
            otherwise
                rotortype
                error('Unrecognized rotortype')
        end      
    end
end

