function [omega_u, omega_l] = mt_find_omega_a_thrust(T_u, T_l, motornum, nblades, thrust_lvl)

    if motornum == 111
        % Get omega_u and omega_l from 
        [a_thrust_m1, a_thrust_m2] = ...
            mt_get_a_thrust(motornum, nblades, thrust_lvl);
        % T_u = a_thrust_m1 * omega_u^2;
        % T_l = a_thrust_m2 * omega_l^2;
        
        omega_u = sqrt( T_u / a_thrust_m1 );
        omega_l = sqrt( T_l / a_thrust_m2 );
        
        return
    end
    if motornum == 12
        % Get omega_u and omega_l from 
        [a_thrust_m1, a_thrust_m2] = ...
            mt_get_a_thrust(motornum, nblades, thrust_lvl);
        % T_u = a_thrust_m1 * omega_u^2;
        % T_l = a_thrust_m2 * (omega_l + 1750*pi/30)^2;
        
        % omega_u = sqrt( T_u / a_thrust_m1 );
        % omega_l = sqrt( T_l / a_thrust_m2 ) - 1750*pi/30 ;
        
        omega_u = sqrt( T_u / a_thrust_m1 );
        omega_l = sqrt( T_l / a_thrust_m1 );

        return
    end
end
