function rotor_xyz = blade_model_rotor_xyz(psi_arr, y_arr, blade_st)

    num_dr      = length(y_arr);
    num_dpsi    = length(psi_arr);   
    
    rotor_xyz   = zeros(3, num_dpsi, num_dr); 
    for i = 1:num_dpsi
        psi     = psi_arr(i);
        beta    =  blade_st.b0 + blade_st.b1c*cos(psi) + blade_st.b1s*sin(psi);        
        for j = 1:num_dr
            y   = y_arr(j);
            x   = y .* cos(beta) .* cos(psi);
            y   = y .* cos(beta) .* sin(psi);
            z   = y .* sin(beta);
            
            rotor_xyz(:, i, j) = [x;y;z];
        end
    end
    % rotor_xyz(:, 1, :)
end
