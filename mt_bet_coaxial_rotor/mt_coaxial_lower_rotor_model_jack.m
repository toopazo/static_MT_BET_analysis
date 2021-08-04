 function [thrust, a_thrust] = mt_coaxial_lower_rotor_model_jack(...
        omega_l, csvp, m1_thrust, omega_u, m2_thrust)
        
        % Check lower_rotor_aT.m file
        % in folder experiment_anecha/toopazo_anecha for more info

        % rho = db.blade_h_st.rho; 
        % rotArea = db.blade_h_st.rotArea;
        % R = db.blade_h_st.R;
        db = csvp;
        
        radius = db.blade_h_st.R; % csvp.rotor_data['radius']
        density =  db.blade_h_st.rho; % csvp.rotor_data['density']
        area = db.blade_h_st.rotArea; % csvp.rotor_data['area']
        % solidity = csvp.rotor_data['solidity']
        nblades = db.blade_h_st.Nb; % csvp.rotor_data['nblades']
        % chord = csvp.rotor_data['chord']

        theta = 12 * pi / 180

        a_thrust_m1 = 0.000336; % RotorDataTools.get_single_a_thrust(nblades)
        omega_at_m2_thrust = sqrt(m2_thrust / a_thrust_m1)

        % alpha (at the tip) of a isolated rotor at a given thrust
        vc_u = 0
        vh_u = sqrt(m2_thrust / (2 * density * area))
        vi_u = vh_u * (
                - (0.5 * vc_u / vh_u) + sqrt((0.5 * vc_u / vh_u)^2 + 1)
        )
        up_u = vc_u + vi_u
        phi_u = atan2(up_u, (omega_at_m2_thrust * radius))
        alpha_u = theta - phi_u

        % phih = atan2(dwl, rpmlhzz R)
        % phil = atan2(dwu + dwl, rpml R)

        % climb velocity experienced by a rotor downstream of m1
        up_u = vc_u + sqrt(m1_thrust / (2 * density * area))

        % Check if omega makes both aoa equal
        vc_l = up_u
        vh_l = sqrt(m2_thrust / (2 * density * area))
        vi_l = vh_l * (
                - (0.5 * vc_l / vh_l) + sqrt((0.5 * vc_l / vh_l)^2 + 1)
        )
        up_l = vc_l + vi_l
        phi_l = atan2(up_l, (omega_l * radius))
        alpha_l = theta - phi_l

        alpha_err = alpha_l - alpha_u

        thrust = m2_thrust + (alpha_err * 180 / pi) * 1
        a_thrust = thrust / (omega_l * omega_l)
        % return [thrust, a_thrust]
