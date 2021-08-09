import numpy as np
import math
from unit_conversion import UnitConversion
from rotor_model import BladeSt
from matmath import MatMath as Mat
from mt_loads import MT


class Bet:
    @staticmethod
    def bet_forces_inflow_bemt(sigma, Cla, Cl0, theta, r, lambda_c):
        # function [lambda, lambda_i] =
        #     bet_forces_inflow_bemt(sigma, Cla, Cl0, theta, r, lambda_c)

        #    @book{leishman2006principles,
        #      title={Principles of helicopter aerodynamics with CD extra},
        #      author={Leishman, Gordon J},
        #      year={2006},
        #      publisher={Cambridge university press}
        #    }

        # lambda = lambda_c + lambda_i
        val = (sigma * Cla / 16 - lambda_c / 2) ** 2 + \
              (sigma * Cla * theta * r / 8 + sigma * Cl0 * r / 8)
        lambda_t = np.sqrt(val) - (sigma * Cla / 16 - lambda_c / 2)

        lambda_i = lambda_t - lambda_c

        if np.isnan(lambda_t):
            print(lambda_t)
            raise RuntimeError('lambda is NaN')
        if np.isreal(lambda_t) is False:
            print(lambda_t)
            raise RuntimeError('lambda is complex')
        return [lambda_t, lambda_i]

    @staticmethod
    def calculate_vel_ut_up_ur(
            r, lambda_t, mu, beta, betadot, psi, omega, radius):
        # function [UT UP UR] =
        #     bet_UT_UP_UR(r, lambda, mu, beta, betadot, psi, omega, R)
        # %    @book{leishman2006principles,
        # %      title={Principles of helicopter aerodynamics with CD extra},
        # %      author={Leishman, Gordon J},
        # %      year={2006},
        # %      publisher={Cambridge university press}
        # %    }
        #

        # % UT = vel along the TPP
        # % UP = vel perpendicular to the TPP
        # % UR = vel radial to blade section
        vel_ut = (omega * radius) * (r + mu * np.sin(psi))
        vel_up = (omega * radius) * (
                    lambda_t + r * betadot / omega + mu * beta * np.cos(psi))
        vel_ur = (omega * radius) * (mu * np.cos(psi))

        if omega == 0:
            vel_up = (omega * radius) * (lambda_t + 0 + mu * beta * np.cos(psi))

        if np.isnan(vel_ut):
            raise RuntimeError('UT is NaN')
        if np.isnan(vel_up):
            raise RuntimeError('UP is NaN')
        if np.isnan(vel_ur):
            raise RuntimeError('UR is NaN')
        if np.isfinite(vel_ut) is False:
            raise RuntimeError('UT is not finite')
        if np.isfinite(vel_up) is False:
            raise RuntimeError('UP is not finite')
        if np.isfinite(vel_ur) is False:
            raise RuntimeError('UR is not finite')
        if np.isreal(vel_ut) is False:
            raise RuntimeError('UT is not real')
        if np.isreal(vel_up) is False:
            raise RuntimeError('UP is not real')
        if np.isreal(vel_ur) is False:
            raise RuntimeError('UR is not real')

        return [vel_ut, vel_up, vel_ur]

    @staticmethod
    def bet_forces_inflow_bemt_c81table(
            r, mu, beta, betadot, psi, omega, R, vsound,
            theta, sigma, lambda_c, blade_st):
        # function [lambda_t, lambda_i] = bet_forces_inflow_bemt_c81table(...
        #     r, mu, beta, betadot, psi, omega, R, Vsound,
        #     theta, sigma, lambda_c, blade_st)

        assert isinstance(blade_st, BladeSt)

        # rotor_radius = blade_st.R
        # rotor_angvel = blade_st.omega

        if r <= blade_st.r_cutoff:
            lambda_i = 0
            lambda_t = lambda_c + lambda_i
            return [lambda_t, lambda_i]

        # options = optimoptions(
        #           @fsolve, 'Display', 'off', 'MaxFunctionEvaluations', 1500)
        # xguess = lambda_c + 0.1
        # [x0, fval, exitflag, output] = fsolve(@funzero, xguess, options)
        xguess = lambda_c + 0.1
        [x0, fval, exitflag, output] = Mat.fsolve(
            Bet.bemt_funzero, xguess,
            # lambda_t,
            r, mu, beta, betadot, psi, omega, R, vsound,
            theta, sigma, lambda_c, blade_st)

        lambda_t = x0
        lambda_i = lambda_t - lambda_c

        # err = funzero(lambda_t)
        # if norm(err) > 10**(-4)
        if exitflag < 1:
            # options = optimoptions(
            #           @fsolve,
            #           'Display', 'iter', 'MaxFunctionEvaluations', 1500)
            # xguess = 0.1
            # [x0, fval, exitflag, output] = fsolve(@funzero, xguess, options)
            options = None
            xguess = 0.1
            [x0, fval, exitflag, output] = Mat.fsolve(xguess, options)

            print(r)
            print(x0)
            print(fval)
            print(exitflag)
            print(output)
            print(lambda_t)
            raise RuntimeError('exitflag >= 1')
        return [lambda_t, lambda_i]

    @staticmethod
    def bemt_funzero(lambda_t,
                     r, mu, beta, betadot, psi, omega, radius, vsound,
                     theta, sigma, lambda_c, blade_st):
        assert isinstance(blade_st, BladeSt)
        # lambda_t = sqrt(x*x)

        [vel_ut, vel_up, vel_ur] = Bet.calculate_vel_ut_up_ur(
            r, lambda_t, mu, beta, betadot, psi, omega, radius)
        vel_u = math.sqrt(vel_ut**2 + vel_up**2 + vel_ur**2)
        mach = vel_u / vsound
        phi = math.atan2(vel_up, vel_ut)    # atan(UP / UT)
        alpha = theta - phi

        [Cl, Cd, Cm] = blade_st.get_ClCdCm(mach, alpha, blade_st.rotortype)

        # leishman2006principles says that Prandtl tip loss factor is a good
        # approximation for helicopter blades, but not so much
        # for propellers
        # So, I assume it means that is not good for "smaller" blades
        # Or most likely, not good for blades with high pitch angles
        if phi > 0:
            f = (blade_st.Nb/2) * (1-r)/(r*phi)
            # ptlf = Prandtl tip loss factor = F
            ptlf = (2/math.pi) * math.acos(math.exp(-f))
        else:
            ptlf = 1
            if phi != 0:
                phi_deg = phi * 180 / np.pi
                print('phi is negative')
        # F = 1

        # BEMT 0.5*Cl*r**2 dr = 4*F*lambda_t*(lambda_t - lambda_c) r dr
        err = (0.5 * sigma * Cl * r) - (
                    4 * ptlf * lambda_t * (lambda_t - lambda_c))

        if (np.isnan(err) is True) or (np.isfinite(err) is False) or (
                np.isreal(err) is False):
            print(omega)
            print(lambda_t)
            print([vel_ut, vel_up, vel_ur])
            print(phi)
            print(alpha)
            print([Cl, Cd, Cm])
            print(ptlf)
            print(err)
            raise RuntimeError('err is NaN, Inf or Complex')

        if phi < 0:
            print(omega)
            print(lambda_t)
            print([vel_ut, vel_up, vel_ur])
            print(phi)
            print(alpha)
            print([Cl, Cd, Cm])
            print(ptlf)
            print(err)
            print(phi * 180 / np.pi)
            raise RuntimeError('err is negative')

        return err

    @staticmethod
    def bet_forces_coax_pos(blade_st, r, r_arr):
        # function lambda_c = bet_forces_coax_pos(blade_st, r, r_arr)

        assert isinstance(blade_st, BladeSt)

        r_cutoff = blade_st.r_cutoff
        # lambda_c = blade_st.lambda_c    # Normalized axial flow velocity
        veltip = blade_st.omega * blade_st.R

        if blade_st.coax_pos == 'upper':   # lower on upper influence
            # Radius of lower rotor downwash at upper rotor location
            r0 = blade_st.coaxu_r0
            if r <= 1:
                # inflow 1
                # lambda_l = np.mean(blade_st.coax_lambda[r_arr > r_cutoff])
                downwash_l = np.mean(blade_st.coax_downwash[r_arr > r_cutoff])
                # inflow 2
                # lambda_l = interp1([0 r_arr(:)],  ...
                #     [0 blade_st.coax_lambda(:)], (r/r0)*r_arr(), 'linear')
                # inflow 3
                # lambda_l = interp1(r_arr, ...
                #     blade_st.coax_lambda, r, 'linear')

                # influence 1
                # lambda_infl = (1 / r0) * (1 / r0) * lambda_l
                downwash_infl = (1 / r0) * (1 / r0) * downwash_l
                lambda_infl = downwash_infl / veltip
                # influence 2
                # lambda_infl = blade_st.coaxu_l0
                if np.isnan(lambda_infl):
                    raise RuntimeError('lambda_u is NaN')

                # Carefull here, make sure all elements in this summation are
                # noramlized wrt the same tip veolocity (omega R)
                lambda_c = blade_st.lambda_c + lambda_infl
                return lambda_c
            else:
                lambda_c = blade_st.lambda_c + 0
                return lambda_c

        if blade_st.coax_pos == 'lower':   # upper on lower influence
            # Radius of upper rotor downwash at lower rotor location
            r0 = blade_st.coaxl_r0
            if r <= r0:
                # inflow 1
                # lambda_u = np.mean(blade_st.coax_lambda[r_arr > r_cutoff])
                downwash_u = np.mean(blade_st.coax_downwash[r_arr > r_cutoff])
                # inflow 2
                # lambda_u = interp1(r_arr,  ...
                #     blade_st.coax_lambda, (r/r0)*r_arr(), 'linear')
                # inflow 3
                # lambda_u = interp1(r_arr, ...
                #     blade_st.coax_lambda, r, 'linear')

                # influence 1
                # lambda_infl = (1 / r0) * (1 / r0) * lambda_u
                downwash_infl = (1 / r0) * (1 / r0) * downwash_u
                lambda_infl = downwash_infl / veltip
                # influence 2
                # lambda_infl = blade_st.coaxl_l0
                if np.isnan(lambda_infl):
                    raise RuntimeError('lambda_l is NaN')

                # Carefull here, make sure all elements in this summation are
                # noramlized wrt the same tip veolocity (omega R)
                lambda_c = blade_st.lambda_c + lambda_infl
                return lambda_c
            else:
                lambda_c = blade_st.lambda_c + 0
                return lambda_c

        raise RuntimeError('blade_st.coax_pos is neither upper nor lower')

    @staticmethod
    def bet_forces(blade_st):

        assert isinstance(blade_st, BladeSt)
    #    @book{leishman2006principles,
    #      title={Principles of helicopter aerodynamics with CD extra},
    #      author={Leishman, Gordon J},
    #      year={2006},
    #      publisher={Cambridge university press}
    #    }

        # blade_st parameters
        nblades = blade_st.Nb          # Number of blades
        rotor_radius = blade_st.R           # Blade length
        omega = blade_st.omega       # Blade angular speed
        y_arr = blade_st.y_arr       # Blade sectional location distribution
        theta_arr = blade_st.theta_arr   # Blade pitch distribution
        sigma_arr = blade_st.sigma_arr   # Blade solidity distribution
        chord_arr = blade_st.chord_arr   # Blade chord distribution
        rho = blade_st.rho         # Fluid density
        vsound = blade_st.Vsound      # Speed of sound
        Cla = blade_st.Cla         # Coeff of drag term associated with alpha**1
        Cl0 = blade_st.Cl0         # Coeff of drag term associated with alpha**0
        Cd0 = blade_st.Cd0         # Coeff of drag term associated with alpha**0
        d1 = blade_st.d1          # Coeff of drag term associated with alpha**1
        d2 = blade_st.d2          # Coeff of drag term associated with alpha**2
        b0 = blade_st.b0          # Blade conning angle
        b1c = blade_st.b1c         # Blade flapping associated with cos(psi)
        b1s = blade_st.b1s         # Blade flapping associated with sin(psi)
        lambda_c = blade_st.lambda_c    # Normalized axial flow velocity
        mu = blade_st.mu          # Normalized TPP flow velocity
        CT_target = blade_st.CT_target
        nsections = blade_st.nsections
        r_cutoff = blade_st.r_cutoff

        dr = 1.0 / nsections
        num_dr = nsections
        # blade_st.r_arr = 0     1     2     3
        # r_arr =   0.5   1.5   2.5
        r_arr = np.linspace(dr, 1, num_dr) - dr/2

        dpsi = UnitConversion.deg2rad(360) / nsections
        num_dpsi = nsections
        psi_arr = np.linspace(0, 2*np.pi-dpsi, num_dpsi)
        # psi_arr = [ 0, 22.5, 45.0, .. 315.0, 337.5 ]
        # psi_arr = deg2rad( linspace(0, 360-22.5, 16) )
        # num_dpsi = size(psi_arr, 2)

        L_arr = np.zeros((num_dpsi, num_dr))
        D_arr = np.zeros((num_dpsi, num_dr))
        T_arr = np.zeros((num_dpsi, num_dr))
        Q_arr = np.zeros((num_dpsi, num_dr))
        P_arr = np.zeros((num_dpsi, num_dr))

        lambda_c_arr = np.zeros((num_dpsi, num_dr))
        lambda_i_arr = np.zeros((num_dpsi, num_dr))
        lambda_arr = np.zeros((num_dpsi, num_dr))
        alpha_arr = np.zeros((num_dpsi, num_dr))

        chord_fnct = Mat.interp1_fnct(y_arr, chord_arr, 'linear')
        sigma_fnct = Mat.interp1_fnct(y_arr, sigma_arr, 'linear')
        theta_fnct = Mat.interp1_fnct(y_arr, theta_arr, 'linear')

        for i in range(0, num_dpsi):
            psi = psi_arr[i]
            # fprintf('psi %.2f/%.2f \n', rad2deg(psi), rad2deg(psi_arr()))

            if (mu == 0) and (i >= 1):
                # Copy BET along r around psi to avoid redundant calculations
                for j in range(0, num_dr):
                    L_arr[i, j] = L_arr[0, j]
                    D_arr[i, j] = D_arr[0, j]
                    T_arr[i, j] = T_arr[0, j]
                    Q_arr[i, j] = Q_arr[0, j]
                    P_arr[i, j] = P_arr[0, j]

                    lambda_c_arr[i, j] = lambda_c_arr[0, j]
                    lambda_i_arr[i, j] = lambda_i_arr[0, j]
                    lambda_arr[i, j] = lambda_arr[0, j]
                    alpha_arr[i, j] = alpha_arr[0, j]
                continue

            # print('Proceed to calculate BET along r')
            for j in range(0, num_dr):

                # Normalized sectional location
                r = r_arr[j]

                # Sectional location
                y = r * rotor_radius
                dy = dr * rotor_radius

                # 1) Calculate sectional blade geometry
                # chord = Mat.interp1(y_arr, chord_arr, y, 'linear')
                # sigma = Mat.interp1(y_arr, sigma_arr, y, 'linear')
                # theta = Mat.interp1(y_arr, theta_arr, y, 'linear')
                chord = chord_fnct(y)
                sigma = sigma_fnct(y)
                theta = theta_fnct(y)

                # 2) Calculate blade motion
                beta = b0 + b1c * math.cos(psi) + b1s * math.sin(psi)
                betadot = omega * (-b1c * math.sin(psi) + b1s * math.cos(psi))

                # Modification of lambda_c for coax rotors
                lambda_c = blade_st.lambda_c
                if (blade_st.coax_pos == 'upper') or (
                        blade_st.coax_pos == 'lower'):
                    lambda_c = Bet.bet_forces_coax_pos(blade_st, r, r_arr)

                # 3) Calculate inflow
                if mu <= 0.01:
                    # inflow according to BEMT = (BET + MT + mu=0 )
                    if blade_st.c81table:
                        # It is necessary to simultaneously solve for [lambda_t, UT, UP, UR, Cl]
                        #   [UT, UP, UR] = bet_UT_UP_UR(...
                        #       r, lambda_t, mu, beta, betadot, psi, omega, R)
                        #   U = sqrt(UT**2 + UP**2 + UR**2)
                        #   mach = U / Vsound
                        #   phi = atan2(UP, UT) # UP / UT
                        #   alpha = theta - phi
                        #   [Cl, Cd, Cm] = blade_st.c81table.get_ClCdCm(...
                        #       mach, alpha, blade_st.rotortype)
                        # Such that BEMT assumption holds True and
                        #   err = ( 0.5*Cl*r ) - ( 4*F*lambda_t*(lambda_t - lambda_c) )
                        #   is zero or very close to zero
                        [lambda_t, lambda_i] = \
                            Bet.bet_forces_inflow_bemt_c81table(
                                r, mu, beta, betadot, psi, omega, rotor_radius,
                                vsound, theta, sigma, lambda_c, blade_st)
                    else:
                        # It is assumed that Cl = Cl0 + Cla*alpha
                        [lambda_t, lambda_i] = Bet.bet_forces_inflow_bemt(
                            sigma, Cla, Cl0, theta, r, lambda_c)

                # inflow according to BET and Glauert
                if mu > 0.01:
                    kx = 1.2
                    # tan(TPP_alpha) = lambda_c / mu
                    TPP_alpha = math.atan(lambda_c / (mu+10**-6))
                    # lambda_i from MT
                    Vtip = np.nan  # Only used when mu == 0
                    lambda_MT = MT.mt_inflow(
                        CT_target, mu, TPP_alpha, Vtip, lambda_c)
                    lambda_0 = lambda_MT - lambda_c
                    # Glauert approximation
                    lambda_i = lambda_0 * ( 1 + kx*r*math.cos(psi) )
                    # total inflow
                    lambda_t = lambda_c + lambda_i

                # 4) Calculate velocities UT UP UR
                # UT = vel along the TPP
                # UP = vel perpendicular to the TPP
                # UR = vel radial to blade section
                [UT, UP, UR] = Bet.calculate_vel_ut_up_ur(r, lambda_t, mu, beta, betadot, psi, omega, rotor_radius)

                # 5) Sectional angle of attack = theta - phi
                # phi = atan(UP / UT)
                # At small angles => phi =  UP / UT
                #   UT = omega y = omega R r
                #   UP = Vc + vi
                #   phi = ( Vc + vi ) / ( omega R r ) = lambda_t / r
                phi = math.atan2(UP, UT) # lambda_t / r
                alpha = theta - phi

                # 6) Calculate sectional [Cl, Cd, Cm]
                if blade_st.c81table:
                    U = math.sqrt(UT**2 + UP**2 + UR**2)
                    mach = U / vsound

                    [Cl, Cd, Cm] = blade_st.get_ClCdCm(
                        mach, alpha, blade_st.rotortype)
                else:
                    Cl = Cl0 + Cla*alpha
                    Cd = Cd0 + d1*alpha + d2*alpha**2
                    Cm = 0

                # 7) Calculate sectional dT dQ dP
                # [dT, dQ, dP] = bet_dT_dQ_dP_using_Cl(...
                #     UT, UP, UR, chord, Cl, Cd, Nb, rho, omega, y, dy)
                # 7.1) Calculate sectional dL dD
                dL = 0.5 * rho * chord * Cl * (UT ** 2 + UP ** 2) * dy
                dD = 0.5 * rho * chord * Cd * (UT ** 2 + UP ** 2) * dy
                # 7.2) Calculate sectional dT dQ dP
                dT = nblades * dL
                dQ = nblades * (phi*dL + dD) * y
                dP = dQ * omega

                # 8) Save local section results
                # NOTE: T_arr(psi, r) is a partition of the continuous thrust
                # and represents the contribution of the ENTIRE section dr
                # whose center is localted at (psi, r)
                # A True distribution is therefore T_arr(psi, r)/dr
                # because it is "dr" -or "nsections"- independent
                # The same is True for Q_arr and P_arr, but NOT for
                # lambda_c_arr, lambda_i_arr or lambda_arr
                if (r <= r_cutoff) or (omega == 0):
                    dL = 0
                    dD = 0
                    dT = 0
                    dQ = 0
                    dP = 0

                    lambda_c = 0  # it is an extrinsic variable
                    lambda_i = 0
                    lambda_t = 0
                    alpha = 0

                L_arr[i, j] = dL
                D_arr[i, j] = dD
                T_arr[i, j] = dT
                Q_arr[i, j] = dQ
                P_arr[i, j] = dP

                lambda_c_arr[i, j] = lambda_c
                lambda_i_arr[i, j] = lambda_i
                lambda_arr[i, j] = lambda_t
                alpha_arr[i, j] = alpha

        # Put results into bet_st
        bet_st = BetResults()
        bet_st.blade_st = blade_st
        bet_st.psi_arr = psi_arr
        bet_st.r_arr = r_arr
        bet_st.dpsi = psi_arr
        bet_st.dr = dr

        bet_st.L_arr = L_arr
        bet_st.D_arr = D_arr
        bet_st.T_arr = T_arr
        bet_st.Q_arr = Q_arr
        bet_st.P_arr = P_arr

        bet_st.lambda_c_arr = lambda_c_arr
        bet_st.lambda_i_arr = lambda_i_arr
        bet_st.lambda_arr = lambda_arr
        bet_st.alpha_arr = alpha_arr

        # bet_st = bet_forces_add_total(bet_st, False)
        return bet_st

    @staticmethod
    def bet_forces_add_total(bet_st, verbose):
        # function bet_st = bet_forces_add_total(bet_st, verbose)
        assert isinstance(bet_st, BetResults)

        # Unpack bet_forces results
        blade_st = bet_st.blade_st
        assert isinstance(blade_st, BladeSt)

        rotor_type = blade_st.rotortype
        rotor_angvel = blade_st.omega
        rotor_radius = blade_st.R
        rotor_veltip = blade_st.omega * blade_st.R
        rotor_area = blade_st.rotArea
        fluid_density = blade_st.rho

        psi_arr = bet_st.psi_arr
        r_arr = bet_st.r_arr
        T_arr = bet_st.T_arr
        Q_arr = bet_st.Q_arr
        P_arr = bet_st.P_arr

        lambda_c_arr = bet_st.lambda_c_arr
        lambda_i_arr = bet_st.lambda_i_arr
        lambda_arr = bet_st.lambda_arr
        alpha_arr = bet_st.alpha_arr

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Net produced thrust, torque and power
        # Net produced loads along azimuth
        Tr_arr = Bet.bet_forces_mean_along_dpsi(T_arr)
        Qr_arr = Bet.bet_forces_mean_along_dpsi(Q_arr)
        Pr_arr = Bet.bet_forces_mean_along_dpsi(P_arr)

        # A blade produces lift at all r, but only during one psi at a time
        # Therefore, the following is wrong
        #   T = sum(T_arr, 'all')
        #   Q = sum(Q_arr, 'all')
        # And the this is right
        T = sum(Tr_arr)    # Sum of avgerage of T(psi, r) along psi
        Q = sum(Qr_arr)    # Sum of avgerage of Q(psi, r) along psi
        P = Q * rotor_angvel

        # percT = T / blade_st.Thover * 100
        # percQ = Q / ( blade_st.Pavail / blade_st.omega ) * 100
        # percP = P / blade_st.Pavail * 100

        CT = T / (fluid_density * rotor_area * rotor_veltip ** 2)
        CQ = Q / (fluid_density * rotor_area * rotor_veltip ** 2 * rotor_radius)
        CP = CQ    # CQ == CP

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Inflow
        lcr_arr = Bet.bet_forces_mean_along_dpsi(lambda_c_arr)
        lir_arr = Bet.bet_forces_mean_along_dpsi(lambda_i_arr)
        lr_arr = Bet.bet_forces_mean_along_dpsi(lambda_arr)
        aoa_arr = Bet.bet_forces_mean_along_dpsi(alpha_arr)

        # lambda_c = bet_st.blade_st.lambda_c
        # Avgerage of { avgerage of li(psi, r) along psi } along r
        lambda_c = np.mean(lcr_arr)
        lambda_i = np.mean(lir_arr)
        lambda_t = np.mean(lr_arr)
        # lambda_t = lambda_c + lambda_i

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Net dissymmetry of lift moment
        [dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr] = Bet.bet_forces_dolM(bet_st)

        mean_dolMx = np.mean(dolMxpsi_arr)
        mean_dolMy = np.mean(dolMypsi_arr)
        mean_dolMz = np.mean(dolMzpsi_arr)

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Add to total bet_st
        bet_st.total_T = T
        bet_st.total_Q = Q
        bet_st.total_P = P
        # print('[bet_forces_add_total] bet_st.total_P %s' % bet_st.total_P)

        # bet_st.total_percT = percT
        # bet_st.total_percQ = percQ
        # bet_st.total_percP = percP

        bet_st.total_CT = CT
        bet_st.total_CQ = CQ
        bet_st.total_CP = CP

        # Average along r
        bet_st.total_lambda_c = lambda_c
        bet_st.total_lambda_i = lambda_i
        bet_st.total_lambda_t = lambda_t

        # distribution along r
        bet_st.total_lcr = lcr_arr
        bet_st.total_lir = lir_arr
        bet_st.total_lr = lr_arr
        # bet_st.total_lcr = lambda_c .* lir_arr ./ lir_arr
        # bet_st.total_lcr(isnan(bet_st.total_lcr)) = 0
        # bet_st.total_lr = bet_st.total_lcr + bet_st.total_lir
        bet_st.total_aoar = aoa_arr

        # dissymmetry of lift
        bet_st.total_mean_dolMx = mean_dolMx
        bet_st.total_mean_dolMy = mean_dolMy
        bet_st.total_mean_dolMz = mean_dolMz

        ####################################

        # if blade_st.omega < 0
        #     blade_st.omega
        #     bet_st.total
        #     error('blade_st.omega < 0')
        # 
        # if bet_st.total_lambda_i < 0
        #     blade_st.omega
        #     bet_st.total
        #     error('bet_st.total_lambda_i < 0')
        # 
        # if bet_st.total_T < 0
        #     blade_st.omega
        #     bet_st.total
        #     error('bet_st.total_T < 0')
        # 

        # 4) Print
        if verbose:
            # Print blade type
            print('  %s ' % rotor_type)

            # Angular speed
            print('  Omega = %.2f rad/s, Vtip = %.2f m/s' %
                  (rotor_angvel, rotor_veltip))

            # Print net produced Thrust and Torque
            print('  T %.2f N, Q %.2f Nm, P = %.2f' % (T, Q, P))
            # print('  percThover = %.2f, percPhover = %.2f \n', percT, percP)
            print('  CT %.4f, CQ %.4f' % (CT, CQ))

            # Print net dissymmetry of lift moment
            print('  Mean dolM about psi [0, 2 pi] is [%.2f %.2f %.2f]' %
                  (mean_dolMx, mean_dolMy, mean_dolMz))

        return bet_st

    @staticmethod
    def bet_forces_mean_along_dpsi(T_arr):
        # function Tr = bet_forces_mean_along_dpsi(T_arr)
        # print('bet_forces_mean_along_dpsi')
        # print(T_arr.shape)
        # print(np.size(T_arr, 2))

        num_dr = np.size(T_arr, 1)
        # num_dpsi = np.size(T_arr, 1)

        #  For each T(psi, r) take average along psi
        Tr = np.zeros(num_dr)
        for j in range(0, num_dr):
            Tr[j] = np.mean(T_arr[:, j])
        return Tr

    @staticmethod
    def bet_forces_dolM(bet_st):
        # function [dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr] = bet_forces_dolM(bet_st)
        assert isinstance(bet_st, BetResults)

        # Unpack bet_st results
        blade_st = bet_st.blade_st
        psi_arr = bet_st.psi_arr
        r_arr = bet_st.r_arr
        T_arr = bet_st.T_arr
        Q_arr = bet_st.Q_arr
        lambda_i_arr = bet_st.lambda_i_arr

        # Get sectional position [x y z]
        rotor_xyz = BladeSt.blade_model_rotor_xyz(psi_arr, r_arr, blade_st)
        num_dpsi = len(psi_arr)
        num_dr = len(r_arr)

        # Calculate moment
        dolM_arr = np.zeros((3, num_dpsi, num_dr))
        for i in range(0, num_dpsi):
            x_arr = np.squeeze( rotor_xyz[0, i, :] )
            y_arr = np.squeeze( rotor_xyz[1, i, :] )
            z_arr = np.squeeze( rotor_xyz[2, i, :] )
            for j in range(0, num_dr):
                x = x_arr[j]
                y = y_arr[j]
                z = z_arr[j]
                Tij = T_arr[i, j]

                dolM = np.cross( [x, y, z], [0, 0, Tij] )
                dolM_arr[:, i, j] = dolM

        # Calculate net mean moment
        dolMx_arr = np.squeeze(dolM_arr[0, :, :])
        dolMy_arr = np.squeeze(dolM_arr[1, :, :])
        dolMz_arr = np.squeeze(dolM_arr[2, :, :])

        dolMxpsi_arr = Bet.bet_forces_sum_along_dr(dolMx_arr)
        dolMypsi_arr = Bet.bet_forces_sum_along_dr(dolMy_arr)
        dolMzpsi_arr = Bet.bet_forces_sum_along_dr(dolMz_arr)

    #    mean_dolMx = mean(dolMxpsi_arr)
    #    mean_dolMy = mean(dolMypsi_arr)
    #    mean_dolMz = mean(dolMzpsi_arr)
    #    dolM = [mean_dolMx mean_dolMy mean_dolMz]
        return [dolMxpsi_arr, dolMypsi_arr, dolMzpsi_arr]

    @staticmethod
    def bet_forces_sum_along_dr(T_arr):
        # function Tpsi = bet_forces_sum_along_dr(T_arr)

        num_dr = np.size(T_arr, 1)
        num_dpsi = np.size(T_arr, 0)

        # For each T(psi, r) add along r
        Tpsi = np.zeros(num_dpsi)
        for i in range(0, num_dpsi):
            Tpsi[i] = np.sum(T_arr[i, :])
        return Tpsi


class BetResults:
    __isfrozen = False

    def __setattr__(self, key, value):
        # https://stackoverflow.com/questions/3603502/
        # prevent-creating-new-attributes-outside-init
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError("%r is a frozen class. No new attr %s can be set" %
                            (self, key))
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True

    def __init__(self):
        self.blade_st = None
        self.psi_arr = None
        self.r_arr = None
        self.dpsi = None
        self.dr = None

        self.L_arr = None
        self.D_arr = None
        self.T_arr = None
        self.Q_arr = None
        self.P_arr = None

        self.lambda_c_arr = None
        self.lambda_i_arr = None
        self.lambda_arr = None
        self.alpha_arr = None

        # self.total = None

        # Add to total bet_st
        self.total_T = None
        self.total_Q = None
        self.total_P = None

        # self.total_percT = percT
        # self.total_percQ = percQ
        # self.total_percP = percP

        self.total_CT = None
        self.total_CQ = None
        self.total_CP = None

        # Average along r
        self.total_lambda_c = None
        self.total_lambda_i = None
        self.total_lambda_t = None

        # distribution along r
        self.total_lcr = None
        self.total_lir = None
        self.total_lr = None
        # self.total_lcr = lambda_c .* lir_arr ./ lir_arr
        # self.total_lcr(isnan(self.total_lcr)) = 0
        # self.total_lr = self.total_lcr + self.total_lir
        self.total_aoar = None

        # dissymmetry of lift
        self.total_mean_dolMx = None
        self.total_mean_dolMy = None
        self.total_mean_dolMz = None

        self._freeze()  # no new attributes after this point.


class BetFindOmega:
    @staticmethod
    def funzero(x, blade_st, thrust_target):
        # function err = funzero(x)
        # print('[funzero]')
        # print(x)
        blade_st.omega = float(x)

        bet_st = Bet.bet_forces(blade_st)
        bet_st = Bet.bet_forces_add_total(bet_st, False)
        thrust_bet = bet_st.total_T
        err = thrust_bet - thrust_target
        return err

    @staticmethod
    def bet_find_omega(blade_st, thrust_target):
        # function blade_st = bet_find_omega(blade_st, T_target)

        assert isinstance(blade_st, BladeSt)

        # options = optimset('Display','off');
        xguess = blade_st.omega
        [x0, fval, exitflag, output] = \
            Mat.fsolve(BetFindOmega.funzero, xguess, blade_st, thrust_target)
        # err = funzero(x0);

        # print('[bet_find_omega]')
        # print(x0)
        blade_st.omega = float(x0)
        bet_st = Bet.bet_forces(blade_st)
        bet_st = Bet.bet_forces_add_total(bet_st, True)
        _ = bet_st

        return blade_st

