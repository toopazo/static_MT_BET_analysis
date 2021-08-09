# from abc import ABC, abstractmethod
import numpy as np
import math
from unit_conversion import UnitConversion as UnitCnv
from gnrl_config import GnrlConfig
import scipy.io
from matmath import MatMath as Mat

# class BladeSt(ABC):
#     @abstractmethod
#     def noofsides(self):
#         pass


class BladeSt:
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
        # Static blade properties (no time dynamics)
        self.Nb = None          # Number of blades
        self.R = None           # Blade length
        self.rotArea = None     # Rotor disk area
        self.rotortype = None

        self.y_arr = None       # Blade sectional location distribution
        self.theta_arr = None   # Blade pitch distribution
        self.sigma_arr = None   # Blade solidity distribution
        self.chord_arr = None   # Blade chord distribution

        self.nsections = None
        self.r_cutoff = None

        self.Cla = None         # Coeff of lift rate
        self.Cl0 = None         # Coeff of lift at zero aoa
        self.Cd0 = None         # Equivalent coeff of drag, assoc. with alpha**0
        self.d1 = None          # Coeff of drag term associated with alpha**1
        self.d2 = None          # Coeff of drag term associated with alpha**2
        self.c81table = None
        self.get_ClCdCm = None

        self.TPP_alpha = None
        self.TPP_alphad = None
        self.b0 = None          # Blade conning angle
        self.b1c = None         # Blade flapping associated with cos(psi)
        self.b1s = None         # Blade flapping associated with sin(psi)

        self.omega = None       # Blade angular speed
        # self.Vtip = None        # tip velocity

        # Static coaxial properties (no time dynamics)
        self.coax_pos = None
        self.coaxu_r0 = None
        self.coaxl_r0 = None
        # self.coax_lambda = None
        self.coax_downwash = None

        # Static fluid flow properties (no time dynamics)
        self.rho = None         # Fluid density
        self.Vsound = None      # Speed of sound
        self.lambda_c = None    # Normalized axial flow velocity
        self.mu = None          # Normalized TPP flow velocity
        # self.normf = None
        # self.mach = None

        # Other properties
        self.casename = None
        self.CT_target = None

        self._freeze()  # no new attributes after this point.

    @staticmethod
    def blade_model_rotor_xyz(psi_arr, y_arr, blade_st):
        # function rotor_xyz = blade_model_rotor_xyz(psi_arr, y_arr, blade_st)

        num_dr = len(y_arr)
        num_dpsi = len(psi_arr)

        rotor_xyz = np.zeros((3, num_dpsi, num_dr))
        for i in range(0, num_dpsi):
            psi = psi_arr[i]
            beta = blade_st.b0 + blade_st.b1c * np.cos(
                psi) + blade_st.b1s * np.sin(psi)
            for j in range(0, num_dr):
                y = y_arr[j]
                x = y * np.cos(beta) * np.cos(psi)
                y = y * np.cos(beta) * np.sin(psi)
                z = y * np.sin(beta)

                rotor_xyz[:, i, j] = [x, y, z]
        # rotor_xyz(:, 1, :)
        return rotor_xyz

    @staticmethod
    def blade_model(rotortype, rho, lambda_c, mu, collpitch, omega, CT_target):
        # function blade_st = blade_model(
        #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)

        Vsound = 340.2941

        # In the case of hover => lambda_c = mu = 0
        #   lambda_c = mu*tan(TPP_alpha)
        #   TPP_alpha = atan(lambda_c/mu)
        #   atan(0/0) = NaN
        #   atan(1/0) = pi/2 = 90 deg
        #   tan(atan(1/0)) = 1.633123935319537e+16
        #   0*tan(atan(1/0)) = 0
        #   blade_st.TPP_alpha = atan((blade_st.lambda_c + 10**-6)/(blade_st.mu))

        blade_st = BladeSt()
        blade_st.rotortype = rotortype
        known_rotortype = False

        # if strfind(rotortype,'constant_twist_blade') == 1
        #     known_rotortype = True
        #
        #     # UH-60A
        #     [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry()
        #
        #     # 1) General rotor properties
        #     blade_st.Nb = Nb
        #     blade_st.R = R
        #     blade_st.rotArea = pi*blade_st.R**2
        #     blade_st.omega = omega
        #     blade_st.Vtip = blade_st.omega*blade_st.R
        #     blade_st.Vsound = Vsound
        #     blade_st.mach = blade_st.Vtip / blade_st.Vsound
        #
        #     # 2) Airflow properties
        #     blade_st.rho = rho            # density in kg/m3
        #     blade_st.normf = blade_st.rho * blade_st.rotArea * blade_st.Vtip**2
        #     blade_st.lambda_c = lambda_c
        #     blade_st.mu = mu
        #     blade_st.CT_target = CT_target
        #     blade_st.TPP_alpha = atan( (blade_st.lambda_c + 10**-6) / (blade_st.mu) )
        #     blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha)
        #
        #     # 3) Blade motion
        #     blade_st.b0 = 0
        #     blade_st.b1c = 0
        #     blade_st.b1s = 0
        #
        #     # 4) Blade lift and drag
        #     blade_st.Cla = Cla
        #     blade_st.Cl0 = 0.0
        #     blade_st.Cd0 = Cd0
        #     blade_st.d1 = 0
        #     blade_st.d2 = 0
        #     # blade_st.c81table.get_ClCdCm = NaN
        #
        #     # 5) Pitch, chord and solidity distribution
        #     # nsections = 16 =>
        #     #                   len(psi_arr) = 17
        #     #                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000
        #     nsections = 16
        #     blade_st.ones_arr = ones(1, nsections + 1)
        #     blade_st.nsections = nsections
        #     blade_st.y_arr = linspace(0, R, nsections + 1)
        #     theta0str = extractBetween(rotortype,'constant_twist_blade_','rad')
        #     if isempty(theta0str)
        #         theta0 = deg2rad(22.0335)    # deg2rad(22.0)
        #     else
        #         theta0 = str2double(theta0str)
        #
        #     blade_st.theta_arr = theta0 .* blade_st.ones_arr + collpitch
        #     blade_st.chord_arr = chord .* blade_st.ones_arr
        #     blade_st.sigma_arr = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R)
        #
        #
        # if rotortype, 'ideally_twisted_blade')
        #     known_rotortype = True
        #
        #     # UH-60A
        #     [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry()
        #
        #     # 1) General rotor properties
        #     blade_st.Nb = Nb
        #     blade_st.R = R
        #     blade_st.rotArea = pi*blade_st.R**2
        #     blade_st.omega = omega
        #     blade_st.Vtip = blade_st.omega*blade_st.R
        #     blade_st.Vsound = Vsound
        #     blade_st.mach = blade_st.Vtip / blade_st.Vsound
        #
        #     # 2) Airflow properties
        #     blade_st.rho = rho            # density in kg/m3
        #     blade_st.normf = blade_st.rho * blade_st.rotArea * blade_st.Vtip**2
        #     blade_st.lambda_c = lambda_c
        #     blade_st.mu = mu
        #     blade_st.CT_target = CT_target
        #     blade_st.TPP_alpha = atan( (blade_st.lambda_c + 10**-6) / (blade_st.mu) )
        #     blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha)
        #
        #     # 3) Blade motion
        #     blade_st.b0 = 0
        #     blade_st.b1c = 0
        #     blade_st.b1s = 0
        #
        #     # 4) Blade lift and drag
        #     blade_st.Cla = Cla
        #     blade_st.Cl0 = 0.0
        #     blade_st.Cd0 = Cd0
        #     blade_st.d1 = 0
        #     blade_st.d2 = 0
        #     # blade_st.c81table.get_ClCdCm = NaN
        #
        #     # 5) Pitch, chord and solidity distribution
        #     # nsections = 16 =>
        #     #                   len(psi_arr) = 17
        #     #                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000
        #     nsections = 16
        #     blade_st.ones_arr = ones(1, nsections + 1)
        #     blade_st.nsections = nsections
        #     blade_st.y_arr = linspace(0, R, nsections + 1)
        #     itb_theta0 = deg2rad(15.020)     # deg2rad(15.2)
        #     blade_st.theta_arr = itb_theta0 .* blade_st.ones_arr ./ linspace(0, R, nsections+1)
        #     blade_st.theta_arr = blade_st.theta_arr + collpitch
        #     theta_max = deg2rad(70)
        #     ind_arr = find(blade_st.theta_arr >= theta_max)
        #     for ind = ind_arr
        #         blade_st.theta_arr(ind) = theta_max
        #
        #     blade_st.chord_arr = chord .* blade_st.ones_arr
        #     blade_st.sigma_arr = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R)
        #
        # if rotortype, 'linearly_twisted_blade')
        #     known_rotortype = True
        #
        #     # UH-60A
        #     [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry()
        #
        #     ltb_theta0 = deg2rad(37.560)     # deg2rad(35.0)
        #     ltb_theta1 = deg2rad(17.000)     # deg2rad(18.0)
        #
        #     # 1) General rotor properties
        #     blade_st.Nb = Nb
        #     blade_st.R = R
        #     blade_st.rotArea = pi*blade_st.R**2
        #     blade_st.omega = omega
        #     blade_st.Vtip = blade_st.omega*blade_st.R
        #     blade_st.Vsound = Vsound
        #     blade_st.mach = blade_st.Vtip / blade_st.Vsound
        #
        #     # 2) Airflow properties
        #     blade_st.rho = rho            # density in kg/m3
        #     blade_st.normf = blade_st.rho * blade_st.rotArea * blade_st.Vtip**2
        #     blade_st.lambda_c = lambda_c
        #     blade_st.mu = mu
        #     blade_st.CT_target = CT_target
        #     blade_st.TPP_alpha = atan( (blade_st.lambda_c + 10**-6) / (blade_st.mu) )
        #     blade_st.TPP_alphad = rad2deg(blade_st.TPP_alpha)
        #
        #     # 3) Blade motion
        #     blade_st.b0 = 0
        #     blade_st.b1c = 0
        #     blade_st.b1s = 0
        #
        #     # 4) Blade lift and drag
        #     blade_st.Cla = Cla
        #     blade_st.Cl0 = 0.0
        #     blade_st.Cd0 = Cd0
        #     blade_st.d1 = 0
        #     blade_st.d2 = 0
        #     # blade_st.c81table.get_ClCdCm = NaN
        #
        #     # 5) Pitch, chord and solidity distribution
        #     # nsections = 16 =>
        #     #                   len(psi_arr) = 17
        #     #                   psi_arr = 0   22.5000   45.0000  ..  337.5000  360.0000
        #     nsections = 16
        #     blade_st.ones_arr = ones(1, nsections + 1)
        #     blade_st.nsections = nsections
        #     blade_st.y_arr = linspace(0, R, nsections + 1)
        #     blade_st.theta_arr = linspace(ltb_theta0, ltb_theta1, nsections + 1) + collpitch
        #     blade_st.chord_arr = chord .* blade_st.ones_arr
        #     blade_st.sigma_arr = (blade_st.Nb * blade_st.chord_arr) / (pi*blade_st.R)

        if (rotortype == 'KDECF245DP') or (rotortype == 'KDECF245TP') or (rotortype == 'KDECF245HP'):
            known_rotortype = True

            [y_arr, chord_arr, theta_arr] = \
                KdeGeometry.kde_rotor_geometry(rotortype)

            nblades = None
            if rotortype == 'KDECF245DP':
                nblades = 2
            if rotortype == 'KDECF245TP':
                nblades = 3
            if rotortype == 'KDECF245HP':
                nblades = 6


            # 1) General rotor properties
            blade_st.Nb = nblades
            blade_st.R = y_arr[-1]
            blade_st.rotArea = np.pi*blade_st.R**2
            blade_st.omega = omega
            # blade_st.Vtip = blade_st.omega*blade_st.R
            blade_st.Vsound = Vsound
            # blade_st.mach = blade_st.Vtip / blade_st.Vsound

            # 2) Airflow properties
            blade_st.rho = rho            # density in kg/m3
            # blade_st.normf = \
            #     blade_st.rho * blade_st.rotArea * blade_st.Vtip ** 2
            blade_st.lambda_c = lambda_c
            blade_st.mu = mu
            blade_st.CT_target = CT_target
            if (blade_st.mu == 0) and (blade_st.lambda_c > 0):
                blade_st.TPP_alpha = + np.pi / 2
            elif (blade_st.mu == 0) and (blade_st.lambda_c < 0):
                blade_st.TPP_alpha = - np.pi / 2
            elif (blade_st.mu == 0) and (blade_st.lambda_c == 0):
                blade_st.TPP_alpha = + np.pi / 2
            else:
                blade_st.TPP_alpha = math.atan(blade_st.lambda_c / blade_st.mu)
            blade_st.TPP_alphad = UnitCnv.rad2deg(blade_st.TPP_alpha)

            # 3) Blade motion
            blade_st.b0 = 0
            blade_st.b1c = 0
            blade_st.b1s = 0

            # 4) Blade lift and drag
            use_c81table = True
            if use_c81table:
                blade_st.Cla = np.nan # 2*pi*0.9488  # clbar = 0.5987281447 => theta approx 5.5deg
                blade_st.Cl0 = np.nan # 0.6158       # Cl at zero angle of attack
                blade_st.Cd0 = np.nan # blade_st.Cla * deg2rad(5) / 29
                blade_st.d1 = np.nan # 0
                blade_st.d2 = np.nan # 0
                # blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm
                blade_st.c81table = True
                blade_st.get_ClCdCm = KdeAero.kde_rotor_ClCdCm
            else:
                blade_st.Cla = 2*np.pi*0.9488  # clbar = 0.5987281447 => theta approx 5.5deg
                blade_st.Cl0 = 0.6158       # Cl at zero angle of attack
                blade_st.Cd0 = blade_st.Cla * UnitCnv.deg2rad(5 / 29)
                blade_st.d1 = 0
                blade_st.d2 = 0
                # blade_st.c81table.get_ClCdCm = @kde_rotor_ClCdCm

            # 5) Pitch, chord and solidity distribution
            nsections = len(y_arr) - 1
            ones_arr = np.ones(nsections + 1)
            blade_st.nsections = nsections
            blade_st.y_arr = y_arr
            blade_st.r_cutoff = 0.15
            blade_st.theta_arr = theta_arr + collpitch
            blade_st.chord_arr = chord_arr * ones_arr
            blade_st.sigma_arr = (blade_st.Nb * blade_st.chord_arr) / (
                        np.pi * blade_st.R)

        if rotortype == 'KDECF305DP':
            known_rotortype = True

            [y_arr, chord_arr, theta_arr] = KdeGeometry.kde_rotor_geometry(rotortype)

            # 1) General rotor properties
            blade_st.Nb = 2
            blade_st.R = y_arr[-1]
            blade_st.rotArea = np.pi*blade_st.R**2
            blade_st.omega = omega
            # blade_st.Vtip = blade_st.omega*blade_st.R
            blade_st.Vsound = Vsound
            # blade_st.mach = blade_st.Vtip / blade_st.Vsound

            # 2) Airflow properties
            blade_st.rho = rho            # density in kg/m3
            # blade_st.normf = \
            #     blade_st.rho * blade_st.rotArea * blade_st.Vtip**2
            blade_st.lambda_c = lambda_c
            blade_st.mu = mu
            blade_st.CT_target = CT_target
            blade_st.TPP_alpha = math.atan( (blade_st.lambda_c + 10**-6) / (blade_st.mu) )
            blade_st.TPP_alphad = UnitCnv.rad2deg(blade_st.TPP_alpha)

            # 3) Blade motion
            blade_st.b0 = 0
            blade_st.b1c = 0
            blade_st.b1s = 0

            # 4) Blade lift and drag
            blade_st.Cla = np.nan # 2*pi*0.9488  # clbar = 0.5987281447 => theta approx 5.5deg
            blade_st.Cl0 = np.nan # 0.6158       # Cl at zero angle of attack
            blade_st.Cd0 = np.nan # blade_st.Cla * deg2rad(5) / 29
            blade_st.d1 = np.nan # 0
            blade_st.d2 = np.nan # 0
            blade_st.c81table = True
            blade_st.get_ClCdCm = KdeAero.kde_rotor_ClCdCm

            # 5) Pitch, chord and solidity distribution
            nsections = len(y_arr) - 1
            ones_arr = np.ones(nsections + 1)
            blade_st.nsections = nsections
            blade_st.y_arr = y_arr
            blade_st.theta_arr = theta_arr + collpitch
            blade_st.chord_arr = chord_arr * ones_arr
            blade_st.sigma_arr = \
                (blade_st.Nb * blade_st.chord_arr) / (np.pi*blade_st.R)

        if rotortype == 'NACA5407':
            known_rotortype = True

            [y_arr, chord_arr, theta_arr] = \
                KdeGeometry.kde_rotor_geometry(rotortype)

            # 1) General rotor properties
            blade_st.Nb = 2
            blade_st.R = y_arr[-1]
            blade_st.rotArea = np.pi * blade_st.R**2
            blade_st.omega = omega
            # blade_st.Vtip = blade_st.omega*blade_st.R
            blade_st.Vsound = Vsound
            # blade_st.mach = blade_st.Vtip / blade_st.Vsound

            # 2) Airflow properties
            blade_st.rho = rho            # density in kg/m3
            # blade_st.normf = \
            #     blade_st.rho * blade_st.rotArea * blade_st.Vtip**2
            blade_st.lambda_c = lambda_c
            blade_st.mu = mu
            blade_st.CT_target = CT_target
            blade_st.TPP_alpha = \
                math.atan((blade_st.lambda_c + 10**-6)/(blade_st.mu))
            blade_st.TPP_alphad = UnitCnv.rad2deg(blade_st.TPP_alpha)

            # 3) Blade motion
            blade_st.b0 = 0
            blade_st.b1c = 0
            blade_st.b1s = 0

            # 4) Blade lift and drag
            blade_st.Cla = np.nan # 2*pi*0.9488 # clbar = 0.5987281447
            blade_st.Cl0 = np.nan # 0.6158       # Cl at zero angle of attack
            blade_st.Cd0 = np.nan # blade_st.Cla * deg2rad(5) / 29
            blade_st.d1 = np.nan # 0
            blade_st.d2 = np.nan # 0
            blade_st.c81table = True
            blade_st.get_ClCdCm = KdeAero.kde_rotor_ClCdCm

            # 5) Pitch, chord and solidity distribution
            nsections = len(y_arr) - 1
            ones_arr = np.ones(nsections + 1)
            blade_st.nsections = nsections
            blade_st.y_arr = y_arr
            blade_st.theta_arr = theta_arr + collpitch
            blade_st.chord_arr = chord_arr * ones_arr
            blade_st.sigma_arr = \
                (blade_st.Nb * blade_st.chord_arr) / (np.pi*blade_st.R)

        if not known_rotortype:
            print(rotortype)
            raise RuntimeError('Unknown rotortype')

        # Case name
        blade_st.casename = \
            'mu' + str(blade_st.mu) + '_lc' + str(blade_st.lambda_c)
        blade_st.casename = str(blade_st.casename).replace('.', 'p')

        return blade_st


class KdeAero:
    NACA5407_Station6 = scipy.io.loadmat('../mt_bet_data/NACA5407_Station6.mat')
    # NACA5407_Station6 = NACA5407_Station6['NACA5407_Station6']
    # print('kde_rotor_ClCdCm')
    # print(NACA5407_Station6)
    # dtype=[('ANGLE_Cl', 'O'), ('ANGLE_Cd', 'O'), ('ANGLE_Cm', 'O'),
    # ('MACH_Cl', 'O'), ('MACH_Cd', 'O'), ('MACH_Cm', 'O'), ('Cl', 'O'),
    # ('Cd', 'O'), ('Cm', 'O')]
    NACA5407_Station6 = NACA5407_Station6['NACA5407_Station6']
    # print(NACA5407_Station6)
    # print(len(NACA5407_Station6))
    NACA5407_Station6 = NACA5407_Station6[0]
    NACA5407_Station6 = NACA5407_Station6[0]
    # print(NACA5407_Station6)
    # print(len(NACA5407_Station6))
    NACA5407_Station6_ANGLE_Cl = NACA5407_Station6[0]
    NACA5407_Station6_ANGLE_Cd = NACA5407_Station6[1]
    NACA5407_Station6_ANGLE_Cm = NACA5407_Station6[2]
    NACA5407_Station6_MACH_Cl = NACA5407_Station6[3]
    NACA5407_Station6_MACH_Cd = NACA5407_Station6[4]
    NACA5407_Station6_MACH_Cm = NACA5407_Station6[5]
    NACA5407_Station6_Cl = NACA5407_Station6[6]
    NACA5407_Station6_Cd = NACA5407_Station6[7]
    NACA5407_Station6_Cm = NACA5407_Station6[8]

    ##########################
    # [Cl, Cd, Cm] @ mach = 0.28
    col = 3
    X = NACA5407_Station6_ANGLE_Cl[:, col]
    V = NACA5407_Station6_Cl[:, col]
    Cl_fnct = Mat.interp1_fnct(X, V, 'linear')

    X = NACA5407_Station6_ANGLE_Cd[:, col]
    V = NACA5407_Station6_Cd[:, col]
    Cd_fnct = Mat.interp1_fnct(X, V, 'linear')

    X = NACA5407_Station6_ANGLE_Cm[:, col]
    V = NACA5407_Station6_Cm[:, col]
    Cm_fnct = Mat.interp1_fnct(X, V, 'linear')
    ##########################

    @staticmethod
    def kde_rotor_ClCdCm(mach, aoa, rotortype):
        # function [Cl, Cd, Cm] = kde_rotor_ClCdCm(mach, aoa, rotortype)   
    
        # persistent NACA5407_Station6
        # if isempty(NACA5407_Station6)
        #     # filename = [bpath '/' 'NACA5407_Station6']
        #     filename = 'mt_bet_data/NACA5407_Station6'
        #     st = load(filename)
        #     NACA5407_Station6 = st.NACA5407_Station6
        naca5407_station6 = KdeAero.NACA5407_Station6
        # NACA5407_Station6 = NACA5407_Station6['NACA5407_Station6']
        # print('kde_rotor_ClCdCm')
        # print(NACA5407_Station6)
    
        # mach_arr = [0.0000 0.2849 0.2850 0.3649 0.3650 0.6000]
        # mach = 0.3649
        if mach > 0.6000:
            print('[kde_rotor_ClCdCm] Saturating mach ..')
            mach = 0.6000
        
        aoa = UnitCnv.rad2deg(aoa)
    #    if aoa > 180
    #        aoa = 180
    #    
    #    if aoa < -180
    #        aoa = -180
    #
        [Cl, Cd, Cm] = [np.nan, np.nan, np.nan]
    
        if rotortype == 'NACA5407':
            rotortype = 'NACA5407_Station6'
        
        if rotortype == 'KDECF245DP':
            rotortype = 'NACA5407_Station6'
        
        if rotortype == 'KDECF245TP':
            rotortype = 'NACA5407_Station6'
        
        if rotortype == 'KDECF245HP':
            rotortype = 'NACA5407_Station6'
        
        if rotortype == 'KDECF305DP':
            rotortype = 'NACA5407_Station6'
        
        if rotortype == 'NACA5407_Station6':
            use_fixed_mach = True
            if use_fixed_mach:
                # [Cl, Cd, Cm] @ mach = 0.28
                col = 3
    
                Xq = aoa
                # X = NACA5407_Station6.ANGLE_Cl[:, col]
                # V = NACA5407_Station6.Cl[:, col]
                # Cl = Mat.interp1(X, V, Xq, 'spline')
                Cl = KdeAero.Cl_fnct(Xq)
    
                Xq = aoa
                # X = NACA5407_Station6.ANGLE_Cd[:, col]
                # V = NACA5407_Station6.Cd[:, col]
                # Cd = Mat.interp1(X, V, Xq, 'spline')
                Cd = KdeAero.Cd_fnct(Xq)
    
                Xq = aoa
                # X = NACA5407_Station6.ANGLE_Cm[:, col]
                # V = NACA5407_Station6.Cm[:, col]
                # Cm = Mat.interp1(X, V, Xq, 'spline')
                Cm = KdeAero.Cm_fnct(Xq)

            else:
                # Vq = interp2(X,Y,V,Xq,Yq)        
                Xq = mach
                Yq = aoa
                X = naca5407_station6.MACH_Cl
                Y = naca5407_station6.ANGLE_Cl
                V = naca5407_station6.Cl
                Cl = Mat.interp2(X, Y, V, Xq, Yq, 'linear')
    
                Xq = mach
                Yq = aoa
                X = naca5407_station6.MACH_Cd
                Y = naca5407_station6.ANGLE_Cd
                V = naca5407_station6.Cd
                Cd = Mat.interp2(X, Y, V, Xq, Yq, 'linear')
    
                Xq = mach
                Yq = aoa
                X = naca5407_station6.MACH_Cm
                Y = naca5407_station6.ANGLE_Cm
                V = naca5407_station6.Cm
                Cm = Mat.interp2(X, Y, V, Xq, Yq, 'linear')

            modify_coeffs = GnrlConfig.tunning_for_kde245['modify_ClCdCm']
            if sum(modify_coeffs) != 3:
                Cl = Cl * modify_coeffs[0]
                Cd = Cd * modify_coeffs[1]
                Cm = Cm * modify_coeffs[2]

        if (np.isnan(Cl) is True) or \
                (np.isfinite(Cl) is False) or (np.isreal(Cl) is False):
            print(Cl)
            print(mach)
            print(aoa)
            print(rotortype)
            raise RuntimeError('Cl is NaN, Inf or Complex')
        
        if (np.isnan(Cd) is True) or \
                (np.isfinite(Cd) is False) or (np.isreal(Cd) is False):
            print(Cd)
            print(mach)
            print(aoa)
            print(rotortype)
            raise RuntimeError('Cd is NaN, Inf or Complex')
        
        if (np.isnan(Cm) is True) or \
                (np.isfinite(Cm) is False) or (np.isreal(Cm) is False):
            print(Cm) 
            print(mach) 
            print(aoa) 
            print(rotortype)
            raise RuntimeError('Cm is NaN, Inf or Complex')

        return [Cl, Cd, Cm]


class KdeGeometry:
    kde245_geometry = scipy.io.loadmat('../mt_bet_data/kde245_geometry.mat')
    kde245_geometry = kde245_geometry['kde245_geometry']

    kde305_geometry = scipy.io.loadmat('../mt_bet_data/kde305_geometry.mat')
    kde305_geometry = kde305_geometry['kde305_geometry']

    naca5407_geometry = scipy.io.loadmat('../mt_bet_data/naca5407_geometry.mat')
    naca5407_geometry = naca5407_geometry['naca5407_geometry']

    @staticmethod
    def kde_rotor_defaults(rotortype):
        # function [rho, lambda_c, mu, collpitch] = \
        #     kde_rotor_defaults(rotortype)

        # #    rho = 1.225            # density at NASA Langley WT
        # #    rho = 1.2115900277     # density at KDE
        # #    rho = 1.1849           # density at AneCha
        # #    rho = 1.1388           # density at AneCha @(400 msl, 25 deg)
        #
        # # collpitch = -deg2rad(0.0)
        # # collpitch = -deg2rad(1.2)
        # # collpitch = -deg2rad(2.5)

        if rotortype == 'KDECF305DP':
            # rotortype = 'KDECF305DP'
            rho = 1.225              # density at NASA Langley WT
            lambda_c = 0.0
            mu = 0.0
            collpitch = -UnitCnv.deg2rad(2.5)
            # omega = omega
            # CT_target = NaN       # CT_target is only needed for mu != 0

            # blade_st = blade_model(...
            #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
            return [rho, lambda_c, mu, collpitch]

        if rotortype == 'KDECF245DP':
            # rotortype = 'KDECF305DP'
            # rotortype = 'KDECF245DP'
            # rho = 1.2115900277     # density at KDE
            # rho = 1.1849           # density at AneCha
            rho = 1.1388           # density at AneCha @(400 msl, 25 deg)
            lambda_c = 0.0
            mu = 0.0
            # collpitch = -deg2rad(1.2)
            collpitch = -UnitCnv.deg2rad(0.0)
            # omega = omega
            # CT_target = NaN       # CT_target is only needed for mu != 0

            # blade_st = blade_model(...
            #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
            return [rho, lambda_c, mu, collpitch]

        if rotortype == 'KDECF245TP':
            # rotortype = 'KDECF305DP'
            rotortype = 'KDECF245DP'
            # rho = 1.2115900277     # density at KDE
            # rho = 1.1849           # density at AneCha
            rho = 1.1388           # density at AneCha @(400 msl, 25 deg)
            lambda_c = 0.0
            mu = 0.0
            # collpitch = -deg2rad(1.2)
            collpitch = -UnitCnv.deg2rad(0.0)
            # omega = omega
            # CT_target = NaN     # CT_target is only needed for mu != 0

            # blade_st = blade_model(...
            #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
            return [rho, lambda_c, mu, collpitch]

        if rotortype == 'KDECF245HP':
            # rotortype = 'KDECF305DP'
            rotortype = 'KDECF245DP'
            # rho = 1.2115900277     # density at KDE
            # rho = 1.1849           # density at AneCha
            rho = 1.1388           # density at AneCha @(400 msl, 25 deg)
            lambda_c = 0.0
            mu = 0.0
            # collpitch = -deg2rad(1.2)
            collpitch = -UnitCnv.deg2rad(0.0)
            # omega = omega
            # CT_target = NaN      # CT_target is only needed for mu != 0

            # blade_st = blade_model(...
            #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)

            return [rho, lambda_c, mu, collpitch]

        raise RuntimeError('No rotortype found')

    @staticmethod
    def kde_rotor_geometry_cf245():
        kde245_geometry = KdeGeometry.kde245_geometry

        # Use KDECF305 scaled down as geoemtry for KDECF245
        # 27.5 x 8.1
        # sizefactor = 24.5 / 30.5
        # pitchfactor = 8.1 / 9.7
        # kde305_geometry = KdeGeometry.kde305_geometry
        # # print(kde305_geometry)
        # # print(len(kde305_geometry))
        # y_arr = kde305_geometry[0]
        # y_arr = y_arr[0]
        # chord_arr = kde305_geometry[1]
        # chord_arr = chord_arr[0]
        # theta_arr = kde305_geometry[2]
        # theta_arr = theta_arr[0]
        # # print(y_arr)
        # # print(len(y_arr))
        # # exit(0)
        # y_arr = y_arr * sizefactor
        # chord_arr = chord_arr * sizefactor
        # theta_arr = theta_arr * pitchfactor

        # Use KDECF245 directly
        y_arr = kde245_geometry['y_arr']
        chord_arr = kde245_geometry['chord_arr']
        theta_arr = kde245_geometry['theta_arr']
        # print('y_arr %s' # y_arr)
        # print('y_arr[0] %s' # y_arr[0])
        # print('type(y_arr[0]) %s' # type(y_arr[0]))
        y_arr = np.array(y_arr).flatten()
        chord_arr = np.array(chord_arr).flatten()
        theta_arr = np.array(theta_arr).flatten()
        y_arr = y_arr[0]
        chord_arr = chord_arr[0]
        theta_arr = theta_arr[0]
        y_arr = y_arr.flatten()
        chord_arr = chord_arr.flatten()
        theta_arr = theta_arr.flatten()
        # print('y_arr %s' # y_arr)
        # print('y_arr[0] %s' # y_arr[0])
        # print('type(y_arr[0]) %s' # type(y_arr[0]))
        # exit(0)
        # bet_cf245_pitch

        pitch0 = GnrlConfig.tunning_for_kde245['modify_collective']
        if pitch0 != 0:
            theta_arr = theta_arr + UnitCnv.deg2rad(pitch0)

        return [y_arr, chord_arr, theta_arr]

    @staticmethod
    def kde_rotor_geometry(rotortype):
        # function [y_arr, chord_arr, theta_arr] = kde_rotor_geometry(rotortype)

        cm2m = UnitCnv.cm2m
        in2cm = UnitCnv.in2cm

        y_arr = None
        chord_arr = None
        theta_arr = None

        if (rotortype == 'KDECF245DP') or \
                (rotortype == 'KDECF245TP') or (rotortype == 'KDECF245HP'):

            [y_arr, chord_arr, theta_arr] = \
                KdeGeometry.kde_rotor_geometry_cf245()

        elif rotortype == 'KDECF275DP':
            # 27.5 x 8.9
            sizefactor = 27.5 / 30.5
            pitchfactor = 8.9 / 9.7

            kde305_geometry = KdeGeometry.kde305_geometry
            y_arr = kde305_geometry.y_arr * sizefactor
            chord_arr = kde305_geometry.chord_arr * sizefactor
            theta_arr = kde305_geometry.theta_arr * pitchfactor

        elif rotortype == 'KDECF305DP':
            # 30.5 x 9.7
            kde305_geometry = KdeGeometry.kde305_geometry
            y_arr = kde305_geometry.y_arr
            chord_arr = kde305_geometry.chord_arr
            theta_arr = kde305_geometry.theta_arr

        elif rotortype == 'NACA5407shape':
            R = 15.25 * in2cm * cm2m
            nsections = 11
            naca5407_geometry = KdeGeometry.naca5407_geometry
            y_arr = np.linspace(0, R, nsections)
            chord_arr = naca5407_geometry.shape.x_arr
            theta_arr = naca5407_geometry.shape.y_arr

        elif rotortype == 'NACA5407':
            naca5407_geometry = KdeGeometry.naca5407_geometry
            y_arr = naca5407_geometry.y_arr
            chord_arr = naca5407_geometry.chord_arr
            theta_arr = naca5407_geometry.theta_arr

        return [y_arr, chord_arr, theta_arr]

