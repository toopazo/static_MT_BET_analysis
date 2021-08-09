import numpy as np
# import matplotlib as plt
from unit_conversion import UnitConversion as UnitCnv
from bet_loads import Bet, BetResults, BetFindOmega
from mt_loads import MT
from rotor_model import BladeSt, KdeGeometry
from matmath import MatMath as Mat
# from anecha_data import AnechaData
# from toopazo_tools.matplotlib import FigureTools, PlotTools
from gnrl_config import GnrlConfig
# from scipy import io
# import collections
import pickle
import copy


class CoaxBetResults:
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
        self.bet_u_st = None    # BetResults()
        self.bet_l_st = None    # BetResults()

        self._freeze()  # no new attributes after this point.


class Dcollpitch:
    __isfrozen = False  # name mangling

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
        self.dcollpitch_arr = None
        self.collpitch_u_arr = None
        self.collpitch_l_arr = None
        self.P_h = None
        self.W = None
        self.coaxu_r0 = None
        self.coaxl_r0 = None

        self.eta_T_arr = None

        self.blade_h_st = BladeSt()
        self.bet_h_st = BetResults()
        self.T_h = None
        self.Q_h = None
        self.P_h = None
        self.P_0 = None
        self.coax = CoaxBetResults()
        self.sbs = CoaxBetResults()

        self._freeze()  # no new attributes after this point.

    @staticmethod
    def writetable(filename, data2d):
        mdict = {}
        cnt = 1
        # nvariables = len(data2d)
        nsamples = None
        for v_arr in data2d:
            mkey = 'Var%s' % cnt
            mdict[mkey] = v_arr
            nsamples = len(v_arr)
            cnt = cnt + 1

        # table = [
        #     Var1 = table_m1_thrust,
        #     Var2 = table_m2_thrust,
        #     Var3 = table_eta_trust,
        #     Var4 = table_m1_omega,
        #     Var5 = table_m2_omega,
        #     Var6 = table_eta_omega,
        #     Var7 = table_m1_power,
        #     Var8 = table_m2_power,
        #     Var9 = table_coax_power
        # ]

        fd = open(filename, 'w')

        line = ''
        for key, val in sorted(mdict.items()):
            # print('[writetable] key %s, val %s' % (key, val))
            line = line + key + ','
        line = line[0:-1]
        line = line + '\n'
        fd.write(line)

        for i in range(0, nsamples):
            line = ''
            for key, val in sorted(mdict.items()):
                # print('[writetable] key %s, val %s' % (key, val))
                line = line + str(val[i]) + ','
            line = line[0:-1]
            line = line + '\n'
            fd.write(line)

        fd.close()

    @staticmethod
    def bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1):
        # function bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1)

        assert isinstance(db, Dcollpitch)

        dcollpitch_arr = db.dcollpitch_arr
        collpitch_u_arr = db.collpitch_u_arr
        collpitch_l_arr = db.collpitch_l_arr
        P_h = db.P_h
        thrust_tot = db.W
        coaxu_r0 = db.coaxu_r0
        coaxl_r0 = db.coaxl_r0
        rho = db.blade_h_st.rho
        rotArea = db.blade_h_st.rotArea
        R = db.blade_h_st.R

        # i = dcollpitch
        # j = eta_thrust

        num_j = len(db.eta_T_arr)

        Tcoax_u_arr = np.zeros(num_j)
        Tcoax_l_arr = np.zeros(num_j)
        eta_Tcoax_arr = np.zeros(num_j)
        lambdacoax_u_arr = np.zeros(num_j)
        lambdacoax_l_arr = np.zeros(num_j)
        eta_lambdacoax_arr = np.zeros(num_j)
        vu_arr = np.zeros(num_j)
        vl_arr = np.zeros(num_j)

        bet_u_st = db.coax.bet_u_st[0, 0]
        len_total_lcr = len(bet_u_st.total_lcr)
        lrcoax_u_arr = np.zeros((num_j, len_total_lcr))
        lrcoax_l_arr = np.zeros((num_j, len_total_lcr))
        omegacoax_u_arr = np.zeros(num_j)
        omegacoax_l_arr = np.zeros(num_j)
        eta_omegacoax_arr = np.zeros(num_j)
        Qcoax_u_arr = np.zeros(num_j)
        Qcoax_l_arr = np.zeros(num_j)
        eta_Qcoax_arr = np.zeros(num_j)
        Pcoax_u_arr = np.zeros(num_j)
        Pcoax_l_arr = np.zeros(num_j)
        eta_Pcoax_arr = np.zeros(num_j)
        Pcoax_arr = np.zeros(num_j)

        Tsbs_u_arr = np.zeros(num_j)
        Tsbs_l_arr = np.zeros(num_j)
        eta_Tsbs_arr = np.zeros(num_j)
        lambdasbs_u_arr = np.zeros(num_j)
        lambdasbs_l_arr = np.zeros(num_j)
        eta_lambdasbs_arr = np.zeros(num_j)
        omegasbs_u_arr = np.zeros(num_j)
        omegasbs_l_arr = np.zeros(num_j)
        eta_omegasbs_arr = np.zeros(num_j)
        Qsbs_u_arr = np.zeros(num_j)
        Qsbs_l_arr = np.zeros(num_j)
        eta_Qsbs_arr = np.zeros(num_j)
        Psbs_u_arr = np.zeros(num_j)
        Psbs_l_arr = np.zeros(num_j)
        Psbs_arr = np.zeros(num_j)
        kint_arr = np.zeros(num_j)
        mt_vind_arr = np.zeros(num_j)

        mt_lambdacoax_u_arr = np.zeros(num_j)
        mt_vu_arr = np.zeros(num_j)
        mt_lambdacoax_l_arr = np.zeros(num_j)
        mt_vl_arr = np.zeros(num_j)

        i = np.argwhere(db.dcollpitch_arr == dcollpitch)
        i = i.squeeze()
        for j in range(0, num_j):   # = 1:length(db.eta_T_arr)
            bet_u_st = db.coax.bet_u_st[i, j]
            bet_l_st = db.coax.bet_l_st[i, j]

            assert isinstance(bet_u_st, BetResults)
            assert isinstance(bet_l_st, BetResults)

            radius_u = bet_u_st.blade_st.R
            radius_l = bet_l_st.blade_st.R
            angvel_u = bet_u_st.blade_st.omega
            angvel_l = bet_l_st.blade_st.omega
            veltip_u = angvel_u * radius_u
            veltip_l = angvel_l * radius_l

            Tcoax_u_arr[j] = bet_u_st.total_T
            Tcoax_l_arr[j] = bet_l_st.total_T
            eta_Tcoax_arr[j] = Tcoax_u_arr[j] / Tcoax_l_arr[j]

            lambdacoax_u_arr[j] = bet_u_st.total_lambda_t
            lambdacoax_l_arr[j] = bet_l_st.total_lambda_t
            eta_lambdacoax_arr[j] = lambdacoax_u_arr[j] / lambdacoax_l_arr[j]

            vu_arr[j] = lambdacoax_u_arr[j] * veltip_u
            vl_arr[j] = lambdacoax_l_arr[j] * veltip_l

            # distribution along r
            lrcoax_u_arr[j, :] = bet_u_st.total_lcr + bet_u_st.total_lir
            lrcoax_l_arr[j, :] = bet_l_st.total_lcr + bet_l_st.total_lir

            omegacoax_u_arr[j] = angvel_u
            omegacoax_l_arr[j] = angvel_l
            eta_omegacoax_arr[j] = omegacoax_u_arr[j] / omegacoax_l_arr[j]

            Qcoax_u_arr[j] = bet_u_st.total_Q
            Qcoax_l_arr[j] = bet_l_st.total_Q
            eta_Qcoax_arr[j] = Qcoax_u_arr[j] / Qcoax_l_arr[j]

            Pcoax_u_arr[j] = bet_u_st.total_P
            Pcoax_l_arr[j] = bet_l_st.total_P
            eta_Pcoax_arr[j] = Pcoax_u_arr[j] / Pcoax_l_arr[j]
            Pcoax_arr[j] = Pcoax_u_arr[j] + Pcoax_l_arr[j]

            ###################################
            bet_u_st = db.sbs.bet_u_st[i, j]
            bet_l_st = db.sbs.bet_l_st[i, j]

            Tsbs_u_arr[j] = bet_u_st.total_T
            Tsbs_l_arr[j] = bet_l_st.total_T
            eta_Tsbs_arr[j] = Tsbs_u_arr[j] / Tsbs_l_arr[j]

            lambdasbs_u_arr[j] = bet_u_st.total_lambda_t
            lambdasbs_l_arr[j] = bet_l_st.total_lambda_t
            eta_lambdasbs_arr[j] = lambdasbs_u_arr[j] / lambdasbs_l_arr[j]

            omegasbs_u_arr[j] = bet_u_st.blade_st.omega
            omegasbs_l_arr[j] = bet_l_st.blade_st.omega
            eta_omegasbs_arr[j] = omegasbs_u_arr[j] / omegasbs_l_arr[j]

            Qsbs_u_arr[j] = bet_u_st.total_Q
            Qsbs_l_arr[j] = bet_l_st.total_Q
            eta_Qsbs_arr[j] = Qsbs_u_arr[j] / Qsbs_l_arr[j]

            Psbs_u_arr[j] = bet_u_st.total_P
            Psbs_l_arr[j] = bet_l_st.total_P
            Psbs_arr[j] = Psbs_u_arr[j] + Psbs_l_arr[j]

            ###################################
            kint_arr[j] = Pcoax_arr[j] / Psbs_arr[j]

            ###################################
            mt_vind_arr[j] = np.sqrt(thrust_tot / (2 * rho * np.pi * R ** 2))

            ###################################
            # Assuming only upper on lower interference and perfect contraction
            # vi_u = np.sqrt(m1_thrust / (2 * density * area))
            # vc_l = vi_u
            # vh_l = np.sqrt(m2_thrust / (2 * density * area))
            # vi_l = vh_l * (
            #         - (0.5 * vc_l / vh_l) + np.sqrt((0.5 * vc_l / vh_l) ** 2 + 1)
            # )

            # MT inflow for upper rotor
            Vtip = omegacoax_u_arr[j] * R
            CT = Tcoax_u_arr[j] / (rho * rotArea * Vtip**2)
            mu = 0
            TPP_alpha = 0
            lambda_c = 0
            mt_lambdacoax_u_arr[j] = MT.mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c)
            mt_vu_arr[j] = mt_lambdacoax_u_arr[j] * Vtip

            ###################################
            # MT inflow for lower rotor
            Vtip = omegacoax_l_arr[j] * R
            CT = Tcoax_l_arr[j] / (rho*rotArea*Vtip**2)
            mu = 0
            TPP_alpha = 0
            lambda_c = mt_lambdacoax_u_arr[j]
            mt_lambdacoax_l_arr[j] = MT.mt_inflow(CT, mu, TPP_alpha, Vtip, lambda_c)
            mt_vl_arr[j] = mt_lambdacoax_l_arr[j] * Vtip

            # ###################################
            # # MT inflow for lower rotor using extra RPM
            # omega_l = omegacoax_u_arr[j]*1.2
            # csvp = db
            # m1_thrust = Tcoax_u_arr[j]
            # omega_u = omegacoax_u_arr[j]
            # m2_thrust = Tcoax_l_arr[j]
            # [thrust, a_thrust] = mt_coaxial_lower_rotor_model_jack(
            #     omega_l, csvp, m1_thrust, omega_u, m2_thrust)

        # print('[bet_plot_coax_dcollpitch] dcollpitch %.4f, i %d, P_h %.4f' %
        #       (dcollpitch, i, P_h))
        # print('[bet_plot_coax_dcollpitch] collpitch_u %.4f, collpitch_l %.4f' %
        #       (UnitCnv.rad2deg(collpitch_u_arr[i]),
        #        UnitCnv.rad2deg(collpitch_l_arr[i])))

        # rpm2rads = np.pi / 30
        # rads2rpm = 30 / np.pi

        # Save data to file
        table_m1_thrust = Tcoax_u_arr
        table_m2_thrust = Tcoax_l_arr
        table_eta_trust = eta_Tcoax_arr
        table_m1_omega = omegacoax_u_arr
        table_m2_omega = omegacoax_l_arr
        table_eta_omega = eta_omegacoax_arr
        table_m1_power = Qcoax_u_arr * table_m1_omega
        table_m2_power = Qcoax_l_arr * table_m2_omega
        table_coax_power = Pcoax_arr
        table = [
            table_m1_thrust,
            table_m2_thrust,
            table_eta_trust,
            table_m1_omega,
            table_m2_omega,
            table_eta_omega,
            table_m1_power,
            table_m2_power,
            table_coax_power
        ]

        # arg1 = ['_' rotortype '_' str(W) '_'
        #     str(coaxu_r0) '_' str(coaxl_r0)]
        # writetable(T, ['img/bet_plot_coax_dcollpitch' arg1 '.txt'])
        str2 = rotortype + '_' + str(round(thrust_tot)) + '_' + str(coaxu_r0) \
               + '_' + str(coaxl_r0)
        filename = 'img/' + str1 + '_' + str2 + '.txt'
        Dcollpitch.writetable(filename, table)

        Dcollpitch.bet_plot_coax_dcollpitch_part2(
            rotortype, dcollpitch, db, str1)

    @staticmethod
    def bet_plot_coax_dcollpitch_part2(rotortype, dcollpitch, db, str1):
        pass
    #     nfig = 1
    #     figure(nfig)
    #     subplot(2, 1, 1)
    #     hold on
    #     grid on
    #     plot(eta_Tcoax_arr, Tcoax_u_arr, '-r*')
    #     plot(eta_Tcoax_arr, Tcoax_l_arr, '-b*')
    #     W_arr = Tcoax_u_arr + Tcoax_l_arr
    #     plot(eta_Tcoax_arr, W_arr, '-k*')
    #     legend('coax u', 'coax l', 'coax W', 'Location','NorthEast')
    #     xlabel('\eta_T')
    #     ylabel('Thrust N')
    #     if rotortype == 'KDECF305DP':
    #         ylim([0, 300])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         xlim([0, 5])
    #         ylim([0, 70])
    #         yticks([0, 10, 20, 30, 40, 50, 60, 70])
    #     end
    #
    #     subplot(2, 1, 2)
    #     hold on
    #     grid on
    #     plot(eta_Tsbs_arr, Tsbs_u_arr, 'r-*')
    #     plot(eta_Tsbs_arr, Tsbs_l_arr, 'b-*')
    #     W_arr = Tsbs_u_arr + Tsbs_l_arr
    #     plot(eta_Tsbs_arr, W_arr, '-k*')
    #     legend('sbs u', 'sbs l', 'sbs W', 'Location','NorthEast')
    #     xlabel('\eta_T')
    #     ylabel('Thrust N')
    #     if rotortype == 'KDECF305DP':
    #         ylim([0, 300])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         xlim([0, 5])
    #         ylim([0, 70])  yticks([0, 10, 20, 30, 40, 50, 60, 70])
    #     end
    #
    #     ##################################################
    #     nfig = 2
    #     figure(nfig)
    #     subplot(2, 1, 1)
    #     hold on
    #     grid on
    #     yyaxis left
    #     plot(eta_Tcoax_arr, omegacoax_u_arr*rads2rpm, 'r-*')
    #     plot(eta_Tcoax_arr, omegacoax_l_arr*rads2rpm, 'b-*')
    #     xlabel('\eta_T')
    #     ylabel('\Omega RPM')
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tcoax_arr(1) - 0.3, omegacoax_u_arr(1)*rads2rpm,
    #             [str(dcollpitch) char(176)])
    #         text(eta_Tcoax_arr(1) - 0.3, omegacoax_l_arr(1)*rads2rpm,
    #             [str(dcollpitch) char(176)])
    #         ylim([2000, 7000])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 4000])
    #         xlim([0, 5])
    #     end
    #     yyaxis right
    #     plot(eta_Tcoax_arr, eta_omegacoax_arr, 'k-*')
    #     legend('coax u', 'coax l', '\eta_{\Omega}', 'Location','SouthEast')
    #     ylabel('\eta_{\Omega}')
    #     ax = gca
    #     ax.YAxis(1).Color = 'k'
    #     ax.YAxis(2).Color = 'k'
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tcoax_arr(end) + 0.1, eta_omegacoax_arr(end),
    #             [str(dcollpitch) char(176)])
    #         ylim([0, 3])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 2])
    #         xlim([0, 5])
    #     end
    #
    #     subplot(2, 1, 2)
    #     hold on
    #     grid on
    #     yyaxis left
    #     plot(eta_Tsbs_arr, omegasbs_u_arr*rads2rpm, 'r-*')
    #     plot(eta_Tsbs_arr, omegasbs_l_arr*rads2rpm, 'b-*')
    #     legend('sbs u', 'sbs l', 'Location','southeast')
    #     xlabel('\eta_T')
    #     ylabel('\Omega RPM')
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tsbs_arr(1) - 0.3, omegasbs_u_arr(1)*rads2rpm,
    #             [str(dcollpitch) char(176)])
    #         text(eta_Tsbs_arr(1) - 0.3, omegasbs_l_arr(1)*rads2rpm,
    #             [str(dcollpitch) char(176)])
    #         ylim([2000, 7000])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 4000])
    #         xlim([0, 5])
    #     end
    #     yyaxis right
    #     plot(eta_Tsbs_arr, eta_omegasbs_arr, 'k-*')
    #     legend('sbs u', 'sbs l', '\eta_{\Omega}', 'Location','southeast')
    #     ylabel('\eta_{\Omega}')
    #     ax = gca
    #     ax.YAxis(1).Color = 'k'
    #     ax.YAxis(2).Color = 'k'
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tsbs_arr(end) + 0.1, eta_omegasbs_arr(end),
    #             [str(dcollpitch) char(176)])
    #         ylim([0, 3])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 2])
    #         xlim([0, 5])
    #     end
    #
    #     ##################################################
    #     nfig = 3
    #     figure(nfig)
    #     subplot(2, 1, 1)
    #     hold on
    #     grid on
    #     plot(omegacoax_u_arr*rads2rpm, Tcoax_u_arr, 'r-*')
    #     plot(omegacoax_l_arr*rads2rpm, Tcoax_l_arr, 'b-*')
    #     legend('coax u', 'coax l', 'Location','northwest')
    #     xlabel('\Omega RPM')
    #     ylabel('Thrust N')
    #     if rotortype == 'KDECF305DP':
    #         text(omegacoax_u_arr(end)*rads2rpm + 10, Tcoax_u_arr(end) + 20,
    #             [str(dcollpitch) char(176)])
    #         text(omegacoax_l_arr(end)*rads2rpm - 10, Tcoax_l_arr(end) - 20,
    #             [str(dcollpitch) char(176)])
    #         xlim([2000, 7000])
    #         ylim([0, 300])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         xlim([0, 4000])
    #         ylim([0, 70]) yticks([0, 10, 20, 30, 40, 50, 60, 70])
    #     end
    #
    #     subplot(2, 1, 2)
    #     hold on
    #     grid on
    #     plot(omegasbs_u_arr*rads2rpm, Tsbs_u_arr, 'r-*')
    #     plot(omegasbs_l_arr*rads2rpm, Tsbs_l_arr, 'b-*')
    #     legend('sbs u', 'sbs l', 'Location','northwest')
    #     xlabel('\Omega RPM')
    #     ylabel('Thrust N')
    #     if rotortype == 'KDECF305DP':
    #         text(omegasbs_u_arr(end)*rads2rpm + 10, Tsbs_u_arr(end) + 20,
    #             [str(dcollpitch) char(176)])
    #         text(omegasbs_l_arr(end)*rads2rpm - 10, Tsbs_l_arr(end) - 20,
    #             [str(dcollpitch) char(176)])
    #         xlim([2000, 7000])
    #         ylim([0, 300])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         xlim([0, 4000])
    #         ylim([0, 70]) yticks([0, 10, 20, 30, 40, 50, 60, 70])
    #     end
    #
    #     ##################################################
    #     nfig = 4
    #     figure(nfig)
    #     hold on
    #     grid on
    #     yyaxis left
    #     # plot(eta_Tcoax_arr, Pcoax_arr / P_h, 'r-*')
    #     # plot(eta_Tsbs_arr, Psbs_arr / P_h, 'b-*')
    #     plot(eta_Tcoax_arr, Pcoax_arr, 'r-*')
    #     plot(eta_Tsbs_arr, Psbs_arr, 'b-*')
    #     xlabel('\eta_T')
    #     ylabel('Power W')
    #     if rotortype == 'KDECF305DP':
    #         ylim([0, 2000]) # ylim([0, 2.0])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 750]) yticks([0, 250, 500, 750])
    #         xlim([0, 5])
    #     end
    #     yyaxis right
    #     plot(eta_Tcoax_arr, kint_arr, 'k-*')
    #     legend('coax', 'sbs', 'k_{int}', 'Location','northeast')
    #     ylabel('k_{int}')
    #     ax = gca
    #     ax.YAxis(1).Color = 'k'
    #     ax.YAxis(2).Color = 'k'
    #     if rotortype == 'KDECF305DP':
    #         text(0.1, kint_arr(1), [str(dcollpitch) char(176)])
    #         ylim([0, 3.0])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 3.0])
    #         xlim([0, 5])
    #     end
    #
    #     ##################################################
    #     nfig = 5
    #     figure(nfig)
    #     subplot(2, 1, 1)
    #     hold on
    #     grid on
    #     yyaxis left
    #     plot(eta_Tcoax_arr, Qcoax_u_arr, 'r-*')
    #     plot(eta_Tcoax_arr, Qcoax_l_arr, 'b-*')
    #     xlabel('\eta_T')
    #     ylabel('Torque Nm')
    #     if rotortype == 'KDECF305DP':
    # #        text(eta_Tcoax_arr(end) + 0.1, Qcoax_l_arr(end) + 0.1,
    # #            [str(dcollpitch) char(176)])
    #         ylim([0, 10])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 2.0])
    #         xlim([0, 5])
    #     end
    #     yyaxis right
    #     plot(eta_Tcoax_arr, eta_Qcoax_arr, 'k-*')
    #     legend('coax u', 'coax l', '\eta_Q', 'Location','north')
    #     ylabel('\eta_Q')
    #     ax = gca
    #     ax.YAxis(1).Color = 'k'
    #     ax.YAxis(2).Color = 'k'
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tcoax_arr(end) + 0.1, eta_Qcoax_arr(end) + 0.1,
    #             [str(dcollpitch) char(176)])
    #         ylim([0, 5])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 5.0])
    #         xlim([0, 5])
    #     end
    #
    #     subplot(2, 1, 2)
    #     hold on
    #     grid on
    #     yyaxis left
    #     plot(eta_Tsbs_arr, Qsbs_u_arr, 'r-*')
    #     plot(eta_Tsbs_arr, Qsbs_l_arr, 'b-*')
    #     xlabel('\eta_T')
    #     ylabel('Torque Nm')
    #     if rotortype == 'KDECF305DP':
    # #        text(eta_Tsbs_arr(end) + 0.1, Qsbs_l_arr(end) + 0.1,
    # #            [str(dcollpitch) char(176)])
    #         ylim([0, 10])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 2.0])
    #         xlim([0, 5])
    #     end
    #     yyaxis right
    #     plot(eta_Tsbs_arr, eta_Qsbs_arr, 'k-*')
    #     legend('sbs u', 'sbs l', '\eta_Q', 'Location','north')
    #     ylabel('\eta_Q')
    #     ax = gca
    #     ax.YAxis(1).Color = 'k'
    #     ax.YAxis(2).Color = 'k'
    #     if rotortype == 'KDECF305DP':
    #         text(eta_Tsbs_arr(end) + 0.1, eta_Qsbs_arr(end) + 0.1,
    #             [str(dcollpitch) char(176)])
    #         ylim([0, 5])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         ylim([0, 5])
    #         xlim([0, 5])
    #     end
    #
    #     pause(30*0)
    #
    #     ##################################################
    #     nfig = 6
    #     figure(nfig)
    #     subplot(2, 1, 1)
    #     hold on
    #     grid on
    #     # yyaxis left
    #     # plot(eta_Tcoax_arr, lambdacoax_u_arr, 'r-*')
    #     # plot(eta_Tcoax_arr, lambdacoax_l_arr, 'b-*')
    #     # plot(eta_Tcoax_arr, mt_lambdacoax_u_arr, 'r-*')
    #     # plot(eta_Tcoax_arr, mt_lambdacoax_l_arr, 'b-*')
    #     # plot(eta_Tcoax_arr, 0.0708 * ones(size(mt_lambdacoax_l_arr)), 'k--*')
    #     plot(eta_Tcoax_arr, vu_arr, 'r-*')
    #     plot(eta_Tcoax_arr, vl_arr, 'b-*')
    #     plot(eta_Tcoax_arr, mt_vu_arr, 'r-*')
    #     plot(eta_Tcoax_arr, mt_vl_arr, 'b-*')
    #     plot(eta_Tcoax_arr, mt_vind_arr, 'k--*')
    #     xlabel('\eta_T')
    #     # ylabel('\lambda = \lambda_c + \lambda_i')
    #     ylabel('v = v_c + v_i')
    #     if rotortype == 'KDECF305DP':
    # #        text(eta_Tcoax_arr(end) + 0.1, lambdacoax_l_arr(end) + 0.1,
    # #            [str(dcollpitch) char(176)])
    #         ylim([0, 10])
    #         xlim([0, 4.5])
    #     end
    #     if rotortype == 'KDECF245DP':
    #         # ylim([0, 0.15])
    #         ylim([0, 15])
    #         xlim([0, 5])
    #     end
    #     legend('coax u', 'coax l', 'mt coax u', 'mt coax l', 'mt single',
    #         'Location','northeast')
    #
    #     # lambda_mtu = mt_lambdacoax_u_arr(1)
    #     # lambda_mtl = mt_lambdacoax_l_arr(1)
    #     # lambda_mtu =
    #     #     0.0659
    #     # lambda_mtl =
    #     #     0.1016
    #
    #     subplot(2, 1, 2)
    #     hold on
    #     grid on
    #     for j = 1:length(db.eta_T_arr)
    #         plot(linspace(0, 1, length(lrcoax_u_arr[j, :])), lrcoax_u_arr[j, :], '-r')
    #         plot(linspace(0, 1, length(lrcoax_l_arr[j, :])), lrcoax_l_arr[j, :], '-b')
    #     end
    #     xlabel('Normalized radius')
    #     ylabel('BET \lambda = \lambda_c + \lambda_i')
    #     # ylim([0, 2])
    #     # xlim([0, 5])
    #
    #     ##################################################
    #     # str1 = 'bet_plot_coax_dcollpitch'
    #     # str2 = rotortype
    #     str2 = [rotortype '_' str(W)
    #         '_' str(coaxu_r0) '_' str(coaxl_r0)]
    #     nfig_arr = [1, 2, 3, 4, 5, 6]
    #     savefig = False
    #     plot_save_nfig_arr(str1, str2, nfig_arr, savefig)
    #
    #     ##################################################
    #     if (length(dcollpitch_arr) > 1) &and (length(eta_T_arr) > 1):
    #         bet_plot_coax_Zsurf(rotortype, db)


class BetMainCoaxial:
    @staticmethod
    def load_from_mat(filename):
        # mdict = io.loadmat(filename)
        filename = filename.replace('.mat', '.pickle')
        mdict = BetMainCoaxial.load_from_pickle(filename)
        return mdict

    @staticmethod
    def save_to_mat(filename, mdict):
        # scipy.io.savemat(
        #     file_name, mdict, appendmat=True, format='5',
        #     long_field_names=False, do_compression=False, oned_as='row')
        # io.savemat(filename, mdict)
        filename = filename.replace('.mat', '.pickle')
        BetMainCoaxial.save_to_pickle(filename, mdict)

    @staticmethod
    def save_to_pickle(filename, obj):
        try:
            print('[save_to_pickle] filename %s' % filename)
            with open(filename, "wb") as f:
                pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        except Exception as ex:
            print("Error during pickling object (Possibly unsupported):", ex)

    @staticmethod
    def load_from_pickle(filename):
        try:
            with open(filename, "rb") as f:
                return pickle.load(f)
        except Exception as ex:
            print("Error during unpickling object (Possibly unsupported):", ex)

    @staticmethod
    def bet_sbs_forces(bladesbs_u_st, bladesbs_l_st):

        bet_u_st = Bet.bet_forces(bladesbs_u_st)
        bet_u_st = Bet.bet_forces_add_total(bet_u_st, False)

        bet_l_st = Bet.bet_forces(bladesbs_l_st)
        bet_l_st = Bet.bet_forces_add_total(bet_l_st, False)

        return [bet_u_st, bet_l_st]

    @staticmethod
    def bet_coax_forces(blade_u_st, blade_l_st):
        # function [bet_u_st, bet_l_st] = \
        #     bet_coax_forces(blade_u_st, blade_l_st)

        assert isinstance(blade_u_st, BladeSt)
        assert isinstance(blade_l_st, BladeSt)

        [bet_u_st, bet_l_st] = [None, None]

        # # Non-iterative claculation of coax rotor
        # ##################################################
        # blade_u_st.coax_pos = 'upper'
        # # Lower on upper rotor inflow
        # blade_u_st.coax_lambda = 0 * blade_u_st.y_arr
        # bet_u_st = bet_forces(blade_u_st)
        # bet_u_st = bet_forces_add_total(bet_u_st, False)
        # ##################################################
        # blade_l_st.coax_pos = 'lower'
        # # Upper on lower rotor inflow
        # blade_l_st.coax_lambda = 0 * (bet_u_st.total_lcr + bet_u_st.total_lir)
        # bet_l_st = bet_forces(blade_l_st)
        # bet_l_st = bet_forces_add_total(bet_l_st, False)
        # ##################################################
        # return

        # Iterative claculation of coax rotor
        veltip_u = blade_u_st.R * blade_u_st.omega
        veltip_l = blade_l_st.R * blade_l_st.omega
        # nsections_u = blade_u_st.nsections
        nsections_l = blade_l_st.nsections
        # lambda_l = 0.0 * np.ones((len(blade_l_st.y_arr)-1, 1))
        downwash_l = 0.0 * np.ones(nsections_l)
        thrust_u_prev = 0
        thrust_l_prev = 0

        niter = 100
        for i in range(0, niter):   # = 1:niter
            # Simulate coax rotor
            ##################################################
            blade_u_st.coax_pos = 'upper'
            # Lower on upper rotor inflow
            # Lower on upper interference
            # blade_u_st.coax_lambda = lambda_l
            blade_u_st.coax_downwash = downwash_l
            bet_u_st = Bet.bet_forces(blade_u_st)
            bet_u_st = Bet.bet_forces_add_total(bet_u_st, False)
            lambda_u = bet_u_st.total_lr
            downwash_u = lambda_u * veltip_u
            ##################################################
            blade_l_st.coax_pos = 'lower'
            # Upper on lower rotor inflow
            # Upper on lower interference
            # blade_l_st.coax_lambda = lambda_u
            blade_l_st.coax_downwash = downwash_u
            bet_l_st = Bet.bet_forces(blade_l_st)
            bet_l_st = Bet.bet_forces_add_total(bet_l_st, False)
            lambda_l = bet_l_st.total_lr
            downwash_l = lambda_l * veltip_l
            ##################################################

            # Check for terminating condition
            delta_thrust_u = abs(bet_u_st.total_T - thrust_u_prev)
            delta_thrust_l = abs(bet_l_st.total_T - thrust_l_prev)
            delta_thrust = delta_thrust_u + delta_thrust_l
            thrust_u_prev = bet_u_st.total_T
            thrust_l_prev = bet_l_st.total_T

            # print('[bet_coax_forces] i %d, delta_T %.4f N, lambda_u %.4f, '
            #       'lambda_l %.4f' %
            #       (i, delta_T, bet_u_st.total_lambda, bet_l_st.total_lambda))
            if delta_thrust < 10**(-2):
                # arg = '[bet_coax_forces] i ' + str(i) + ' delta ' + str(delta)
                # disp(arg)
                # print('[bet_coax_forces] i %d, delta %.4f N, lambda_u %.4f, '
                #       'lambda_l %.4f' %
                #       (i, delta, bet_u_st.total_lambda,
                #        bet_l_st.total_lambda))
                break

            if i >= (niter-1):
                raise RuntimeError('Too many iterations without convergence')

        # if i < (2)
        #     raise RuntimeError('Too few iterations for convergence')
        # end
        return [bet_u_st, bet_l_st]

    @staticmethod
    def bet_coax_eta_thrust(coax_thrust, rotortype, coaxu_r0, coaxl_r0):
        # function 
        # [
        #     db            
        # ] = bet_coax_eta_thrust(
        #     coax_thrust ,  
        #     rotortype   ,   
        #     coaxu_r0    ,  
        #     coaxl_r0      
        # )

        print('######################################################')
        print('Given:')
        print('  constant total thrust')
        print('  theta0_u = collpitch_u')
        print('  theta0_l = collpitch_l')
        print('  T_u = W*(eta_T)/(1+eta_T)')
        print('  T_l = W - T_u')
        # print('  eta_T = eta_T_arr[i]')
        print('Find:')
        print('  omega_u and omega_l')     

        rpm2rads = np.pi / 30
        # rads2rpm = 30 / np.pi

        [rho, lambda_c, mu, collpitch] = \
            KdeGeometry.kde_rotor_defaults(rotortype)
        omega = 3000 * rpm2rads
        CT_target = np.nan      # CT_target is only needed for mu != 0

        collpitch0 = collpitch

        # filename = ['img/bet_coax_eta_T_arr_' rotortype '.mat']
        # filename = 'img/bet_coax_eta_thrust' + '_' + \
        #            rotortype + '_' + str(coax_thrust) + '_' + \
        #            str(coaxu_r0) + '_' + str(coaxl_r0) + '.mat'

        # Hover case for KDECF305DP
        # omega = 3890 * rpm2rads
        # db.blade_h_st = blade_model(
        #     rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)

        db = Dcollpitch()

        # Hover case for KDECF245DP
        db.blade_h_st = BladeSt.blade_model(
            rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
        T_target = coax_thrust/2
        db.blade_h_st = BetFindOmega.bet_find_omega(db.blade_h_st, T_target)

        db.bet_h_st = Bet.bet_forces(db.blade_h_st)
        db.bet_h_st = Bet.bet_forces_add_total(db.bet_h_st, False)
        db.T_h = db.bet_h_st.total_T 
        db.Q_h = db.bet_h_st.total_Q
        db.P_h = np.sqrt(2) * 2 * db.bet_h_st.total_P
        db.W = 2 * db.T_h
        db.P_0 = db.W / (2 * db.blade_h_st.rho * db.blade_h_st.rotArea)

        db.coaxu_r0 = coaxu_r0
        db.coaxl_r0 = coaxl_r0

        # db.dcollpitch_arr = -1.0:0.5:+2.0
        db.dcollpitch_arr = np.array([0])
        # db.eta_T_arr = 0.5:0.5:4
        db.eta_T_arr = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
        # db.blade_h_st
        # db
        num_i = len(db.dcollpitch_arr)
        num_j = len(db.eta_T_arr)

        db.collpitch_u_arr = np.zeros(num_i)
        db.collpitch_l_arr = np.zeros(num_i)
        db.coax.bet_u_st = np.zeros((num_i, num_j), dtype=BetResults)
        db.coax.bet_l_st = np.zeros((num_i, num_j), dtype=BetResults)
        db.sbs.bet_u_st = np.zeros((num_i, num_j), dtype=BetResults)
        db.sbs.bet_l_st = np.zeros((num_i, num_j), dtype=BetResults)
        for i in range(0, num_i):  # = 1:length(db.dcollpitch_arr)
            dcollpitch = db.dcollpitch_arr[i]
            collpitch_u = UnitCnv.deg2rad( collpitch0 - dcollpitch )
            collpitch_l = UnitCnv.deg2rad( collpitch0 + dcollpitch )

            db.collpitch_u_arr[i] = collpitch_u
            db.collpitch_l_arr[i] = collpitch_l

            # Get blade model
            blade_u_st = BladeSt.blade_model(
                rotortype, rho, lambda_c, mu, collpitch_u, omega, CT_target)
            blade_l_st = BladeSt.blade_model(
                rotortype, rho, lambda_c, mu, collpitch_l, omega, CT_target)

            # Set blades for coax and sbs cases
            bladecoax_u_st = blade_u_st
            bladecoax_l_st = blade_l_st 
            bladesbs_u_st = blade_u_st
            bladesbs_l_st = blade_l_st

            # Add coaxial r0 parameters to both coax blades
            bladecoax_u_st.coaxu_r0 = coaxu_r0
            bladecoax_u_st.coaxl_r0 = coaxl_r0
            bladecoax_l_st.coaxu_r0 = coaxu_r0
            bladecoax_l_st.coaxl_r0 = coaxl_r0

            for j in range(0, num_j):   # = 1:length(db.eta_T_arr)
                # Calculate desired upper and lower thrust
                eta_T = db.eta_T_arr[j]
                T_u_target = db.W * (eta_T) / (1 + eta_T)
                T_l_target = db.W - T_u_target  

                # Find RPM that matches desired upper and lower thrust
                [bladecoax_u_st, bladecoax_l_st] = \
                    BetCoaxMatchThrust.bet_coax_match_thrust(
                        bladecoax_u_st, bladecoax_l_st, T_u_target, T_l_target)
                # [bladecoax_u_st, bladecoax_l_st] = \
                #     BetCoaxMatchThrust.bet_coax_find_omega(
                #         bladecoax_u_st, bladecoax_l_st, T_u_target, T_l_target)
                # omegacoax_u = bladecoax_u_st.omega
                # omegacoax_l = bladecoax_l_st.omega

                # Calculate coaxial loads
                [bet_u_st, bet_l_st] = BetMainCoaxial.bet_coax_forces(
                    bladecoax_u_st, bladecoax_l_st)

                # https://realpython.com/python-pass-by-reference/
                # https://www.geeksforgeeks.org/
                #   is-python-call-by-reference-or-call-by-value/
                # Making a deep copy here is EXTREMELY important, because
                # otherwise the values stored in db are update with the
                # new calculations in the loop
                db.coax.bet_u_st[i, j] = copy.deepcopy(bet_u_st)
                db.coax.bet_l_st[i, j] = copy.deepcopy(bet_l_st)

                ##################################################
                # Find RPM that matches desired upper and lower thrust
                [bladesbs_u_st, bladesbs_l_st] = \
                    BetCoaxMatchThrust.bet_sbs_match_thrust(
                        bladesbs_u_st, bladesbs_l_st, T_u_target, T_l_target)
                # omegasbs_u = bladesbs_u_st.omega
                # omegasbs_l = bladesbs_l_st.omega

                # Calculate side-by-side loads
                [bet_u_st, bet_l_st] = BetMainCoaxial.bet_sbs_forces(
                    bladesbs_u_st, bladesbs_l_st)

                # https://realpython.com/python-pass-by-reference/
                # https://www.geeksforgeeks.org/
                #   is-python-call-by-reference-or-call-by-value/
                # Making a deep copy here is EXTREMELY important, because
                # otherwise the values stored in db are update with the
                # new calculations in the loop
                db.sbs.bet_u_st[i, j] = copy.deepcopy(bet_u_st)
                db.sbs.bet_l_st[i, j] = copy.deepcopy(bet_l_st)

        # for j in range(0, num_j):
        #     print('[bet_coax_eta_thrust] j %s' % j)
        #     omega_u = db.coax.bet_u_st[0, j].blade_st.omega
        #     omega_l = db.coax.bet_l_st[0, j].blade_st.omega
        #     print('[bet_coax_eta_thrust] omega_u %s' % omega_u)
        #     print('[bet_coax_eta_thrust] omega_l %s' % omega_l)
        #     thrust_u = db.coax.bet_u_st[0, j].total_T
        #     thrust_l = db.coax.bet_l_st[0, j].total_T
        #     print('[bet_coax_eta_thrust] thrust_u %s' % thrust_u)
        #     print('[bet_coax_eta_thrust] thrust_l %s' % thrust_l)
        return db

    @staticmethod
    def bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0):
        # function bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0)

        read_from_file = False
        # read_from_file = True

        # rotortype = 'KDECF245DP'               
        # rotortype = 'KDECF245TP'     
        # rotortype = 'KDECF245HP'  
        # rotortype = 'KDECF305DP'

        coax_thrust = None

        # thrust_lvl
        if rotortype == 'KDECF245DP':
            coax_thrust = 35
        if rotortype == 'KDECF245TP':
            coax_thrust = 40

        filename = 'img/bet_coax_eta_thrust' + '_' + \
                   rotortype + '_' + str(coax_thrust) + '_' + \
                   str(coaxu_r0) + '_' + str(coaxl_r0) + '.mat'

        if read_from_file:
            db = BetMainCoaxial.load_from_mat(filename)
        else:
            db = BetMainCoaxial.bet_coax_eta_thrust(
                coax_thrust, rotortype, coaxu_r0, coaxl_r0)
            BetMainCoaxial.save_to_mat(filename, db)

        dcollpitch = 0
        str1 = 'bet_coax_eta_thrust'
        Dcollpitch.bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1)

        # thrust_lvl
        if rotortype == 'KDECF245DP':
            coax_thrust = 45
        if rotortype == 'KDECF245TP':
            coax_thrust = 50

        filename = 'img/bet_coax_eta_thrust' + '_' + \
                   rotortype + '_' + str(coax_thrust) + '_' + \
                   str(coaxu_r0) + '_' + str(coaxl_r0) + '.mat'

        if read_from_file:
            db = BetMainCoaxial.load_from_mat(filename)
        else:
            db = BetMainCoaxial.bet_coax_eta_thrust(
                coax_thrust, rotortype, coaxu_r0, coaxl_r0)
            BetMainCoaxial.save_to_mat(filename, db)

        dcollpitch = 0
        str1 = 'bet_coax_eta_thrust'
        Dcollpitch.bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1)

        # thrust_lvl
        if rotortype == 'KDECF245DP':
            coax_thrust = 55
        if rotortype == 'KDECF245TP':
            coax_thrust = 60

        filename = 'img/bet_coax_eta_thrust' + '_' + \
                   rotortype + '_' + str(coax_thrust) + '_' + \
                   str(coaxu_r0) + '_' + str(coaxl_r0) + '.mat'

        if read_from_file:
            db = BetMainCoaxial.load_from_mat(filename)
        else:
            db = BetMainCoaxial.bet_coax_eta_thrust(
                coax_thrust, rotortype, coaxu_r0, coaxl_r0)
            BetMainCoaxial.save_to_mat(filename, db)

        dcollpitch = 0
        str1 = 'bet_coax_eta_thrust'
        Dcollpitch.bet_plot_coax_dcollpitch(rotortype, dcollpitch, db, str1)

    @staticmethod
    def bet_main_coaxial_batch():
        # Radius of lower rotor wake at upper rotor location
        coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0]
        # Radius of upper rotor wake at lower rotor location
        coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0]
        # Upper and lower rotor type 
        rotortype = 'KDECF245DP'
    
        for coaxu_r0 in coaxu_r0_arr:
            for coaxl_r0 in coaxl_r0_arr:
                print('[bet_main_coaxial_batch] coaxu_r0 %s' % coaxu_r0)
                print('[bet_main_coaxial_batch] coaxl_r0 %s' % coaxl_r0)
                # Calculate results
                BetMainCoaxial.bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0)
    
        # Radius of lower rotor wake at upper rotor location
        coaxu_r0_arr = [2.0, 3.0, 4.0, 5.0]
        # Radius of upper rotor wake at lower rotor location
        coaxl_r0_arr = [0.7, 0.8, 0.9, 1.0]
        # Upper and lower rotor type 
        rotortype = 'KDECF245TP'
    
        for coaxu_r0 in coaxu_r0_arr:
            for coaxl_r0 in coaxl_r0_arr:
                print('[bet_main_coaxial_batch] coaxu_r0 %s' % coaxu_r0)
                print('[bet_main_coaxial_batch] coaxl_r0 %s' % coaxl_r0)
                # Calculate results
                BetMainCoaxial.bet_main_coaxial(rotortype, coaxu_r0, coaxl_r0)


class BetCoaxMatchThrust:
    funzero_coax_dict = None

    @staticmethod
    def bet_sbs_match_thrust(
            bladesbs_u_st, bladesbs_l_st, T_u_target, T_l_target):
        return [bladesbs_u_st, bladesbs_l_st]

    @staticmethod
    def fsolve_funzero_coax(
            xguess, blade_u_st, blade_l_st, T_u_target, T_l_target):
        # funzero_coax(x, blade_u_st, blade_l_st, T_u_target, T_l_target)
        funzero_coax = BetCoaxMatchThrust.funzero_coax
        # bet_coax_match_thrust_fzero(
        #     omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)
        # funzero = BetCoaxMatchThrust.bet_coax_match_thrust_fzero

        x = xguess
        niter = 10
        for i in range(0, niter):
            err = funzero_coax(
                x, blade_u_st, blade_l_st, T_u_target, T_l_target)
            # err_u = T_u - T_u_target
            # err_l = T_l - T_l_target
            # err = [err_u, err_l]
            omega_u = x[0]
            omega_l = x[1]
            rpm_u = UnitCnv.rads2rpm(omega_u)
            rpm_l = UnitCnv.rads2rpm(omega_l)
            err_u = err[0]
            err_l = err[1]
            t1 = (i, rpm_u, T_u_target, err_u)
            t2 = (i, rpm_l, T_l_target, err_l)
            # print('m1: rpm_u %s, T_u %s, T_u_target %s, err_u %s' % t1)
            # print('m2: rpm_l %s, T_l %s, T_u_target %s, err_u %s' % t2)
            print('m1: i %s, rpm_u %s, T_u_target %s, err_u %s' % t1)
            print('m2: i %s, rpm_l %s, T_u_target %s, err_u %s' % t2)
            if np.linalg.norm(err) < 1:
                ier = 1
                break
            else:
                kp = UnitCnv.rpm2rads(10)
                x = x + np.array([-kp * err_u, -kp * err_l])

        x0 = x
        fval = err
        exitflag = ier
        output = ""
        return [x0, fval, exitflag, output]

    @staticmethod
    def bet_coax_match_thrust_fzero(
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function [err_u, err_l] = bet_coax_match_thrust_fzero(
        #     omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)

        omega_u = np.sqrt(omega_u * omega_u)
        omega_l = np.sqrt(omega_l * omega_l)

        # Apply changes in rotor speed
        blade_u_st.omega = omega_u
        # blade_u_st.Vtip = blade_u_st.omega * blade_u_st.R
        blade_l_st.omega = omega_l
        # blade_l_st.Vtip = blade_l_st.omega * blade_l_st.R

        # Calculate BET loads
        [bet_u_st, bet_l_st] = BetMainCoaxial.bet_coax_forces(
            blade_u_st, blade_l_st)
        T_u = bet_u_st.total_T
        Q_u = bet_u_st.total_Q
        lambda_u = bet_u_st.total_lambda_t
        T_l = bet_l_st.total_T
        Q_l = bet_l_st.total_Q
        lambda_l = bet_l_st.total_lambda_t

        err_u = T_u - T_u_target
        err_l = T_l - T_l_target
        # err = [err_u, err_l]
        return [err_u, err_l]

    @staticmethod
    def funzero_coax(x, blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function err = funzero_coax(x)
        omega_u = x[0]
        omega_l = x[1]

        [err_u, err_l] = BetCoaxMatchThrust.bet_coax_match_thrust_fzero(
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)
        err = [err_u, err_l]
        return err

    @staticmethod
    def funzero_upper(x, blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function err = funzero_upper(x)
        omega_u = x[0]
        omega_l = 0     # x(2)

        [err_u, err_l] = BetCoaxMatchThrust.bet_coax_match_thrust_fzero(
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)
        err = [err_u]
        return err

    @staticmethod
    def funzero_lower(x, blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function err = funzero_lower(x)
        omega_u = 0     # x(1)
        omega_l = x[1]  # x(2)

        [err_u, err_l] = BetCoaxMatchThrust.bet_coax_match_thrust_fzero(
            omega_u, omega_l, blade_u_st, blade_l_st, T_u_target, T_l_target)
        err = [err_l]
        return err

    @staticmethod
    def bet_coax_match_thrust(blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function [blade_u_st, blade_l_st] = bet_coax_match_thrust(
        #     blade_u_st, blade_l_st, T_u_target, T_l_target)

        # rpm2rads = np.pi / 30
        rads2rpm = 30 / np.pi

        print('[bet_coax_match_thrust] Running fsolve ')
        x0 = None

        # Special case, zero upper thrust => zero upper rpm
        if T_u_target == 0:
            xguess = [blade_l_st.omega]
            x0 = Mat.fsolve(BetCoaxMatchThrust.funzero_lower, xguess,
                            blade_u_st, blade_l_st, T_u_target, T_l_target)
            x0 = [0, x0]

        # Special case, zero lower thrust => zero lower rpm
        if T_l_target == 0:
            xguess = [blade_u_st.omega]
            x0 = Mat.fsolve(BetCoaxMatchThrust.funzero_upper, xguess,
                            blade_u_st, blade_l_st, T_u_target, T_l_target)
            x0 = [x0, 0]

        # General case
        if T_u_target != 0 and T_l_target != 0:
            xguess = [blade_u_st.omega, blade_l_st.omega]
            [x0, fval, exitflag, output] = \
                Mat.fsolve(BetCoaxMatchThrust.funzero_coax, xguess,
                           blade_u_st, blade_l_st, T_u_target, T_l_target)
            x0 = [x0[0], x0[1]]

        # Calculate error
        err = BetCoaxMatchThrust.funzero_coax(
            x0, blade_u_st, blade_l_st, T_u_target, T_l_target)
        err_u = err[0]
        err_l = err[1]

        # Check error magnitude
        if np.linalg.norm(err) > 10**(-4):
            print(err)
            print(blade_u_st)
            print(blade_l_st)
            raise RuntimeError('norm(err) > 10**-4')

        # Apply changes
        blade_u_st.omega = float(x0[0])
        # blade_u_st.Vtip = blade_u_st.omega * blade_u_st.R
        blade_l_st.omega = float(x0[1])
        # blade_l_st.Vtip = blade_l_st.omega * blade_l_st.R

        # Calculate BET loads
        [bet_u_st, bet_l_st] = BetMainCoaxial.bet_coax_forces(
            blade_u_st, blade_l_st)
        T_u = bet_u_st.total_T
        Q_u = bet_u_st.total_Q
        lambda_u = bet_u_st.total_lambda_t
        T_l = bet_l_st.total_T
        Q_l = bet_l_st.total_Q
        lambda_l = bet_l_st.total_lambda_t

        print('[bet_coax_match_thrust] omega_u %.4f, omega_l %.4f' %
              (blade_u_st.omega * rads2rpm, blade_l_st.omega * rads2rpm))
        print('[bet_coax_match_thrust] lambda_u %.4f, lambda_l %.4f' %
              (lambda_u, lambda_l))
        print('[bet_coax_match_thrust] T_u %.4f, T_l %.4f' %
              (T_u, T_l))
        print('[bet_coax_match_thrust] err_u %.4f, err_l %.4f' %
              (err_u, err_l))

        return [blade_u_st, blade_l_st]

    @staticmethod
    def bet_coax_find_omega(blade_u_st, blade_l_st, T_u_target, T_l_target):
        # function [blade_u_st, blade_l_st] = bet_coax_match_thrust(
        #     blade_u_st, blade_l_st, T_u_target, T_l_target)

        # rpm2rads = np.pi / 30
        rads2rpm = 30 / np.pi
        print('[bet_coax_find_omega] Running fsolve ')

        xguess = [blade_u_st.omega, blade_l_st.omega]
        [x0, fval, exitflag, output] = \
            BetCoaxMatchThrust.fsolve_funzero_coax(
                xguess, blade_u_st, blade_l_st, T_u_target, T_l_target)

        # Calculate error
        err = BetCoaxMatchThrust.funzero_coax(
            x0, blade_u_st, blade_l_st, T_u_target, T_l_target)
        err_u = err[0]
        err_l = err[1]

        # Check error magnitude
        if np.linalg.norm(err) > 10 ** (-4):
            print(err)
            print(blade_u_st)
            print(blade_l_st)
            raise RuntimeError('norm(err) > 10**-4')

        # Apply changes
        blade_u_st.omega = float(x0[0])
        # blade_u_st.Vtip = blade_u_st.omega * blade_u_st.R
        blade_l_st.omega = float(x0[1])
        # blade_l_st.Vtip = blade_l_st.omega * blade_l_st.R

        # Calculate BET loads
        [bet_u_st, bet_l_st] = BetMainCoaxial.bet_coax_forces(
            blade_u_st, blade_l_st)
        T_u = bet_u_st.total_T
        Q_u = bet_u_st.total_Q
        lambda_u = bet_u_st.total_lambda_t
        T_l = bet_l_st.total_T
        Q_l = bet_l_st.total_Q
        lambda_l = bet_l_st.total_lambda_t

        print('[bet_coax_find_omega] omega_u %.4f, omega_l %.4f' %
              (blade_u_st.omega * rads2rpm, blade_l_st.omega * rads2rpm))
        print('[bet_coax_find_omega] lambda_u %.4f, lambda_l %.4f' %
              (lambda_u, lambda_l))
        print('[bet_coax_find_omega] T_u %.4f, T_l %.4f' %
              (T_u, T_l))
        print('[bet_coax_find_omega] err_u %.4f, err_l %.4f' %
              (err_u, err_l))

        return [blade_u_st, blade_l_st]


if __name__ == "__main__":
    ufilename = 'img/' + BetMainCoaxial.bet_main_coaxial_batch.__name__ + '.txt'
    GnrlConfig.write_config(ufilename)
    BetMainCoaxial.bet_main_coaxial_batch()
