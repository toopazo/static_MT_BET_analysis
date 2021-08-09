import numpy as np
import matplotlib as plt
from unit_conversion import UnitConversion as uc
from bet_loads import Bet
from rotor_model import BladeSt, KdeGeometry
from matmath import MatMath as Mat
from anecha_data import AnechaData
from toopazo_tools.matplotlib import FigureTools, PlotTools
from gnrl_config import GnrlConfig


class BetVsAnecha:
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
        self.omegarpm_arr = None
        self.ane_T_arr =  None
        self.ane_Q_arr =  None
        self.bet_T_arr =  None
        self.bet_Q_arr =  None
        self.err_T_arr =  None
        self.err_Q_arr =  None
        self.mean_err_T =  None
        self.mean_err_Q =  None
        self.std_err_T =  None
        self.std_err_Q =  None
        self.perr_T_arr =  None
        self.perr_Q_arr =  None
        # self.str1 =  None
        # self.str2 =  None

        self._freeze()  # no new attributes after this point.

    def bet_main_omega_plot(self, str1, str2):
        # function bet_main_omega_plot(self, str1, str2)
    
        # Loadcel precision
        # According to https://www.ati-ia.com/products/ft/ft_models.aspx?id=Delta
        # and model US-150-600
        # https://www.convertunits.com/from/lbf/to/N
        # the max error in thrust is 1/16 = 0.0625 lbf = 0.2780 N
        # http://convert-to.com/conversion/torque/convert-lbf-in-to-n-m.html
        # the max error in torque is 1/16 = 0.0625 lbf-in =  0.01 Nm
    
        # Experimentally determined error are about
        # Error in thrust is 1.00Nm
        error_thrust_arr = 1.0 * np.ones(np.size(self.ane_T_arr))
        # Error in torque is 0.03Nm
        error_torque_arr = 0.03 * np.ones(np.size(self.ane_Q_arr))
        # Error in RPM is 30 RPM
        error_rpm_arr = 30 * np.ones(np.size(self.omegarpm_arr))
        # Error in power is 0.03Nm * 30rpm = 0.03 * 3.14 rad/s = 0.1 W
        error_power_arr = 0.1
    
        plot_error = True
    
        # nfig = 1
        # fig = figure(nfig)
        # subplot(2, 1, 1)
        # hold on
        # grid on
        # set(gca,'FontSize', 14)
        # # plot(omegarpm_arr, kde_T_arr, '-*')
        # # plot(self.omegarpm_arr, self.ane_T_arr, '-*')
        # plot(self.omegarpm_arr, self.bet_T_arr, '-b')
        # errorbar(self.omegarpm_arr, self.ane_T_arr, error_thrust_arr, '-r')
        # errorbar(self.omegarpm_arr, self.ane_T_arr, error_rpm_arr, '-r', 'horizontal')
        # xlabel('Rotor speed RPM')
        # ylabel('Thrust N')
        # xlim([2000, 4500])
        # ylim([0, 80])
    
        # if plot_error:
        #     yyaxis right
        #     plot(self.omegarpm_arr, self.perr_T_arr, '-k')
        #     # leg('KDE', 'AneCha', 'BET', 'error', 'Location', 'northwest')
        #     # leg('AneCha', 'BET', 'error', 'Location', 'northwest')
        #     ylabel('Error %')
        #     ax = gca
        #     ax.YAxis(1).Color = 'k'
        #     ax.YAxis(2).Color = 'k'
        #     ylim([-10, 10])
        #
        # #    ind = find(omegarpm_arr==3890)
        # #    omegaind = omegarpm_arr(ind)
        # #    val = round(kde_T_arr(ind), 1)
        # #    text(2600, 17, ['KDE: ' num2str(omegaind) ' RPM, ' num2str(val) ' N'])
        # #    val = round(bet_T_arr(ind), 1)
        # #    text(2600, 12, ['BET: ' num2str(omegaind) ' RPM, ' num2str(val) ' N'])
        #     # 3890 rpm 126.5090 N 4.7044 Nm
        #     text(2800, -5, ['std dev error: ' num2str(self.std_err_T) ' N'])
        #     text(2800, -8, ['   mean error: ' num2str(self.mean_err_T) ' N'])
        
    
        # subplot(2, 1, 2)
        # hold on
        # grid on
        # set(gca,'FontSize', 14)
        # # plot(omegarpm_arr, kde_Q_arr, '-*')
        # # plot(self.omegarpm_arr, self.ane_Q_arr, '-*')
        # plot(self.omegarpm_arr, self.bet_Q_arr, '-b')
        # errorbar(self.omegarpm_arr, self.ane_Q_arr, error_torque_arr, '-r')
        # errorbar(self.omegarpm_arr, self.ane_Q_arr, error_rpm_arr, '-r', 'horizontal')
        # xlabel('Rotor speed RPM')
        # ylabel('Torque Nm')
        # xlim([2000, 4500])
        # ylim([0, 3])
        #
        # if plot_error:
        #     yyaxis right
        #     plot(self.omegarpm_arr, self.perr_Q_arr, '-k')
        #     # leg('KDE', 'AneCha', 'BET', 'error', 'Location', 'northwest')
        #     # leg('AneCha', 'BET', 'error', 'Location', 'northwest')
        #     leg('BET', 'Exp. uncertainty', 'Exp. uncertainty', 'Perc. error', 'Location', 'northwest')
        #     ylabel('Error %')
        #     ax = gca
        #     ax.YAxis(1).Color = 'k'
        #     ax.YAxis(2).Color = 'k'
        #     ylim([-10, 10])
        #
        # #    ind = find(omegarpm_arr==3890)
        # #    omegaind = omegarpm_arr(ind)
        # #    val = round(kde_Q_arr(ind), 1)
        # #    text(2600, 17, ['KDE: ' num2str(omegaind) ' RPM, ' num2str(val) ' Nm'])
        # #    val = round(bet_Q_arr(ind), 1)
        # #    text(2600, 12, ['BET: ' num2str(omegaind) ' RPM, ' num2str(val) ' Nm'])
        #     # text(2800, 8, ['std dev error: ' num2str(self.std_err_Q) ' Nm'])
        #     # text(2800, 5, ['   mean error: ' num2str(self.mean_err_Q) ' Nm'])
        #     text(3300, -2, ['std dev error: ' num2str(self.std_err_Q) ' Nm'])
        #     text(3300, -5, ['  mean error: ' num2str(self.mean_err_Q) ' Nm'])
        # else:
        #     legend('BET', 'Exp. uncertainty', 'Exp. uncertainty', 'Location', 'northwest')

        # fig = figure(nfig)
        # filename = ['img/' str1 '_' str2 '.jpg']
        # saveas(fig, filename)
        #
        # close
        # all

        [fig, ax_arr] = FigureTools.create_fig_axes(2, 1)

        x_arr = [self.omegarpm_arr, self.omegarpm_arr]
        y_arr = [self.ane_T_arr, self.bet_T_arr]
        xlabel_arr = ['']
        ylabel_arr = ['Thrust N']
        PlotTools.ax1_x2_y2([ax_arr[0]], x_arr, xlabel_arr, y_arr, ylabel_arr)
        ax_arr[0].set_xlim([2000, 4500])
        ax_arr[0].set_ylim([0, 80])
        ax_arr[0].set_xticks([2000, 2500, 3000, 3500, 4000, 4500])
        ax_arr[0].set_yticks([0, 20, 40, 60, 80])
        ax_arr[0].legend(['ane', 'bet'])

        x_arr = [self.omegarpm_arr, self.omegarpm_arr]
        y_arr = [self.ane_Q_arr, self.bet_Q_arr]
        xlabel_arr = ['Rotor speed RPM']
        ylabel_arr = ['Torque Nm']
        PlotTools.ax1_x2_y2([ax_arr[1]], x_arr, xlabel_arr, y_arr, ylabel_arr)
        ax_arr[1].set_xlim([2000, 4500])
        ax_arr[1].set_ylim([0, 3])
        ax_arr[1].set_xticks([2000, 2500, 3000, 3500, 4000, 4500])
        ax_arr[1].set_yticks([0, 1, 2, 3])

        filename = 'img/' + str1 + '_' + str2 + '.jpg'
        FigureTools.savefig(filename=filename, closefig=True)


class BetMainOmega:
    @staticmethod
    def bet_main_omega():
        # Add functions to the path
        # addpath('mt_bet')
        # addpath('mt_bet_data')
        # addpath('mt_bet_coaxial_rotor')
        # addpath('mt_bet_plot')

        # rpm2rads = UnitConversion.rpm2rads
        # rads2rpm = UnitConversion.rads2rpm

        rotortype = 'KDECF245DP'
        # rotortype = 'KDECF245TP'
        # rotortype = 'KDECF245HP'

        # omegarpm_arr = [2000, 2500, 3000, 3500, 3890, 4000, 4038, 4500]
        # omegarpm_arr = [2000, 2500, 3000, 3500, 4000, 4500]
        omegarpm_arr = np.linspace(2000, 4500, 10)
        npoints = len(omegarpm_arr)

        # T and Q vs RPM for anecha and BET
        ################################################
        ane_T_arr = np.zeros(npoints)
        ane_Q_arr = np.zeros(npoints)
        bet_T_arr = np.zeros(npoints)
        bet_Q_arr = np.zeros(npoints)
        rhoAvtip2_arr = np.zeros(npoints)
        rhoAvtip3_arr = np.zeros(npoints)
        for i in range(0, npoints):
            omega = uc.rpm2rads(omegarpm_arr[i])
            print('[bet_main_omega] omega #d RPM \n', omegarpm_arr[i])

            # rotortype = 'KDE6213XF185_KDECF305DP'
            # rotortype = 'KDE6213XF185_KDECF245DP'
            # rotortype = 'KDE6213XF185_KDECF245TP'
            # modeltype = 'interp1'
            # modeltype = 'quadratic'
            # [kde_T, kde_Q, kde_P] = kde_rotor_TQP(omega, rotortype, modeltype)
            # kde_T_arr[i] = kde_T
            # kde_Q_arr[i] = kde_Q

            # Anecha experimental results
            modeltype = 'quadratic'
            [ane_T, ane_Q, ane_P] = AnechaData.anecha_rotor_TQP(omega,
                                                                rotortype,
                                                                modeltype)
            ane_T_arr[i] = ane_T
            ane_Q_arr[i] = ane_Q

            # BET results
            [rho, lambda_c, mu, collpitch] = \
                KdeGeometry.kde_rotor_defaults(rotortype)
            omega = omega
            CT_target = np.nan  # CT_target is only needed for mu != 0
            blade_st = BladeSt.blade_model(
                rotortype, rho, lambda_c, mu, collpitch, omega, CT_target)
            bet_st = Bet.bet_forces(blade_st)
            bet_st = Bet.bet_forces_add_total(bet_st, False)

            fluid_density = blade_st.rho
            rotor_radius = bet_st.blade_st.R
            rotor_angvel = bet_st.blade_st.omega
            rotor_area = bet_st.blade_st.rotArea
            rotor_veltip = rotor_angvel * rotor_radius

            # Save results
            bet_T_arr[i] = bet_st.total_T
            bet_Q_arr[i] = bet_st.total_Q
            rhoAvtip2_arr[i] = fluid_density * rotor_area * rotor_veltip ** 2
            rhoAvtip3_arr[i] = fluid_density * rotor_area * rotor_veltip ** 3

        # RMS error
        # err_T_arr = bet_T_arr - kde_T_arr
        # err_Q_arr = bet_Q_arr - kde_Q_arr
        err_T_arr = bet_T_arr - ane_T_arr
        err_Q_arr = bet_Q_arr - ane_Q_arr

        mean_err_T = np.mean(err_T_arr)
        mean_err_Q = np.mean(err_Q_arr)
        std_err_T = np.std(err_T_arr)  # sqrt(np.mean( err_T_arr*err_T_arr ) )
        std_err_Q = np.std(err_Q_arr)  # sqrt(np.mean( err_Q_arr*err_Q_arr ) )

        print('stdev err_T %.4f \n', std_err_T)
        print('stdev err_Q %.4f \n', std_err_Q)
        print('mean err_T %.4f \n', mean_err_T)
        print('mean err_Q %.4f \n', mean_err_Q)

        perr_T_arr = err_T_arr / ane_T_arr * 100
        perr_Q_arr = err_Q_arr / ane_Q_arr * 100

        ################################################
        # Plot
        st = BetVsAnecha()
        st.omegarpm_arr = omegarpm_arr
        st.ane_T_arr = ane_T_arr
        st.ane_Q_arr = ane_Q_arr
        st.bet_T_arr = bet_T_arr
        st.bet_Q_arr = bet_Q_arr
        st.err_T_arr = err_T_arr
        st.err_Q_arr = err_Q_arr
        st.mean_err_T = mean_err_T
        st.mean_err_Q = mean_err_Q
        st.std_err_T = std_err_T
        st.std_err_Q = std_err_Q
        st.perr_T_arr = perr_T_arr
        st.perr_Q_arr = perr_Q_arr
        str1 = 'bet_main_omega'
        str2 = rotortype
        BetVsAnecha.bet_main_omega_plot(st, str1, str2)

        # # Save data to file
        # table_thrust = bet_T_arr
        # table_omega = omegarpm_arr * rpm2rads
        # table_power = bet_Q_arr * table_omega
        # table_torque = bet_Q_arr
        # table_rhoAvtip2 = rhoAvtip2_arr
        # table_rhoAvtip3 = rhoAvtip3_arr

        # T = table(
        #     table_thrust,
        #     table_omega,
        #     table_power,
        #     table_torque,
        #     table_rhoAvtip2,
        #     table_rhoAvtip3
        # )
        #
        # str1 = 'bet_main_omega_'
        # str2 = rotortype
        # filename = ['img/' str1 '_' str2 '.txt']
        # writetable(T, filename)


if __name__ == "__main__":
    ufilename = 'img/' + BetMainOmega.bet_main_omega.__name__ + '.txt'
    GnrlConfig.write_config(ufilename)
    BetMainOmega.bet_main_omega()
