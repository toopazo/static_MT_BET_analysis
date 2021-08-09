from toopazo_tools.file_folder import FileFolderTools as FFTool
# from unit_conversion import UnitConversion as UnitCnv


class GnrlConfig:
    # bpath_bet_results = '/home/tzo4/Dropbox/tomas/pennState/avia/software/' \
    #                      'static_MT_BET_coaxRotor/img/' \
    #                      'bet_cf245_pitch'

    bpath_bet_results = '/home/tzo4/Dropbox/tomas/pennState/avia/software/' \
                         'static_MT_BET_analysis/python_code/img'

    tunning_for_kde245 = {
        # 'modify_collective': -1.7,            # img_collective
        # 'modify_ClCdCm': [1, 1, 1]            # img_collective
        # 'modify_collective': 0,               # img_coeffs_v1
        # 'modify_ClCdCm': [0.8, 1.5, 1.0]      # img_coeffs_v1
        'modify_collective': -2.0,
        'modify_ClCdCm': [1.0, 2.0, 1.0]
    }

    @staticmethod
    def write_config(filename):
        mdict = GnrlConfig.tunning_for_kde245
        fd = open(filename, 'w')
        fd.write(str(mdict))
        fd.close()


    @staticmethod
    def bet_filename(nblades, thrust_lvl, r0_upper, r0_lower):
        bpath_bet = GnrlConfig.bpath_bet_results

        # Parse file
        betfile_v1 = ''
        betfile_v2 = ''
        if nblades == 2 and thrust_lvl in [35, 45, 55]:
            betfile_v1 = 'bet_coax_eta_thrust_KDECF245DP_%i_%.01f_%.01f.txt' % \
                         (thrust_lvl, r0_upper, r0_lower)
            betfile_v2 = 'bet_coax_eta_thrust_KDECF245DP_%i_%i_%.01f.txt' % \
                         (thrust_lvl, r0_upper, r0_lower)

        if nblades == 3 and thrust_lvl in [40, 50, 60]:
            betfile_v1 = 'bet_coax_eta_thrust_KDECF245TP_%i_%.01f_%.01f.txt' % \
                         (thrust_lvl, r0_upper, r0_lower)
            betfile_v2 = 'bet_coax_eta_thrust_KDECF245TP_%i_%i_%.01f.txt' % \
                         (thrust_lvl, r0_upper, r0_lower)

        # Check for file in both notations
        # betfile_v1 = bet_coax_eta_thrust_KDECF245DP_35_2.0_0.7.txt
        # betfile_v2 = bet_coax_eta_thrust_KDECF245DP_35_2_0.7.txt
        filename_v1 = bpath_bet + '/%s' % betfile_v1
        filename_v2 = bpath_bet + '/%s' % betfile_v2
        if FFTool.is_file(filename_v1):
            filename = filename_v1
        elif FFTool.is_file(filename_v2):
            filename = filename_v2
        else:
            print(filename_v1)
            print(filename_v2)
            raise RuntimeError('No such file')
        return filename
