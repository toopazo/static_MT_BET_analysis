function mt_main_FM_single_coax()
    clear
    clc
    close all

    % Add functions to the path
    addpath('mt_bet');
    
    density = 1.1849;
    area = 0.304151164328273;
    
    % [operating_point] motornum 111 nblades 2 omegarpm 3086.75
    % [operating_point] power 319.481949 cpower 0.000871 thrust 35.002975 cthrust 0.009601
    single_nblades2_thrust35 = [3086.75, 319.481949, 0.000871, 35.002975, 0.009601];
    % [operating_point] motornum 111 nblades 2 omegarpm 3299.75
    % [operating_point] power 390.287883 cpower 0.000871 thrust 40.000380 cthrust 0.009601
    single_nblades2_thrust40 = [3299.75, 390.287883, 0.000871 40.000380, 0.009601];
    % [operating_point] motornum 111 nblades 2 omegarpm 3500.0
    % [operating_point] power 465.742734 cpower 0.000871 thrust 45.002655 cthrust 0.009601
    single_nblades2_thrust45 = [3500.00, 465.742734, 0.000871, 45.002655, 0.009601];
    % [operating_point] motornum 111 nblades 2 omegarpm 3689.25
    % [operating_point] power 545.451591 cpower 0.000871 thrust 50.000945 cthrust 0.009601
    single_nblades2_thrust50 = [3689.25, 545.451591, 0.000871, 50.000945, 0.009601];
    % [operating_point] motornum 111 nblades 2 omegarpm 3869.25
    % [operating_point] power 629.248703 cpower 0.000871 thrust 54.999105 cthrust 0.009601
    single_nblades2_thrust55 = [3869.25, 629.248703, 0.000871, 54.999105, 0.009601];
    % [operating_point] motornum 111 nblades 2 omegarpm 4041.25
    % [operating_point] power 716.950395 cpower 0.000871 thrust 59.997545 cthrust 0.009601
    single_nblades2_thrust60 = [4041.25, 716.950395, 0.000871, 59.997545, 0.009601];

    % [operating_point] motornum 111 nblades 3 omegarpm 2750.0
    % [operating_point] power 322.521908 cpower 0.001244 thrust 34.997343 cthrust 0.012095
    single_nblades3_thrust35 = [2750.00, 322.521908, 0.001244, 34.997343, 0.012095];
    % [operating_point] motornum 111 nblades 3 omegarpm 2940.0
    % [operating_point] power 394.097001 cpower 0.001244 thrust 40.000401 cthrust 0.012095
    single_nblades3_thrust40 = [2940.00, 394.097001, 0.001244, 40.000401, 0.012095];
    % [operating_point] motornum 111 nblades 3 omegarpm 3118.25
    % [operating_point] power 470.212242 cpower 0.001244 thrust 44.997828 cthrust 0.012095
    single_nblades3_thrust45 = [3118.25, 470.212242, 0.001244, 44.997828, 0.012095];
    % [operating_point] motornum 111 nblades 3 omegarpm 3287.0
    % [operating_point] power 550.757285 cpower 0.001244 thrust 49.999895 cthrust 0.012095
    single_nblades3_thrust50 = [3287.00, 550.757285, 0.001244, 49.999895, 0.012095];
    % [operating_point] motornum 111 nblades 3 omegarpm 3447.5
    % [operating_point] power 635.439142 cpower 0.001244 thrust 55.001969 cthrust 0.012095
    single_nblades3_thrust55 = [3447.50, 550.757285, 0.001244, 49.999895, 0.012095];
    % [operating_point] motornum 111 nblades 3 omegarpm 3600.75
    % [operating_point] power 724.002468 cpower 0.001244 thrust 60.000604 cthrust 0.012095
    single_nblades3_thrust60 = [3600.75, 724.002468, 0.001244, 60.000604, 0.012095];

    % [operating_point] motornum 111 nblades 6 omegarpm 2343.75
    % [operating_point] power 335.464202 cpower 0.002090 thrust 34.999024 cthrust 0.016652
    single_nblades6_thrust35 = [2343.75, 335.464202, 0.002090, 34.999024, 0.016652];
    % [operating_point] motornum 111 nblades 6 omegarpm 2505.5
    % [operating_point] power 409.822265 cpower 0.002090 thrust 39.996518 cthrust 0.016652
    single_nblades6_thrust40 = [2505.50, 409.822265, 0.002090, 39.996518, 0.016652];
    % [operating_point] motornum 111 nblades 6 omegarpm 2657.5
    % [operating_point] power 489.026223 cpower 0.002090 thrust 44.996622 cthrust 0.016652
    single_nblades6_thrust45 = [2657.50, 489.026223, 0.002090, 44.996622, 0.016652];
    % [operating_point] motornum 111 nblades 6 omegarpm 2801.25
    % [operating_point] power 572.753744 cpower 0.002090 thrust 49.996213 cthrust 0.016652
    single_nblades6_thrust50 = [2801.25, 572.753744, 0.002090, 49.996213, 0.016652];
    % [operating_point] motornum 111 nblades 6 omegarpm 2938.0
    % [operating_point] power 660.796452 cpower 0.002090 thrust 54.996741 cthrust 0.016652
    single_nblades6_thrust55 = [2938.00, 660.796452, 0.002090, 54.996741, 0.016652];
    % [operating_point] motornum 111 nblades 6 omegarpm 3068.75
    % [operating_point] power 753.003260 cpower 0.002090 thrust 60.000710 cthrust 0.016652
    single_nblades6_thrust60 = [3068.75, 753.003260, 0.002090, 60.000710, 0.016652];
    
    
    % [operating_point] motornum 12 nblades 2 eta_thrust [0.92114186]
    % [operating_point] m1 power 89.303767 thrust 16.844000
    % [operating_point] m2 power 230.698063 thrust 18.286000
    % [operating_point] coax power 320.001830 thrust 35.130000
    coax_nblades2_thrust35 = [0.92114186, 89.303767, 16.844000, 230.698063, 18.286000, 320.001830, 35.130000];
    % [operating_point] motornum 12 nblades 2 eta_thrust [0.91870058]
    % [operating_point] m1 power 147.385774 thrust 21.380000
    % [operating_point] m2 power 317.624478 thrust 23.272000
    % [operating_point] coax power 465.010253 thrust 44.652000
    coax_nblades2_thrust45 = [0.91870058, 147.385774, 21.380000, 317.624478, 23.272000, 465.010253, 44.652000];
    % [operating_point] motornum 12 nblades 2 eta_thrust [0.90263819]
    % [operating_point] m1 power 191.820923 thrust 25.866000
    % [operating_point] m2 power 425.529950 thrust 28.656000
    % [operating_point] coax power 617.350873 thrust 54.522000
    coax_nblades2_thrust55 = [0.90263819, 191.820923, 25.866000, 425.529950, 28.656000, 617.350873, 54.522000];

    % [operating_point] motornum 12 nblades 3 eta_thrust [1.03171689]
    % [operating_point] m1 power 125.862083 thrust 20.233000
    % [operating_point] m2 power 246.733827 thrust 19.611000
    % [operating_point] coax power 372.595910 thrust 39.844000
    coax_nblades3_thrust40 = [1.03171689, 125.862083, 20.233000, 246.733827, 19.611000, 372.595910, 39.844000];
    % [operating_point] motornum 12 nblades 3 eta_thrust [0.9157594]
    % [operating_point] m1 power 173.893318 thrust 23.720000
    % [operating_point] m2 power 373.061219 thrust 25.902000
    % [operating_point] coax power 546.954536 thrust 49.622000
    coax_nblades3_thrust50 = [0.9157594, 173.893318, 23.720000, 373.061219, 25.902000, 546.954536, 49.622000];
    % [operating_point] motornum 12 nblades 3 eta_thrust [0.90027314]
    % [operating_point] m1 power 228.889262 thrust 28.346000
    % [operating_point] m2 power 473.836782 thrust 31.486000
    % [operating_point] coax power 702.726044 thrust 59.832000
    coax_nblades3_thrust60 = [0.90027314, 228.889262, 28.346000, 473.836782, 31.486000, 702.726044, 59.832000];
    
    single_nblades2_FM_arr(1) = get_FM_from_single(single_nblades2_thrust35);
    single_nblades2_FM_arr(2) = get_FM_from_single(single_nblades2_thrust40);
    single_nblades2_FM_arr(3) = get_FM_from_single(single_nblades2_thrust45);
    single_nblades2_FM_arr(4) = get_FM_from_single(single_nblades2_thrust50);
    single_nblades2_FM_arr(5) = get_FM_from_single(single_nblades2_thrust55);
    single_nblades2_FM_arr(6) = get_FM_from_single(single_nblades2_thrust60);
    
    single_nblades3_FM_arr(1) = get_FM_from_single(single_nblades3_thrust35);
    single_nblades3_FM_arr(2) = get_FM_from_single(single_nblades3_thrust40);
    single_nblades3_FM_arr(3) = get_FM_from_single(single_nblades3_thrust45);
    single_nblades3_FM_arr(4) = get_FM_from_single(single_nblades3_thrust50);
    single_nblades3_FM_arr(5) = get_FM_from_single(single_nblades3_thrust55);
    single_nblades3_FM_arr(6) = get_FM_from_single(single_nblades3_thrust60);
    
    single_nblades6_FM_arr(1) = get_FM_from_single(single_nblades6_thrust35);
    single_nblades6_FM_arr(2) = get_FM_from_single(single_nblades6_thrust40);
    single_nblades6_FM_arr(3) = get_FM_from_single(single_nblades6_thrust45);
    single_nblades6_FM_arr(4) = get_FM_from_single(single_nblades6_thrust50);
    single_nblades6_FM_arr(5) = get_FM_from_single(single_nblades6_thrust55);
    single_nblades6_FM_arr(6) = get_FM_from_single(single_nblades6_thrust60);
    
    coax_nblades2_FM_arr(1) = get_FM_from_coax(coax_nblades2_thrust35);
    coax_nblades2_FM_arr(2) = get_FM_from_coax(coax_nblades2_thrust45);
    coax_nblades2_FM_arr(3) = get_FM_from_coax(coax_nblades2_thrust55);
    
    coax_nblades3_FM_arr(1) = get_FM_from_coax(coax_nblades3_thrust40);
    coax_nblades3_FM_arr(2) = get_FM_from_coax(coax_nblades3_thrust50);
    coax_nblades3_FM_arr(3) = get_FM_from_coax(coax_nblades3_thrust60);
    
    fig = figure(1);
    single_thrust_arr = [35, 40, 45, 50, 55, 60];
    hold on;
    plot(single_thrust_arr, single_nblades2_FM_arr, '-*', 'LineWidth', 2); %, 'Color', 'red')
    plot(single_thrust_arr, single_nblades3_FM_arr, '-*', 'LineWidth', 2); %, 'Color', 'green')
    plot(single_thrust_arr, single_nblades6_FM_arr, '-*', 'LineWidth', 2); %, 'Color', 'blue')
    plot([35, 45, 55], coax_nblades2_FM_arr, '-*', 'LineWidth', 2); %, 'Color', 'black')
    plot([40, 50, 60], coax_nblades3_FM_arr, '-*', 'LineWidth', 2); %, 'Color', 'magenta')
    ylim([0.7, 1])
    xlim([30, 65])
    grid on
    xlabel('Thrust N')
    ylabel('FM')
    legend('N_b=2', 'N_b=3', 'N_b=6', 'N_b=2 coax', 'N_b=3 coax', 'Location', 'NorthEast')
    arg = ['\rho ' num2str( round( density , 6 ) ) ' kg/m3'];
    text(35, 0.94, arg)
    arg = ['R ' num2str( round( sqrt(area/pi)*100 , 2 ) ) ' cm'];
    text(35, 0.92, arg)
    arg = ['A ' num2str( round( area , 6 ) ) ' m2'];
    text(35, 0.90, arg)
    set(gca,'FontSize', 14)
    saveas(fig, 'img/mt_main_FM_single_coax.jpg')
    close all;
    
    function FM = get_FM_from_single(arr)
        disp('FM for conditions')
        thrust = arr(4)
        Pideal  = (thrust)^(3/2) / sqrt(2*density*area)
        Preal   = arr(2)
        FM      = Pideal / Preal
        disp('------')
    end
    function FM = get_FM_from_coax(arr)   
        disp('FM for conditions')
        thrust = arr(7)
        Pideal  = (thrust)^(3/2) / sqrt(2*density*area)
        Preal   = arr(6)
        FM      = Pideal / Preal
        disp('------')
    end    
end





