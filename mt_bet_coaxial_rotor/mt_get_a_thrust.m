function [a_thrust_m1, a_thrust_m2] = mt_get_a_thrust(motornum, nblades, thrust_lvl)

    if motornum == 111
        % [thrust_vs_rpm] m1, motornum 111, nblades 2, afit 0.0003353000000000073, lsqe 18.5527698185157
        % [thrust_vs_rpm] m1, motornum 111, nblades 3, afit 0.0004221000000000094, lsqe 19.409315284780362
        % [thrust_vs_rpm] m1, motornum 111, nblades 6, afit 0.0005811999999999628, lsqe 37.44608385206588
        if nblades == 2
            a_thrust_m1 = 0.000335;
            a_thrust_m2 = 0.000000;
            return
        end
        if nblades == 3
            a_thrust_m1 = 0.000422;
            a_thrust_m2 = 0.000000;
            return
        end
        if nblades == 6
            a_thrust_m1 = 0.000582;
            a_thrust_m2 = 0.000000;
            return
        end
    end
    if motornum == 12 && nblades == 2
        % [thrust_vs_rpm] m1, motornum 12, nblades 2, afit 0.00034960000000000763, lsqe 5.4404910218461655
        % [thrust_vs_rpm] m2, motornum 12, nblades 2, afit 0.001323899999999955, lsqe 10.12220052463191
        % [thrust_vs_rpm] m1, motornum 12, nblades 2, afit 0.00033460000000000727, lsqe 4.048238291570132
        % [thrust_vs_rpm] m2, motornum 12, nblades 2, afit 0.0010240999999997851, lsqe 6.082302428829693
        % [thrust_vs_rpm] m1, motornum 12, nblades 2, afit 0.00034310000000000747, lsqe 11.037161519012397
        % [thrust_vs_rpm] m2, motornum 12, nblades 2, afit 0.0008027999999998481, lsqe 13.06405943361954
        if thrust_lvl == 35
            a_thrust_m1 = 0.000350;
            a_thrust_m2 = 0.000206;
            return 
        end
        if thrust_lvl == 45
            a_thrust_m1 = 0.000335;
            a_thrust_m2 = 0.000218;
            return 
        end
        if thrust_lvl == 55
            a_thrust_m1 = 0.000343;
            a_thrust_m2 = 0.000207;
            return 
        end
    end
    if motornum == 12 && nblades == 3
        % [thrust_vs_rpm] m1, motornum 12, nblades 3, afit 0.0004389000000000098, lsqe 4.5850608370652886
        % [thrust_vs_rpm] m2, motornum 12, nblades 3, afit 0.001963300000000295, lsqe 15.476399652523318
        % [thrust_vs_rpm] m1, motornum 12, nblades 3, afit 0.00041970000000000933, lsqe 4.09240641714146
        % [thrust_vs_rpm] m2, motornum 12, nblades 3, afit 0.0014551000000000292, lsqe 15.253076805235015
        % [thrust_vs_rpm] m1, motornum 12, nblades 3, afit 0.0004227000000000094, lsqe 8.359580816971995
        % [thrust_vs_rpm] m2, motornum 12, nblades 3, afit 0.0012105999999998908, lsqe 23.03952352645525
        if thrust_lvl == 40
            a_thrust_m1 = 0.000439;
            a_thrust_m2 = 0.001963;
            return 
        end
        if thrust_lvl == 50
            a_thrust_m1 = 0.000420;
            a_thrust_m2 = 0.001455;
            return 
        end
        if thrust_lvl == 60
            a_thrust_m1 = 0.000423;
            a_thrust_m2 = 0.001211;
            return 
        end
    end
    % If we get here report a problem
    % raise RuntimeErrorend
end
