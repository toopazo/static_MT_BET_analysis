function fig = blade_plot(nfig1, nfig2, blade_st)

    y_arr       = blade_st.y_arr;
    theta_arr   = blade_st.theta_arr; 
    sigma_arr   = blade_st.sigma_arr;
    chord_arr   = blade_st.chord_arr;
    R           = blade_st.R;    
    nsections   = blade_st.nsections;
    
    psi_arr     = linspace(0, 2*pi, nsections+1);

    fig = figure(nfig1);
    
    subplot(3, 1, 1);
    hold on;
    grid on;
    plot(y_arr, chord_arr, '-*'); 
    xlabel('radius m');
    ylabel('chord m');
    ylim([0, 0.2]);
    
    subplot(3, 1, 2);
    hold on;
    grid on;
    plot(y_arr, sigma_arr, '-*');
    xlabel('radius m');
    ylabel('solidity');
    ylim([0, 0.2]);
    
    subplot(3, 1, 3);
    hold on;
    grid on;
    plot(y_arr, rad2deg(theta_arr), '-*');
    xlabel('radius m');
    ylabel('pitch deg');
    ylim([-10, 40]);
    
    %%%%%%%%%%%%%%
    
    nfig = nfig2;
    fig = figure(nfig);
    hold on;
    grid on;     
    num_dpsi    = length(psi_arr);
    num_dr      = length(y_arr);      
    rotor_xyz   = blade_model_rotor_xyz(psi_arr, y_arr, blade_st);
    for i = 1:num_dpsi
        x_arr = squeeze( rotor_xyz(1, i, :) );
        y_arr = squeeze( rotor_xyz(2, i, :) );
        z_arr = squeeze( rotor_xyz(3, i, :) );
        plot3(x_arr, y_arr, z_arr, '-*b');
    end                
    arg = [ 'blade flapping \beta(\psi)' ];
    title(arg);
    xlabel('x axis m');
    ylabel('y axis m');
    zlabel('z axis m'); 
    plot3D_rotorAxes(nfig)

end
