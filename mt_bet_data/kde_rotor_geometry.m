function [y_arr, chord_arr, theta_arr] = kde_rotor_geometry(rotortype)

    rpm2rads   = pi/30;
    rads2rpm   = 30/pi;
    in2cm      = 2.54;     % 1in = 2.54cm
    ft2cm      = 30.48;    % 1ft = 30.48cm
    cm2m       = 1/100;

    persistent kde_geometry
    if isempty(kde_geometry)
        % st = load([bpath '/' 'kde_geometry']);
        st = load('mt_bet_data/kde_geometry');
        kde_geometry = st.kde_geometry;
    end

    if strcmp(rotortype, 'KDECF245DP') || strcmp(rotortype, 'KDECF245TP') || strcmp(rotortype, 'KDECF245HP')
        % 27.5 x 8.1
        sizefactor  = 24.5 ./ 30.5;
        pitchfactor = 8.1 ./ 9.7;
        y_arr       = kde_geometry.CF305.y_arr .* sizefactor;
        chord_arr   = kde_geometry.CF305.chord_arr .* sizefactor;
        theta_arr   = kde_geometry.CF305.theta_arr .* pitchfactor;        
    end
    if strcmp(rotortype, 'KDECF275DP')
        % 27.5 x 8.9
        sizefactor  = 27.5 ./ 30.5;
        pitchfactor = 8.9 ./ 9.7;
        y_arr       = kde_geometry.CF305.y_arr .* sizefactor;
        chord_arr   = kde_geometry.CF305.chord_arr .* sizefactor;
        theta_arr   = kde_geometry.CF305.theta_arr .* pitchfactor;        
    end
    if strcmp(rotortype, 'KDECF305DP')
        % 30.5 x 9.7
        y_arr       = kde_geometry.CF305.y_arr;
        chord_arr   = kde_geometry.CF305.chord_arr;
        theta_arr   = kde_geometry.CF305.theta_arr;        
    end
    if strcmp(rotortype, 'NACA5407shape')        
        R = 15.25 .* in2cm .* cm2m;
        nsections = 11;
        y_arr = linspace(0, R, nsections);
        chord_arr = kde_geometry.NACA5407.shape.x_arr;
        theta_arr = kde_geometry.NACA5407.shape.y_arr;
    end    
    if strcmp(rotortype, 'NACA5407')        
        y_arr       = kde_geometry.NACA5407.y_arr;
        chord_arr   = kde_geometry.NACA5407.chord_arr;
        theta_arr   = kde_geometry.NACA5407.theta_arr;  
        
    end
end

