function [R, chord, Nb, Cla, Cd0, Thover, omega_h, Pavail] = uh60_rotor_geometry()
                    
    rpm2rads        = (pi/30);  % 60 RPM = 1 RPsecond = 1 Hz = 2*pi rad/s
    in2m            = (0.0254); 
    ft2m            = (0.3048);
    lb2kg           = (0.453592);                    
                    
    % UH-60A %%%%%%%%%%%%%%%%%%%%%%%%                       
    % CT is between 0.0017 to 0.01, average being 0.008
    R         = 26.83 * ft2m;             % m
    chord     = 15 * in2m; % 1.75 * ft2m; % m
    omega_h   = 27;                       % rad/s
    % solidity  = 0.083;
    Nb        = 4;
    Cla       = 2*pi*0.9104; % 5.73
    Cd0       = Cla * deg2rad(20) / 25;
    
    % https://en.wikipedia.org/wiki/Sikorsky_UH-60_Black_Hawk
    vehicleM  = 16260 * lb2kg;    %  kg
    Thover    = vehicleM * 9.8;
    Pavail    = 2*1410*1000*0.80;  %  W                 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
end
