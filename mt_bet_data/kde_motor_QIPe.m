function [...
    Qmotor  , ... % Motor torque
    Imotor  , ... % current allowed into the motor
    Vmotor  , ... % effective voltage applied to the motor
    Pmotor  , ... % shaft power
    effcy     ... % Pelec / Pshaft
    ] = kde_motor_QIPe(...
        omega       , ... % wrel of rotors
        uinput      , ... % 
        motortype   , ... % 
        modeltype     ... % 
    )

    % Km    = Q / sqrt(Ploss)   = Kt / sqrt(R0) = motor size constant, in units of Nm/sqrt(W)
    % Kv    = w_noload / Vpeak  = 1 / Ke        = motor velocity constant, in units of rad/Vs
    % Kt    = Q / I             = 1 / Kv = Ke   = motor torque constant, in units of Nm/A
    % Vemf  = Ke * w            = w / Kv        = back-emf voltage produced by the winding, in units of V    
    % mRPM  = eRPM * 60 * 2 / Np                = mechanical RPM, in units of RPM
    % eRPM  =                                   = electrical RPM, in units of RPM
    % Np    =                                   = number of magnetic motor poles (KDECAN_Bus_Protocol_1.0.3.pdf)  
    rpm2rads = pi/30;    
        
    switch motortype
        case 'KDE6213XF'
            % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde6213xf-185
            motor_Kv    = 185 * rpm2rads;
            motor_Kt    = 0.0516;
            motor_Km    = 0.1924;
            motor_Ke    = motor_Kt; % = 1/Kv
            motor_Imax  = 62;
            motor_Pmax  = 2755;
            motor_Vmax  = 52.2;
            motor_I0    = 0.6 * (motor_Vmax-5)/10;    % KDE reported I0 = 0.6A @ 10V
            motor_R     = 0.072;
            motor_Np    = 28;               % number of magnetic motor poles
        case 'KDE8218XF'
            % https://www.kdedirect.com/collections/uas-multi-rotor-brushless-motors/products/kde8218xf-120 
            motor_Kv    = 120 * rpm2rads;
            motor_Kt    = 0.0796;
            motor_Km    = 0.4137;
            motor_Ke    = motor_Kt; % = 1/Kv
            motor_Imax  = 110;
            motor_Pmax  = 5695;
            motor_Vmax  = 60.9;
            motor_I0    = 0.8 * (motor_Vmax-5)/10;    % KDE reported I0 = 0.6A @ 10V
            motor_R     = 0.037;
            motor_Np    = 28;               % number of magnetic motor poles                   
        otherwise
            motor_Kv    = NaN;
            motor_Kt    = NaN;
            motor_Km    = NaN;
            motor_Ke    = NaN;
            motor_Imax  = NaN;
            motor_Pmax  = NaN;
            motor_Vmax  = NaN;
            motor_I0    = NaN;
            motor_R     = NaN;
            motor_Np    = NaN;
            fprintf('[kde_motor_QIPe] Unrecognized rotortype %s \n', rotortype);
    end   
    
    if strcmp(modeltype, 'throttle')    
        throttle = uinput;   
        % Veff should be Vmax whenever throttle is 0.5 and omega is < omega_max/3 
%        omega
%        omega_max = (motor_Kv*motor_Vmax)
%        omega_ref = omega_max / 3
%        
%    %    err_omega   = esc_p1*0.5 - omega_ref
%    %    motor_Vmax  = esc_p2*err_omega  = esc_p2*esc_p1*0.5 - esc_p2*omega_ref
%    %    motor_Vmax  = 1*esc_p1*0.5 - 1*omega_ref
%    %    esc_p1      = motor_Vmax + omega_ref
%        esc_p1      = 2 * ( motor_Vmax + omega_ref )
%        esc_p2      = 1
        esc_p1 = 1200;
        esc_p2 = 1;

        err_omega   = esc_p1*throttle - omega;
        Veff        = esc_p2*err_omega;
    end
    if strcmp(modeltype, 'voltage')
        Veff = uinput;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Saturate effective voltage (equivalent continuous voltage comming from ESC)
    if Veff <= 0
        Veff = 0;
    end 
    if Veff >= motor_Vmax
        Veff = motor_Vmax;
    end
        
    % Current allowed into the motor
    % Veff = R*I + Ke*omega
    allowed_I = (Veff - motor_Kt * omega) / motor_R;
    
    % Saturate allowed current
    if allowed_I <= 0
        allowed_I = 0;
    end
    if allowed_I >= motor_Imax
        allowed_I = motor_Imax;
    end
    % Saturate according to max_Pelec
    % Max Pelec for all (V, I) combinations
    % max_Pelec = motor_Vmax * motor_Imax;
    % Imax @ given V
    % Imax = max_Pelec / Veff;
    % if allowed_I >= Imax
    %     allowed_I = Imax;
    % end
    
    % Torque
    Qmotor = (allowed_I - motor_I0) * motor_Kt;
    if Qmotor < 0
        Qmotor = 0;
    end
    
    % Motor output
    Qmotor  = Qmotor;               
    Imotor  = allowed_I;
    Vmotor  = Veff;
    Pmotor  = Qmotor * omega;
    Pelec   = Veff * allowed_I;
    effcy   = Pmotor / (Pelec + 10^-4); % Avoid 0/0 = NaN
    
end

