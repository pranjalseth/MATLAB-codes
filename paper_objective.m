function y = paper_objective(x)
   
    % x(1) = Velocity in [m/s] 
    % x(2) = Altitude in [m]  
    % x(3) = Propeller 
    % x(4) = Maximum TOW in [kg]  
    % x(5) = Motor / Propeller RPM
    % x(6) = Battery 
    % x(7) = Motor
    
    %% For camera Canon SX 230HS
    
    pixel_size = 1.5 * 10 ^(-6);                      % Pixel Size in [m] 
    focal_length = 0.005;                             % Equivalent focal length in [m] 
    time_betweenimages = 3;                           % Time between consecutive images in [s]
    pixel_width = 3000;                               % Number of pixels in image width
    pixel_length = 4000;                              % Number of pixels in image length
    frontal_overlap = 70;                             % Frontal Overlap in [%]
    side_overlap = 60;                                % Side Overlap in [%]
    
    %% Propeller
    
    d_propeller = [0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.2032 0.2032 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.254 0.254 0.254 0.254 0.254 0.254 0.254 0.254 0.2667 0.2667 0.2794 0.2794 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.381];
    p_propeller = [0.254 0.0762 0.1016 0.127 0.1524 0.1651 0.1778 0.2032 0.2286 0.127 0.2032 0.254 0.0762 0.1016 0.1143 0.127 0.1524 0.1778 0.2032 0.2286 0.0762 0.1016 0.127 0.1524 0.1778 0.2032 0.2286 0.254 0.1143 0.1524 0.0762 0.1778 0.1524 0.254 0.2794 0.3048 0.2921 0.3175 0.3302 0.3556 0.1016 0.127 0.1524 0.254 0.2794 0.3302 0.1016 0.1397 0.1778 0.2032 0.2286 0.2540 0.2794 0.3048 0.3302 0.3556 0.1016 0.127 0.1524 0.1778 0.2032 0.2159 0.1778];
       
    
    %% Battery
    
    battery_capacity = [1 2.2 3.7 4 4 5 5 5 5.8 5.8 6.2 8 1.4 3 4 4 5.2 5.2 5.2 6.6 8 5.2 1.6 1.4 5.2 3.3 3.3];                                                         % Battery Capacity in [Ah]
    battery_voltage = [11.1 11.1 18.5 11.1 14.8 11.1 14.8 18.5 18.5 22.2 14.8 11.1 11.1 11.1 11.1 14.8 11.1 14.8 22.2 22.2 14.8 11.1 11.1 14.8 11.1 14.8 11.1];            % Battery Voltage in [V]
    battery_crating = [25 25 25 25 40 25 60 25 25 40 40 25 40 10 10 10 10 10 10 10 10 15 40 65 20 25 25];

    %% Motor
    
    motor_kv = [2600 750 950 1100 1450 1630 2608 2640 3200 3900 4100 1250 790 1100 1000 1250 1450 910 1000];
    motor_resistance = [0.07 0.16 0.070 0.107 0.023 0.079 0.048 0.063 0.040 0.064 0.11 0.034 0.04 0.023 0.031 0.021 0.019 0.063 0.052];
    motor_noloadcurrent = [0.9 0.8 1 1 3.5 0.8 0.8 0.8 0.8 0.8 1.3 3 1.8 3.1 2.4 3 4 1.5 1.7]; 
    
    %% Flight Time
       
    number_rotors = 4;
    density_air = 1.2041;
    acceleration_gravity = 9.81;
    drone_efficiency = 0.3;     
    k_v = drone_efficiency * (d_propeller(x(3)) / 2) * (( 2 * pi * density_air * number_rotors)^0.5)/(acceleration_gravity ^ 1.5);       % Obtained from a research paper titled 'Endurance Optimisation of Battery-Powered Rotorcraft'
    k_u = abs(-17.2 * (x(4)^2) + 16.7 * x(4) - 3);          % x(4) is take off mass of drone

    k_b = 0.9;                                         % Battery Efficiency
    
    flight_time = battery_capacity(x(6)) * k_u * k_b * k_v * battery_voltage(x(6))/(x(4) ^ 1.5)
    
    %% Quadrotor Productivity
 
    quadrotor_productivity = (1 - 0.01 * side_overlap) * pixel_size * pixel_width * x(1) * x(2) / focal_length;  
       
    %% Combined Optimization Function
    
    y =  - flight_time - x(4) - quadrotor_productivity;
    
end
