function [c,ceq] = hydrogen_constraintfunc(x)
%% For camera Canon SX 230HS
    
    pixel_size = 1.5 * 10 ^(-6);                      % Pixel Size in [m] 
    focal_length = 0.005;                             % Equivalent focal length in [m] 
    time_betweenimages = 3;                           % Time between consecutive images in [s]
    pixel_width = 3000;                               % Number of pixels in image width
    pixel_length = 4000;                              % Number of pixels in image length
    frontal_overlap = 70;                             % Frontal Overlap in [%]
    side_overlap = 60;                                % Side Overlap in [%]
    
    
    number_rotors = 6;
    density_air = 1.2041;
    acceleration_gravity = 9.81;
    drone_efficiency = 0.3;  
    
    %% Propeller
    
    d_propeller = [0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.1778 0.2032 0.2032 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.2286 0.254 0.254 0.254 0.254 0.254 0.254 0.254 0.254 0.2667 0.2667 0.2794 0.2794 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3048 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3302 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.3556 0.381];
    p_propeller = [0.254 0.0762 0.1016 0.127 0.1524 0.1651 0.1778 0.2032 0.2286 0.127 0.2032 0.254 0.0762 0.1016 0.1143 0.127 0.1524 0.1778 0.2032 0.2286 0.0762 0.1016 0.127 0.1524 0.1778 0.2032 0.2286 0.254 0.1143 0.1524 0.0762 0.1778 0.1524 0.254 0.2794 0.3048 0.2921 0.3175 0.3302 0.3556 0.1016 0.127 0.1524 0.254 0.2794 0.3302 0.1016 0.1397 0.1778 0.2032 0.2286 0.2540 0.2794 0.3048 0.3302 0.3556 0.1016 0.127 0.1524 0.1778 0.2032 0.2159 0.1778];
       
    
    %% Motor
    
    motor_kv = [2600 750 950 1100 1450 1630 2608 2640 3200 3900 4100 1250 790 1100 1000 1250 1450 910 1000];
    motor_resistance = [0.07 0.16 0.070 0.107 0.023 0.079 0.048 0.063 0.040 0.064 0.11 0.034 0.04 0.023 0.031 0.021 0.019 0.063 0.052];
    motor_noloadcurrent = [0.9 0.8 1 1 3.5 0.8 0.8 0.8 0.8 0.8 1.3 3 1.8 3.1 2.4 3 4 1.5 1.7];  
    
    
    %% Hydrogen Fuel Cell (Palcan PC6-1200) Parameters
    
    number_cells = 25;
    active_area = 96;                              % [A] - cm^2
    mass_stack = 1.5;                              % [m_s] - kg
    maximum_currentdensity = 22;                   % [I_l] - Ampere / cm^2
    gross_power = 1200;                            % [W_gross] - Watt
    
    %% Universal Constants
    
    faraday_constant = 96485.34;                   % [F] - Coulombs/Equiv
    gas_constant = 8.3145;                         % [R] - J/mol/K
    lambda = 1.4;                                  % lambda = C_p / C_v
    
    density_water = 997;                           % [rho] - kg/m^3 (at 298 K)
    specificheat_water = 4128;                     % [C_pw] - J/K/kg
    
    specificheat_air = 1006;                       % [C_pa] - J/K/kg
    molarmass_air = 0.02897;                       % [M_a] - kg/mol
    pressure_atmospheric = 1.01325;                % [P_atm] - Bars
    
    %% Model 
       
    hydrogen_stoichiometry = 1.1;
    hydrogen_feed = x(6) * hydrogen_stoichiometry * number_cells / (2 * faraday_constant);       % [N_Hydrogen] - mol/s
    oxygen_feed = x(6) * x(8) * number_cells / (4 * faraday_constant);                           % [N_Oxygen] - mol/s
    
    %% Nernst Voltage
    
    current_density = x(6) * number_cells / active_area;                                                     % [i] - Ampere / cm^2
    pressure_satwater = exp(70.434643 - (7362.6981 / x(10)) + 0.006952085 * x(10) - 9 * log(x(10)));         % [P_sat] - Bars
    
    x_satwater = pressure_satwater / x(9);
    x_inhumother = (1 - x_satwater) * 0.79;
    x_outhumother = (1 - x_satwater) / (1 + ((x(8) - 1 / x(8)) * 0.21 / 0.79));
    x_channelother = (x_inhumother - x_outhumother) / log(x_inhumother / x_outhumother);    
    
    pressure_oxygeninterface = x(9) * (1 - x_satwater - x_channelother * exp(0.291 * current_density / (x(10) ^ 0.832)));                 % [P_O2] - Bars
    pressure_hydrogeninterface = 0.5 * pressure_satwater * ( -1 + 1/(x_satwater * exp(1.653 * current_density / (x(10) ^ 1.334))));       % [P_H2] - Bars
    
    e_nernst = 1.229 - 0.00085 * (x(10) - 298.15) + 0.000043085 * x(10) * (log(pressure_hydrogeninterface) + 0.5 * log(pressure_oxygeninterface))  % [E_Nernst] - Volt
    
    %% Activation Overvoltage
    
    concentration_oxygeninterface = pressure_oxygeninterface / (5080000 * exp(-498/x(10)));
    activation_overvoltage = -0.944957 + 0.00301801 * x(10) + 0.00007401 * x(10) * log(concentration_oxygeninterface) - 0.0001880 * x(10) * log(x(6))   % [Eta_act] - Volt
    
    %% Ohmic Overvoltage
    
    lambda_membrane= 14.6;                                 
    resistivity_membrane = 181.6 * (1 + 0.03 * current_density + (0.062 * (x(10) / 303) ^ 2) * (current_density ^ 2.5)) / ((lambda_membrane - 0.634 - 3 * current_density) * exp(4.18 * (x(10) - 303) / x(10)));   % Ohm-cm
    thickness_membrane = 0.01275;                          % [t] - cm
    resistance_membrane = resistivity_membrane * thickness_membrane * number_cells / active_area;    % [R_ohmic] - Ohm
    ohmic_overvoltage = - x(6) * resistance_membrane      % [Eta_ohmic] - V
    
    %% Concentration Overvoltage
       
    concentration_overvoltage = (gas_constant * x(10)* log(1 - current_density/maximum_currentdensity) / ( 2 * faraday_constant))      % [Eta_conc] - Volt
    
    %% Cell Voltage
    
    voltage_cell = e_nernst + activation_overvoltage + ohmic_overvoltage + concentration_overvoltage;   % Volt
    voltage_stack = voltage_cell * number_cells;                                                     % Volt
    
    %% Work Cooling Pump
    
    maximum_voltage = 1.48;                                                    % [E_max] - Volt
    heat_total = (maximum_voltage - voltage_cell) * x(6) * number_cells;       % [Q_total] - Watt
    temperature_waterin = 298;                                                 % [T_waterin] - Kelvin
    factor_heat = 0.8;   
    mass_coolwater = factor_heat * heat_total / (specificheat_water * (x(10) - temperature_waterin));     % [m_coolwater] - kg /s
      
    s_factor = 1.5;
    efficiency_pump = 0.6;
    efficiency_pumpmotor = 0.8;
    pressure_dropcoolant = 1 * 10^4;                                           % [P_dropcoolant] - Pascals
    
    work_coolingpump = mass_coolwater * pressure_dropcoolant * s_factor / (density_water * efficiency_pump * efficiency_pumpmotor);    % [W_coolpump] - Watt
    
    %% Work Compressor
    
    mass_air = oxygen_feed * molarmass_air;                                    % [m_air] - kg/s
    temperature_airin = 298;                                                   % [T_airin] - Kelvin
    efficiency_compressor = 0.6;
        
    work_compressor = specificheat_air * mass_air * temperature_airin * (((x(9) / pressure_atmospheric) ^ ((lambda-1)/lambda)) - 1) / efficiency_compressor;    % [W_comp] - Watt
    
    %% Parasitic Work
    
    parasitic_work = work_coolingpump + work_compressor;                       % [W_parasitic] - Watt
    net_power = gross_power - parasitic_work;                                  % [W_net] - Watt
        
    
    %% Weight Minimization
    
    sigma_maxcomp = 1.9 * 10^4;                                                % Composite Overwrap Maximum Stress - Bars
    factor_safety = 2.5;
    thickness_composite = factor_safety * 1.5 * x(11) * (x(9) - pressure_atmospheric) / sigma_maxcomp;      % [t_composite] - m
    
    length_radiusratio = 2;
    tank_length = length_radiusratio * x(11);                                           % - m
    
    cylinder_compositevolume = 2 * pi * thickness_composite * tank_length * x(11);      % - m^3
    cylinder_tankvolume = pi * tank_length * (x(11) - thickness_composite)^2;           % - m^3
    
    hemisphere_compositevolume = 4 * pi * (x(11)^3 - (x(11) - thickness_composite)^3) / 3;    % - m^3
    hemisphere_tankvolume = 4 * pi * (x(11) - thickness_composite)^3 / 3;                     % - m^3
    
    f_mount = 0.1;
    density_composite = 1530;                                                  % [rho_composite] - kg/m^3
    regulator_mass = 0.35;                                                     % [m_reg] - kg
    
    mass_tank = (1 + f_mount) * (cylinder_compositevolume + hemisphere_compositevolume) * density_composite + regulator_mass + x(12);
    
    %% Feature Normalization and Regression

    mu = [0.2492 0.1574 12.1387 6817];
    sigma = [0.0329 0.0534 11.31086 3811.7];
    
    X = [d_propeller(x(3)) p_propeller(x(3)) x(1) x(5)];
    X = bsxfun(@minus, X, mu);
    X = bsxfun(@rdivide, X, sigma);
    X = [1 X];
    
    theta_ct = [0.0717; -0.0024; 0.0250; -0.0369; 0.0230];
    theta_cp = [0.0497; -0.0048; 0.0269; -0.0198; 0.0152];
    
    ct = X * theta_ct;
    cp = X * theta_cp;   
    
    %% Velocity Constraint
    
    upperbound_velocity = 0.01 * frontal_overlap * pixel_size * x(2) * pixel_length/(time_betweenimages * focal_length);
    c(1) = x(1) - upperbound_velocity;
    
    %% Volume Constraint
    
    required_tankvolume = x(12) * 4124.0079 * x(10) / (x(9) * 10 ^ 5)
    c(2) = required_tankvolume - (cylinder_tankvolume + hemisphere_tankvolume);
    
    %% Payload Capacity Constraint
    
    c(3) = (mass_tank + mass_stack) - x(4);

    %% Thrust Constraint
    
    thrust = number_rotors * density_air * ct * ((x(5)/60) ^ 2) * (( d_propeller(x(3)))^4); 
    c(4) =  x(4) * acceleration_gravity - thrust;  
    
    %% Motor/Propeller Matching Constraint
  
    motor_torque = (((voltage_stack - x(5) / motor_kv(x(7)))/motor_resistance(x(7))) - 0.8)/motor_kv(x(7));
    prop_torque = (density_air * cp * ((x(5)/60) ^ 2) * (( d_propeller(x(3)))^5))/(2*pi);

    c(5) = motor_torque - 2.5 * prop_torque;
    c(6) = 0.2 * prop_torque - motor_torque;
    
    %% Net Power Constraint
    
    motor_current = number_rotors * (voltage_stack - x(5) / motor_kv(x(7)))/motor_resistance(x(7))   
    c(7) = x(6) * voltage_stack - net_power;
    c(8) = motor_current - x(6);
    c(9) = current_density - maximum_currentdensity;
        
    %% GSD Constraint
    c(9) = pixel_size * x(2) / focal_length - 0.05;
    
    %% Equality Constraint
    ceq = [];
      
end
    
    
    
    