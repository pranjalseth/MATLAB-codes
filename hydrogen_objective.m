function y = hydrogen_objective(x)
   
    % x(1) = Velocity in [m/s] 
    % x(2) = Altitude in [m]  
    % x(3) = Propeller 
    % x(4) = Maximum TOW in [kg]  
    % x(5) = Motor / Propeller RPM
    % x(6) = Current 
    % x(7) = Motor
    % x(8) = Oxygen Stoichiometry
    % x(9) = Stack Pressure  (Bars)
    % x(10) = Stack Temperature (Kelvin)
    % x(11) = Tank Radius
    % x(12) = Mass Hydrogen
    
    %% For camera Canon SX 230HS
    
    pixel_size = 1.5 * 10 ^(-6);                      % Pixel Size in [m] 
    focal_length = 0.005;                             % Equivalent focal length in [m] 
    time_betweenimages = 3;                           % Time between consecutive images in [s]
    pixel_width = 3000;                               % Number of pixels in image width
    pixel_length = 4000;                              % Number of pixels in image length
    frontal_overlap = 70;                             % Frontal Overlap in [%]
    side_overlap = 60;                                % Side Overlap in [%]
       
    %% Maximize Quadrotor Productivity
 
    quadrotor_productivity = (1 - 0.01 * side_overlap) * pixel_size * pixel_width * x(1) * x(2) / focal_length;  
    
        
    %% Hydrogen Fuel Cell (Palcan PC6-1200) Parameters
    
    number_cells = 25;
    faraday_constant = 96485.34;                   % [F] - Coulombs/Equiv
    
    %% Model 
       
    hydrogen_stoichiometry = 1.1;
    hydrogen_feed = x(6) * hydrogen_stoichiometry * number_cells / (2 * faraday_constant)       % [N_Hydrogen] - mol/s
          
    %% Maximize Flight Time
    
    molarmass_hydrogen = 0.002016;                                             % [M_H2] - kg/mol
    massflow_hydrogen = hydrogen_feed * molarmass_hydrogen;                    % kg/s
    
    flight_time = x(12) / (60 * massflow_hydrogen)                            % [t_f] - minutes               
    
    %% Optimization Function
        
    y = - flight_time - x(4) - quadrotor_productivity;
    
end
