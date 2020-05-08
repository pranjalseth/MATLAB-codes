function E = equations( Var )
% EQUATIONS FOR FSOLVE as per: AIAA 2005-7357
%   Simultaneous Non Linear equations for FSOLVE function

P_o = 2511; %Static Pressure: Pascals
T_o = 225; %Freestream Temperature: Kelvin
M_o = 2.2; %Freestream Mach Number
Gamma = 1.4; %Ratio of Specific Heats: taken constant

M_f = 0.8; 
P_f = 26000; %Final Pressure

M0 = M_o;
M4 = M_f;





M1 = Var(1,1); %Mach Numbers
M2 = Var(2,1);
M3 = Var(3,1);
theta1 = Var(1,2); %Shock Angles
theta2 = Var(2,2);
theta3 = Var(3,2);
delta1 = Var(1,3); %Flow Deflection Angles
delta2 = Var(2,3);
delta3 = Var(3,3);


E(1) = M1^2 - (((Gamma+1)^2*M0^4*sind(theta1)^2 - 4*(M0^2*sind(theta1)^2-1)*(Gamma*M0^2*sind(theta1)^2+1))/...
    ((2*Gamma*M0^2*sind(theta1)^2-(Gamma-1))*((Gamma-1)*M0^2*sind(theta1)^2+2)));
E(2) = M2^2 - (((Gamma+1)^2*M1^4*sind(theta2)^2 - 4*(M1^2*sind(theta2)^2-1)*(Gamma*M1^2*sind(theta2)^2+1))/...
    ((2*Gamma*M1^2*sind(theta2)^2-(Gamma-1))*((Gamma-1)*M1^2*sind(theta2)^2+2)));
E(3) = M3^2 - (((Gamma+1)^2*M2^4*sind(theta3)^2 - 4*(M2^2*sind(theta3)^2-1)*(Gamma*M2^2*sind(theta3)^2+1))/...
    ((2*Gamma*M2^2*sind(theta3)^2-(Gamma-1))*((Gamma-1)*M2^2*sind(theta3)^2+2)));
E(4) = tand(delta1) - ((2*cotd(theta1)*(M0^2*sind(theta1)^2 - 1))/...
    (2 + M0^2*(Gamma+1-2*sind(theta1)^2)));
E(5) = tand(delta2) - ((2*cotd(theta2)*(M1^2*sind(theta2)^2 - 1))/...
    (2 + M1^2*(Gamma+1-2*sind(theta2)^2)));
E(6) = tand(delta3) - ((2*cotd(theta3)*(M2^2*sind(theta3)^2 - 1))/...
    (2 + M2^2*(Gamma+1-2*sind(theta3)^2)));
E(7) = M0*sind(theta1) - M1*sind(theta2);
E(8) = M1*sind(theta2) - M2*sind(theta3);
E(9) = M4^2 - (((Gamma-1)*M3^2+2)/(2*Gamma*M3^2 - (Gamma-1)));

end

