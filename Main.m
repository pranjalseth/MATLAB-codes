
%MAIN PROGRAM FOR SOLUTIONS

%% Given conditions

% Reference: AIAA 2005-7357

P_o = 2511; %Static Pressure: Pascals
T_o = 225; %Freestream Temperature: Kelvin
M_o = 2.2; %Freestream Mach Number - reference
Gamma = 1.4; %Ratio of Specific Heats: taken constant

M_f = 0.8; %final Mach Number - from LoopDeLoop.m
P_f = 26000; %Final Pressure

M0 = M_o;
M4 = M_f;

%%

%options = optimoptions(@fsolve,'MaxIter', 1000, 'MaxFunEvals', 5000, 'Algorithm', 'levenberg-marquardt');
%options = optimoptions(@fsolve,'MaxIter', 5000, 'MaxFunEvals', 10000);

varGuess = [1.96 32.2 6.2; 1.73 36.3 6.4; 1.25 50.6 13.1]; %Guess Values from LoopDeLoop.m
%Var = fsolve(@equations, varGuess, options);

Var = fsolve(@equations, varGuess);

Ma1 = Var(1,1); %Mach Numbers
Ma2 = Var(2,1);
Ma3 = Var(3,1);
theta1 = Var(1,2); %Shock Angles
theta2 = Var(2,2);
theta3 = Var(3,2);
delta1 = Var(1,3); %FLow Deflection Angles
delta2 = Var(2,3);
delta3 = Var(3,3);

beta1 = delta1; %Ramp Angles
beta2 = beta1 + delta2;
beta3 = beta2 - delta3;



PR1 = PRec(M0, theta1); %Pressure Recoveries
PR2 = PRec(Ma1, theta2);
PR3 = PRec(Ma2, theta3);

PR4 = PRec(Ma3, 90);

PR = PR1*PR2*PR3*PR4;

