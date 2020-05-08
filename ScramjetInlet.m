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

function [ expressionValue ] = evaluate( MachNo, ShockAngle, FlowDeflectionAngle )
%Diiference between Left and Right sides of theta-beta-M equation

Gamma = 1.4;

R = 2*cotd(ShockAngle)*((MachNo^2*sind(ShockAngle)^2) - 1)/(MachNo^2*(Gamma+cosd(2*ShockAngle)) + 2);
L = tand(FlowDeflectionAngle);
diff = L - R;

expressionValue = abs(diff);

end

%%
%ITERATIVE PROGRAM TO GET INITIAL GUESS VALUES

% take 20-30s to return values
%%
Gamma = 1.4;


% Reference: AIAA 2005-7357

M0 = 2.2; %reference 

for RA1 = 0:0.1:89 % Ramp Angle
    for SA1 = 0:0.1:89 % Shock Angle
        e1 = evaluate(M0, SA1, RA1);
        if (e1<=0.01)
            break
        end
    end
    M1sqr = Msqr(M0, SA1); % Mach Number
    M1 = sqrt(M1sqr);
    d = tanDelta(M0, SA1);
    FDA1 = atan(d)/pi*180; % Flow Deflection Angle

     for RA2 = 0:0.1:89
        for SA2 = 0:0.1:89
            e2 = evaluate(M1, SA2, (RA2-RA1));
            if (e2<=0.01)
             break
            end
        end
     if (M0*sind(SA1)-M1*sind(SA2)<=0.01)
        break
     end
    end
    M2sqr = Msqr(M1, SA2);
    M2 = sqrt(M2sqr);
    d = tanDelta(M1, SA2);
    FDA2 = atan(d)/pi*180;
    
    for RA3 = 0:0.1:89
     for SA3 = 0:0.1:89
         e3 = evaluate(M2, SA3, RA2-RA3);
         if (e3<=0.01)
             break
         end
     end
        if (M1*sind(SA2)-M2*sind(SA3)<=0.01)
         break
        end
    end
    M3sqr = Msqr(M2, SA3);
    M3 = sqrt(M3sqr);
    d = tanDelta(M2, SA3);
    FDA3 = atan(d)/pi*180;
    
    
    
    if (M3-1.25<=0.01) % reference
        break
    end
    
end

M4sqr = (((Gamma-1)*M3^2+2)/(2*Gamma*M3^2 - (Gamma-1)));
M4 = sqrt(M4sqr);


%values obtained from here are passsed as guess values in
%PranjalInletProblem.m



% M0 = 6;
% %for RA1 = 1:0.1:89 % Ramp Angle
% RA1 = 0.5;
%     for SA1 = 0:0.1:89 % Shock Angle
%         e1 = evaluate(M0, SA1, RA1);
%         if (e1<=0.01)
%             break
%         end
%     end
%     M1sqr = Msqr(M0, SA1);
%     M1 = sqrt(M1sqr);
%     d = tanDelta(M0, SA1);
%     FDA1 = atan(d)/pi*180;
% 
%      for RA2 = 0:0.1:89
%         for SA2 = 0:0.1:89
%             e2 = evaluate(M1, SA2, RA2);
%             if (e2<=0.01)
%              break
%             end
%         end
%      if (M0*sind(SA1)-M1*sind(SA2)<=0.01)
%         break
%      end
%     end
%     M2sqr = Msqr(M1, SA2);
%     M2 = sqrt(M2sqr);
%     d = tanDelta(M1, SA2);
%     FDA2 = atan(d)/pi*180;
%     
%     for RA3 = 0:0.1:89
%      for SA3 = 0:0.1:89
%          e3 = evaluate(M2, SA3, RA3);
%          if (e3<=0.01)
%              break
%          end
%      end
%         if (M1*sind(SA2)-M2*sind(SA3)<=0.01)
%          break
%         end
%     end
%     M3sqr = Msqr(M2, SA3);
%     M3 = sqrt(M3sqr);
%     d = tanDelta(M2, SA3);
%     FDA3 = atan(d)/pi*180;
%     
%     M4sqr = (((Gamma-1)*M3^2+2)/(2*Gamma*M3^2 - (Gamma-1)));
%     M4 = sqrt(M4sqr);
% %     
% %     if (M4-0.3<=0.01)
% %         break
% %     end
% %     
% % end

function [ Msq ] = Msqr( M, SA )
%MACH NUMBER Calculation
%   calculates mach number

Gamma = 1.4;

Msq = (((Gamma+1)^2*M^4*sind(SA)^2 - 4*(M^2*sind(SA)^2-1)*(Gamma*M^2*sind(SA)^2+1))/...
    ((2*Gamma*M^2*sind(SA)^2-(Gamma-1))*((Gamma-1)*M^2*sind(SA)^2+2)));

end

function [ PR ] = PRec( M, SA )
%PRESSURE RECOVERY CALCULATION FUNCTION
%   Caltulates pressure recovery across a shock for a given Mach Number and Shock Angle

G = 1.4;

a = (((G+1)*M^2*sind(SA)^2)/((G-1)*M^2*sind(SA)^2+2))^(G/(G-1));
b = ((G+1)/(2*G*M^2*sind(SA)^2-(G-1)))^(1/(G-1));

PR = a*b;

end

function [ TDvalue ] = tanDelta( M, SA )
%Tan(delta) as per theta-beta-M equation

Gamma = 1.4;

TDvalue = ((2*cotd(SA)*(M^2*sind(SA)^2 - 1))/(2 + M^2*(Gamma+1-2*sind(SA)^2)));

end

%%MAIN


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

