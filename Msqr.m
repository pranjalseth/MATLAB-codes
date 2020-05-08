function [ Msq ] = Msqr( M, SA )
%MACH NUMBER Calculation
%   calculates mach number

Gamma = 1.4;

Msq = (((Gamma+1)^2*M^4*sind(SA)^2 - 4*(M^2*sind(SA)^2-1)*(Gamma*M^2*sind(SA)^2+1))/...
    ((2*Gamma*M^2*sind(SA)^2-(Gamma-1))*((Gamma-1)*M^2*sind(SA)^2+2)));

end

