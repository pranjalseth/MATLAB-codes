function [ PR ] = PRec( M, SA )
%PRESSURE RECOVERY CALCULATION FUNCTION
%   Caltulates pressure recovery across a shock for a given Mach Number and Shock Angle

G = 1.4;

a = (((G+1)*M^2*sind(SA)^2)/((G-1)*M^2*sind(SA)^2+2))^(G/(G-1));
b = ((G+1)/(2*G*M^2*sind(SA)^2-(G-1)))^(1/(G-1));

PR = a*b;

end

