function [ TDvalue ] = tanDelta( M, SA )
%Tan(delta) as per theta-beta-M equation

Gamma = 1.4;

TDvalue = ((2*cotd(SA)*(M^2*sind(SA)^2 - 1))/(2 + M^2*(Gamma+1-2*sind(SA)^2)));

end

