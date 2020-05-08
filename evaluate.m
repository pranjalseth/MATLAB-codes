function [ expressionValue ] = evaluate( MachNo, ShockAngle, FlowDeflectionAngle )
%Diiference between Left and Right sides of theta-beta-M equation

Gamma = 1.4;

R = 2*cotd(ShockAngle)*((MachNo^2*sind(ShockAngle)^2) - 1)/(MachNo^2*(Gamma+cosd(2*ShockAngle)) + 2);
L = tand(FlowDeflectionAngle);
diff = L - R;

expressionValue = abs(diff);

end

