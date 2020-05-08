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

