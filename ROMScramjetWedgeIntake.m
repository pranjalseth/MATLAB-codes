%% Pranjal Seth

%%%

function s=cpinteg(t,k)  
    %%function to integrat c_p*dT given the limits 
    %t=temperature which is the VARIABLE
    %k=parameter matrix to calculate c_p (value specified later in fsolve
    %in the main code)
     %Gas constat
   b=4:-1:0;
     l=t.^b;
   s=k*l';
end



function hy = eachReflection (th1, t_pre, v_pre, r_pre, h_tot, t_ref, a1) 

%function to calculate properties after each foreshock

%t_pre = temperature after previous shock
%v_pre = velocity after previous shock
%r_pre = density after previous shock

f = @(t) cpinteg(t,a1);
R = 287.035; %gas constant

            b_j = th1 : 0.1 : 89;
            z = length(b_j);
            it(1) = 1;
                    
            for j = 1:z
                             b = b_j(j);
                             s_b = ((R * t_pre + (sind(b) * v_pre)^2) * (tand(b - th1)/tand(b)) - (v_pre * cosd(b) * tand(b - th1))^2)/R;
                             
                             c = integral(f ,t_ref,s_b,'ArrayValued',true);
                             it(j+1) = 2 * (h_tot - c) * ((cosd(b - th1))^2) - (v_pre * cosd(b))^2;
                             
                            if it(j+1) * it(j) <0
                                
                                b_ob = b;
                                break;
                                
                            end
            end
                    
                 
            hy(1) = b_ob; %beta of the current shock
            
            %calculating properties if there wasnt unstart
            
            %temperature after shock reflection
            hy(2) = ((R * t_pre + (sind(b_ob) * v_pre)^2) * (tand(b_ob - th1)/tand(b_ob)) - (v_pre * cosd(b_ob) * tand(b_ob - th1))^2)/R; 
        
            
            hy(3) = v_pre * cosd(b_ob)/cosd(b_ob - th1);                 %velocity after each shock reflection
            hy(4) = r_pre * tand(b_ob)/tand(b_ob - th1);                 %density after shock 
                                        
            hy(5) = cpinteg(hy(2),a1);                                   %cp after each shock reflection
            hy(6) = hy(5)/(hy(5)- R);                                    %gamma after each shock reflection
                        
            hy(7) = hy(3)/(hy(6) * hy(2) * R)^0.5;                       %mach number after each shock reflection
            hy(8) = hy(4) * R * hy(2);                                   %pressure after each shock reflection  
            
end

%%%

function [tf, vf, rf, pf, bf, cpf, gf, mf]  = forebody(theta_n, aoa, t_0, v_0, r_0, h_tot, t_ref, e,a1)

    %function to return propertirs after the forebody
    %all the variabes mean the same as said in the main code
      
    R = 287.035;      % Gas constant
    f = @(t) cpinteg(t,a1); 
    
        th1 = (theta_n+aoa);
        
        t = t_0;
        v = v_0;
        r = r_0;
        
  for i=1:e        %controlling loop for the number of ramps  
        
        if e==1
            th1 = theta_n+aoa;
        else
            th1 = theta_n(i)-theta_n(i-1);
        end
        
        b_j = th1 : 0.1 : 89;
        z = length(b_j);
        it(1) = 1;
        
    for j = 1 : z  %loop to calculate the value of beta
        
           b = b_j(j);
           s_b = ((R * t + (sind(b) * v)^2)*(tand(b - th1) / tand(b)) - (v * cosd(b) * tand(b - th1))^2) /R;
           c = integral(f ,t_ref,s_b,'ArrayValued',true);
           it(j+1) = 2 * (h_tot - c) * ((cosd(b - th1))^2) - (v_0 * cosd(b))^2;
           
                            if it(j+1) * it(j)<0
                                
                                b_fb(i)=b;
                                break;
                                
                            end
    end
  
           
            t_fb(i) = ((R * t + (sind(b_fb(i)) * v)^2)*(tand(b_fb(i) - th1)/tand(b_fb(i))) - (v * cosd(b_fb) * tand(b_fb - th1))^2)/R;  %after foreshock
            r_fb(i) = r*tand(b_fb(i)) / tand(b_fb(i) - th1);     %density after foreshock
            
            v_fb(i) = v*cosd(b_fb(i)) / cosd(b_fb(i) - th1);     %velocity after foreshock
            p_fb(i) = r_fb(i) * R * t_fb(i);                           %pressure after foreshock
            
            t = t_fb(i);
            v = v_fb(i);
            p = p_fb(i);
            
            cp_fb(i) = cpinteg(t_fb(i), a1);
            g_fb(i) = cp_fb(i) /(cp_fb(i) - R);           %gamma after weach foreeshock 
            m_fb(i) = v_fb(i)/(g_fb(i) * t_fb(i) * R)^0.5;   %mach number after each foreshock 
  end
  
  tf = t_fb;
  vf = v_fb;
  rf = r_fb;
  pf = p_fb;
  bf = b_fb;
  cpf = cp_fb;
  gf = g_fb;
  mf = m_fb;
end

%%%

function [x,y]=intersection(m,xb,yb,l,t1,t2, x_c, y_c)
    %%function to find the intersection of shock with the ramp OR cowl
    %m=slope of oblique shock
%     xb,yb=previous co-ordinate before shock
%     l= used to command whether intersection is with cowl or ramp. 
%     l is user defined. In this code,l=0 means intersection with ramp 
%     and l=1 means intersection with cowl. Input by user.
%     t1=ramp angle th_n
%     t2=cowl angle th_c
%     x_c, y_c= Co-ordinates of cowl
        

        if l==0 %odd shock i.e. intersecting with ramp
            x=(xb*m-yb)/(m+tand(t1));
            y=-x*tand(t1);
        end
        
        if l==1 %even shock i.e. intersecting cowl
            x=(xb*m-yb+y_c+x_c*tand(t2))/(m+tand(t2));
            y=-x*tand(t2)+y_c+x_c*tand(t2);
        end 
end
      
%%%

function st = machplotter(l_r, x_ob, m_ob, m_fb,q)


            %function to plot pressure after each shock. 
            
            v = cumsum(l_r);
            
                for i = 1:length(v)
                    
                    if(any(v(i) > x_ob(1)))
                        
                        v(i) = x_ob(1);
                    end
                    
                    
                end
                k = [0,v];
                    
                    
                for j = 1 : length(k)-1   
                    plot([k(j:j+1)], m_fb(j)*ones(2,1), 'r-o');
                    hold on;
                end
                
                plot(x_ob(1:q+1), m_ob(1:q+1),'-o');   
                
                hold off;
                
                xlabel('Length from nose to the isolator inlet');
                ylabel('Mach Number Variation');
                legend('ramp','cowl');
              st=0;  
end
                                      
%%%

function st = presplotter(l_r, x_ob, p_ob, p_fb,q)


            %function to plot pressure after each shock. 
            
            v = cumsum(l_r);
            
                for i = 1:length(v)
                    
                    if(any(v(i) > x_ob(1)))
                        
                        v(i) = x_ob(1);
                    end
                    
                    
                end
                k = [0,v];
                    
                    
                for j = 1 : length(k)-1   
                    plot([k(j:j+1)], p_fb(j)*ones(2,1),'r-o');
                    hold on;
                end
                
                plot(x_ob(1:q+1), p_ob(1:q+1),'-o');  
                hold off;
                
                xlabel('Length from nose to the isolator inlet');
                ylabel('Pressure Variation');
                legend('ramp','cowl');
                


st=0;
end
                                      
%%%

%%All the inpur variables here

m_0 = 6;             %freestream mach hitting the nose
t_ref = 298.15;      %reference temperature

p_0 = 1867;       %ambient pressure
R = 287.035;      % Gas constant

t_0 = 223.5;      %ambient temperature
aoa = 0;          %angle of attack


l_c = 6.083;        %length of cowl along its inclination
h_t = 0.152;        %isolator entrance height
l_r = [9.687];      % horizontal length of each ramp
l_n = sum(l_r);

theta_c = 0;        %cowl angle
theta_n = [3.65];   %ramp angle (array in the case of multiple wedges of the size of variable e)

e = 1;              %number of ramps (excluding the first cone if dealing with conical forebody)
n = 25;             %expected number of reflecions in cowl. (assume a higher number initially to avoid confusion)


%parameter array for calculating c_p
a1 = [1.9327e-10, -7.9999e-7, 1.1407e-3, -4.4890e-1, 1.0575e+3];


%% Dont edit below this line.
%%Free stream values calculation
cp_0 = cpinteg(t_0,a1);
g_0 = cp_0 /(cp_0 - R);             %ambient gamme
v_0 = m_0 * (g_0 * R * t_0)^0.5;    %ambient velocity
r_0 = p_0 / (R * t_0);              %ambient density

f = @(t) cpinteg(t,a1);
h_tot = integral(f, t_ref, t_0,'ArrayValued', true) + 0.5 * (v_0)^2;%total enthalphy of free stream conditions


x_c = l_n - l_c * cosd(theta_c);
y_c = -(l_n * tand(theta_n)+h_t - l_c * sind(theta_c));



%% Forebody yo

            %function to return properties after forebody
            [tf, vf, rf, pf, bf, cpf, gf, mf] = forebody(theta_n, aoa, t_0, v_0, r_0, h_tot, t_ref, e,a1);
            
            t_fb = tf;       %temperature post each  foreshock
            v_fb = vf;       %velocity post shock each foreshock
            
            r_fb = rf;       %density post each foreshock
            p_fb = pf;       %pressure post each foreshock    
            
            b_fb = bf;              %shock angle of each foreshock   
            cp_fb = cpf;            %c_p after each foreshock
   
            g_fb = gf;           %gamma after weach foreeshock 
            m_fb = mf;           %mach number after each foreshock        
            
           
         if e == 1
            m_fbs = -tand(b_fb(end) - aoa);                                %slope of final foreshock 
         else
            m_fbs = -tand(b_fb(end) + theta_n(end-1));
         end
            x_in= (y_c+x_c*tand(theta_c))/(m_fbs+tand(theta_c));    %intersection of the final foreshock and cowl
            
            
%%            
%Cowl           

 

%where xyz_ob(k,1)= properties after (k-1)th oblique shock inside cowl and k!=1 

v_tob = zeros(n,1);   %total velocity after each oblique shock. 
t_ob = zeros(n,1);    %temperatures after each oblique shock

cp_ob = zeros(n,1);   %cp after each oblique shock
r_ob = zeros(n,1);    %densities after each oblique shock

g_ob = zeros(n,1);    %gamma after each oblique shock
m_ob = zeros(n,1);    %mach numbers after each oblique shock

p_ob = zeros(n,1);    %pressure after each oblique shock
b_ob = zeros(n,1);    %shock angle wrt the initial flow direction inside cowl

v_tob(1,1) = v_fb(end);    %total velocity after final foreshock
t_ob(1,1) = t_fb(end);     %temperature after the final fore shock.

r_ob(1,1) = r_fb(end);     %density after final shock
cp_ob(1,1) = cp_fb(end);   %cp after final foreshock

g_ob(1,1) = g_fb(end);     %gamma after final foreshock
m_ob(1,1) = m_fb;     %mach after final foreshock

M = zeros(n,1);     %slope matrixfor shocks
M(1,1)= m_fbs;      %forebody shock slope

p_ob(1,1) = p_fb(end);       %pressure after final ramp
b_ob(1,1) = b_fb(end);      %shock angle wrt the initial flow direction inside cowl

%co-ordinates for shock intersection with cowl and ramp
x_ob = zeros(n,1);
y_ob = zeros(n,1);

th1 = theta_n(end) - theta_c;    %flow turn angle

x_ob(1,1) = x_c;        %cowl lip x-cordinate x_c
y_ob(1,1) = y_c;        %cowl lip y-coordinate y_c

x_t = l_n;              % the x-coord of isolator entrance is also known

b_ob(1,1) = b_fb;


if x_in <= x_c      %on design and mo<md off design
    
    for i = 1:n
    
        r=mod(i,2);
        
        if (r~=0) %odd numbered shoch which intersect ramp
        
            %function to calculate properties after each reflection
            hy = eachReflection (th1, t_ob(i), v_tob(i), r_ob(i), h_tot, t_ref, a1); 
            b_ob(i+1) = hy(1);    %shock angle after each reflection
            
                if b_ob == 90     %checking for unstart i.e. normal shock occurance..2
                     
                     disp(unstart);
                     break
                     
                 end
                 
            %properties after foreshock if there was no 
            t_ob(i+1) = hy(2);      %temperature
            v_tob(i+1) = hy(3);     %velocity
            
            r_ob(i+1) = hy(4);      %density
            cp_ob(i+1) = hy(5);     %cp
            
            g_ob(i+1) = hy(6);      %gamma
            m_ob(i+1) = hy(7);      %mach number
            p_ob(i+1) = hy(8);      %pressure
            
            
            %calculating the intersection points with ramp
            M(i+1) = tand(abs(b_ob(i+1))-theta_n);                                                            % slope of oblique shock wave
           [x_ob(i+1), y_ob(i+1)] = intersection(M(i+1),x_ob(i), y_ob(i),0, theta_n, theta_c, x_c, y_c);      %x and y coord of intersection
           
                if(x_ob(i+1)>x_t)
                    
                    break
                    
                end
        end
        
        if r == 0 %even numbered shock i.e. intersecting cowl
            
            %function to calculate properties after each reflection
            hy = eachReflection (th1, t_ob(i), v_tob(i), r_ob(i), h_tot, t_ref, a1); 
            b_ob(i+1) = hy(1);    %shock angle after each reflection
            
                if b_ob == 90     %checking for unstart i.e. normal shock occurance..2
                     
                     disp(unstart);
                     break
                     
                 end
                 
            %properties after foreshock if there was no 
            t_ob(i+1) = hy(2);      %temperature
            v_tob(i+1) = hy(3);     %velocity
            
            r_ob(i+1) = hy(4);      %density
            cp_ob(i+1) = hy(5);     %cp
            
            g_ob(i+1) = hy(6);      %gamma
            m_ob(i+1) = hy(7);      %mach number
            p_ob(i+1) = hy(8);      %pressure
            
            %calculating the intersection points with ramp
            M(i+1) = -tand(abs(b_ob(i+1)) + theta_c);                                                       % slope of oblique shock wave
            [x_ob(i+1), y_ob(i+1)] = intersection(M(i+1),x_ob(i), y_ob(i),1, theta_n, theta_c, x_c, y_c);     %x and y coord of intersection
           
                if(x_ob(i+1) > x_t)
                    
                    break
                    
                end
        end
        

    end
    
    q=i; %to be used for plotting


%Plotting

%Plotting

figure(1)
st = presplotter(l_r, x_ob, p_ob, p_fb,q);
figure(2)
st2 =  machplotter(l_r,x_ob, m_ob, m_fb, q);
    

                
end

    
    if x_in>x_c
        
        w =(tand(theta_c) - tand(theta_n));  %area gradient wrt length
        gr = g_fb * R;
        df = @(u,t) [w * (u/(l_n * tand(theta_n) + h_t - (l_n - x) * tand(theta_c)- x * tand(theta_n)) * (gr * t/(u^2 - gr * t))); -w *(u^2/(cp_fb(l_n * tand(theta_n) + h_t -(l_n - x) * tand(theta_c) - x * tand(theta_n))) * (gr * t/(u^2 - gr * t)))]; %funcgtion for ode45
        
        xspan = x_c : 0.001 : x_t;      %boundaries of control volume with a step size of 0.001  
        I = [v_fb(end) ; t_fb(end)];               %initial values at the start of control volume i.e. conditions at cowl
        
        [v_ise, t_ise] = ode45(df, xspan, I);        %x_se= properties at isolator entrance
        
        cp_ise = cpinteg(t_ise, a1);         %Cp at isolator entrance
        g_ise = cp_ise/(cp_ise-R);           %gamma at isolator entrance
        
        m_ise = v_ise/(g_ise * R * t_ise)^0.5;
    
    end
    




    
    


















































































