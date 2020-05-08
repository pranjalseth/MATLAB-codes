function [x,y]=intersection(m,xb,yb,l,t1,t2, x_c, y_c)
    %%fnction to find the intersection of shock with the ramp OR cowl
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
      