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














































