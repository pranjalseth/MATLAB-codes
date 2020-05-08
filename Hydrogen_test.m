    %% Final Optimization
    
    LB = [0;100;1;2;1000;0;1;1.5;1;338;0;0];
    UB = [30;500;63;10;10000;500;19;5;15;373;5;1000];
    IntCon = [3 7];
    %opts = optimoptions(@ga,'PlotFcn', {@gaplotrange,@gaplotbestindiv,@gaplotdistance});
    
    [x2, fval2, exitFlag, Output] = ga(@hydrogen_objective,12,[],[],[],[],LB,UB,@hydrogen_constraintfunc, IntCon)

    
    
   