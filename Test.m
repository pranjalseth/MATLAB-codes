    %% Final Optimization
    
    LB = [0;100;1;2;1000;1;1];
    UB = [30;500;63;4;9000;27;19];
    IntCon = [3 6 7];
    %opts = optimoptions(@ga,'PlotFcn', {@gaplotrange,@gaplotbestindiv,@gaplotdistance});
    
    %[x1, fval1] = fmincon(@objfun1, [30 15 160 0.16 3.3 150], [],[],[],[],LB,UB,@constraintfunc);
    [x2, fval2, exitFlag, Output] = ga(@paper_objective,7,[],[],[],[],LB,UB,@paper_constraintfunc, IntCon)
    %[x3, fval3] = patternsearch(@discrete, [30 15 140 3 3.5 8000 5 2], [],[],[],[],LB,UB,@constraintfunc, IntCon);
    %[x4, fval4] = simulannealbnd(@objfun1,[30 15 160 0.16 3.3 150], LB, UB);  
    
    
    
   