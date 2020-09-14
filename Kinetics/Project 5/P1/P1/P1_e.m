
%1e
%The conversion of A, i.e. XAmm = 0.8673


function [] = P1_e
    
    V = 1000 ; %L
    Q = 10   ; %L/s
    tau  = V/Q;  %s
    t_1e = load('P1_time.txt');
    t = 60.*t_1e; %convert to sec
    
    E_t0 = load('P1_E(t).txt');
    E_t1 = E_t0./60; %convert to 1/s
    
    F_t1 = cumtrapz(t,E_t1);
    
    ppE  = spline(t,E_t1); %we're making the RTD a continuos function of time
    ppF  = spline(t,F_t1);
    E_t  = @(t)(ppval(ppE,t));
    F_t  = @(t)(ppval(ppF,t));
    
    %Double check that your regression or interpolation fit the E & F curves well enough.
    tplot = linspace(0,max(t),500)';
    figure(1)
    plot(t,E_t1,'ro',tplot,E_t(tplot),'k-')
    xlabel('time')
    ylabel('E(t)')
    
    figure(2)
    plot(t,F_t1,'ro',tplot,F_t(tplot),'k-')
    xlabel('time')
    ylabel('F(t)')
    
    %Based on initial warning (below), if F(t) was maxed at 1 as usual, I
    %figured that the max time limit should be just before t=1.083184e+04 s.
    %Therefore tmax = 10831.83 s. Also see on F_t plot, as this time
    %should give F_t of 0.999999
    
    %"Warning: Failure at t=1.083184e+04.  Unable to meet integration tolerances without 
    %reducing the step size below the smallest value allowed (2.910383e-11) at time t. 
    % > In ode15s (line 730)
    %   In P1_e (line 53) "
    
    tmax = 10831.83; %s

    CAf  = 0.0313; %mol/L
    CBf  = 0.0313; %mol/L
    CCf  = 0;      %mol/L
    CDf  = 0;      %mol/L

    param.E  = E_t;
    param.F  = F_t;
    param.Cf = [CAf;CBf;CCf;CDf];

    
    %Set up ODE solver to integrate Maximum Mixedness model
    tspan = [0 tmax];
    lspan = tmax - tspan;
    C0    = [0.0313;0.0313;0;0];
    f     = @(l,C)P1_e0(l,C,param);
    %options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [l,C] = ode15s(f,lspan,C0);
        
    CA   = C(end,1);
    
    XAmm   = (CAf - CA)/CAf
  

end



function [D] = P1_e0(l,C,param)
    
    %This function returns the differential equations required to solve the
    %Maximum Mixedness model for a given RTD.
    
    %Assign concentrations as a function of lambda (time remaning in reactor)
    CA = C(1);
    CB = C(2);
    CC = C(3);
    CD = C(4);
    
    k1   = 175; %L^2/mol^2/s
    
    %Extract parameter values from structure
    E_t = param.E;
    F_t = param.F;
    CAf = param.Cf(1);
    CBf = param.Cf(2);
    CCf = param.Cf(3);
    CDf = param.Cf(4);
    
    %Calculate rates in each packet
    r1 = k1*CA*CB^2;
    
    %Calculate production rates Rj in each packet
    RA = -r1;
    RB = -r1;
    RC =  r1;
    RD =  r1;
    
    %Allocate space for ODEs
    D = zeros(size(C));    
    
    %Write balances for each species in space lambda + dlambda (i.e., MM model) for DE 1 to 4
    D(1) = E_t(l)/(1-F_t(l))*(CA - CAf) - RA;
    D(2) = E_t(l)/(1-F_t(l))*(CB - CBf) - RB;
    D(3) = E_t(l)/(1-F_t(l))*(CC - CCf) - RC;
    D(4) = E_t(l)/(1-F_t(l))*(CD - CDf) - RD;
    
end

