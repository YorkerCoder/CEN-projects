
%1c
%The conversion of A, i.e. XAcs = 0.9908

function [] = P1_c

    V = 1000 ; %L
    Q = 10   ; %L/s
    tau  = V/Q;  %s
    t_1c = load('P1_time.txt');
    t = 60.*t_1c; %convert to sec

    E_t1 = zeros(length(t),1);
    F_t1 = zeros(length(t),1);


for i = 1:length(t)
    
    if t(i) < tau/2
        E_t1(i) = 0; %1/s
        F_t1(i) = 0;
    else
        E_t1(i) = (tau^2)/(2*t(i).^3); %1/s
        F_t1(i) = 1 - 1/(4 * (t(i)./tau).^2);
    end
    
end
    
    ppE  = spline(t,E_t1); %we're making the RTD a continuos function of time
    ppF  = spline(t,F_t1);
    E_t  = @(t)(ppval(ppE,t));
    F_t  = @(t)(ppval(ppF,t));
    
    %Double check that your regression or interpolation fit the E curve well enough.
    tplot = linspace(0,max(t),500)';
    figure(1)
    plot(t,E_t1,'ro',tplot,E_t(tplot),'k-')
    xlabel('time')
    ylabel('E(t)')
    
    figure(2)
    plot(t,F_t1,'ro',tplot,F_t(tplot),'k-')
    xlabel('time')
    ylabel('F(t)')
    
    
    tmax = max(t);

    CAf  = 0.0313; %mol/L
    CBf  = 0.0313; %mol/L
    CCf  = 0;      %mol/L
    CDf  = 0;      %mol/L

    param.E  = E_t;
    param.Cf = [CAf;CBf;CCf;CDf];


    %Set up ODE solver to integrate Complete Segregation model
    tspan   = [0 tmax];
    C0      = zeros(8,1);
    C0(1:2) = 0.0313;
    f       = @(t,C)P1_c0(t,C,param);
    [t,C]   = ode15s(f,tspan,C0);

    CAbar   = C(end,5);

    XAcs    = (CAf - CAbar)/CAf


end



function [D] = P1_c0(t,C,param)

    %This function returns the differential equations required to solve the
    %segregation model for a given RTD.
    
    %Assign "packet" concentrations for individual batch reactors
    CA = C(1);
    CB = C(2);
    CC = C(3);
    CD = C(4);
    
    %Assign "average" concentrations at reactor exit
    CAbar = C(5);
    CBbar = C(6);
    CCbar = C(7);
    CDbar = C(8);
    
    %Extract parameter values from structure
    E_t = param.E;
    
    k1   = 175;    %L^2/mol^2/s
    
    %Calculate rates in each packet
    r1 = k1*CA*CB^2;
    
    %Calculate production rates Rj in each packet
    RA = -r1;
    RB = -r1;
    RC =  r1;
    RD =  r1;
    
    %Allocate space for ODEs
    D = zeros(size(C));    
    
    %Write balances for each packet (as batch reactor)
    D(1) = RA;
    D(2) = RB;
    D(3) = RC;
    D(4) = RD;
    %Calculate average exit concentrations using RTD
    D(5) = CA*E_t(t);
    D(6) = CB*E_t(t);
    D(7) = CC*E_t(t);
    D(8) = CD*E_t(t);
    
end

