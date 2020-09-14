
%1a
%Ideal PFR 
%The conversion of A, i.e. XApfr = 0.9642

function [] = P1_a

    V = 1000 ; %L
    Q = 10   ; %L/s
    tau = (V/Q);  %s
    t_1a = load('P1_time.txt');
    t = 60.*t_1a; %convert to sec
    
    E_t0 = load('P1_E(t).txt');
    E_t = E_t0./60; %convert to 1/s

    
    tm_res = trapz(t,t.*E_t);


    CAf  = 0.0313; %mol/L
    CBf  = 0.0313; %mol/L
    CCf  = 0;      %mol/L
    CDf  = 0;      %mol/L


    %Set up ODE solver to integrate the Ideal CSTR
    tspan   = [0 tm_res];
    C0      = zeros(4,1);
    C0(1:2) = 0.0313;
    f       = @(t,C)P1_a0(t,C);
    [t,C]   = ode15s(f,tspan,C0);

    CA   = C(end,1);

    XApfr    = (CAf - CA)/CAf


end



function [D] = P1_a0(t,C)

    %This function returns the differential equations required to solve the
    %segregation model for a given RTD.
    
    %Assign "packet" concentrations for individual batch reactors
    CA = C(1);
    CB = C(2);
    CC = C(3);
    CD = C(4);
    
    
    k1   = 175;    %L^2/mol^2/s
    
    r1 = k1*CA*CB^2;
    
    RA = -r1;
    RB = -r1;
    RC =  r1;
    RD =  r1;
    
    %Allocate space for ODEs
    D = zeros(4,1);    
    
    %Calculate average exit concentrations using RTD
    D(1) = RA; %dCj/dtau
    D(2) = RB;
    D(3) = RC;
    D(4) = RD;
    
end

