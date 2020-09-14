
%1b
%Ideal CSTR
%The conversion of A, i.e. XAcstr = 0.6620

function [] = P1_b

    V = 1000 ; %L
    Q = 10   ; %L/s
    tau = (V/Q);  %s
    t_1b = load('P1_time.txt');
    t = 60.*t_1b; %convert to sec
    
    E_t0 = load('P1_E(t).txt');
    E_t = E_t0./60; %convert to 1/s


    tm_res = trapz(t,t.*E_t);


    CAf  = 0.0313; %mol/L
    CBf  = 0.0313; %mol/L
    CCf  = 0;      %mol/L
    CDf  = 0;      %mol/L

    param.tau = tau;
    param.Cf  = [CAf;CBf;CCf;CDf];


    %Set up ODE solver to integrate the Ideal CSTR
    tspan   = [0 tm_res];
    C0      = zeros(4,1);
    C0(1:2) = 0.0313;
    f       = @(t,C)P1_b0(t,C,param);
    [t,C]   = ode15s(f,tspan,C0);

    CA   = C(end,1);

    XAcstr    = (CAf - CA)/CAf


end



function [D] = P1_b0(t,C,param)

    %This function returns the differential equations required to solve the
    %segregation model for a given RTD.
    
    %Assign "packet" concentrations for individual batch reactors
    CA = C(1);
    CB = C(2);
    CC = C(3);
    CD = C(4);
    
    %Extract parameter values from structure
    tau = param.tau;
    CAf = param.Cf(1);
    CBf = param.Cf(2);
    CCf = param.Cf(3);
    CDf = param.Cf(4);
    
    k1   = 175;    %L^2/mol^2/s
    
    r1 = k1*CA*CB^2;
    
    RA = -r1;
    RB = -r1;
    RC =  r1;
    RD =  r1;
    
    %Allocate space for ODEs
    D = zeros(4,1);    
    
    %Calculate average exit concentrations using RTD
    D(1) = RA + CAf/tau - CA/tau; %dCj/dt
    D(2) = RB + CBf/tau - CB/tau;
    D(3) = RC + CCf/tau - CC/tau;
    D(4) = RD + CDf/tau - CD/tau;
    
end

