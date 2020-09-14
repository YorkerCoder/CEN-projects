clear
%constants for water and Acetone Antoine Equation
%Note only work between 273 and 350K
A = [16.3872 14.3145];
B = [3885.70 2756.22];
C = [-42.980 -45.090 ];

P = 101; %kPa
x1 = 0.0:0.01:1.0;
Temps = zeros(length(x1),1);

%Solve for each composition
for i = 1:length(x1)
    fun = @(T) eqs(T,x1(i),A,B,C,P);
    Temps(i) = fsolve(fun, [300]);
end

%calculate the gas phase compositions
y1 = zeros(length(x1),1);
y2 = zeros(length(x1),1);
x2 = 1 - x1;
for i = 1:length(x1)
    %get the antoine equation values 
    ps1 = Antoine(A(1),B(1),C(1),Temps(i));
    ps2 = Antoine(A(2),B(2),C(2),Temps(i));
    
    gam = activ([x1(i), 1-x1(i)], Temps(i)) ;
    %calculate the y value using modified raoults law
    y1(i) = x1(i) *gam(1) * ps1 / P;
    y2(i) = x2(i) *gam(2) * ps2 / P;
end

%water
figure
plot(x1,Temps , y1,Temps)
title('phase diagram of water')
xlabel('x1/y1')
ylabel('T(K)')
legend('liquid phase','vapor phase')

%acetone
figure 
plot(x2,Temps , y2,Temps)
title('phase diagram of acetone')
xlabel('x1/y1')
ylabel('T(K)')
legend('liquid phase','vapor phase')


%%
%Functions are defined below here
%The equation that gets solved for bubble point
function F = eqs(T, x, A, B, C, P)
    gams = activ([x(1), 1-x(1)], T);
    F(1) = x(1) * Antoine(A(1),B(1),C(1),T) * gams(1) + (1-x(1)) * Antoine(A(2),B(2),C(2),T) * gams(2) - P;
end

%A function to calculate vapor pressure using the antoine equation 
function Ps = Antoine(A,B,C,T)
    Ps = exp(A - ( B / ( T+C ) ) );
end

%%
%Unifaq function
function gam = activ(x, T)
   % UNIFAC method for estimating activity coefficients
    %% R and Q values from Table H.1 (Smith, Van Ness, Abbott)
    R(1:42)=0; Q(1:42)=0;
    
    R(1)=0.90110; R(2)=0.67440; R(3)=0.44690; R(4)=0.21950;
    Q(1)=0.84800; Q(2)=0.54000; Q(3)=0.22800; Q(4)=0.21950;
    R(10)=0.5313; R(12)=1.2663; R(13)=1.0396; R(15)=1.0000;
    Q(10)=0.4000; Q(12)=0.9680; Q(13)=0.6600; Q(15)=1.2000;
    R(17)=0.9200; R(19)=1.6724; R(20)=1.4457; R(25)=1.0880;
    Q(17)=1.4000; Q(19)=1.4880; Q(20)=1.1800; Q(25)=1.0880;
    R(26)=0.9183; R(27)=0.6908; R(32)=1.4337; R(33)=1.2070;
    Q(26)=0.7800; Q(27)=0.4680; Q(32)=1.2440; Q(33)=0.9360;
    R(34)=0.9795; R(41)=1.8701; R(42)=1.6434;
    Q(34)=0.6240; Q(41)=1.7240; Q(42)=1.4160;
    
    %%a_mk (in K) UNIFAC-VLE Interaction parameters from Table H.2
    amk=zeros(42,42);
    amk(1:4,1:4)=0; amk(1:4,10)=61.13; amk(1:4,12:13)=76.50; amk(1:4,15)=986.5; amk(1:4,17)=1318; amk(1:4,19:20)=476.4; amk(1:4,25:27)=251.5; amk(1:4,32:34)=255.7; amk(1:4,41:42)=597;
    amk(10,1:4)=-11.12;
    amk(15,1:4)=156.4; amk(15,10)=89.6; amk(15,12:13)=28.82; amk(15,15)=0; amk(15,17)=353.5; amk(15,19:20)=84; amk(15,25:27)=28.06; amk(15,32:34)=42.7; amk(15,41:42)=6.712;
    amk(17,1:4)=300; amk(17,10)=362.3; amk(17,12:13)=377.6; amk(17,15)=-229.1; amk(17,17)=0; amk(17,19:20)=-195.4; amk(17,25:27)=540.5; amk(17,32:34)=168; amk(17,41:42)=112.6;
    amk(19:20,1:4)=26.76; amk(19:20,10)=140.1; amk(19:20,12:13)=365.8; amk(19:20,15)=164.5; amk(19:20,17)=472.5; amk(19:20,19:20)=0; amk(19:20,25:27)=-103.6; amk(19:20,32:34)=-174.2; amk(19:20,41:42)=481.7;
    amk(25:27,1:4)=83.36; amk(25:27,10)=52.13; amk(25:27,12:13)=65.69; amk(25:27,15)=237.7; amk(25:27,17)=-314.7; amk(25:27,19:20)=191.1; amk(25:27,25:27)=0; amk(25:27,32:34)=251.5; amk(25:27,41:42)=-18.51;
    amk(32:34,1:4)=65.33; amk(32:34,10)=-22.31; amk(32:34,12:13)=223; amk(32:34,15)=-150; amk(32:34,17)=-448.2; amk(32:34,19:20)=394.6; amk(32:34,25:27)=-56.08; amk(32:34,32:34)=0; amk(32:34,41:42)=147.1;
    amk(41:42,1:4)=24.82; amk(41:42,10)=-22.97; amk(41:42,12:13)=-138.4; amk(41:42,15)=185.4; amk(41:42,17)=242.8; amk(41:42,19:20)=-287.5; amk(41:42,25:27)=38.81; amk(41:42,32:34)=-108.5; amk(41:42,41:42)=0;

    %% Specify mixture
    N=2; % number of species i=1: Acetone (CH3CO CH3)(1); i=2: Water (H2O(2))
    nu(1:N,1:42)=0;
    nu(1,1)=1; nu(1,19)=1;
    nu(2,17)=1;
    %% Calculate activity coefs

    %fixed all of these by transposing nu, q, and r
    r=R*nu';
    q=Q*nu';
    tau=exp(-amk./T);
    xjqj=x*q';
    rjxj=x*r';

    %calculate eki
    eki = nu.*Q./q';

    betaik=eki*tau;

    %calculate theta using broadcasting and matrix multiplication 
    theta = x.*q*eki/xjqj;

    %calculate sk using matrix multiplication
    sk = theta*tau;

    J=r./rjxj;
    L=q./xjqj;
    
    %Calculate ln(gammaC)
    gamlnc = 1 - J + log(J) - 5*q.*(1- (J./L) + log(J./L));
    
    %Calculate ln(gammaR)
    A = betaik./sk;
    B = eki.*log(A);
    gamlnr = q'.*(1- (A*theta' - sum(B,2)));
    
    %Calculate the activity coeficcients 
    gam =exp(gamlnc + gamlnr');
end
