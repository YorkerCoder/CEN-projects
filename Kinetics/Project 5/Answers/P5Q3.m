clear
syms t t1

%get the RTD as a function
E(t) = 0.1 -(0.1*t)/t1;

%a)
%get mean residence time
tau = int(E*t, t, 0, t1)


%b)
%get standard deviation
stdDev = int(E*(t-tau)^2,0,t1);

%get the ratio
ratio = stdDev/(tau^2)

%From here I would solve for it using fsolve, but since i'm not given
%t1 I can't, and I can't solve for it by hand because of the exponential
%mixed with the algabraic terms ie:
% ratio = 2/Pe + 2*(1-exp(-Pe))/Pe^2

%Again how many tanks would be given by t1, but will follow this formula
% the equation for this is just the inverse of the ratio
% that is: ntanks = tau^2/stdDev

ntanks = tau^2/stdDev %note you would round all of these up