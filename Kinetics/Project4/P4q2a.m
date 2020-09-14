clear
%make a mesh of points
xmesh = linspace(0,0.25,52);
%put in some initial
guess = [1.2e-28; 0.0];
%initialize the bvp problem
solinit = bvpinit(xmesh, guess);
%https://www.mathworks.com/help/matlab/math/solve-bvp-with-singular-term.html
S = [0 0; 0 -2]; %need this matrix to hand the singularity

options = bvpset('SingularTerm',S);
sol = bvp4c(@bvpfcn, @bcfcn, solinit, options);
%sol.y(1,:) are the values of Ca
%sol.y(2,:) are the values of dCa/dr
plot(sol.x, sol.y(1,:), '-o')

%get the concentration at 0
c0 = sol.y(0,0);

%get the conversion of a




function residual = bcfcn(ya,yb) %boundary condition
    %note ya is the value at x = 0 and yb is the value at x = 0.25
    %the # in () denotes if it's the first defined or second defined in 
    residual = [ya(2)%BC at x = 0 | dc/dr = 0
                yb(1)-5.83e-5];%BC at x = 0.25 | c = 5.83e-5
       
end

function dcdr = bvpfcn(~,c)%Equation to solve
    
    CA = c(1);%concentration
    CA1 = c(2);%derivative of the concentration
    
    K1 = 90100; %cm^3/mol
    K2 = 6500; %cm^3/mol
    K4 = 64400; %cm^3/mol
    T  = 523; %K
    Da = 0.045; %cm^2/s
    CAs = 5.83e-5; %mol/cm^3
    CBs = 1.40e-4; %mol/cm^3
    CCs = 1.17e-5; %mok/cm^3
    k3 = 7.41e8; %g^2/mol/cm^3/s
    Cm = 1.8*10^-5; %mol/g
    
    %Calculate conversion
    xa = (CAs - CA)/CAs;

    %calculate the concentration of b and c
    CB = CBs - CAs*xa;

    CC = CCs + CAs*xa;
    
    %calculate the rate
    r = ((k3* K1*CA * ((K2*CB)^(1/2)) *Cm^2) / (1 + K1*CA + ((K2*CB)^(1/2)) + K4*CC));

    %gotta reverse the signs because it's going
    Ra = -r;

    %derivatives
    dcdr = zeros(2,1);
    dcdr(1) = CA1; 
    %The S vector will lead to -2CA1/r - RA/Da
    dcdr(2) = - Ra/Da; %note look at the documentation 
end