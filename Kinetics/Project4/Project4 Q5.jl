#=======================================================
# The answer is 201 grams 
# run function p4a0() for this
=======================================================#

using DifferentialEquations
using Interpolations


#I gotta put this in a file somewhere to call to make it easier to read
function interp1(Xpts, Ypts, Xvalue)
    itp = LinearInterpolation(Ypts, Xpts)
    itp(Xvalue)
end

function p4a0()
    #declare all of the constants
    FA0 = 15.0  #mol/s
    P   = 3.0 #atm
    T   = 550.0 #K
    Da  = 0.0008 # cm^2/s
    R   = 8.205736e-5 * (100^3) #cm^3*atm/K/mol
    #Q = Qf
    Qf = FA0*R*T/P #cm^3/s
    CA0 = FA0/Qf #mol/cm^3
    k1 = 1.21e5 #cm^3/mol/s
    Rp = 1.5#cm
    vbar = Rp/2
    rho_c = 0.73 #g/cm^3

    #calculate the fraction of volume that is in the hollow cylinder compared to a full one
    BigV = pi*1.3^2*0.7
    SmallV = pi*0.5^2*0.7
    realV =  BigV-SmallV
    FractionalV = realV/BigV

    #solve the ODE 
    params = [Da, k1, Qf, vbar, FractionalV]
    C0 = [FA0]
    span = (0.0, 325.0)#span of the catalyst
    prob = ODEProblem(p4a1, C0, span, params)
    sol = solve(prob, saveat = 0.5)

    #set everything up for the interpolation function 
    #had to reverse it because it only takes things going in ascending order for the y axis
    fa = convert(Array,sol)
    fa = reverse(fa[1,:])
    vcat = interp1(reverse(sol.t), fa, 1.5)

    #calculate the mass of catalyst needed
    mcat = vcat*rho_c
end

function p4a1(du,u,p,t)
    FA = u[1]
    Da, k1, Qf, vbar, FractionalV = p

    #calculate concentration of the bulk 
    CA = FA/Qf

    #no Boundary layer over the particle so bulk = surface concentration
    CAs = CA

    #calculate phi 
    phi = sqrt(1.5*k1*CAs*vbar^2/Da)

    
    #calculate eta
    eta = 1/phi * FractionalV

    #calculate the rate
    r = k1*CAs*eta 

    #calculate the production rate of A 
    RA = -r

    #solve the DifferentialEquations
    du[1] = RA
end

answwer = p4a0()

#This was if i wanted to solve into the particle without and approximation
#but it leads to severely increased times
#= 
#function to parameterize to solve for
function p4b0()
    FA0 = 15  #mol/s
    P   = 3.0 #atm
    T   = 550 #K
    Da  = 0.0008 # cm^2/s
    Router = 1.5 #cm
    Rinner = 0.7 #cm
    R = 8.205736e-5 * (100^3) #cm^3*atm/K/mol
    #Q = Qf
    Qf = FA0*R*T/P #cm^3/s
    CA0 = FA0/Qf #mol/cm^3
    k1 = 1.21e5 #cm^3/mol/s

    params = [Da, Rinner, k1]
    C0 = [1.74097911e-7, 0.0]
    span = (0.0, 1.3)
    prob = ODEProblem(p4b1,C0,span, params)
    sol = solve(prob, RKO65(), dt = 0.000001, abstol = 1e-7, reltol = 1e-7)


end

#for solving into the cylinder
function p4b1(du,u,p,t)
    #variables
    CA, C1A = u
    #parameters
    Da, Rinner, k1 = p

    #reaction rate
    r = k1*CA^2

    #production rate of A
    RA = -r

    if t < Rinner
        du[1] = 0
        du[2] = 0
    else
        du[1] = C1A
        du[2] = -RA/Da - C1A/t
    end

end
=#
