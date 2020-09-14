#===============================================================================
#a) 124kg of catalyst p4a0()
#b) 124kg of catalyst -> the mass transfer coefficient is super large p4b0()
===============================================================================#

using DifferentialEquations
using Plots; gr()
using Interpolations
using SpecialFunctions

function interp1(Xpts, Ypts, Xvalue)
    itp = LinearInterpolation(Ypts, Xpts)
    itp(Xvalue)
end

function p4a0()
    #set up the constants
    R   = 1.987204 #cal/K/mol
    R1  = 8.205736e-5 * (100^3) #cm^3*atm/K/mol
    P0  = 1.0      #atm
    T   = 630      #K
    FA0 = 0.5      #mol/s
    Rp = 0.35      #cm
    L = 0.5        #cm
    rho_c = 0.84/1000 #kg/cm^3
    rho_b = 0.52 #g/cm^3
    Da = 1.40e-3 #cm^2/s
    k = 5.74e13 * exp(-38000/(R*T)) #1/s
    xa = 0.90
    Qf = FA0*R1*T/P0 #cm^3/s

    #solve for eta
    vbar = Rp/2
    phi = (k*vbar^2/Da)^0.5
    eta = besseli(1, 2.0phi)/besseli(0,2.0phi)/phi

    #Isothermal and no pressure drop

    #solve the ODE
    param = [k,eta,Qf]
    F0 = [FA0]
    span = (0.0, 2.5e5)
    prob = ODEProblem(p4a1, F0, span, param)
    sol = solve(prob, saveat=0.5)

    #set up the span
    #note had to reverse them because they y points need to be in increasing order to begin with
    fa = convert(Array,sol)
    fa = reverse(fa[1,:])
    vcat = interp1(reverse(sol.t),fa,0.05)
    mcat = rho_c*vcat

end

function p4a1(du,u,p,t)
    #declare the constants
    k, eta, Qf = p

    #solve diffeq
    du[1] = -k*eta*u[1]/Qf

end

function p4b0()
    #set up the constants
    R   = 1.987204 #cal/K/mol
    R1  = 8.205736e-5 * (100^3) #cm^3*atm/K/mol
    P0  = 1.0      #atm
    T   = 630      #K
    FA0 = 0.5      #mol/s
    Rp = 0.35    #cm
    L = 0.5     #cm
    rho_c = 0.84/1000 #kg/cm^3
    rho_b = 0.52 #g/cm^3
    Da = 1.40e-3 #cm^2/s
    k = 5.74e13 * exp(-38000/(R*T)) #1/s
    xa = 0.90
    Qf = FA0*R1*T/P0 #cm^3/s
    vbar = Rp/2
    phi = (k*vbar^2/Da)^0.5
    eta = besseli(1, 2.0phi)/besseli(0,2.0phi)/phi
    kc = 375 #cm/s
    a = 1/vbar

    #solve the ODE
    param = [k, eta, Qf, kc, a]
    F0 = [FA0, 0.0]#the second value is a guess for the algebraic part
    span = (0.0, 2.5e5)#in cm

    #mass matrix
    m = [1 0
         0 0]
    f = ODEFunction(p4b1, mass_matrix = m)
    prob = ODEProblem(f, F0, span, param)
    sol = solve(prob, Rosenbrock23(), saveat=1)

    #make a range, and concentration vector to be able to use interpolations
    v = reverse([0.0:1:250000 ...])#making a span of vectors and reversing them
    fa = convert(Array,sol)
    fa = reverse(fa[1,:])
    vcat = interp1(v,fa,0.05)#0.05 translates to 90% conversion of a
    mcat = rho_c*vcat
    #plot(sol, vars =(0,1))
end

function p4b1(du, u, p, t)
    k, eta, Qf, kc, a = p
    FA, CAs = u #we want the surface concentration so we choose it to be the variable we solve for

    #concentration in the bulk
    CA = FA/Qf

    #production rate at the surface(first order)
    Rsurface = -k*CAs

    #average production rate of the particle
    Rparticle = eta*Rsurface

    #production rate
    RA = Rparticle

    #note that this is solving for A as a function of catalyst volume
    du[1] = RA
    du[2] = kc*a*(CA-CAs)-k*CAs*eta
end
