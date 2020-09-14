#===============================================================================
#Answer to question a) Ca Center = 1.2591759e-28
#Answer to question b) eta = 0.02571107
===============================================================================#
using DifferentialEquations
using Plots; gr()
using Interpolations
using Optim
using QuadGK

#function I made to interpolate easily
function interp1(Xpts, Ypts, Xvalue)
    itp = LinearInterpolation(Ypts, Xpts)
    itp(Xvalue)
end

#function i made to numerically integrate easily
function NumInt(Xpts, Ypts, lbound, ubound ; reltol = 1e-15)
    itpc = CubicSplineInterpolation(Xpts, Ypts)
    integral, error = quadgk(x -> itpc(x), lbound, ubound, rtol = reltol)
    return (integral, error)
end


function p4a0()
    K1 = 90100 #cm^3/mol
    K2 = 6500 #cm^3/mol
    K4 = 64400 #cm^3/mol
    T  = 523 #K
    Da = 0.045 #cm^2/s
    CAs = 5.83e-5 #mol/cm^3
    CBs = 1.40e-4 #mol/cm^3
    CCs = 1.17e-5 #mok/cm^3
    k3 = 7.41e8 #g^2/mol/cm^3/s
    Cm = 1.8*10^-5 #mol/g
    Rp = 0.25 #cm


    #guessed concentrations at r = 0
    CAcenter = 1.2591759e-28 #5.1621e-19
    dCAcenter = 0.0


    param = [K1, K2, K4, T, Da, k3, Cm, CAs, CBs, CCs]
    C0 = [CAcenter, dCAcenter]
    span = (0.0, Rp)
    prob = ODEProblem(p4a1, C0, span, param)
    #prob = TwoPointBVProblem(p4a1, bc!, C0, span, param)
    sol = solve(prob, RKO65(), dt = 0.001, abstol = 1e-30, reltol = 1e-30)

    #get the production rate of A at the surface
    Ras = (k3* K1*CAs * ((K2*CBs)^(1/2)) *Cm^2) / (1 + K1*CAs + ((K2*CBs)^(1/2)) + K4*CCs)

    #get the particles average production rate
    #generate a radius
    r = 0:Rp/1000:Rp
    #convert Ca to Ca bar
    C = sol(r)
    CA = C[1,:]
    CAbar = CA/CAs
    #generate r bar
    rbar = r/(Rp/3)

    #generate the part that gets integrated
    foo = CAbar.^2 .* rbar.^2
    #do the numerical integration
    foo1 = NumInt(rbar, foo, 0.0, 3.0, reltol = 1e-30)#note returns a tuple(integral value and error)
    #get eta
    eta = foo1[1]/9

    #return the effectiveness and center concentration
    return(CAcenter, eta)
end

function p4a1(du, u, p, t)
    CA, CA1 = u
    K1, K2, K4, T, Da, k3, Cm, CAs, CBs, CCs = p

    #Calculate conversion
    xa = (CAs - CA)/CAs

    #calculate the concentration of b and c
    CB = CBs - CAs*xa

    CC = CCs + CAs*xa

    r = ((k3* K1*CA * ((K2*CB)^(1/2)) *Cm^2) / (1 + K1*CA + ((K2*CB)^(1/2)) + K4*CC))

    #gotta reverse the signs because it's going
    Ra = -r

    #diffeqs
    du[1] = CA1
    du[2] = -Ra/Da - 2CA1/t
end
