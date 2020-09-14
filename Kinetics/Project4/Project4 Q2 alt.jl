#p4ab0 returns the answers to both a and b which are 0.049, and 1.259 respectively
using DifferentialEquations
using Plots; gr()
using Interpolations
using Optim

function interp1(Ypts, Xpts, Xvalue)
    knots = (Ypts,)
    itp = interpolate(knots, Xpts, Gridded(Linear()))
    itp(Xvalue)
end

function p4ab0()
    K1 = 90100/(100^3) #cm^3/mol
    K2 = 6500/(100^3) #cm^3/mol
    K4 = 64400/(100^3) #cm^3/mol
    T  = 523 #K
    Da = 0.045/(100^2) #cm^2/s
    CAs = 5.83e-5*(100^3) #mol/cm^3
    CBs = 1.40e-4*(100^3) #mol/cm^3
    CCs = 1.17e-5*(100^3) #mok/cm^3
    k3 = 7.41e8*(100^3) #g^2/mol/cm^3/s
    Cm = 1.8*10^-5 #mol/g
    Rp = 0.25*100 #cm



    #guessed concentrations at r = 0
    CAcenter = 1.259176e-15 #5.1621e-19
    dCAcenter = 0.0


    param = [K1, K2, K4, T, Da, k3, Cm, CAs, CBs, CCs]
    C0 = [CAcenter, dCAcenter]
    span = (0.0, Rp)
    prob = ODEProblem(p4a1, C0, span, param)
    #prob = TwoPointBVProblem(p4a1, bc!, C0, span, param)
    sol = solve(prob, RKO65(), dt = 0.001, abstol = 1e-16, reltol = 1e-16)

    #get the production rate of A at the surface
    Ras = (k3* K1*CAs * ((K2*CBs)^(1/2)) *Cm^2) / (1 + K1*CAs + ((K2*CBs)^(1/2)) + K4*CCs)

    #get the particles average production rate
    vbar = Rp/3
    a = 1/vbar
    Rap = a*Da*sol[end][2]

    eta = Rap/Ras

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
