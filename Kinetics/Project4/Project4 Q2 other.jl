using DifferentialEquations
using Plots; gr()
using Interpolations

function interp1(Ypts, Xpts, Xvalue)
    knots = (Ypts,)
    itp = interpolate(knots, Xpts, Gridded(Linear()))
    itp(Xvalue)
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
    CAcenter = 1.259176e-28 #5.1621e-19
    dCAcenter = 0.0


    param = [K1, K2, K4, T, Da, k3, Cm, CAs, CBs, CCs]
    C0 = [CAcenter, dCAcenter, 8.17e-5, 0.0, 7e-5,0.0]
    span = (0.0, Rp)
    prob = ODEProblem(p4a1, C0, span, param)
    #prob = TwoPointBVProblem(p4a1, bc!, C0, span, param)
    sol = solve(prob, RKO65(), dt = 0.001, abstol = 1e-32, reltol = 1e-32)
    plot(sol, vars = (0,1))
end

function p4a1(du, u, p, t)
    CA, CA1, CB, CB1, CC, CC1= u
    K1, K2, K4, T, Da, k3, Cm, CAs, CBs, CCs = p

    #rate
    r = ((k3* K1*CA * ((K2*CB)^(1/2)) *Cm^2) / (1 + K1*CA + ((K2*CB)^(1/2)) + K4*CC))

    #gotta reverse the signs because it's going
    Ra = -r
    Rb = -r
    Rc = r

    #diffeqs
    du[1] = CA1
    du[2] = -Ra/Da - 2CA1/t
    du[3] = CB1
    du[4] = -Rb/Da - 2CB1/t
    du[5] = CC1
    du[6] = -Rc/Da - 2CC1/t
end
