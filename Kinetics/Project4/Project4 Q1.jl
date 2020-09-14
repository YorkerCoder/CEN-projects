#===============================================================================
a) 314 kg of catalyst
b) Yes the mass of the catalyst needs to be changed
c) 553 kg of catalyst
#when your run the questions they show the value of phi to show
# that the assumptions for the equations are valid
p4a0 answers a, and p4b0 answer both b and c
===============================================================================#



using DifferentialEquations
using Plots; gr()
using Interpolations

function interp1(Xpts, Ypts, Xvalue)
    itp = LinearInterpolation(Ypts, Xpts)
    itp(Xvalue)
end

function p4a0()
    T = 650 #K
    P = 2.0 #atm
    R = 8.206e-5 #m3*atm/mol/K
    k = 4.54e7/(100^3) #m^3/mol/s
    n = 2
    Rp = 0.635e-2/2 #m
    vbar = Rp/3
    a = 1/vbar
    DA = 9.5e-7 #m^2/sec
    rho_c = 1.75 #g/cm^3
    rho_b = 0.84 #g/cm^3
    Q = 1000*(1/1000) #L/s * (m^3/1000L) = m^3/s

    Fa0 = P*Q/R/T

    parameters = [T, P, R, k, vbar, a, DA, n, Rp]
    F0 = [Fa0, 0.0]
    Vspan = (0.0, 0.25)

    prob = ODEProblem(p4a1, F0, Vspan, parameters, saveat = 0.001)
    sol = solve(prob, Rodas4(), saveat = 0.00001)

    solu = convert(Array,sol)



    xa = (F0[1] .- solu[1,:])./F0[1]
    FT = solu[1,:] .+ solu[2,:]
    Q = FT .* R .* T ./P
    CA = solu[1,:] ./ Q
    PHI = ((n+1)./2 .* k .* CA .* vbar^2 ./ DA).^(1/2)
    g1 = plot(sol.t, xa)
    g2 = plot(sol.t, PHI)
    vcat = interp1(sol.t, xa, 0.93)
    mcat = vcat*rho_c*100^3/1000 #in kg

    return (mcat, vcat, g1, g2, PHI[end])#the PHI is to show that the assumption is safe
end

function p4a1(du, u, p, t)
    FA, FB = u
    T, P, R, k, vbar, a, DA, n, Rp = p

    FTOT = FA + FB
    Q   = FTOT*R*T/P
    CA  = FA/Q
    CAs = CA

    PHI = ((n+1)/2 * k * CAs * vbar^2 / DA)^(1/2)

    eta = 1/PHI


    RS  = k*CAs^2 #rate at the surfaceo f the particle
    RP  = eta*RS  #average production rate of the particle
    RA  = -RP     #production rate of A
    RB  =  RP     #production rate of B


    du[1] = RA
    du[2] = RB
end

#note the only differences here are some constants, and that vbar = Rp/2 for a cylinder
function p4b0()
    T = 650 #K
    P = 2.0 #atm
    R = 8.206e-5 #m3*atm/mol/K
    k = 4.54e7/(100^3) #m^3/mol/s
    n = 2
    Rp = 0.635e-2/2 #m
    L = 8.90e-3 #m
    vbar = Rp/2
    a = 1/vbar
    DA = 7.2e-7 #m^2/sec
    rho_c = 1.79 #g/cm^3
    rho_b = 0.96 #g/cm^3
    Q = 1000*(1/1000) #L/s * (m^3/1000L) = m^3/s

    Fa0 = P*Q/R/T

    parameters = [T, P, R, k, vbar, a, DA, n, Rp, Q]
    F0 = [Fa0, 0.0]
    Vspan = (0.0, 0.50)

    prob = ODEProblem(p4b1, F0, Vspan, parameters, saveat = 0.001)
    sol = solve(prob, Rodas4(), saveat = 0.0001)

    solu = convert(Array,sol)



    xa = (F0[1] .- solu[1,:])./F0[1]
    FT = solu[1,:] .+ solu[2,:]
    Q = FT .* R .* T ./P
    CA = solu[1,:] ./ Q
    PHI = ((n+1)./2 .* k .* CA .* vbar^2 ./ DA).^(1/2)
    g1 = plot(sol.t, xa)
    g2 = plot(sol.t, PHI)
    vcat = interp1(sol.t, xa, 0.93)
    mcat = vcat*rho_c*(100^3)/1000 #in kg

    return (mcat, vcat, g1, g2, PHI[end])#the PHI is to show that the assumption is safe

end

function p4b1(du, u, p, t)
    FA, FB = u
    T, P, R, k, vbar, a, DA, n, Rp, Q = p



    CAs = FA/Q
    PHI = ((n+1)/2 * k * CAs * vbar^2 / DA)^(1/2)

    eta = 1/PHI


    RS  = k*CAs^2 #rate at the surfaceo f the particle
    RP  = eta*RS  #average production rate of the particle
    RA  = -RP     #production rate of A
    RB  =  RP     #production rate of B


    du[1] = RA
    du[2] = RB
end
