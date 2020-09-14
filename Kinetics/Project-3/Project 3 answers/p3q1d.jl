using DifferentialEquations
using Plots; gr()

#around .6hours it hits steady state
function p32()
    ΔH = - 36000 #BTU/lb-mol
    cpa, cpb, cpc, cpm = 35, 18, 46, 19.5 #BTU/lb-mol/R
    va, vb, vm = 1.083, 0.290, 0.649 #ft^3 / mol
    ua = 16000 #BTU/h/R
    Vr = 500/7.48052 #ft^3
    cb0 = 3.45 #lb-mol/ft^3
    T0 = 75 + 459.67 #R
    faf, fbf, fmf = 80, 1000, 100 #lb-mol/h
    Ta = 60 + 459.67 #R
    R = 1.9859 #BTU/lb-mol/R
    Qf = faf*va + fbf*vb + fmf*vm #ft^3/h
    tau = Vr/Qf #h

    #calculate feed concentrations
    caf, cbf, ccf, cmf = [faf, fbf, 0.0, fmf]/Qf

    parameters = [ΔH, cpa, cpb, cpc, cpm, va, vb, vm, ua, Vr, T0,
                  faf, fbf, fmf, Ta, R, Qf, tau, caf, cbf, ccf, cmf]

    initial_conditions = [caf, cbf+cb0, 0.0, T0]

    tspan = (0.0, 1.0)
    prob = ODEProblem(p32!, initial_conditions, tspan, parameters)
    sol = solve(prob, Tsit5())
    plot(sol, vars = [(0,1), (0,3)])
end


function p32!(du, u, p, t)
    ca, cb, cc, T = u
    ΔH, cpa, cpb, cpc, cpm, va, vb, vm, ua, Vr, T0, faf, fbf, fmf, Ta, R, Qf, tau, caf, cbf, ccf, cmf = p

    #set the molar flow rate of methanol
    fm = fmf

    #caluclate k
    k = 16.96e12 * exp(-32400.0/R/T)

    #calculate ΔH of reaction
    ΔHr = ΔH + (cpc - cpa - cpb)*(T - 536.73)

    #calculate reaction rate
    r = k*ca

    #calculate production rates
    Ra, Rb, Rc = [-r, -r, r]

    #differential equations
    du[1] = dCa = ((caf - ca)/tau) + Ra
    du[2] = dCb = ((cbf - cb)/tau) + Rb
    du[3] = dCc = ((ccf - cc)/tau) + Rc
    du[4] = dT  = (-ΔHr*r*Vr +(faf*cpa + fbf*cpb + fmf*cpm)*(T0 - T)+ua*(Ta-T))/(ca*cpa + cb*cpb + cc*cpc + cmf*cpm)/Vr

end
