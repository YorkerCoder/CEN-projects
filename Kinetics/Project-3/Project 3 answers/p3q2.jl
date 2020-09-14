using DifferentialEquations
using Plots; gr()
#This function will return a plot of temperature, and then a plot of concentrations, and finally the time where T > 300C
function p3q2()
    ΔHr = -5.9e5 #kcal/kmol
    cpa, cpb, cpw = 40.0, 8.38, 18.0 #kcal/kmol/K
    ua = 35.85 #kcal/min/K
    k0 = 1.4e-4 #m^3/kmol/min
    T0 = 461.0 #K
    Ea = 11273 #cal/mol
    Vl = 5.119 #m^3 volume of liquid
    na0, nb0, nw0 = 9.044, 33.0, 103.7 #kmol
    Ti = 175+273.15 #K initial temmperature of the reactor
    R = 1.9872 #cal/mol/K
    Ta = 25+273.15 #K
    #Initial Conditions
    IC = [na0, nb0, 0.0, 0.0, Ti]
    #Parameters
    params = [ΔHr, cpa, cpb, cpw, ua, k0, T0, Ea, Vl, na0, nb0, nw0, R, Ta]
    #time span
    tspan = (0.0, 200.0) #in minutes

    #set it up and solve it
    prob = ODEProblem(p3q2!, IC, tspan, params)
    sol  = solve(prob, Tsit5(), saveat = 0.1)
    Temp = plot(sol, vars = [(0,5)])
    Conc = plot(sol, vars = [(0,1), (0,2), (0,3), (0,4)])
    data = hcat(sol.u...)
    t = sol.t

    answer = Interp1d(data[5,:], t, 573.15)

    return(Temp, Conc, answer)
end

function p3q2!(du, u, p, t)
    na, nb, nc, nd, T = u
    ΔHr, cpa, cpb, cpw, ua, k0, T0, Ea, Vl, na0, nb0, nw0, R, Ta = p

    #calculate change in enthalpy
    ΔH = ΔHr #assume ΔCp is 0

    #calculate the new k
    k = k0*exp(-(Ea/R) * ((1/T) - (1/T0)) )

    #get the concentration of a and b
    ca, cb = [na, nb]/Vl

    #get the rate
    r = k*ca*cb

    #get production rate of A
    Ra, Rb, Rc, Rd = -r, -2r, r, r

    #write the differential equations
    du[1] = dna = Ra*Vl
    du[2] = dnb = Rb*Vl
    du[3] = dnc = Rc*Vl
    du[4] = dnd = Rd*Vl

    #change the heat balance based on the time
    if t < 45.0
        du[5] = dT = 0
    elseif t>=45.0 && t<=55.0
        du[5] = dT = (-ΔHr*r*(Vl))/(na*cpa+nb*cpb+nw0*cpw)
    else
        #the volume of the liquid is ratiod to the full volume of the reactor
        du[5] = dT = (-ΔHr*r*Vl + ua*(Vl/(3*3.26))*(Ta-T))/(na*cpa+nb*cpb+nw0*cpw)
    end


end

#Made this to simulate the results of Interp1d in Matlab
function Interp1d(Data, range, solution)
    count = 1;
    δ = 0.0
    for point in Data
        δ = solution - point
        if δ < 1e-9
            break
        end
        count += 1
    end
    return range[count]
end
