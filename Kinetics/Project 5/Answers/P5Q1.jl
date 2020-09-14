#================================================================
#p5a0 is the answer to a = 0.9844
#p5b0 is the answer to b = 0.8689
#p5c0 is the answer to c = 0.9933
#p5d0 is the answer to d = 0.9379 
#p5e0 is the answer to e = 0.8673 -> note due to the solvers not being able to handle it I had to cut this on off early
================================================================#
#Just run the whole file at once and the answers will print out
using DifferentialEquations
using Plots
using DelimitedFiles
using Interpolations, QuadGK

#function i made to numerically integrate data
function numint(Xpts, Ypts, lbound, ubound ; reltol = 1e-15)
    itpc = LinearInterpolation(Xpts, Ypts)
    integral, error = quadgk(x -> itpc(x), lbound, ubound, rtol = reltol)
    return integral
end


#Note had to edit the data file by adding commas next to them 
Data = readdlm("P2DATA", ',' , Float64, '\n', skipstart = 0)
t = Data[:,1]
E = Data[:,2]

#calculate the residence time from the RTD data given
E_t = E.*t #this is just a setup for the data 
tau = numint(t, E_t, t[1], t[end])

#calculate normalized RTD
F = numint(t, E, t[1], t[end])

#Question 1a=====================================================
function p5a0(tau)#We're passing in the residence time here
    #declare all my constant variables up here
    k = 175.0 #L^2/mol^2/s
    caf = cbf = 0.0313 #mol/L
    ccf = cdf = 0.0#mol/L
    V = 1000.0 #L
    Qf = 10.0 #L/s
    T = 320.0 #k

    parameters = [k,Qf]
    tauSpan = (0.0, 200*60)#seconds #200 is from the time data
    u0 = [caf, cbf, ccf, cdf]
    prob = ODEProblem(p5a1, u0, tauSpan, parameters)
    sol = solve(prob)
    tempsol = convert(Array, sol)
    conversion = (caf.-tempsol[1,:])./caf
    return conversion[end] #

end

function p5a1(du, u, p, tau)
    #define the variable 
    ca, cb, cc, cd = u 

    #assign the parameters
    k, Qf = p 

    #calculate rate constant 
    r = k*ca*cb^2

    #calculate production rates
    Ra = Rb = -r
    Rc = Rd = r

    #calculate the change in concentration 
    #this is dc/d(tau) = R
    du[1] = dca = Ra
    du[2] = dcb = Rb
    du[3] = dcc = Rc
    du[4] = dcd = Rd
end

#Question 1b=====================================================
function p5b0(tau)
    #declare all my constant variables up here
    k = 175.0 #L^2/mol^2/s
    caf = cbf = 0.0313 #mol/L
    ccf = cdf = 0.0#mol/L
    V = 1000.0 #L
    Qf = 10.0 #L/s
    T = 320.0 #k
    
    parameters = [k,tau,caf,cbf,ccf,cdf]
    tspan = (0.0, tau*60)
    u0 = [caf, cbf, ccf, cdf]
    prob = ODEProblem(p5b1, u0, tspan, parameters)
    sol = solve(prob)
    tempsol = convert(Array, sol)
    conversion = (caf.-tempsol[1,:])./caf
    return conversion[end] #

end

function p5b1(du, u, p, t)
    #define the variable 
    ca, cb, cc, cd = u 

    #assign the parameters
    k, tau, caf, cbf, ccf, cdf = p 

    #calculate rate constant 
    r = k*ca*cb^2

    #calculate production rates
    Ra = Rb = -r
    Rc = Rd = r

    #calculate the change in concentration
    du[1] = dca = caf - ca + Ra*tau*60
    du[2] = dcb = cbf - cb + Rb*tau*60
    du[3] = dcc = ccf - cc + Rc*tau*60
    du[4] = dcd = cdf - cd + Rd*tau*60
end

#Question 1c=====================================================
function p5c0()#note maxTime is 200 in this case
    #declare all my constant variables up here
    k = 175.0 #L^2/mol^2/s
    caf = cbf = 0.0313 #mol/L
    ccf = cdf = 0.0#mol/L
    V = 1000.0 #L
    Qf = 10.0 #L/s
    T = 320.0 #k
    
    #calculate the residence time
    tau = V/Qf #seconds 

    #get the amount of time steps
    Data = readdlm("P2DATA", ',' , Float64, '\n', skipstart = 0)
    #store time in seconds
    t = Data[:,1].*60 #seconds 
    #store max time in second 
    tmax = maximum(t)#Seconds 

    #preallocate the RTD 
    E = []

    #calculate RTD
    for i in t
        if i < tau/2
            push!(E, 0.0)
        else
            #calculate a temp variable and then push it
            temp = (tau^2)/(2*i^3) #1/seconds 
            push!(E,temp)
        end
    end

 
    #plot(t,E) #A check to see if I did it correctly

    #make an interpolating function of RTD
    Et = LinearInterpolation(t,E)

    
    #set up the ODE solver
    parameters = [Et, k, caf, cbf, ccf, cdf]
    tspan = (0.0, tmax)
    u0 = [caf, cbf, ccf, cdf, 0.0, 0.0, 0.0, 0.0]
    prob = ODEProblem(p5c1, u0, tspan, parameters)
    sol = solve(prob)

    #extract the average ending concentration of A 
    cend = sol.u[end]
    cabar = cend[5]

    #conversion of A 
    Xa = (caf - cabar)/caf
    
end

function p5c1(du, u, p, t)
    #first set are packet concentration second set are average concentrations
    ca, cb, cc, cd, cabar, cbbar, ccbar, cdbar = u
    
    #get the parameters 
    Et, k, caf, cbf, ccf, cdf = p

    #calculate the rate 
    r1 = k*ca*cb^2

    #calculate the production rates 
    Ra, Rb = -r1, -r1
    Rc, Rd = r1, r1

    #solve the balance on each packet
    du[1] = Ra
    du[2] = Rb
    du[3] = Rc
    du[4] = Rd

    #solve for the average exit concentration
    du[5] = ca*Et(t)
    du[6] = cb*Et(t)
    du[7] = cc*Et(t)
    du[8] = cd*Et(t)
end

#Question 1d=====================================================
function p5d0()#note maxTime is 200 in this case
    #declare all my constant variables up here
    k = 175.0 #L^2/mol^2/s
    caf = cbf = 0.0313 #mol/L
    ccf = cdf = 0.0#mol/L
    V = 1000.0 #L
    Qf = 10.0 #L/s
    T = 320.0 #k
    
    #calculate the residence time
    tau = V/Qf #seconds 

    #get the amount of time steps
    Data = readdlm("P2DATA", ',' , Float64, '\n', skipstart = 0)
    #store time in seconds
    t = Data[:,1].*60 #seconds 
    #store max time in second 
    tmax = maximum(t)#Seconds 

    #preallocate the RTD 
    E = Data[:,2]./60 #1/seconds 

    #plot(t,E) #A check to see if I did it correctly

    #make an interpolating function of RTD
    Et = LinearInterpolation(t,E)

    
    #set up the ODE solver
    parameters = [Et, k, caf, cbf, ccf, cdf]
    tspan = (0.0, tmax)
    u0 = [caf, cbf, ccf, cdf, 0.0, 0.0, 0.0, 0.0]
    prob = ODEProblem(p5d1, u0, tspan, parameters)
    sol = solve(prob)

    #extract the average ending concentration of A 
    cend = sol.u[end]
    cabar = cend[5]

    #conversion of A 
    Xa = (caf - cabar)/caf
    
end

function p5d1(du, u, p, t)
    #first set are packet concentration second set are average concentrations
    ca, cb, cc, cd, cabar, cbbar, ccbar, cdbar = u
    
    #get the parameters 
    Et, k, caf, cbf, ccf, cdf = p

    #calculate the rate 
    r1 = k*ca*cb^2

    #calculate the production rates 
    Ra, Rb = -r1, -r1
    Rc, Rd = r1, r1

    #solve the balance on each packet
    du[1] = Ra
    du[2] = Rb
    du[3] = Rc
    du[4] = Rd

    #solve for the average exit concentration
    du[5] = ca*Et(t)
    du[6] = cb*Et(t)
    du[7] = cc*Et(t)
    du[8] = cd*Et(t)
end


#Question 1e=====================================================
function p5e0()#note maxTime is 200 in this case
    #declare all my constant variables up here
    k = 175.0 #L^2/mol^2/s
    caf = cbf = 0.0313 #mol/L
    ccf = cdf = 0.0#mol/L
    V = 1000.0 #L
    Qf = 10.0 #L/s
    T = 320.0 #k

    #get the amount of time steps
    Data = readdlm("P2DATA", ',' , Float64, '\n', skipstart = 0)
    #store time in seconds
    t = Data[:,1].*60 #seconds 
    #store max time in second 
    tmax = maximum(t)#Seconds 

    #preallocate the RTD 
    E = Data[:,2]./60 #1/seconds 

    #plot(t,E) #A check to see if I did it correctly

    #make an interpolating function of RTD
    Et = LinearInterpolation(t,E)

    #preallocate the normalized RTD
    F = []
    
    #make the normalized RTD
    for i in t
        #temporarily store a data value
        temp = numint(t,E,0.0,i)
        push!(F,temp)
    end

    #make an interpolating function for the normalized RTD
    Ft = LinearInterpolation(t,F)
    
    #
    #set up the ODE solver
    parameters = [tmax, Et, Ft, k, caf, cbf, ccf, cdf]
    λspan = (11000, 0.0)#note that past 11000 seconds it's unstable for some reason
    u0 = [caf, cbf, ccf, cdf]
    prob = ODEProblem(p5e1, u0, λspan, parameters)
    sol = solve(prob)

    cend = sol[end]
    caend = cend[1]

    Xa = (caf-caend)/caf
end

function p5e1(du, u, p, λ)
    #first set are packet concentration second set are average concentrations
    ca, cb, cc, cd= u
    
    #get the parameters 
    tmax, Et, Ft, k, caf, cbf, ccf, cdf = p

    #calculate the rate 
    r1 = k*ca*cb^2

    #get the diffeq's
    constant = Et(λ)/(1-Ft(λ))
    du[1] = constant*(ca-caf)+r1
    du[2] = constant*(cb-cbf)+r1
    du[3] = constant*(cc-ccf)-r1
    du[4] = constant*(cd-cdf)-r1
end

Q1a = p5a0(tau)
Q1b = p5b0(tau)
Q1c = p5c0()
Q1d = p5d0()
Q1e = p5e0()

print(Q1a, ", ", Q1b, ", ", Q1c, ", ", Q1d, ", ", Q1e)