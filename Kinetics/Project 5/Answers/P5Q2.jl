#=
# if you run the function p5a0(Et,Ft) you'll get 0.48 for the segregation model and 0.47 for the maximum mixidness model
# run p5b0(Et,Ft) for part b 
# a) Xas =0.478 , Xamm = 0.467 
# b) Ybs = 0.47, Ybmm = 0.46, Ycs = 0.023, Ycmm = 0.0253
=#


using DifferentialEquations
using Plots
using DelimitedFiles
using Interpolations, QuadGK

function numint(Xpts, Ypts, lbound, ubound ; reltol = 1e-15)
    itpc = LinearInterpolation(Xpts, Ypts)
    integral, error = quadgk(x -> itpc(x), lbound, ubound, rtol = reltol)
    return integral
end

#constructing the concentration curve
function curve(x)
    if 0.0 < x <1.0
        return x 
    elseif 1.0 <= x < 2.0
        return 2-x
    else
        return 0
    end
end

#generate the curve
y = curve.(0.0:0.01:2.0)
#get the integral of the curve
intc, error = quadgk(x->curve(x),0.0,2.0,rtol = 1e-15)

#get the RTD
E = y./intc
#make it an interpolating function
Et = LinearInterpolation(0.0:0.01:2.0,E)

t = 0.0:0.01:2.0
F = []
for i in t
    #temporarily store a data value
    temp = numint(t,E,0.0,i)
    push!(F,temp)
end

Ft = LinearInterpolation(t,F)

#this holds the answers to both segregation and max mix 
#ANSWER = (0.477, 0.467)
function p5a0(Et,Ft)
    T = 300 #K 
    k = 0.5 #L/mol

    #inital concentrations
    caf = 1.0 #mol/L 
    cbf = 0.0 #mol/L 

    #parameters
    p1 = [Et, k, caf, cbf]
    c01 = [caf, cbf, 0.0, 0.0] #0.0 are the average concentrations
    tspan = (0.0, 2.0)#going from 0 to 2.0 minutes just like what E(t) goes to 

    #set up the problem for the segregation model
    prob1 = ODEProblem(p5a1s, c01, tspan, p1)
    #solve it for the segregation method
    sol1 = solve(prob1)

    #get the conversion for the problem
    cend = sol1.u[end]
    cabar = cend[3]
    Xa1 = (caf - cabar)/caf

    #set up the maximum mixing problem ========================================
    p2 = [Et, Ft, k, caf, cbf]
    c02 = [caf, cbf]
    λspan = (1.99999999999, 0.0) #note had to do it this way otherwise i divide by 0
    prob2 = ODEProblem(p5a1m, c02, λspan, p2)
    sol2 = solve(prob2)

    #calculate the conversion for maximum mixidness 
    cend = sol2.u[end]
    caend = cend[1]
    Xa2 = (caf - caend)/caf

    return (Xa1, Xa2)
end

#solves for the segregation
function p5a1s(du, u, p, t)
    ca, cb, cabar, cbbar = u
    Et, k, caf, cbf = p
    
    #calculate rate
    r = k*ca^2

    #calculate the production rate 
    Ra = -2r
    Rb = r

    #solve the balance on a packet
    du[1] = Ra
    du[2] = Rb

    #get the average concentration
    du[3] = ca*Et(t)
    du[4] = cb*Et(t)
end

#solves for the max mixing
function p5a1m(du, u, p, λ)
    ca, cb = u
    Et, Ft, k, caf, cbf = p

    #get the rate
    r = k*ca^2

    Ra = -2r
    Rb = r

    #get the DifferentialEquations

    constant = Et(λ)/(1-Ft(λ))
    du[1] = constant*(ca-caf) - Ra
    du[2] = constant*(cb-cbf) - Rb 
end




#############################################################################################
#2b)
function p5b0(Et,Ft)
    T = 300 #K 
    k1 = 0.5 #L/mol 
    k2 = 0.12 #L/mol/min 

    #inital concentrations
    caf = 1.0 #mol/L 
    cbf = ccf = 0.0 #mol/L 

    #set up the segregation mode==============================================
    #parameters
    p1 = [Et, k1, k2]
    c01 = [caf, cbf, ccf, 0.0, 0.0, 0.0] #0.0 are the average concentrations
    tspan = (0.0, 2.0)#going from 0 to 2.0 minutes just like what E(t) goes to 

    #set up the problem for the segregation model
    prob1 = ODEProblem(p5b1s, c01, tspan, p1)
    #solve it for the segregation method
    sol1 = solve(prob1, Rosenbrock23())


    #get the overall selectivity for the problem
    cend = sol1.u[end]
    cabar = cend[4]
    cbbar = cend[5]
    ccbar = cend[6]

    Yb1 = cbbar/(caf-cabar)
    Yc1 = ccbar/(caf-cabar)
    con1 = (caf - cabar)/caf
   
    #set up the maximum mixing problem ========================================
    p2 = [Et, Ft, k1, k2, caf, cbf, ccf]
    c02 = [caf, cbf, ccf]
    λspan = (1.99999999999, 0.0) #note had to do it this way otherwise i divide by 0
    prob2 = ODEProblem(p5b1m, c02, λspan, p2)
    sol2 = solve(prob2)

    #calculate the overall selectivity for maximum mixidness 
    cend = sol2.u[end]
    caend = cend[1]
    cbend = cend[2]
    ccend = cend[3]

    Yb2 = cbend/(caf-caend)
    Yc2 = ccend/(caf-caend)
    con2 = (caf-caend)/caf

    #selectivity of b seg, selectivity of b mm, 
    return (Yb1, Yb2, Yc1, Yc2)
end

#solves for the segregation
function p5b1s(du, u, p, t)
    ca, cb, cc, cabar, cbbar, ccbar = u
    Et, k1, k2 = p
    
    #calculate rate
    r1 = k1*ca^2
    r2 = k2*ca*cb 

    #calculate the production rate 
    Ra = -2r1 - r2 
    Rb = r1 - r2 
    Rc = r2 

    #solve the balance on a packet
    du[1] = Ra
    du[2] = Rb
    du[3] = Rc 

    #get the average concentration
    du[4] = ca*Et(t)
    du[5] = cb*Et(t)
    du[6] = cc*Et(t)
end

#solves for the max mixing
function p5b1m(du, u, p, λ)
    ca, cb, cc = u
    Et, Ft, k1, k2, caf, cbf, ccf = p

    #get the rate
    r1 = k1*ca^2
    r2 = k2*ca*cb

    #get production rates 
    Ra = -2r1-r2
    Rb = r1 -r2
    Rc = r2

    #solve the DifferentialEquations
    constant = Et(λ)/(1-Ft(λ))
    du[1] = constant*(ca-caf) - Ra
    du[2] = constant*(cb-cbf) - Rb
    du[3] = constant*(cc-ccf) - Rc 
end
