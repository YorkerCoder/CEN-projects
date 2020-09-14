# name: Emir Hadzic
# Class: CEN 786 Kinetics
# Assignment: Project 3 Question 1 a,b,c


# ==============================================================================
# (a) Tss = 541.4 °R
# (b) function p31() plot rate of heat removal and generation
# (c) function p30() will get the eigen values and steady state temperature
#       λ1 and λ2 are negative so the system is stable
# ==============================================================================

using Plots; gr()
using NLsolve
using LinearAlgebra#in this case i'm using it to get eigen values and vectors


#p30 and p30! solve using the non linear solver for T
function p30()
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

    parameters = [ΔH, cpa, cpb, cpc, cpm, va, vb, vm, ua, Vr, T0,
                  faf, fbf, fmf, Ta, R, Qf, tau]
    Guess = [548.0]

    #tell it to solve autodiff forward skips us feeding a jacobian matrix
    answer = nlsolve((F,x)->p30!(F,x,parameters), Guess, autodiff = :forward)
    Tss = answer.zero[1]
    #rate constant at steady state
    k = 16.96e12*exp(-32000/R/Tss)

    #conversion at steady state
    xa = k*tau/(1+k*tau)

    #concentrations at feed
    caf, cbf, cmf = [faf, fbf, fmf]/Qf

    #concentrations at steady state
    ca, cb, cc, cm = [caf*(1-xa), cbf-caf*xa, caf*xa, cmf]

    #calculate ΔHr
    ΔHr = ΔH + (cpc - cpa - cpb)*(Tss - 536.73)

    #define the long terms in the bottom and top here
    top    = faf*cpa + fbf*cpb + fmf*cpm
    bottom = ca*cpa + cb*cpb + cc*cpc + cm*cpm

    df1dCa = (-1.0/tau) - k
    df1dT  = (32400.0/R/(Tss^2))*k*ca
    df2dCa = -ΔHr*k/bottom
    df2dT  = (-top - ua + ΔHr*(32400.0/R/Tss^2)*k*Vr - 7*k*Vr)/bottom #the -7 is from the derivative of deltaHr with respect to T

    #make a 2x2 array for these
    J = [df1dCa df1dT; df2dCa df2dT]

    #get the eigenvalues of these
    EigenValues = eigvals(J)

    return(Tss, EigenValues)
end

#this is the function that the non linear solver solves for
function p30!(F,x,p)
    #pass all of the parameters
    ΔH, cpa, cpb, cpc, cpm, va, vb, vm, ua, Vr, T0, faf, fbf, fmf, Ta, R, Qf, tau = p

    #caluclate k
    k = 16.96e12 * exp(-32400/R/x[1])

    #calculate ΔH of reaction
    ΔHr = ΔH + (cpc - cpa - cpb)*(x[1] - 536.73)

    F[1] = -ΔHr*faf*(k*tau/(1+k*tau)) +((faf*cpa + fbf*cpb + fmf*cpm)*(T0 - x[1])+ua*(Ta-x[1]))
end




#the following 2 functions are used to graph generation and removal
function p31()
    T = 450:1:700
    answer1 = p31!.(T)#vecotorized my function
    answer2 = hcat(answer1...)
    g1 = plot(T,answer2[1,:])
    g1 = plot!(T, answer2[2,:])

    #g2 = plot(T, answer2[3,:])

    #return(g1,g2)

    #plot it
    g1
end

function p31!(var)
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

    #calculate initial feed concentrations
    caf, cbf, cdf = [faf, fbf, fmf]/Qf

    #caluclate k
    k = 16.96e12 * exp(-32400/R/var)

    #calculate ΔH of reaction
    ΔHr = ΔH + (cpc - cpa - cpb)*(var - 536.73)

    G = -ΔHr*faf*(k*tau/(1+k*tau))
    R = -((faf*cpa + fbf*cpb + fmf*cpm)*(T0 - var)+ua*(Ta-var))
    F = G - R

    return [G,R,F]
end
