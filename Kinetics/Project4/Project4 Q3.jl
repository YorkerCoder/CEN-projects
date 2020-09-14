

using NLsolve

function p4a0()
    Rp = 0.20e-2 #m
    Da = 1.5e-2 * (100^2) #m^2/s
    T0 = 398 #k
    r0 = 2.63*(100^3) #mol/m^3/s
    Ca0 = 3.25e-5 * (100^3) #mol/m^3
    Ea = 20 #kcal/mol
    T1 = 448 #K
    R = 1.9872e-3 #kcal/K/mol

    vbar = Rp/3

    parameters = [vbar, Da, Ca0, r0]
    guesses = [3.0, 100000.0, 0.5]

    #solve for phi0, k0, and eta0
    ans1 = nlsolve((F,x)->p4a1(F,x,parameters), guesses, autodiff = :forward, iterations = 1000000)

    #store them
    phi0, k0, eta0 = ans1.zero

    #calculate k1
    k1 = k0*exp(Ea/R * ( (1/T0) - (1/T1) ) )

    #calculate phi
    phi = (k1*vbar/Da)^(1/2)

    #calculate eta
    eta =  (1/tanh(3phi) - 1/3phi )/phi

    #get the new rate
    r = k1*Ca0*eta

    #return the new rate
    return r
end

function p4a1(F,x,p)
    vbar, Da, Ca0, r0 = p
    phi, k, eta = x

    F[1] = k*vbar/Da - phi^2
    F[2] = k*Ca0*eta - r0
    F[3] = ((1/tanh(3phi)) - (1/3phi))/phi - eta
end
