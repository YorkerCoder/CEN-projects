using Interpolations
using QuadGK


function NumInt(Xpts, Ypts, lbound, ubound ; reltol = 1e-15)
    itpc = LinearInterpolation(Xpts, Ypts)
    integral, error = quadgk(x -> itpc(x), lbound, ubound, rtol = reltol)
    return (integral, error)
end



x = [0.0, 0.9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.2, 10.0]
f(x) = x^3
y = f.(x)


answer = NumInt(x, y, 0.0, 10)
realvalue = 10^4/4
