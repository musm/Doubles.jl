using Doubles
using Base.Test

# write your own tests here
x = sqrt(2.0)
bx = big(x)
sx = Double(x)
dx = Double(x,0.0)

y = 0.1
by = big(y)
sy = Double(y)
dy = Double(y,0.0)

@test x == sx == dx
@test y == sy == dy

dxy = dx*dy
bxy = bx*by
@test sx*sy == dxy
@test x*y == Float64(dxy)
@test dxy == Double(bxy)