using Doubles
using Base.Test


x = sqrt(2.0)
bx = big(x)
sx = x
dx = Double(x)

y = 0.1
by = big(y)
sy = y
dy = Double(y)

@test x == sx == dx
@test y == sy == dy

dxy = dx*dy
bxy = bx*by
# @test sx*sy == dxy
@test x*y == Float64(dxy)
# @test dxy == Double(bxy)

@test x+y == Float64(dx+dy)
# @test dx+dy == Double(bx+by)

@test x-y == Float64(dx-dy)
# @test dx-dy == Double(bx-by)

@test x/y == Float64(dx/dy)
# @test dx/dy == Double(bx/by)

@test sqrt(y) == Float64(sqrt(dy))
# @test sqrt(dy) == Double(sqrt(by))