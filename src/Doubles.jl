module Doubles

export Double

import Base: +, -, *, /, sqrt, abs, convert, promote_rule, show

FloatTypes = Union{Float32,Float64}
abstract AbstractDouble{T} <: Real

immutable Single{T<:FloatTypes} <: AbstractDouble{T}
    hi::T
end

immutable Double{T<:FloatTypes} <: AbstractDouble{T}
    hi::T
    lo::T
    # function Double(u::T, v::T)
    #     r = u + v
    #     new(r, v + (u - r))
    # end
end
Double(x) = Single(x)



### promotions and conversions ###

# The following hack promotes the float types to AbstractDouble so that 
# float types get properly converted to a Single type.
# We need this since we do not want to promote floats to Double since
# we want to dispatch to methods with the Single type for effiency reasons.
# Similar for the other types.

convert{T}(::Type{AbstractDouble{T}}, z::T)               = Single(z)
convert{T}(::Type{AbstractDouble{T}}, z::Type{Single{T}}) = z
convert{T}(::Type{AbstractDouble{T}}, z::Type{Double{T}}) = z

promote_rule{T}(::Type{Single{T}}, ::Type{T})         = AbstractDouble{T}
promote_rule{T}(::Type{Double{T}}, ::Type{T})         = AbstractDouble{T}
promote_rule{T}(::Type{Double{T}}, ::Type{Single{T}}) = AbstractDouble{T}



### utility functions


@inline function normalize(x::Double)
    r = x.hi + x.lo
    Double(r, (x.hi - r) + x.lo)
end

@inline function normalize{T<:FloatTypes}(x::T, y::T) # same as fast two-sum
    r = x + y
    r, y + (x - r)
end


# the following are only used for non fma systems
# clear lower 27 bits (leave upper 26 bits)
@inline trunclo(x::Float64) = reinterpret(Float64, reinterpret(UInt64, x) & 0xffff_ffff_f800_0000)
# clear lowest 12 bits (leave upper 12 bits)
@inline trunclo(x::Float32) = reinterpret(Float32, reinterpret(UInt32, x) & 0xffff_f000)

@inline function splitprec(x::FloatTypes)
    hx = trunclo(x)
    hx, x-hx
end



### basic error free float arithmetic ###


# fast two-sum addition, if |x| ≥ |y|
@inline function _tadd{T<:FloatTypes}(x::T, y::T)
    r = x + y
    r, y + (x - r)
end

# two-sum addition
@inline function tadd{T<:FloatTypes}(x::T, y::T)
    r = x + y
    v = r - x
    r, (y - v) + (x - (r - v))
end

# two-product fma
@inline function tmul{T<:FloatTypes}(x::T, y::T)
    r = x*y
    r, fma(x,y,-r)
end

# two-product non fma
@inline function tmul_{T<:FloatTypes}(x::T, y::T)
    hx, lx = splitprec(x)
    hy, ly = splitprec(y)
    z = x*y
    z, ((hx*hy-z) + lx*hy + hx*ly) + lx*ly
end


### Double arithmetic ###


## negation

@inline -(x::Double)  = Double(-x.hi, -x.lo)
@inline -(x::Single)  = Single(-x.hi)


## addition

@inline function +{T}(x::Double{T}, y::Double{T})
    r, s = tadd(x.hi, y.hi)
    Double(r, s + x.lo + y.lo)
end

@inline function +{T}(x::Double{T}, y::Single{T})
    r, s = tadd(x.hi, y.hi)
    Double(r, s + x.lo)
end
@inline +{T}(x::Single{T}, y::Double{T}{T}) = y + x

@inline function +{T}(x::Single{T}, y::Single{T})
    r, s = tadd(x.hi, y.hi)
    Double(r,s)
end


## subtraction

@inline function -{T}(x::Double{T}, y::Double{T})
    r, s = tadd(x.hi, -y.hi)
    Double(r, s + (x.lo - y.lo))
end

@inline function -{T}(x::Double{T}, y::Single{T})
    r, s = tadd(x.hi, -y.hi)
    Double(r, s + x.lo)
end

@inline function -{T}(x::Single{T}, y::Double{T})
    r, s = tadd(x.hi, -y.hi)
    Double(r, s + y.lo)
end

@inline function -{T}(x::Single{T}, y::Single{T})
    Double(tsub(x.hi, -y.hi)...)
end


## multiplication


@inline function *{T}(x::Double{T}, y::Double{T})
    r, s = tmul(x.hi, y.hi)
    Double(r, s + x.hi*y.lo + x.Lo*y.hi)
end

@inline function *{T}(x::Double{T}, y::Single{T})
    z0, z1 = tmul(x.hi, y.hi)
    Double(z0, z1 + x.lo)
end
@inline *{T}(x::Single{T}, y::Double{T}) = y*x

@inline function *{T}(x::Single{T}, y::Single{T})
    Double(tmul(x.hi, y.lo)...)
end


# ## division

# # both x and y are twofold
# @inline function /{T}(x::TwoFold{T}, y::TwoFold{T})
#     q0 = x.value/y.value
#     TwoFold(q0, (fma(-q0, y.value, x.value) + fma(-q0, y.error, x.error))/(y.value + y.error))
# end

# # x is twofold and y is dotted
# @inline function /{T}(x::TwoFold{T}, y::Dotted{T})
#     q0 = x.value/y.value
#     TwoFold(q0, (x.error + fma(-q0, y.value, x.value))/y.value)
# end

# # x is dotted, y is twofold
# @inline function /{T<:Dotted}(x::Dotted{T}, y::TwoFold{T})
#     q0 = x.value/y.value
#     TwoFold(q0, (fma(-q0, y.value, x.value) + -q0*y.error)/(y.value + y.error))
# end

@inline function /{T}(x::Single{T}, y::Single{T}) # fma only
    ry = 1/y.hi
    r = x.hi*ry
    Double(r, fma(-r, y.hi, x.hi)*ry)
end


@inline function div_nonfma{T}(x::Single{T}, y::Single{T})
    ry = 1/y
    r = x*ry
    hx, lx = splitprec(r)
    hy, ly = splitprec(y)
    z = r*y
    Double(r, (((-hx*hy+z) - lx*hy - hx*ly) - lx*ly)*ry)
end

## square root


# fast sqrt (no domain checking) make sure to handle errors in calling method
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T, x)))

# x is single, z is double
@inline function sqrt(x::Single)
    r = _sqrt(x.hi)
    Double(r, fma(-r, r, x.hi)/(r+r))
end

# x is double, z is double
@inline function sqrt(x::Double)
    r = _sqrt(x.hi)
    Double(r, (x.lo + fma(-r, r, x.hi))/(r+r))
end



### auxiliary


function abs(x::Single)
    Single(abs(x.hi))
end

function abs(x::Double)
    Double(abs(x.hi), abs(x.lo))
end

function show{T}(io::IO, x::Double{T})
    println(io, "Double{$T}")
    print(io, x.hi, ", ", x.lo)
end

function show{T}(io::IO, x::Single{T})
    println(io, "Double{$T}")
    print(io, x.hi)
end


end