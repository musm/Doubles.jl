module Doubles

export Double, Single, TwoFold, Dotted
import Base: +, -, *, /, sqrt,
             convert, promote_rule, show

FloatTypes = Union{Float32,Float64}
abstract AbstractDouble{T} <: Real

immutable Double{T<:FloatTypes} <: AbstractTwoFold{T}
    hi::T
    lo::T
end

immutable Single{T<:FloatTypes} <: AbstractTwoFold{T}
    hi::T
end

abstract AbstractTwoFold{T} <: Real

immutable TwoFold{T<:FloatTypes} <: AbstractTwoFold{T}
    value::T
    error::T
end

immutable Dotted{T<:FloatTypes} <: AbstractTwoFold{T}
    value::T
end

# The following hack promotes the float types to AbstractTwoFold so that 
# float types get properly converted to a Dotted type.
# We need this since we do not want to promote floats to TwoFold since
# we want to dispatch to methods with the Dotted type for effiency reaons
convert{T<:FloatTypes}(::Type{AbstractTwoFold{T}}, z::T) = Dotted(z)
promote_rule{T<:FloatTypes}(::Type{Dotted{T}}, ::Type{T}) = AbstractTwoFold{T}
promote_rule{T}(::Type{TwoFold{T}}, ::Type{T}) = AbstractTwoFold{T}

# see comment above
convert{T<:FloatTypes}(::Type{AbstractDouble{T}}, z::T) = Single(z)
promote_rule{T<:FloatTypes}(::Type{Single{T}}, ::Type{T}) = AbstractDouble{T}
promote_rule{T}(::Type{Double{T}}, ::Type{T}) = AbstractDouble{T}


# fast sqrt (no domain checking) make sure to handle errors in calling method
_sqrt{T<:FloatTypes}(x::T) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T, x)))



### basic error free floating point arithmetic ###


# fast two-sum addition, if |x| ≥ |y|
@inline function _tadd{T<:FloatTypes}(x::T, y::T)
    r0 = x + y
    r0, y - (r0 - x)
end

# two-sum addition
@inline function tadd{T<:FloatTypes}(x::T, y::T)
    r0 = x + y
    yt = r0 - x
    xt = r0 - yt
    r0, (y - yt) + (x - xt)
end

# fast two-sum subtraction, if |x| ≥ |y|
@inline function _tsub{T<:FloatTypes}(x::T, y::T)
    r0 = x - y
    r0, (x - r0) - y
end

# two-sum subtraction
@inline function tsub{T<:FloatTypes}(x::T, y::T)
    r0 = x - y
    yt = x - r0
    xt = yt - r0
    r0, (y - yt) + (xt - x)
end

# two-product fma
@inline function tmul{T<:FloatTypes}(x::T, y::T)
    r0 = x*y
    r0, fma(x,y,-r0)
end



### Double & TwoFold arithmetic ###



## addition


# both x and y are twofold
@inline function +{T}(x::TwoFold{T}, y::TwoFold{T})
    z0, z1 = tadd(x.value, y.value)
    TwoFold(z0, z1 + (x.error + y.error))
end

# x is twofold and y is dotted
@inline function +{T}(x::TwoFold{T}, y::Dotted{T})
    z0, z1 = tadd(x.value, y.value)
    TwoFold(z0, z1 + x.error)
end

# x is dotted and y is twofold
@inline +{T}(x::Dotted{T}, y::TwoFold{T}{T}) = y + x

# both x and y are dotted
@inline function +{T}(x::Dotted{T}, y::Dotted{T})
    TwoFold(tadd(x.value, y.value)...)
end


## subtraction


# both x and y are twofold
@inline function -{T}(x::TwoFold{T}, y::TwoFold{T})
    z0, z1 = psub(x.value, y.value)
    TwoFold(z0, z1 + (x.error - y.error))
end

# x is twofold and y is dotted
@inline function -{T}(x::TwoFold{T}, y::Dotted{T})
    z0, z1 = psub(x.value, y.value)
    TwoFold(z0, z1 + x.error)
end

# x is dotted and y is twofold
@inline function -{T}(x::Dotted{T}, y::TwoFold{T})
    z0, z1 = psub(x.value, y.value)
    TwoFold(z0, z1 + y.error)
end

# both x and y are dotted
@inline function -{T}(x::Dotted{T}, y::Dotted{T})
    TwoFold(psub(x.value, y.value)...)
end



## multiplication

# both x and y are twofold
@inline function *{T}(x::TwoFold{T}, y::TwoFold{T})
    z0, z1 = tmul(x.value, y.value)
    TwoFold(z0, (z1 + x.error*y.error) + (x.value*y.error + x.error*y.value))
end

# x is twofold and y is dotted
@inline function *{T}(x::TwoFold{T}, y::Dotted{T})
    z0, z1 = tmul(x.value, y.value)
    TwoFold(z0, z1 + x.error)
end

# x is dotted and y is twofold
@inline *{T}(x::Dotted{T}, y::TwoFold{T}) = y*x

# both x and y are dotted
@inline function *{T}(x::Dotted{T}, y::Dotted{T})
    TwoFold(tmul(x.value, y.value)...)
end



## division


# both x and y are twofold
@inline function /{T}(x::TwoFold{T}, y::TwoFold{T})
    q0 = x.value/y.value
    TwoFold(q0, (fma(-q0, y.value, x.value) + fma(-q0, y.error, x.error))/(y.value + y.error))
end

# x is twofold and y is dotted
@inline function /{T}(x::TwoFold{T}, y::Dotted{T})
    q0 = x.value/y.value
    TwoFold(q0, (x.error + fma(-q0, y.value, x.value))/y.value)
end

# x is dotted, y is twofold
@inline function /{T<:Dotted}(x::Dotted{T}, y::TwoFold{T})
    q0 = x.value/y.value
    TwoFold(q0, (fma(-q0, y.value, x.value) + -q0*y.error)/(y.value + y.error))
end

# both x and y are dotted
@inline function /{T}(x::Dotted{T}, y::Dotted{T})
    q0 = x.value/y.value
    TwoFold(q0, (fma(-q0, y.value, x.value))/y.value)
end


## square root (see warning for _sqrt above)


# x is dotted, z is twfold
@inline function sqrt(x::Dotted)
    z0 = _sqrt(x.value)
    TwoFold(z0, fma(-z0, z0, x.value)/(z0+z0)) # actually as good as Double here
end

# x is twofold, z is twofold
@inline function sqrt(x::TwoFold)
    z0 = _sqrt(x.value)
    u0, u1 = tadd(x.value, x.error)
    v0 = _sqrt(u0)
    v1 = (u1 + fma(-v0, v0, u0)) / (v0 + v0)
    w  = tsub(v, Dotted(z0))
    TwoFold(z0, w.value + w.error)
end

# x is single, z is double
@inline function sqrt(x::Single)
    z0 = _sqrt(x.hi)
    Double(z0, fma(-z0, z0, x.hi)/(z0+z0))
end

# x is double, z is double
@inline function sqrt(x::Double)
    z0 = _sqrt(x.hi)
    z1 = x.lo + fma(-z0, z0, x.hi)
    Double(_normalize(z0, (x.lo + fma(-z0, z0, x.hi))/(z0+z0))...)
end

### auxiliary


@inline -(x::TwoFold) = TwoFold(-x.value, -x.error)
@inline -(x::Dotted)  = Dotted(-x.value)

@inline -(x::Double) = Double(-x.value, -x.error)
@inline -(x::Single)  = Single(-x.value)

@inline scale{T}(x::TwoFold{T}, s::Dotted{T}) = TwoFold(s.value*x.value, s.value*x.error)

function show{T<:TwoFold}(io::IO, x::T)
    println(io, T)
    print(io, x.value, " [", x.error, "]")
end

function show{T<:Dotted}(io::IO, x::T)
    println(io, T)
    print(io, x.value)
end

function show{T<:Double}(io::IO, x::T)
    println(io, T)
    print(io, x.hi, " hi, ", x.lo, " lo")
end

function show{T<:Single}(io::IO, x::T)
    println(io, T)
    print(io, x.hi)
end

end