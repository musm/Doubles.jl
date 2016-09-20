module TwoFolds

export TwoFold, Dotted
import Base: show, +, -, *, /, sqrt, convert, promote_rule

FloatTypes = Union{Float32,Float64}
abstract ATwoFold{T}  <: Real

immutable TwoFold{T<:FloatTypes} <: ATwoFold{T}
    value::T
    error::T
end

immutable Dotted{T<:FloatTypes} <: ATwoFold{T}
    value::T
end

# The following hack promotes the float types to ATwoFold and then 
# this gets properly converted to a Dotted type.
# We need this since we do not want to promote floats to a TwoFold since
# we want to dispatch to the proper methods with the Dotted type for effiency reaons
convert{T<:FloatTypes}(::Type{ATwoFold{T}}, z::T) = Dotted(z)
promote_rule{T<:FloatTypes}(::Type{Dotted{T}}, ::Type{T}) = ATwoFold{T}
promote_rule{T}(::Type{TwoFold{T}}, ::Type{T}) = ATwoFold{T}


# fast sqrt (no domain checking) make sure to handle errors in calling method
_sqrt{T}(x::Dotted{T}) = Base.box(T, Base.sqrt_llvm_fast(Base.unbox(T, x)))


### basic error free floating point arithmetic ###


# private fast addition, if |x| ≥ |y|
@inline function _padd{T<:FloatTypes}(x::T, y::T)
    r0 = x + y
    r0, y - (r0 - x)
end

# private addition
@inline function padd{T<:FloatTypes}(x::T, y::T)
    r0 = x + y
    yt = r0 - x
    xt = r0 - yt
    r0, (y - yt) + (x - xt)
end

# private fast subtraction, if |x| ≥ |y|
@inline function _psub{T<:FloatTypes}(x::T, y::T)
    r0 = x - y
    r0, (x - r0) - y
end

# private subtraction
@inline function psub{T<:FloatTypes}(x::T, y::T)
    r0 = x - y
    yt = x - r0
    xt = yt - r0
    r0, (y - yt) + (xt - x)
end

# private multiply, using fma
@inline function pmul{T<:FloatTypes}(x::T, y::T)
    r0 = x*y
    r0, fma(x,y,-r0)
end



### TwoFold arithmetic ###



## addition


# both x and y are twofold
@inline function +{T}(x::TwoFold{T}, y::TwoFold{T})
    z0, z1 = padd(x.value, y.value)
    TwoFold(z0, z1 + (x.error + y.error))
end

# x is twofold and y is dotted
@inline function +{T}(x::TwoFold{T}, y::Dotted{T})
    z0, z1 = padd(x.value, y.value)
    TwoFold(z0, z1 + x.error)
end

# x is dotted and y is twofold
@inline +{T}(x::Dotted{T}, y::TwoFold{T}{T}) = y + x

# both x and y are dotted
@inline function +{T}(x::Dotted{T}, y::Dotted{T})
    TwoFold(padd(x.value, y.value)...)
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
    z0, z1 = pmul(x.value, y.value)
    TwoFold(z0, (z1 + x.error*y.error) + (x.value*y.error + x.error*y.value))
end

# x is twofold and y is dotted
@inline function *{T}(x::TwoFold{T}, y::Dotted{T})
    z0, z1 = pmul(x.value, y.value)
    TwoFold(z0, z1 + x.error)
end

# x is dotted and y is twofold
@inline *{T}(x::Dotted{T}, y::TwoFold{T}) = y*x

# both x and y are dotted
@inline function *{T}(x::Dotted{T}, y::Dotted{T})
    TwoFold(pmul(x.value, y.value)...)
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
    z0 = _sqrt(x)
    TwoFold(z0, fma(-z0, z0, x)/(z0+z0)) # actually as good as Double here
end

# x is twofold, z is twofold
@inline function sqrt(x::TwoFold)
    z0 = _sqrt(x.value)
    u0, u1 = padd(x.value, x.error)
    v0 = _sqrt(u0)
    v1 = add(u1, fma(-v0, v0, u0)) / (v0 + v0)
    w  = sub(v, Dotted(z0))
    TwoFold(z0, w.value + w.error)
end



### auxiliary


@inline -(x::TwoFold) = Double(-x.value, -x.error)
@inline -(x::Dotted)  = Dotted(-x.value)

@inline scale{T}(x::TwoFold{T}, s::Dotted{T}) = TwoFold(s.value*x.value, s.value*x.error)

function show{T<:TwoFold}(io::IO, x::T)
    println(io, T)
    print(io, x.value, " [", x.error, "]")
end

function show{T<:Dotted}(io::IO, x::T)
    println(io, T)
    print(io, x.value)
end

end