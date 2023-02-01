# A basic implementation of Forward Mode AD using dual numbers
# Why Dual Numbers:
# f(x+ϵ)=f(x)+f′(x)ϵ

using Base

# Define Dual data type with value and derivative components
struct Dual{T}
    val::T
    der::T
end

# Overload basic operators to work with dual numbers
# See DiffRules.jl for more 
Base.:+(f::Dual, g::Dual) = Dual(f.val + g.val, f.der + g.der)
Base.:+(f::Dual, α::Number) = Dual(f.val + α, f.der)
Base.:+(α::Number, f::Dual) = f + α

Base.:-(f::Dual, g::Dual) = Dual(f.val - g.val, f.der - g.der)

Base.:*(f::Dual, g::Dual) = Dual(f.val * g.val, f.val * g.der + g.val * f.der)
Base.:*(f::Dual, α::Number) = Dual(f.val * α, f.der * α)
Base.:*(α::Number, f::Dual) = f * α

Base.:/(f::Dual, g::Dual) = Dual(f.val / g.val, (f.der * g.val - g.der * f.val) / (g.val^2))
Base.:/(f::Dual, α::Number) = Dual(f.val / α, f.der / α)
Base.:/(α::Number, f::Dual) = Dual(α / f.val, 1 / f.der)

Base.:^(f::Dual, n::Integer) = Base.power_by_squaring(f,n)

# Define exponential function for Dual data type
exp(f::Dual) = Dual(exp(f.val), exp(f.val) * f.der)

# Example 1 
h(x) = x^2 + 2
a = 3
x̂ = Dual(a,1)

# [3,1]->[3,1]^2+2=[9,6]+2=[11,6]
h(x̂)

# Define derivative function 
derivative(f, x) = f(Dual(x, one(x))).der

# Example 2
r(x) = 3x^3+2x+3
b = 2
derivative(r, b)
# Inline 
derivative(x -> 3x^3+2x+3, 2)

# TODO define the Dual operators to work with vectors
# to find gradient with one pass

# Example 3
# ∂d/∂a
a=Dual(1,1)
b=Dual(2,0)
c=Dual(3,0)
d=a*(b+c)^2
# ∂d/∂b
a=Dual(1,0)
b=Dual(2,1)
c=Dual(3,0)
d=a*(b+c)^2
# ∂d/∂c
a=Dual(1,0)
b=Dual(2,0)
c=Dual(3,1)
d=a*(b+c)^2