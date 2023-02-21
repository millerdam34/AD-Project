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
# why?
Base.:/(α::Number, f::Dual) = Dual(α / f.val, (-α * f.der) / (f.val^2))

Base.:^(f::Dual, n::Integer) = Base.power_by_squaring(f,n)
Base.:^(f::Dual, n::Integer) = Dual(f.val^n, n * (f.val ^ (n - 1)) * f.der)

Base.:log(f::Dual) = Dual(log(f.val), f.der/f.val)

Base.:sin(f::Dual) = Dual(sin(f.val), cos(f.val) * f.der)

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

# Quotient Example
a=Dual(1,0)
b=Dual(2,1)

d=a/b
d1=1/b
d2=a/2

# power rule test

f(x,y)=(x+y)^10

a=Dual(1,1)
b=Dual(1,0)

Base.:^(f::Dual, n::Integer) = Dual(f.val^n, n * (f.val ^ (n - 1)) * f.der)

@time f(a,b)

Base.:^(f::Dual, n::Integer) = Base.power_by_squaring(f,n)

@time f(a,b)

# Example 

@time begin
x_1 = Dual(1.0,1.0)
x_2 = Dual(2.0,0.0)

y=log(x_1)+x_1*x_2-sin(x_2)

# ∂y/∂x_1
y.der

x_1 = Dual(1.0,0.0)
x_2 = Dual(2.0,1.0)

y=log(x_1)+x_1*x_2-sin(x_2)

# ∂y/∂x_2
y.der

end

@time begin
    x_1 = Dual(1.0,1.0)
    x_2 = Dual(2.0,0.0)
    y1=log(x_1)+x_1*x_2-sin(x_2)
    x_1 = Dual(1.0,0.0)
    x_2 = Dual(2.0,1.0)
    y2=log(x_1)+x_1*x_2-sin(x_2)
    return y1.der, y2.der
end

@time test1()