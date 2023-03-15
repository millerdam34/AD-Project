using Base

struct Dual{T}
    val::T
    der::T
end

Base.:+(f::Dual, g::Dual) = Dual(f.val + g.val, f.der + g.der)
Base.:+(f::Dual, α::Number) = Dual(f.val + α, f.der)
Base.:+(α::Number, f::Dual) = f + α

Base.:-(f::Dual, g::Dual) = Dual(f.val - g.val, f.der - g.der)
Base.:-(f::Dual) = Dual(-f.val, -f.der)

Base.:*(f::Dual, g::Dual) = Dual(f.val * g.val, g.der * f.val + f.der * g.val)
Base.:*(f::Dual, α::Number) = Dual(f.val * α, f.der * α)
Base.:*(α::Number, f::Dual) = f * α

Base.:/(f::Dual, g::Dual) = Dual(f.val / g.val, (f.der * g.val - g.der * f.val) / (g.val * g.val))
Base.:/(f::Dual, α::Number) = Dual(f.val / α, f.der / α)
Base.:/(α::Number, f::Dual) = Dual(α / f.val, (-α * f.der) / (f.val^2))

# Base.:^(f::Dual, n::Integer) = Base.power_by_squaring(f,n)
Base.:^(f::Dual, n::Integer) = Dual(f.val^n, n * (f.val ^ (n - 1)) * f.der)

Base.:log(f::Dual) = Dual(log(f.val), f.der/f.val)

Base.:sin(f::Dual) = Dual(sin(f.val), f.der * cos(f.val))
Base.:cos(f::Dual) = Dual(cos(f.val), f.der * -sin(f.val))

Base.:exp(f::Dual) = Dual(exp(f.val), f.der * exp(f.val))

# how to generic?
Base.:one(::Dual{Float64}) = Dual(one(Float64), zero(Float64)) 
Base.:zero(::Dual{Float64}) = Dual(zero(Float64), zero(Float64))
Base.:one(f::Dual) = Dual(one(f.val), zero(f.val))
Base.:zero(f::Dual) = Dual(zero(f.val), zero(f.val))

# returns the derivative of an arbitraty direction vector
function derivative_f(f, input::Vector, direction::Vector, normalize::Bool=true)::typeof(input[1])
    input_dual = to_dual(input, direction, normalize)
    return f(input_dual).der
end

# returns the derivative and value of an arbitrary direction vector
function derivative_value_f(f, input::Vector, direction::Vector, normalize::Bool=true)::typeof(input[1])
    input_dual = to_dual(input, direction, normalize)
    return f(input_dual)
end

# calculates the second derivative using Dual{Dual} numbers
# input: x⃗
# direction1: x⃗₁
# direction2: x⃗₂
# output: ∂²x⃗/∂x⃗₁∂x⃗₂
function second_derivative_f(f, input::Vector, direction1::Vector, direction2::Vector, normalize1::Bool=true, normalize2::Bool=true)
    # Dual(∂x/∂x⃗₁, ∂x⃗/∂x⃗₁∂x⃗₂)
    if (!normalize1)
        input_dual1 = map(x -> Dual(x, zero(input[1])), direction1)
    else
        norm = norm(direction1)
        input_dual1 = map(x -> Dual(x / norm, zero(input[1])), direction1)
    end
    # Dual(x⃗, ∂x/∂x⃗₂)
    input_dual2 = to_dual(input, direction2, normalize2)
    # Dual(Dual(x⃗, ∂x⃗/∂x⃗₂), Dual(∂x⃗/∂x⃗₁, ∂x⃗/∂x⃗₁∂x⃗₂))
    input_dual_dual = to_dual(input_dual2, input_dual1, false)
    return f(input_dual_dual).der.der
end

# calculates the gradient of a function 
# input: x⃗
# f: f(x⃗)
# output: ∇f
function gradient_f(f, input)::Vector
    n = length(input)
    type = typeof(input[1])
    one_val = one(type)
    zero_val = zero(type)
    grad = Array{type,1}(undef, n)
    input_dual = collect(map(x -> Dual(x, zero_val), input))
    for i in 1:n
        temp = input[i]
        input_dual[i] = Dual(temp, one_val)
        grad[i] = f(input_dual).der
        input_dual[i] = Dual(temp, zero_val)
    end
    return grad
end


# input: x⃗
# direction: x⃗₁
# x⃗ -> Dual(x⃗, ∂x⃗/∂x⃗₁)
# x⃗ᵢ = Dual(x⃗ᵢ, x⃗₁ᵢ)
function to_dual(input::Vector, direction::Vector, normalize::Bool=true)
    n = length(input)
    if (normalize)
        norm = norm(direction)
        if (norm)
            map(x -> x / norm, direction)
        end
    end
    input_dual = Array{Dual{typeof(input[1])}, 1}(undef, n)
    for i in 1:n
        input_dual[i] = Dual(input[i], direction[i])
    end
    return input_dual
end
