using Base

struct Dual{T}
    val::T
    der::T
end

Base.:+(f::Dual, g::Dual) = Dual(f.val + g.val, f.der + g.der)
Base.:+(f::Dual, α::Number) = Dual(f.val + α, f.der)
Base.:+(α::Number, f::Dual) = f + α

Base.:-(f::Dual, g::Dual) = Dual(f.val - g.val, f.der - g.der)

Base.:*(f::Dual, g::Dual) = Dual(f.val * g.val, f.val * g.der + g.val * f.der)
Base.:*(f::Dual, α::Number) = Dual(f.val * α, f.der * α)
Base.:*(α::Number, f::Dual) = f * α

Base.:/(f::Dual, g::Dual) = Dual(f.val / g.val, (f.der * g.val - g.der * f.val) / (g.val^2))
Base.:/(f::Dual, α::Number) = Dual(f.val / α, f.der / α)
Base.:/(α::Number, f::Dual) = Dual(α / f.val, (-α * f.der) / (f.val^2))

# Base.:^(f::Dual, n::Integer) = Base.power_by_squaring(f,n)
Base.:^(f::Dual, n::Integer) = Dual(f.val^n, n * (f.val ^ (n - 1)) * f.der)

Base.:log(f::Dual) = Dual(log(f.val), f.der/f.val)

Base.:sin(f::Dual) = Dual(sin(f.val), cos(f.val) * f.der)
Base.:cos(f::Dual) = Dual(cos(f.val), -sin(f.val) * f.der)

Base.:exp(f::Dual) = Dual(exp(f.val), exp(f.val) * f.der)

# Base.:one(::Dual{T}) = Dual(one(T), zero(T))
Base.:one(f::Dual) = Dual(one(f.val), zero(f.val))
Base.:zero(f::Dual) = Dual(zero(f.val), zero(f.val))

function gradient_f(f, input)::Vector
    n = length(input)
    grad = Array{typeof(input[1]),1}(undef, n)
    input_dual = to_dual(input)
    for i in 1:n
        temp = input_dual[i]
        input_dual[i] = Dual(input[i], one(input[i]))
        grad[i] = f(input_dual).der
        input_dual[i] = temp
    end
    return grad
end

function hessian_ff(f,input)::Matrix
    n = length(input)
    hessian = Matrix{typeof(input[1])}(undef, n, n)
    input_dual = to_dual(input)
    for i in 1:n
        temp = input_dual[i]
        input_dual[i] = Dual(input[i], one(input[i]))
        for j in 1:n
            grad = gradient_f(f, input_dual)
            hessian[i, j] = grad[j].der
        end
        input_dual[i] = temp
    end
    return hessian
end

function to_dual(input)
    n = length(input)
    input_dual = Array{Dual{typeof(input[1])}, 1}(undef, n)
    for i in 1:length(input)
        input_dual[i] = Dual(input[i], zero(input[i]))
    end
    return input_dual
end