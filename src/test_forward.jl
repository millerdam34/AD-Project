include("forward.jl")

f₁((x₁,x₂)) = log(x₁) + x₁*x₂ - sin(x₂)
f₂((x₁, x₂, x₃)) = (x₁ * x₂ * sin(x₃) + exp(x₁ * x₂)) / x₃

x₁ = 1.0
x₂ = 2.0
x₃ = π/2

gradient_f(f₁, (x₁, x₂))
gradient_f(f₂, (x₁, x₂, x₃))

hessian_ff(f₁, (x₁, x₂))
hessian_ff(f₂, (x₁, x₂, x₃))

