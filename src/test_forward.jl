include("forward.jl")

f₁((x₁,x₂)) = log(x₁) + x₁*x₂ - sin(x₂)
f₂((x₁, x₂, x₃)) = (x₁ * x₂ * sin(x₃) + exp(x₁ * x₂)) / x₃


x₁ = 1.0
x₂ = 2.0
x₃ = π/2

# Vector{Dual{Float64}}  0.028514 seconds (66.02 k allocations: 3.454 MiB, 95.82% compilation time)
# 2-element Vector{Float64}:
@time gradient_f(f₁, (x₁, x₂))

# Vector{Dual{Float64}}  0.051784 seconds (98.94 k allocations: 5.120 MiB, 96.51% compilation time)
# 3-element Vector{Float64}:
@time gradient_f(f₂, (x₁, x₂, x₃))

# Vector{Dual{Float64}} 21.498908 seconds (567.11 k allocations: 32.850 MiB, 1.86% compilation time)
# 100000-element Vector{Float64}:
f₃(x) = prod(x)
x = map(x -> x + 0.5, rand(Float64, 100000))
@time gradient_f(f₃, x)
