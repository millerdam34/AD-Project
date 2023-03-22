include("reverse.jl")

f₁((x₁, x₂)) = log(x₁) + x₁ * x₂ - sin(x₂)
f₂((x₁, x₂, x₃)) = (x₁ * x₂ * sin(x₃) + exp(x₁ * x₂)) / x₃

x₁ = 1.0
x₂ = 2.0
x₃ = π/2

f₄((x₁,x₂)) = x₁^3+x₂^2



@time gradient_f(f₄, (1.0,2.0))
@time gradient_r(f₄, (1.0,2.0))

# 0.129204 seconds (243.83 k allocations: 12.681 MiB, 6.21% gc time, 99.42% compilation time)
# 2-element Vector{Float64}:
@time gradient_r(f₁, (x₁, x₂))

# 0.070658 seconds (155.43 k allocations: 8.152 MiB, 99.69% compilation time)
# 3-element Vector{Float64}:
@time gradient_r(f₂, (x₁, x₂, x₃))

# 0.167119 seconds (2.17 M allocations: 56.244 MiB, 5.43% gc time, 51.08% compilation time)
# 100000-element Vector{Float64}:
f₃(x) = prod(x)
x = map(x -> x + 0.5, rand(Float64, 100000))
@time gradient_r(f₃, x)

# 0.991577 seconds (20.00 M allocations: 472.999 MiB, 22.68% gc time)
# 1000000-element Vector{Float64}:
x = map(x -> x + 0.5, rand(Float64, 1000000))
@time gradient_r(f₃, x)

# 9.660034 seconds (200.00 M allocations: 4.619 GiB, 22.23% gc time)
# 10000000-element Vector{Float64}:
x = map(x -> x + 0.5, rand(Float64, 10000000))
@time gradient_r(f₃, x)

f((x₁,x₂)) = 3.0 * x₁ + x₂

gradient_r(f, [1.0,2.0])