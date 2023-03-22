include("hessian.jl")

f₂((x₁, x₂, x₃)) = (x₁ * x₂ * sin(x₃) + exp(x₁ * x₂)) / x₃

v = [1.0,2.0,π/2]

@time hessian_ff(f₂, v)
@time hessian_ff!(f₂, v)
@time hessian_rf(f₂, v)
@time hessian_rf!(f₂, v)
@time hessian_fr(f₂, v)

# compare
# f(x) = prod(rand(1000))

f(x) = prod(x)
x = map(x -> x + 0.5, rand(Float64, 1000))

@time hessian_ff(f, x)
@time hessian_ff!(f, x)
@time hessian_rf(f, x)
@time hessian_rf!(f, x)
@time hessian_fr(f, x)


# hessian_ff
# 1.595040 seconds (194.49 k allocations: 17.863 MiB, 6.37% compilation time)
# 1000×1000 Matrix{Float64}:
# 1505.325857 seconds (6 allocations: 763.550 MiB, 0.00% gc time)
# 10000×10000 Matrix{Float64}:

# hessian_rf
# 1.100445 seconds (20.70 M allocations: 805.963 MiB, 20.28% gc time, 14.34% compilation time)
# 1000×1000 Matrix{Float64}:
# 107.864502 seconds (2.18 G allocations: 78.000 GiB, 9.97% gc time)
# 10000×10000 Matrix{Float64}:
# 72.990674 seconds (2.18 G allocations: 78.012 GiB, 5.94% gc time, 0.20% compilation time)
# 10000×10000 Matrix{Float64}:

# hessian_rf!
# 0.653308 seconds (20.46 M allocations: 777.741 MiB, 5.75% gc time)
# 1000×1000 Matrix{Float64}:
# 159.001224 seconds (3.29 G allocations: 91.488 GiB, 8.10% gc time)
# 10000×10000 Matrix{Float64}:

@time gradient_f(f₂, v)
@time gradient_r(f₂, v)

f′((x₁, x₂, x₃)) = (x₂ * sin(x₃) + x₂ * exp(x₁ * x₂)) / x₃

@time forward_reverse(f₂, v, [1, 0, 0], true)
@time reverse_forward(f₂, v, [1, 0, 0], true)

gradient_f(f′, v)

gradient_r(f′, v)