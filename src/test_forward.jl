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
x = map(x -> x + 0.5, rand(Float64, 10000))
@time gradient_f(f₃, x)

f((x,y)) = x * x * y * y

testin = [2.0, 4.0]
testdirs = [[1/2,sqrt(3)/2] [1/sqrt(2),1/sqrt(2)]]
testnorm = [[1.0,sqrt(3)] [1.0,1.0]]

n_derivative(f, [2.0,4.0], testdirs, true)
n_derivative(f, [2.0,4.0], testnorm)

f((x₁,x₂,x₃)) = x₁*x₂*x₃+exp(x₁)+sin(x₂)+x₃^2

input = [1.0,2.0,3.0]

directions = [[1.0,1.0,1.0] [1.0,0.0,1.0] [1.0,1.0,0.0] [1.0,1.0,1.0]]

@time n_derivative(f, input, directions, false)

f(x)=exp(x[1])

# stack limit
D1 = ones(1,16)

@time n_derivative(f, [2.0], D1, true)

D2 = ones(3,15)

@time n_derivative(f, input, D2, false)

f((x₁,x₂)) = 3 * x₁ + x₂
