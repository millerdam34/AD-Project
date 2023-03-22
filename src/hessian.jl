include("derivatives.jl")

# straightforward implementation of calculating every second order derivative of a function
# using minimal external function calls
# input: x⃗
# f: f(x⃗)
# Dual(Dual(f(x⃗), ∂f/∂e⃗ⱼ), Dual(∂f/∂e⃗ᵢ, ∂²f/∂e⃗ᵢ∂e⃗ⱼ))
# output: Mᵢⱼ =  ∂²f/∂e⃗ᵢ∂e⃗ⱼ 
function hessian_ff!(f, input)::Matrix
    # init values
    n = length(input)
    type = typeof(input[1])
    one_val = one(type)
    zero_val = zero(type)
    one_dual = Dual(one_val, zero_val)
    zero_dual = Dual(zero_val, zero_val)
    hessian = Matrix{type}(undef, n, n)

    # cast input to Dual{Dual} and collect tuple to vector for indexing
    input_dual_dual = collect(map(x -> Dual(Dual(x, zero_val), zero_dual), input))

    # for each direction
    for i in 1:n
        # set direction derivative for second direction
        temp_i = input_dual_dual[i]
        input_dual_dual[i] = Dual(Dual(input[i], one_val), zero_dual)

        # for each direction
        for j in i:n
            # set direction derivative for first direction
            temp_j = input_dual_dual[j]
            input_dual_dual[j] = Dual(temp_j.val, one_dual)

            # update hessian
            hessian[i,j] = f(input_dual_dual).der.der
            if (i != j)
                hessian[j,i] = hessian[i,j]
            end

            # reset first direction
            input_dual_dual[j] = temp_j
        end
        # reset second direction
        input_dual_dual[i] = temp_i
    end
    return hessian
end

# straightforward implementation of calculating every second order derivative of a function
# using second derivative function call
# input: x⃗
# f: f(x⃗)
# Dual(Dual(f(x⃗), ∂f/∂e⃗ⱼ), Dual(∂f/∂e⃗ᵢ, ∂²f/∂e⃗ᵢ∂e⃗ⱼ))
# output: Mᵢⱼ =  ∂²f/∂e⃗ᵢ∂e⃗ⱼ 
function hessian_ff(f, input)::Matrix
    # init values
    n = length(input)
    type = typeof(input[1])
    hessian = Matrix{type}(undef, n, n)
    direction1 = zeros(n)
    direction2 = zeros(n)

    # for each direction
    for i in 1:n
        # set direction
        direction1[i] = 1;
        # for each direction
        for j in i:n
            # set direction
            direction2[j] = 1;
            # call second derivative with normal direction and update hessian
            hessian[i,j] = second_derivative_f(f, input, direction1, direction2, false, false)
            if i != j
                hessian[j,i] = hessian[i,j]
            end
            # reset direction
            direction2[j] = 0
        end
        # reset direction
        direction1[i] = 0
    end
    return hessian
end

# calculates the hessian based on derivatives of gradients
# recreates tape for each direction
# input: input vector
# f: function 
# return: hessian 
function hessian_rf(f, input)::Matrix
    # init length, type and hessian
    n = length(input)
    type = typeof(input[1])
    hessian = Matrix{type}(undef, n, n)

    # cast inputs to duals collect tuple to vector
    input_dual = collect(map(x -> Dual(x, zero(x)), input))

    # for each elementary direction
    for i in 1:n
        # choose direction by setting der component to one
        temp = input_dual[i]
        input_dual[i] = Dual(temp.val, one(type))

        # get gradient of input with derivative components
        grad = gradient_r(f, input_dual)

        # update hessian
        for j in i:n
            hessian[i,j] = grad[j].der
            if i != j
                hessian[j,i] = hessian[i,j]
            end
        end

        # reset direction component
        input_dual[i] = temp
    end
    return hessian
end


# calculates the hessian based on derivatives of gradients
# tries to optimize by reusing tape and updating values
# input: input vector
# f: function 
# return: hessian 
function hessian_rf!(f, input)::Matrix
    # init length, type and hessian
    n = length(input)
    type = typeof(input[1])
    hessian = Matrix{type}(undef, n, n)
    # cast all inputs to dual numbers
    input_dual = collect(map(x -> Dual(x, zero(x)), input))
    # record abstract (independent of direction of second derivative) tape for obtaining gradient of f
    tape_struct = record(f, input_dual)

    for i in 1:n
        # define direction by updating dual of input node
        temp = tape_struct.tape[i].val
        tape_struct.tape[i].val = Dual(temp.val, one(type))

        for node in tape_struct.tape
            # reset adjoints from previous pass
            node.adj = Dual(zero(type), zero(type))
            # update non input nodes with new val and ders based on direction 
            if (node.index > n)
                update_node(node)
            end
        end

        # rewind tape to get gradient with specific direction derivative component
        rewind(tape_struct)

        # update hessian
        for j in i:n
            hessian[i,j] = tape_struct.tape[j].adj.der
            if i != j
                hessian[j,i] = hessian[i,j]
            end
        end

        # reset direction component
        tape_struct.tape[i].val = temp
    end
    return hessian
end



# WIP 
function hessian_fr(f, input)::Matrix
    n = length(input)
    type = typeof(input[1])
    hessian = Matrix{type}(undef, n, n)
    direction = zeros(n)

    for i in 1:n
        direction[i] = 1
        row = forward_reverse(f, input, direction, true)

        # store hessian values
        for j in i:n
            hessian[i, j] = row[j]
            if i != j
                hessian[j, i] = hessian[i, j]
            end
        end

        # reset direction
        direction[i] = 0
    end
    return hessian
end

