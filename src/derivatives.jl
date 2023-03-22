include("forward.jl")
include("reverse.jl")

function forward_reverse(f, input::Vector, direction::Vector, normal::Bool=false)::Vector
    if (!normal) 
        mag = norm(direction)
        direction = map(x -> x / mag, direction)
    end
    n = length(input)
    type = typeof(input[1])
    gradient_of_derivative = Array{type}(undef, n)

    # count length of f and f_prime tapes
    input_counter = fill(Dual(Counter(0), Counter(0)), n)
    input_count = f(input_counter)
    input_count_f = input_count.val.count + input_count.der.count + 2n

    # initialize tape
    tape_struct = Tape(Array{Node,1}(undef, input_count_f), 0)

    # initialize input vector of Dual(Node, Node)
    input_dual = Array{Dual,1}(undef, n)

    # initialize first 2n elements of f_tape
    for i in 1:n
        x = input[i]

        # f inputs 
        node1 = Node(x, -1, -1, zero(x), zero(x), zero(x), i, tape_struct, 0)
        tape_struct.tape[i] = node1

        # directional inputs for f
        # *one() in case input is not same type as direction components
        node2 = Node(direction[i] * one(x), -1, -1, zero(x), zero(x), zero(x), i + n, tape_struct, 0)
        tape_struct.tape[i + n] = node2

        # increment sizes
        tape_struct.size += 2
        
        # add to input_dual for generation of f′_tape
        input_dual[i] = Dual(node1, node2)
    end

    f(input_dual[1:n])

    rewind(tape_struct)

    for i in 1:n
        gradient_of_derivative[i] = tape_struct.tape[i].adj
    end

    return gradient_of_derivative
end

function reverse_forward(f, input, direction, normal::Bool=false)::Vector
    if (!normal)
        mag = norm(direction)
        direction = map(x -> x / mag, direction)
    end
    n = length(input)

    input_dual = Array{Dual}(undef, n)

    for i in 1:n
        input_dual[i] = Dual(input[i], direction[i] * one(input[i]))
    end

    return map(x -> x.der, gradient_r(f, input_dual))
end
