include("computation_size.jl")

mutable struct Tape{T}
    tape::Vector{T}
    size::Int64
end

mutable struct Node{T}
    val::T
    par1::Int64
    par2::Int64
    adj::T
    der1::T
    der2::T
    index::Int64
    tape_struct::Tape{Node}
    op::Int64
end

# unary
Base.:log(f::Node) = log_AD(f)
Base.:sin(f::Node) = sin_AD(f)
Base.:cos(f::Node) = cos_AD(f)
Base.:exp(f::Node) = exp_AD(f)

# binary
Base.:*(f::Node, g::Node) = mult_AD(f,g)
Base.:+(f::Node, g::Node) = add_AD(f,g)
Base.:-(f::Node, g::Node) = sub_AD(f,g)
Base.:/(f::Node, g::Node) = div_AD(f,g)
Base.:^(f::Node, g::Node) = pow_AD(f,g)

# constants
Base.:*(f::Number, g::Node) = insert_constant(f, g) * g
Base.:+(f::Number, g::Node) = insert_constant(f, g) + g
Base.:-(f::Number, g::Node) = insert_constant(f, g) - g
Base.:/(f::Number, g::Node) = insert_constant(f, g) / g
Base.:/(f::Node, g::Number) = f / insert_constant(g, f)
Base.:^(f::Number, g::Node) = insert_constant(f,g)^g
Base.:^(f::Node, g::Number) = f^insert_constant(g,f)

function insert_constant(f::Number, g::Node)
    g.tape_struct.size += 1
    return g.tape_struct.tape[g.tape_struct.size] = Node(f * one(g.val), -1, -1, zero(g.val), zero(g.val), zero(g.val), g.tape_struct.size, g.tape_struct, 0)
end

function pow_AD(f::Node,g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(f.val^g.val, f.index, g.index, zero(f.val), g.val * f.val ^ (g.val - 1), f.val^g.val * log(f.val), f.tape_struct.size, f.tape_struct, 1)
end
function exp_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(exp(f.val), f.index, -1, zero(f.val), exp(f.val), zero(f.val), f.tape_struct.size, f.tape_struct, 1)
end
function div_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(f.val/g.val, f.index, g.index, zero(f.val), 1/g.val, (-1 * f.val)/g.val^2, f.tape_struct.size, f.tape_struct, 2)
end
function log_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(log(f.val), f.index, -1, zero(f.val), 1/f.val, zero(f.val), f.tape_struct.size, f.tape_struct, 3)
end
function mult_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(f.val * g.val, f.index, g.index, zero(f.val), g.val, f.val, f.tape_struct.size, f.tape_struct, 4)
end
function sin_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(sin(f.val), f.index, -1, zero(f.val), cos(f.val), zero(f.val), f.tape_struct.size, f.tape_struct, 5)
end
function add_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(f.val + g.val, f.index, g.index, zero(f.val), one(f.val), one(f.val), f.tape_struct.size, f.tape_struct, 6)
end
function sub_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(f.val - g.val, f.index, g.index, zero(f.val), one(f.val), -1 * one(f.val), f.tape_struct.size, f.tape_struct, 7)
end
function cos_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size] = Node(cos(f.val), f.index, -1, zero(f.val), -sin(f.val), zero(f.val), f.tape_struct.size, f.tape_struct, 8)
end

# updates adjoint values from back to front of tape
function rewind(tape_struct::Tape)
    tape_struct.tape[tape_struct.size].adj = one(tape_struct.tape[1].adj)
    for i in tape_struct.size:-1:1
        if (tape_struct.tape[i].par1 != -1)
            tape_struct.tape[tape_struct.tape[i].par1].adj += tape_struct.tape[i].adj * tape_struct.tape[i].der1
        end
        if (tape_struct.tape[i].par2 != -1)
            tape_struct.tape[tape_struct.tape[i].par2].adj += tape_struct.tape[i].adj * tape_struct.tape[i].der2
        end
    end
end

# calculates gradient by recording and rewinding tape
function gradient_r(f, input)::Vector
    n = length(input)
    tape_struct = record(f, input)
    grad = Array{typeof(input[1]),1}(undef, n)

    rewind(tape_struct)

    for i in 1:n
        grad[i] = tape_struct.tape[i].adj
    end

    return grad
end

# records a tape based on input nodes and function
function record(f, input)::Tape
    TAPE_SIZE = get_computation_size(f, input)
    n = length(input)
    tape_struct = Tape(Array{Node,1}(undef, TAPE_SIZE), 0)

    for x in input
        tape_struct.size += 1
        tape_struct.tape[tape_struct.size] = Node(x * one(x), -1, -1, zero(x), zero(x), zero(x), tape_struct.size, tape_struct, 0)
    end

    f(tape_struct.tape[1:n])

    return tape_struct
end

# used to update reused tape nodes
function update_node(node::Node)
    op = node.op
    if op == 0
        return
    end
    val1 = node.tape_struct.tape[node.par1].val
    if node.par2 != -1
        val2 = node.tape_struct.tape[node.par2].val
    end
    type = typeof(val1)
    if op == 1
        node.val = exp(val1)
        # node.der1 = exp(val1)
        node.der1 = node.val
    elseif op == 2
        node.val = val1 / val2
        node.der1 = one(val1) / val2
        node.der2 = -val1 / (val2 * val2)
    elseif op == 3
        node.val = log(val1)
        node.der1 = one(val1) / val1
    elseif op == 4
        node.val = val1 * val2
        node.der1 = val2
        node.der2 = val1
    elseif op == 5
        node.val =  sin(val1)
        node.der1 = cos(val1)
    elseif op == 6
        node.val = val1 + val2
        node.der1 = one(val1)
        node.der2 = one(val2)
    elseif op == 7
        node.val = val1 - val2
        node.der1 = one(val1)
        node.der2 = -one(val2)
    elseif op == 8
        node.val = cos(val1)
        node.der1 = -sin(val1)
    end
    return
end
