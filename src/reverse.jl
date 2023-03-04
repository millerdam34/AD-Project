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
end

Base.:log(f::Node) = log_AD(f)
Base.:*(f::Node, g::Node) = mult_AD(f,g)
Base.:sin(f::Node) = sin_AD(f)
Base.:+(f::Node, g::Node) = add_AD(f,g)
Base.:-(f::Node,g::Node) = sub_AD(f,g)
Base.:/(f::Node,g::Node) = div_AD(f,g)
Base.:exp(f::Node) = exp_AD(f)

function exp_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(exp(f.val), f.index, -1, zero(f.val), exp(f.val), zero(f.val), f.tape_struct.size - 1, f.tape_struct)
end
function div_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(f.val/g.val, f.index, g.index, zero(f.val), 1/g.val, (-1 * f.val)/g.val^2, f.tape_struct.size - 1, f.tape_struct)
end
function log_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(log(f.val), f.index, -1, zero(f.val), 1/f.val, zero(f.val), f.tape_struct.size - 1, f.tape_struct)
end
function mult_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(f.val * g.val, f.index, g.index, zero(f.val), g.val, f.val, f.tape_struct.size - 1, f.tape_struct)
end
function sin_AD(f::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(sin(f.val), f.index, -1, zero(f.val), cos(f.val), zero(f.val), f.tape_struct.size - 1, f.tape_struct)
end
function add_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] = Node(f.val + g.val, f.index, g.index, zero(f.val), one(f.val), one(f.val), f.tape_struct.size - 1, f.tape_struct)
end
function sub_AD(f::Node, g::Node)::Node
    f.tape_struct.size += 1
    return f.tape_struct.tape[f.tape_struct.size - 1] =  Node(f.val - g.val, f.index, g.index, zero(f.val), one(f.val), -1 * one(f.val), f.tape_struct.size - 1, f.tape_struct)
end

function rewind(tape_struct::Tape)
    tape_struct.tape[tape_struct.size - 1].adj = one(tape_struct.tape[1].adj)
    for i in tape_struct.size-1:-1:1
        if (tape_struct.tape[i].par1 != -1)
            tape_struct.tape[tape_struct.tape[i].par1].adj += tape_struct.tape[i].adj * tape_struct.tape[i].der1
        end
        if (tape_struct.tape[i].par2 != -1)
            tape_struct.tape[tape_struct.tape[i].par2].adj += tape_struct.tape[i].adj * tape_struct.tape[i].der2
        end
    end
end

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

function record(f, input)::Tape
    # TODO figure out way to estimate tape size or maintain dynamically
    TAPE_SIZE = 10
    n = length(input)
    tape_struct = Tape(Array{Node,1}(undef, TAPE_SIZE), 1)

    for x in input
        tape_struct.tape[tape_struct.size] = Node(x * one(x), -1, -1, zero(x), zero(x), zero(x), tape_struct.size, tape_struct)
        tape_struct.size += 1
    end

    f(tape_struct.tape)

    return tape_struct
end