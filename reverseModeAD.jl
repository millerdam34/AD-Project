module reverseAD
mutable struct Node
    val::Float64
    # parent index
    par1::Int64
    par2::Int64
    # adjoint ∂f/∂v
    adj::Float64
    der1::Float64
    der2::Float64
    index::Int64
end

mutable struct Tape
    tape::Vector{Node}
    size::Int64
end

Base.:log(f::Node) = log_AD(f, tape)
Base.:*(f::Node, g::Node) = mult_AD(f,g, tape)
Base.:sin(f::Node) = sin_AD(f, tape)
Base.:+(f::Node, g::Node) = add_AD(f,g, tape)
Base.:-(f::Node,g::Node) = sub_AD(f,g, tape)
Base.:/(f::Node,g::Node) = div_AD(f,g, tape)
Base.:exp(f::Node) = exp_AD(f,tape)

function exp_AD(f::Node, tape)::Node
    ans = Node(exp(f.val), f.index, -1, 0, exp(f.val), 0, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function div_AD(f::Node, g::Node, tape)::Node
    ans = Node(f.val/g.val, f.index, g.index, 0, 1/g.val, -f.val/g.val^2, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function log_AD(f::Node, tape)::Node
    ans = Node(log(f.val), f.index, -1, 0, 1/f.val, 0, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function mult_AD(f::Node, g::Node, tape)::Node
    ans = Node(f.val * g.val, f.index, g.index, 0, g.val, f.val, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function sin_AD(f::Node, tape)::Node
    ans = Node(sin(f.val), f.index, -1, 0, cos(f.val), 0, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function add_AD(f::Node, g::Node, tape)
    ans = Node(f.val + g.val, f.index, g.index, 0, 1, 1, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function sub_AD(f::Node, g::Node, tape)
    ans = Node(f.val - g.val, f.index, g.index, 0, 1, -1, tape.size)
    tape.tape[tape.size] = ans
    tape.size += 1
    return ans
end

function rewind(tape::Tape)
    tape.tape[tape.size - 1].adj = 1
    for i in tape.size-1:-1:1
        if (tape.tape[i].par1 != -1)
            tape.tape[tape.tape[i].par1].adj += tape.tape[i].adj * tape.tape[i].der1
        end
        if (tape.tape[i].par2 != -1)
            tape.tape[tape.tape[i].par2].adj += tape.tape[i].adj * tape.tape[i].der2
        end
    end
end

x_1 = Node(1.0, -1, -1, 0.0, 0.0, 0.0, 1)
x_2 = Node(2.0, -1, -1, 0.0, 0.0, 0.0, 2)
x_3 = Node(π/2, -1, -1, 0.0, 0.0, 0.0, 3)

x = 10

tape = Tape(Array{Node,1}(undef, x), 1)

tape.tape[tape.size] = x_1
tape.size += 1
tape.tape[tape.size] = x_2
tape.size += 1
tape.tape[tape.size] = x_3
tape.size += 1

y=(x_1*x_2*sin(x_3)+exp(x_1*x_2))/x_3


rewind(tape)

# ∂y/∂x_1
tape.tape[1].adj
# ∂y/∂x_2
tape.tape[2].adj
# ∂y/∂x_2
tape.tape[3].adj
end

