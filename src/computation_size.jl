using Base

struct Counter
    count::Int64
end

Base.:+(x::Counter, y::Counter) = Counter(x.count + y.count + 1)
Base.:-(x::Counter, y::Counter) = Counter(x.count + y.count + 1)
Base.:*(x::Counter, y::Counter) = Counter(x.count + y.count + 1)
Base.:/(x::Counter, y::Counter) = Counter(x.count + y.count + 1)
Base.:sin(x::Counter) = Counter(x.count + 1)
Base.:cos(x::Counter) = Counter(x.count + 1)
Base.:exp(x::Counter) = Counter(x.count + 1)
Base.:log(x::Counter) = Counter(x.count + 1)
Base.:^(x::Counter, n::Int) = Counter(x.count + 2)

# constants
Base.:+(x::Number, y::Counter) = Counter(y.count + 2)
Base.:-(x::Number, y::Counter) = Counter(y.count + 2)
Base.:*(x::Number, y::Counter) = Counter(y.count + 2)
Base.:/(x::Number, y::Counter) = Counter(y.count + 2)

# counts the number of basic operations performed in the evaluation of a function
# and returns the number of nodes in a corresponding computation graph
function get_computation_size(f, input)::Int64
    # needed to init input vector
    n = length(input)
    # all initial values get 0 so when inputs are used in multiple operations we don't over count
    # then we add n to account for the inputs to get the total nodes in a computation graph
    return f(fill(Counter(0), n)).count + n
end