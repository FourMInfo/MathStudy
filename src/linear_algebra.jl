using GeometryBasics
using Plots
"""
    param_line(p::Point, q::Point, n::Int64) → Ps::[Point]
Creates n points on a line defined by p and q, using the parametric equation of a line, then plots
"""
function param_line(p::Point, q::Point, n::Int64)
    v = q - p
    Ps = Point[]
    for i = 1:n
        t = i / n
        l = p + (t * v)
        s = scatter!(l, legend=false)
        push!(Ps, l)
    end
    s = plot!(Ps, legend=false)
    display(s)
    Ps
end
