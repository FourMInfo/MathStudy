using GeometryBasics
using Plots
using LinearAlgebra

"""
    distance_2_points(p::Point, q::Point) -> Norm(v)
Distance between two points
"""
function distance_2_points(p::Point, q::Point)
    v = q - p
    norm(v)
end

"""
    center of gravity(p::Point, q::Point, t::Float64) → r
Creates center of gravity of points `p` and `q`using parametric equation of a line
"""
function center_of_gravity(p::Point, q::Point, t::Float64)
    v = q - p
    r = p + (t * v)
end

"""
    barycentric_coord(p::Point, q::Point, r::Point) -> t::Float64
Calculate the barycentric co-ordinate or parater of points p and q with center of gravity point r
"""
function barycentric_coord(p::Point, q::Point, r::Point)
    t = norm(r - p) / norm(q - r)
end

"""
    plot_param_line(p::Point, q::Point, n::Int64) → [Point]
Creates `n` points on a line defined by `p` and `q`, using the parametric equation of a line, then plot
"""

function plot_param_line(p::Point, q::Point, n::Int64)
    Ps = Point[]
    for i = 1:n
        t = i / n
        r = center_of_gravity(p, q, t)
        s = scatter!(r, legend=false)
        push!(Ps, r)
    end
    s = plot!(Ps, legend=false)
    display(s)
    Ps
end
