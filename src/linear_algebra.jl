using GeometryBasics
using Plots
using LinearAlgebra
using RationalRoots

"""
    distance_2_points(p::Point, q::Point) -> Norm(v)
Distance between two points
"""
function distance_2_points(p::Point, q::Point)
    v = q - p
    norm(v)
end

"""
    center of gravity(p::Point, q::Point, t) → r
Creates center of gravity of points `p` and `q`using parametric equation of a line
"""
function center_of_gravity(p::Point, q::Point, t)
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

"""
    vector_angle_cos(p::Point, q::Point) -> cos θ
Calculate cfosine of angle between 2 vectors using cosines
"""
function vector_angle_cos(p::Point, q::Point)
    s = dot(p, q) / (norm(p) * norm(q))
end

"""
    is_orthogonal(p::Point, q::Point) -> bool
Check if two vectors are is_orthogonal
"""
function is_orthogonal(p::Point, q::Point)
    dot(p, q) == 0
end

"""
    function polar_unit(y::Vector) -> Vector{Float64}
Return unit vector in polar form for vector `y`
"""
function polar_unit(y::Vector)
    [(y[1] / norm(y)), (y[2] / norm(y))]
end

"""
    orthproj(v::Vector, w::Vector) -> Vector
Orthogonol projection of vector `w` onto `v`
"""
function orthproj(v::Vector, w::Vector)
    u = (dot(v, w) / norm(v)^2) * v
    [round(Int, x) for x in u]
end

"""
    point_in_implicit_line(p::Point, q::Point, x::Point) -> Float64
The orthogonal vector α is calculated as:
    v = Vector(q - p)
    α = [v[2], -v[1]]
The implicit equation of the line is:
    α * (x -p) = 0
    a = α[1]
    b = α[2]
    c = -(a * p[1]) - (b * p[2])
Then calculate the distance of point `x` from the line using the formula:
    (a * x[1] + b * x[2] + c) / norm([a, b]
If 0 the point is in the line!
"""
function point_in_implicit_line(p::Point, q::Point, x::Point)
    a = q[2] - p[2]
    b = p[1] - q[1]
    c = -(a * p[1]) - (b * p[2])
    (a * x[1] + b * x[2] + c) / norm([a, b])
end

"""
    parametric_to_implicit_line(p::Point, v::Vector) -> Tuple{Int64, Int64, Int64}
The parametric equation is:
    l : l(t) = p + tv
Use p and v to calculate:
    l : ax1 + bx2 + c = 0.
And return co-efficient as tuple
"""
function parametric_to_implicit_line(p::Point, v::Vector)
    a = -v[2]
    b = v[1]
    c = -(a * p[1]) - (b * p[2])
    (a, b, c)
end

"""
   implicit_to_parametric line(a::Int64, b::int64, c::int64) -> Tuple(Vector, Point)
"""
function implicit_to_parametric_line(a::Number, b::Number, c::Number)
    v = [b, -a]
    if abs(a) > abs(b)
        # the linter gives error on this constructor so shut off lint.call setting
        p = Point(- c / a, 0)
    else
        p = Point(0, -c / b)
    end
    (v, p)
end

"""
    explicit_line(p::Point, q::Point) -> Tuple(Float64, Float64)
The orthogonal vector α is calculated as:
    v = Vector(q - p)
    α = [v[2], -v[1]]
The explicit equation of the line requires slope & intercept:
    a = -v[2]
    b = v[1]
    slope = -a / b
    c = -(a * p[1]) - (b * p[2])
    intercept = -c / b
"""
function explicit_line(p::Point, q::Point)
    a = q[2] - p[2]
    b = p[1] - q[1]
    c = -(a * p[1]) - (b * p[2])
    slope = -a / b
    intercept = -c / b
    (slope, intercept)
end

"""
   distance_to_implicit_line(a::Int64, b::int64, c::int64, r::Point) -> Float64
"""
function distance_to_implicit_line(a::Number, b::Number, c::Number, r::Vector)
    v = [a, b]
    d = ((a * r[1]) + (b * r[2]) + c) / norm(v)
end

"""
   implicit_line_point_normal_form(a::Int64, b::int64, c::int64) -> Tuple(RationalRoot, RationalRoot, RationalRoot)
"""
function implicit_line_point_normal_form(a::Number, b::Number, c::Number)
    v = [a, b]
    â = a/RationalRoot(norm(v))
    b̂ = b/RationalRoot(norm(v))
    ĉ = c/RationalRoot(norm(v))
    (â, b̂, ĉ)
end

"""
   distance_to_pnf_implicit_line(â::RationalRoot, b̂::RationalRoot, ĉ::RationalRoot, r::Point) -> Float64
"""
function distance_to_pnf_implicit_line(â::RationalRoot, b̂::RationalRoot, ĉ::RationalRoot, r::Vector)
    d = ((â * r[1]) + (b̂ * r[2]) + ĉ)
end