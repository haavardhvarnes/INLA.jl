"""
    convex_hull_polygon(points) -> Matrix{Float64}

Return the convex-hull polygon of `points` (`n × 2` matrix) as a
`k × 2` matrix in counter-clockwise order, without repeating the first
vertex at the end. Uses `DelaunayTriangulation.convex_hull`, which
emits CCW-ordered vertices.

Collinear points are excluded from the hull vertex set. Throws
`ArgumentError` if the hull is degenerate (fewer than three distinct
vertices).
"""
function convex_hull_polygon(points::AbstractMatrix{<:Real})
    size(points, 2) == 2 ||
        throw(ArgumentError("points must be n × 2; got size $(size(points))"))
    size(points, 1) >= 3 ||
        throw(ArgumentError("need at least 3 points for a convex hull; got $(size(points, 1))"))

    pts_tup = [(Float64(points[i, 1]), Float64(points[i, 2])) for i in axes(points, 1)]
    ch = DelaunayTriangulation.convex_hull(pts_tup)
    verts = DelaunayTriangulation.get_vertices(ch)
    # DT returns a closed loop with verts[1] == verts[end]; strip the duplicate.
    idx = @view verts[1:(end - 1)]
    length(idx) >= 3 ||
        throw(ArgumentError("convex hull is degenerate (collinear or duplicate points)"))

    P = Matrix{Float64}(undef, length(idx), 2)
    for (k, v) in enumerate(idx)
        P[k, 1], P[k, 2] = pts_tup[v]
    end
    return P
end

"""
    expand_polygon(polygon, offset) -> Matrix{Float64}

Expand a simple convex polygon outward by a perpendicular distance
`offset`. `polygon` is a `k × 2` matrix of CCW-ordered vertices (no
closing duplicate). Each vertex is shifted along the outward bisector
such that each edge of the output polygon lies `offset` away from the
corresponding input edge.

Used to construct the outer mesh domain from the convex hull of the
observation locations (R-INLA's `offset` argument to `inla.mesh.2d`).
"""
function expand_polygon(polygon::AbstractMatrix{<:Real}, offset::Real)
    size(polygon, 2) == 2 ||
        throw(ArgumentError("polygon must be k × 2; got size $(size(polygon))"))
    n = size(polygon, 1)
    n >= 3 || throw(ArgumentError("polygon needs ≥ 3 vertices; got $n"))
    offset >= 0 || throw(ArgumentError("offset must be non-negative; got $offset"))

    iszero(offset) && return Matrix{Float64}(polygon)

    P = Matrix{Float64}(undef, n, 2)
    for i in 1:n
        prev = i == 1 ? n : i - 1
        nxt = i == n ? 1 : i + 1
        pc = (Float64(polygon[i, 1]), Float64(polygon[i, 2]))
        pp = (Float64(polygon[prev, 1]), Float64(polygon[prev, 2]))
        pn = (Float64(polygon[nxt, 1]), Float64(polygon[nxt, 2]))

        # Outward normals of edges (pp → pc) and (pc → pn) on a CCW polygon
        # are (dy, -dx) of each edge direction.
        e_in = (pc[1] - pp[1], pc[2] - pp[2])
        e_out = (pn[1] - pc[1], pn[2] - pc[2])
        len_in = hypot(e_in[1], e_in[2])
        len_out = hypot(e_out[1], e_out[2])
        (len_in > 0 && len_out > 0) ||
            throw(ArgumentError("polygon has coincident consecutive vertices"))
        n_in = (e_in[2] / len_in, -e_in[1] / len_in)
        n_out = (e_out[2] / len_out, -e_out[1] / len_out)

        # Shift along the sum-of-normals direction by the magnitude that
        # pushes each edge outward by exactly `offset`. Derivation:
        # |n_in + n_out|² = 2 (1 + n_in · n_out) = 4 cos²(α), with α the
        # half exterior angle. Correct shift is 2 · offset · (n_in + n_out)
        # / |n_in + n_out|².
        bx = n_in[1] + n_out[1]
        by = n_in[2] + n_out[2]
        bn2 = bx^2 + by^2
        bn2 > 0 ||
            throw(ArgumentError("polygon bisector is degenerate (spike or reflex corner)"))
        s = 2 * offset / bn2
        P[i, 1] = pc[1] + s * bx
        P[i, 2] = pc[2] + s * by
    end
    return P
end

"""
    cutoff_dedup(points, cutoff) -> Matrix{Float64}

Remove near-duplicate rows of `points` so that the returned set has
pairwise distance ≥ `cutoff`. First occurrence is kept (greedy). Runs
in `O(n²)` — adequate for mesh vertex counts; for 10⁶-point clouds a
KD-tree version would replace this.

`cutoff ≤ 0` returns `Matrix{Float64}(points)` unchanged.
"""
function cutoff_dedup(points::AbstractMatrix{<:Real}, cutoff::Real)
    size(points, 2) == 2 ||
        throw(ArgumentError("points must be n × 2; got size $(size(points))"))
    cutoff <= 0 && return Matrix{Float64}(points)

    kept = Int[]
    n = size(points, 1)
    for i in 1:n
        xi = Float64(points[i, 1])
        yi = Float64(points[i, 2])
        ok = true
        for j in kept
            if hypot(xi - points[j, 1], yi - points[j, 2]) < cutoff
                ok = false
                break
            end
        end
        ok && push!(kept, i)
    end
    P = Matrix{Float64}(undef, length(kept), 2)
    for (k, i) in enumerate(kept)
        P[k, 1] = points[i, 1]
        P[k, 2] = points[i, 2]
    end
    return P
end

"""
    nonconvex_hull_polygon(points; α) -> Matrix{Float64}

Return an α-shape (Edelsbrunner-Mücke) of `points` (`n × 2`) as a CCW
polygon, without repeating the first vertex. Larger `α` yields the
convex hull in the limit; smaller `α` follows the point cloud more
tightly.

# Arguments

- `points::AbstractMatrix{<:Real}` — `n × 2` point cloud.
- `α::Real` — circumradius bound. A Delaunay triangle is *kept* iff its
  circumradius is `≤ α`. If `α` is omitted it defaults to twice the
  median nearest-neighbour distance, a heuristic that gives a
  "reasonable concavity" on roughly uniform point clouds.

# Errors

- `ArgumentError` if `points` has fewer than 3 rows or is not 2D.
- `ArgumentError` if the alpha-shape is empty (`α` too small) or
  multi-component (the simply-connected boundary tracer fails); raise
  `α` and try again.

See also: [`convex_hull_polygon`](@ref) for the convex-hull case.
"""
function nonconvex_hull_polygon(points::AbstractMatrix{<:Real};
        α::Real=_default_alpha(points))
    size(points, 2) == 2 ||
        throw(ArgumentError("points must be n × 2; got size $(size(points))"))
    size(points, 1) >= 3 ||
        throw(ArgumentError("need at least 3 points for an alpha-shape; got $(size(points, 1))"))
    α > 0 ||
        throw(ArgumentError("α must be positive; got $α"))

    pts = [(Float64(points[i, 1]), Float64(points[i, 2])) for i in axes(points, 1)]
    tri = DelaunayTriangulation.triangulate(pts)

    α2 = Float64(α)^2
    kept_tris = NTuple{3, Int}[]
    for T in DelaunayTriangulation.each_solid_triangle(tri)
        i, j, k = T
        if _circumradius_sq(pts[i], pts[j], pts[k]) ≤ α2
            push!(kept_tris, (i, j, k))
        end
    end
    isempty(kept_tris) &&
        throw(ArgumentError("alpha-shape is empty for α=$α; raise α"))

    return _trace_alpha_boundary(pts, kept_tris)
end

# Squared circumradius of a triangle with vertices a, b, c. Using the
# squared form avoids one sqrt and lets the comparison `r² ≤ α²` stay
# numerically clean. Formula: R = a·b·c / (4·area), so
# R² = (a²·b²·c²) / (16·area²).
function _circumradius_sq(
        a::NTuple{2, Float64}, b::NTuple{2, Float64}, c::NTuple{2, Float64})
    ab2 = (a[1] - b[1])^2 + (a[2] - b[2])^2
    bc2 = (b[1] - c[1])^2 + (b[2] - c[2])^2
    ca2 = (c[1] - a[1])^2 + (c[2] - a[2])^2
    cross = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    area2 = 0.25 * cross^2
    iszero(area2) && return Inf       # collinear → infinite circumradius
    return ab2 * bc2 * ca2 / (16.0 * area2)
end

# Default α heuristic: twice the median nearest-neighbour distance over
# the point cloud. O(n²) — fine for mesh-design point clouds (≤ a few
# thousand). Larger clouds: pass α explicitly.
function _default_alpha(points::AbstractMatrix{<:Real})
    n = size(points, 1)
    n >= 2 || return 1.0
    nn = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        best = Inf
        for j in 1:n
            i == j && continue
            d2 = (points[i, 1] - points[j, 1])^2 + (points[i, 2] - points[j, 2])^2
            d2 < best && (best = d2)
        end
        nn[i] = sqrt(best)
    end
    sort!(nn)
    med = nn[(n + 1) ÷ 2]
    return 2.0 * med
end

# Trace the boundary of the kept-triangle set into a single CCW polygon.
# Boundary edges are those incident to exactly one kept triangle; we
# build a `next[v] = w` table from the directed boundary edges (each
# kept triangle contributes its outward-oriented edges) so the walk
# follows the boundary in order without ambiguity at junction points.
function _trace_alpha_boundary(pts::Vector{NTuple{2, Float64}},
        kept_tris::Vector{NTuple{3, Int}})
    edge_count = Dict{NTuple{2, Int}, Int}()
    for (i, j, k) in kept_tris
        for (a, b) in ((i, j), (j, k), (k, i))
            key = a < b ? (a, b) : (b, a)
            edge_count[key] = get(edge_count, key, 0) + 1
        end
    end

    # Directed boundary edges: each kept triangle's three CCW-oriented
    # edges contribute, but only the ones whose undirected key has
    # count == 1 are on the alpha-shape boundary.
    nxt = Dict{Int, Int}()
    for (i, j, k) in kept_tris
        for (a, b) in ((i, j), (j, k), (k, i))
            key = a < b ? (a, b) : (b, a)
            if edge_count[key] == 1
                haskey(nxt, a) &&
                    throw(ArgumentError("alpha-shape boundary is not simply connected; raise α"))
                nxt[a] = b
            end
        end
    end

    isempty(nxt) &&
        throw(ArgumentError("alpha-shape produced no boundary edges; raise α"))

    start = first(keys(nxt))
    loop = Int[start]
    v = nxt[start]
    while v != start
        push!(loop, v)
        haskey(nxt, v) ||
            throw(ArgumentError("alpha-shape boundary chain broke at vertex $v; raise α"))
        v = nxt[v]
        length(loop) > length(nxt) &&
            throw(ArgumentError("alpha-shape boundary did not close; raise α"))
    end
    length(loop) == length(nxt) ||
        throw(ArgumentError("alpha-shape boundary is not simply connected; raise α"))

    P = Matrix{Float64}(undef, length(loop), 2)
    for (k, v) in enumerate(loop)
        P[k, 1] = pts[v][1]
        P[k, 2] = pts[v][2]
    end

    # DT triangles are CCW, so the outward boundary chain inherits CCW
    # orientation; check via signed area and reverse if a degenerate
    # input slipped through.
    if _signed_area(P) < 0
        reverse!(@view P[:, 1])
        reverse!(@view P[:, 2])
    end
    return P
end

function _signed_area(P::AbstractMatrix{<:Real})
    n = size(P, 1)
    s = 0.0
    @inbounds for i in 1:n
        j = i == n ? 1 : i + 1
        s += P[i, 1] * P[j, 2] - P[j, 1] * P[i, 2]
    end
    return 0.5 * s
end

"""
    subdivide_polygon(polygon, max_edge) -> Matrix{Float64}

Insert Steiner points along each edge of `polygon` (a `k × 2` CCW
matrix) so that every consecutive pair of returned vertices is at most
`max_edge` apart. Returns the augmented polygon (no closing duplicate).

Used by `inla_mesh_2d(...; subdivide_boundary = true)` to enforce a
strict per-edge `max_edge` bound on the input boundary, before Ruppert
refinement runs over the interior. Without this, the boundary edges
survive at their input length and the area-based Ruppert bound only
softly limits interior edges.
"""
function subdivide_polygon(polygon::AbstractMatrix{<:Real}, max_edge::Real)
    size(polygon, 2) == 2 ||
        throw(ArgumentError("polygon must be k × 2; got size $(size(polygon))"))
    n = size(polygon, 1)
    n >= 3 || throw(ArgumentError("polygon needs ≥ 3 vertices; got $n"))
    max_edge > 0 ||
        throw(ArgumentError("max_edge must be positive; got $max_edge"))

    out = Vector{NTuple{2, Float64}}()
    sizehint!(out, n)
    for i in 1:n
        ax = Float64(polygon[i, 1])
        ay = Float64(polygon[i, 2])
        push!(out, (ax, ay))
        j = i == n ? 1 : i + 1
        bx = Float64(polygon[j, 1])
        by = Float64(polygon[j, 2])
        d = hypot(bx - ax, by - ay)
        if d > max_edge
            n_seg = ceil(Int, d / max_edge)
            for s in 1:(n_seg - 1)
                t = s / n_seg
                push!(out, (ax + t * (bx - ax), ay + t * (by - ay)))
            end
        end
    end

    P = Matrix{Float64}(undef, length(out), 2)
    for (k, p) in enumerate(out)
        P[k, 1] = p[1]
        P[k, 2] = p[2]
    end
    return P
end
