# Solid angles for `DebyeWolf`: areas of the Voronoi cells of the ray directions (in the
# plane of direction cosines), clipped to their convex hull, via DelaunayTriangulation.jl.

module OpticsBaseDelaunayTriangulationExt

using OpticsBase: OpticsBase
using DelaunayTriangulation: triangulate, voronoi, get_area

function OpticsBase._delaunay_dual_areas(x::AbstractVector{T},
        y::AbstractVector{T}) where {T <: AbstractFloat}
    points = [(x[j], y[j]) for j in eachindex(x, y)]
    tri = triangulate(points)
    vorn = voronoi(tri; clip = true)
    return T[get_area(vorn, j) for j in eachindex(points)]
end

end # module
