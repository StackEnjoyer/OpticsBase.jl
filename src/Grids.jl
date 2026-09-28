"""
    RegularGrid{D,T} <: AbstractGrid{D}

Regular Cartesian sampling grid with `D` dimensions in the local coordinates of a port.

    RegularGrid(dims::NTuple{D,Integer}, spacing::NTuple{D,Real})

Creates a grid with `dims[d]` samples and sample spacing `spacing[d]` in \\[m\\] along
dimension `d`. The spacings are promoted to a common floating-point type `T`.

# Conventions

  - The grid lives in the port-local coordinates `(ξ, η[, ζ])`: for `D = 2` the dimensions
    run along the port axes `(u, v)`, for `D = 3` the third dimension `ζ` runs along the
    port normal `n`.
  - Sample `i` along dimension `d` sits at `(i − (n÷2 + 1))·Δ` with `n = dims[d]` and
    `Δ = spacing[d]`, i.e. the port origin is at sample `n÷2 + 1` (the center used by
    `fftshift`) for even and odd `n` alike. The same rule applies to `ζ`.
  - The grid has no offset of its own; to move the samples, move the port.

# Fields

  - `dims`: number of samples per dimension, all `≥ 1`
  - `spacing`: sample spacing per dimension in \\[m\\], all finite and `> 0`

# Interface

Implements the [`AbstractGrid`](@ref) interface: `size(grid)`, [`spacing`](@ref) and
[`coordinates`](@ref).
"""
struct RegularGrid{D, T <: Real} <: AbstractGrid{D}
    dims::NTuple{D, Int}
    spacing::NTuple{D, T}

    function RegularGrid{D, T}(dims::NTuple{D, Integer},
            spacing::NTuple{D, Real}) where {D, T <: Real}
        D ≥ 1 || throw(ArgumentError("RegularGrid needs at least one dimension"))
        all(≥(1), dims) ||
            throw(ArgumentError("RegularGrid: all dims must be ≥ 1, got $dims"))
        all(s -> isfinite(s) && s > 0, spacing) ||
            throw(ArgumentError("RegularGrid: all spacings must be finite and > 0, got $spacing"))
        return new{D, T}(Int.(dims), T.(spacing))
    end
end

function RegularGrid(dims::NTuple{D, Integer}, spacing::NTuple{D, Real}) where {D}
    T = float(promote_type(map(typeof, spacing)...))
    return RegularGrid{D, T}(dims, spacing)
end

"""
    size(grid::RegularGrid) -> NTuple{D,Int}
    size(grid::RegularGrid, d)

Number of samples per dimension of `grid` (along `dim d` for the second form).
"""
Base.size(g::RegularGrid) = g.dims
Base.size(g::RegularGrid, d::Integer) = g.dims[d]

"""
    spacing(grid::AbstractGrid{D}) -> NTuple{D,Real}

Returns the sample spacing of `grid` per dimension in \\[m\\], in the port-local
coordinates `(ξ, η[, ζ])`. Part of the [`AbstractGrid`](@ref) interface.
"""
function spacing end

spacing(g::RegularGrid) = g.spacing

"""
    coordinates(grid::AbstractGrid, d::Integer) -> AbstractVector

Returns the port-local coordinates in \\[m\\] of the samples of `grid` along dimension `d`
(`ξ` along `u` for `d = 1`, `η` along `v` for `d = 2`, `ζ` along `n` for `d = 3`). Part of
the [`AbstractGrid`](@ref) interface.

For a [`RegularGrid`](@ref) the result is a range with sample `i` at `(i − (n÷2 + 1))·Δ`,
so sample `n÷2 + 1` is exactly `0` (the port origin).
"""
function coordinates end

function coordinates(g::RegularGrid{D, T}, d::Integer) where {D, T}
    1 ≤ d ≤ D || throw(BoundsError(g, d))
    n = g.dims[d]
    # ref value 0 at index n÷2 + 1, so element i is exactly (i − (n÷2 + 1))·Δ
    return StepRangeLen(zero(T), g.spacing[d], n, n ÷ 2 + 1)
end
