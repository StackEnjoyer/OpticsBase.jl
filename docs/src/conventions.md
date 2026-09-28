```@meta
CurrentModule = OpticsBase
```

# Conventions

The following conventions are binding for every exchange format, converter and solver
that connects to `OpticsBase`. If a format deviates from one of them, its docstring must
say so explicitly.

!!! warning
    Mixing conventions silently corrupts phases, powers or polarization states across a
    solver chain. When in doubt, check this page first.

## Units

All quantities are in SI units internally: meters, seconds, watts. Wavelengths are
vacuum wavelengths in meters; the refractive index of the medium is a property of the
port. Support for [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities is
optional and only provided via a package extension; it converts at the API boundary.

## Time convention

Time-harmonic fields use the ``\exp(-i\omega t)`` convention. A plane wave is written as

```math
\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 \exp\big(i(\mathbf{k}\cdot\mathbf{r} - \omega t)\big).
```

Consequently, propagation over a distance ``d`` in a medium with refractive index ``n``
adds the phase ``+ n k_0 d`` with ``k_0 = 2\pi/\lambda``.

## Ports and frames

Every handover of a field between solvers happens at a port: a surface in global
coordinates with an origin and a local frame. There is no field handover without a port.
`OpticsBase` defines no global optical axis; solvers with a preferred axis map it onto
the port frame in their extension.

- The local frame ``(\mathbf{u}, \mathbf{v}, \mathbf{n})`` is orthonormal and
  right-handed, ``\mathbf{u} \times \mathbf{v} = \mathbf{n}``.
- The normal ``\mathbf{n}`` points downstream: fields at the port propagate into the
  half space ``\mathbf{n}\cdot\mathbf{k} > 0``.
- ``\mathbf{u}`` is always given explicitly. It fixes the polarization basis, so there is
  no default for it.
- Local coordinates ``(\xi, \eta, \zeta)`` are measured along
  ``(\mathbf{u}, \mathbf{v}, \mathbf{n})`` from the port origin.

See [`AbstractPort`](@ref) and [`PlanarPort`](@ref).

## Phase reference

- Stored complex amplitudes exclude only the temporal carrier ``\exp(-i\omega t)``. They
  contain the full spatial phase; no tilt or defocus is removed.
- Optical path lengths are absolute: they are counted from the reference point of the
  coherent component, which the first solver of a chain sets (typically the source), and
  every later solver continues the count. Fields that reach one port on different paths
  therefore keep their relative phase.
- The port origin ``\mathbf{r}_0`` is the reference point only for position-dependent
  phase terms such as ``\exp(i\mathbf{k}\cdot(\mathbf{r} - \mathbf{r}_0))``.
- For a ray, the phase is ``k_0 \cdot \mathrm{OPL} + \arg(\text{phasor})``: the optical
  path length carries the propagation phase, the unit phasor all other phase (Fresnel and
  coating phases, Gouy phase). See [`RayBundle`](@ref).

## Power normalization

Field amplitudes ``\mathbf{E}`` are physical complex amplitudes in V/m (peak values, not
RMS): the real field is ``\mathrm{Re}(\mathbf{E}\,e^{-i\omega t})``. The power through a
port is

```math
P = \kappa \int |\mathbf{E}|^2 \,\mathrm{d}A, \qquad \kappa = \frac{n}{2 Z_0},
```

with ``n`` the refractive index at the port and
``Z_0 = \mu_0 c = 376.730313668\,\Omega`` (CODATA 2018, the same value as in
BeamletOptics.jl). The same factor is used by every format and converter, see
[`power_normalization`](@ref) and [`VACUUM_IMPEDANCE`](@ref).

The formula is exact for fields that propagate along the port normal and a paraxial
approximation otherwise. ``\mathbf{E}`` carries no obliquity weighting; formats that need
one (e.g. angular spectra) carry it in their own sample weights.

## Polarization

- Three-dimensional ``\mathbf{E}`` vectors are given in the global frame.
- Jones vectors are given in a right-handed transverse basis. At a port this is
  ``(\mathbf{u}, \mathbf{v})`` for fields along ``\mathbf{n}``. For a ray with direction
  ``\mathbf{d}`` it is the ray basis ``(R\mathbf{u}, R\mathbf{v})``, where ``R`` is the
  minimal rotation about ``\mathbf{n}\times\mathbf{d}`` that turns ``\mathbf{n}`` into
  ``\mathbf{d}``. This equals the Richards–Wolf mapping of an aplanatic lens. See
  [`ray_basis`](@ref) and [`jones_to_global`](@ref).
- Circular polarization is named by helicity: right-circular means positive helicity,
  Jones vector ``(1, i)/\sqrt{2}``. With ``\exp(-i\omega t)`` its field rotates from the
  first to the second basis vector, i.e. counter-clockwise for an observer facing the
  source. This is the opposite of the traditional naming in Born & Wolf and Hecht. See
  [`circular_jones`](@ref).
- For every sample of a plane-wave spectrum, transversality ``\mathbf{k}\cdot\mathbf{E} = 0``
  must hold.

## Coherence

Each field is one monochromatic component, and all its parts (rays, samples) add
coherently as complex amplitudes. Polychromatic or mutually incoherent light is a
collection of fields. Components of the same coherence group add coherently (field
amplitudes), components of different groups add incoherently (intensities). Coherence
groups are not part of the API yet; until then, every field is its own group.

## Sampling grids

Sampled fields live on regular Cartesian grids in the local coordinates of their port.
Sample ``i`` along an axis with ``N`` samples and spacing ``\Delta`` lies at
``(i - (N \div 2 + 1))\,\Delta``, so the port origin sits at the center of the grid in the
FFT sense (the `fftshift` center). Grids have no offset of their own: to move a grid,
move the port. See [`RegularGrid`](@ref) and [`coordinates`](@ref).
