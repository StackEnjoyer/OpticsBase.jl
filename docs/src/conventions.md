```@meta
CurrentModule = OpticsBase
```

# Conventions

These conventions are binding for every package that produces or consumes a
[`PlaneField`](@ref).

!!! warning
    Mixing conventions silently corrupts phases, powers or polarization states across a
    solver chain. When in doubt, check this page first.

## What is exchanged

A [`PlaneField`](@ref) holds the tangential electric and magnetic field,
``(E_u, E_v, H_u, H_v)``, sampled on a plane. By the surface equivalence theorem these
four components determine the field on both sides of the plane, for any medium and both
directions of travel. The normal components are not stored; they follow from Maxwell's
curl equations on the plane,

```math
E_n = \frac{i Z_0}{k_0 n^2}\left(\partial_u H_v - \partial_v H_u\right), \qquad
H_n = -\frac{i}{k_0 Z_0}\left(\partial_u E_v - \partial_v E_u\right).
```

Solvers that only know ``\mathbf{E}`` use the E-only constructor. It fills in
``\mathbf{H}`` for a local plane wave along the *reference direction* of each sample:
``\mathbf{n}`` without a reference sphere, otherwise the ray of the reference sphere
through the sample. For ``R = \infty`` this is
``\mathbf{H} = (n/Z_0)\,\mathbf{n}\times\mathbf{E}``. [`forward`](@ref) and
[`backward`](@ref) split a field into local plane waves along and against the same
reference directions. Both are exact for waves along the reference direction; a plane wave
at the angle ``\theta`` to it is off by the order of ``1 - \cos\theta``. `forward` is
therefore meant for fields that may contain both directions: a solver that knows its
input travels one way takes `f.E` directly.

There is no scalar field type. A scalar solver hands over ``E_u = \psi``, ``E_v = 0`` and
reads ``E_u``.

## Medium

A plane lies in a homogeneous, isotropic, lossless and non-magnetic medium of refractive
index `n`, at least across the sampled region. The tangential fields are exact in any
medium, but `n` enters the E-only constructor, the forward/backward split, the reference
phase and the normal components. A plane through an inhomogeneous structure, such as a
fiber cross-section, has no single `n`: place it in the surrounding medium.

## Units

SI throughout: meters, seconds, watts, V/m, A/m. Wavelengths are vacuum wavelengths.
``Z_0 = \mu_0 c = 376.730313668\,\Omega`` ([`VACUUM_IMPEDANCE`](@ref), the same literal as
BeamletOptics.jl).

## Time convention

Time-harmonic fields use ``\exp(-i\omega t)``. A plane wave is

```math
\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 \exp\big(i(\mathbf{k}\cdot\mathbf{r} - \omega t)\big),
```

so propagation over a distance ``d`` in a medium with index ``n`` adds the phase
``+n k_0 d`` with ``k_0 = 2\pi/\lambda``.

## Amplitudes and power

``\mathbf{E}`` and ``\mathbf{H}`` are complex peak amplitudes: the real field is
``\mathrm{Re}(\mathbf{E}\,e^{-i\omega t})``. The power through the plane is the
time-averaged Poynting flux along ``\mathbf{n}``,

```math
P = \tfrac{1}{2}\,\mathrm{Re}\int (\mathbf{E}\times\mathbf{H}^*)\cdot\mathbf{n}\,\mathrm{d}A
  = \tfrac{1}{2}\,\mathrm{Re}\int (E_u H_v^* - E_v H_u^*)\,\mathrm{d}A,
```

see [`power`](@ref). It is exact and signed; no paraxial factor is involved.

## Frame

- The plane has an `origin` and `axes` ``(\mathbf{u}, \mathbf{v}, \mathbf{n})`` in global
  coordinates, orthonormal (to within ``10^{-6}``) and right-handed,
  ``\mathbf{u}\times\mathbf{v} = \mathbf{n}``. There is no global optical axis.
- ``\mathbf{n}`` is the positive normal: [`power`](@ref) counts the flux along
  ``+\mathbf{n}`` as positive, and [`forward`](@ref) is the part travelling towards the
  ``+\mathbf{n}`` side. The *reference direction* of a sample is a different thing: it
  equals ``\mathbf{n}`` only without a reference sphere, see [What is exchanged](@ref).
- ``\mathbf{u}`` fixes the polarization basis and is always given explicitly.
- Field components are stored in the local frame. Moving or rotating a plane changes
  `origin` and `axes` only, never the arrays.

## Sampling

Sample ``i`` along an axis with ``N`` samples and spacing ``\Delta`` sits at
``(i - 1 - \lfloor N/2 \rfloor)\,\Delta`` from the origin (the fftshift center), see
[`coordinates`](@ref). There are no grid offsets: move the origin instead.

The stored arrays must resolve their phase: the phase difference between neighboring
samples stays below ``\pi``. A strongly curved wavefront is made cheap to sample with the
reference sphere.

## Phase

- The physical fields are `f.E .* reference_phase(f)` and `f.H .* reference_phase(f)`
  (see [`reference_phase`](@ref)). With `R = Inf` the stored arrays are the physical
  fields.
- Physical fields carry the full spatial phase; only the carrier ``\exp(-i\omega t)`` is
  excluded. No tilt or defocus is removed, apart from the reference sphere.
- Optical path lengths are absolute: they are counted from the reference point of the
  chain, which the first solver sets (typically the source). Every later solver continues
  the count, so fields that reach one plane on different paths keep their relative phase.

## Polarization

- Jones vectors are given in the basis ``(\mathbf{u}, \mathbf{v})``.
- Circular polarization is named by helicity: right-circular means positive helicity,
  Jones vector ``(1, i)/\sqrt{2}``. With ``\exp(-i\omega t)`` the field rotates from
  ``\mathbf{u}`` to ``\mathbf{v}``, counter-clockwise for an observer facing the source.
  This is the opposite of the traditional naming in Born & Wolf and Hecht.

## More than one plane

A `Vector{PlaneField}` has no meaning of its own: the function that returns or takes it
states in its docstring which of the following it is. The rules for each:

- **Closed Huygens surface** (e.g. the faces of a box around a scatterer): one coherent,
  monochromatic field; all faces share `λ`, the phase reference and the time origin, and
  every face's `n` points outward, so the summed [`power`](@ref) is the net power leaving
  the enclosed volume.
- **Spectrum** (discrete lines, e.g. a sampled pulse): one field per wavelength, all with
  the common time origin ``t = 0``. Each field holds the peak amplitude of its line, so
  its `power` is the power of that line. The total power is the sum, since different
  frequencies do not interfere in the time average.
- **Mutually incoherent components** (e.g. two polarizations of unpolarized light, or
  independent sources): powers add, phases between the components have no meaning.
