# Conventions

The following conventions are binding for every exchange format, converter and solver
that connects to `OpticsBase`. If a format deviates from one of them, its docstring must
say so explicitly.

!!! warning
    Mixing conventions silently corrupts phases, powers or polarization states across a
    solver chain. When in doubt, check this page first.

## Units

All quantities are in SI units internally: meters, seconds, watts. Wavelengths are
vacuum wavelengths in meters. Support for
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities is optional and
only provided via a package extension; it converts at the API boundary.

## Time convention

Time-harmonic fields use the ``\exp(-i\omega t)`` convention. A plane wave is written as

```math
\mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0 \exp\big(i(\mathbf{k}\cdot\mathbf{r} - \omega t)\big).
```

Consequently, propagation over a distance ``d`` in a medium with refractive index ``n``
adds the phase ``+ n k_0 d``.

## Ports and frames

Every handover of a field between solvers happens at a port: a surface in global
coordinates with an origin, a normal and local transverse axes. There is no field
handover without a port. `OpticsBase` defines no global optical axis; solvers with a
preferred axis map it onto the port frame in their extension.

The port also defines the polarization basis (see below).

## Phase reference

The carrier wave is not contained in the stored complex amplitudes. Optical path lengths
are referenced to the port origin.

!!! warning "Not yet fixed"
    Whether "carrier" means only the temporal factor ``\exp(-i\omega t)`` or also a
    spatial carrier along the port normal has not been decided yet.

## Power normalization

Field amplitudes are normalized such that the power through the port is

```math
P = \kappa \int |\mathbf{E}|^2 \,\mathrm{d}A,
```

with one documented factor ``\kappa`` that accounts for refractive index and impedance.
The same factor is used by every format and converter.

!!! warning "Not yet fixed"
    The explicit form of ``\kappa`` has not been decided yet.

## Polarization

- Jones vectors are given in the local basis of the port.
- Three-dimensional ``\mathbf{E}`` vectors are given in the global frame.
- For every sample of an angular spectrum, transversality ``\mathbf{k}\cdot\mathbf{E} = 0``
  must hold.

!!! warning "Not yet fixed"
    The handedness convention for circular polarization (viewed along or against the
    propagation direction) has not been decided yet.

## Coherence

Polychromatic fields are lists of monochromatic components. Each component belongs to a
coherence group. Components within a group add coherently (field amplitudes),
components of different groups add incoherently (intensities).
