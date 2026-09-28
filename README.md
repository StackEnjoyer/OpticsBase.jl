<p align="center">
  <img src="docs/src/assets/logo.svg" width="300" alt="OpticsBase.jl logo">
</p>

# OpticsBase.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://StackEnjoyer.github.io/OpticsBase.jl/dev/)
[![CI](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/StackEnjoyer/OpticsBase.jl/actions/workflows/CI.yml?query=branch%3Amain)

A lightweight base package for an ecosystem of interoperable optical solvers in Julia.
Solvers do not know each other, only `OpticsBase`: optical propagation data is handed
from one solver to the next through well-defined exchange formats (ray bundles, sampled
fields, angular spectra) at ports with explicit frames and conventions.

`OpticsBase` contains no solver code. It provides abstract types, exchange formats,
traits, ports, conventions and generic converters, and uses the
[CommonSolve.jl](https://github.com/SciML/CommonSolve.jl) interface for solvers.

**Status:** early development, the API is not stable yet.
