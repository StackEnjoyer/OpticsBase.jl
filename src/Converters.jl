"""
    convert_field(conv::AbstractFieldConverter, data::AbstractOpticalData) -> AbstractOpticalData

Converts `data` into another exchange format with the converter `conv`, e.g. a
[`PlaneWaveSpectrum`](@ref) into a [`SampledField`](@ref) with [`PlaneWaveSummation`](@ref).

Steps: [`check_compatibility`](@ref)`(data, conv)` (throws a
[`MissingConverterError`](@ref) before any converter code runs), then
[`__convert_field`](@ref)`(conv, data)`, then a check that the result is an instance of
`output_representation(conv)` (`ArgumentError` otherwise).

A conversion is an approximation; what is assumed and what is conserved is documented by
each converter. Converter packages must not add methods to `convert_field`; they implement
`__convert_field`.
"""
function convert_field(conv::AbstractFieldConverter, field::AbstractOpticalData)
    check_compatibility(field, conv)
    out = __convert_field(conv, field)
    out isa output_representation(conv) ||
        throw(ArgumentError("convert_field: the converter returned a $(typeof(out)), " *
                            "but output_representation(conv) is $(output_representation(conv))"))
    return out
end

"""
    __convert_field(conv::AbstractFieldConverter, data::AbstractOpticalData) -> AbstractOpticalData

Converter-side entry point of [`convert_field`](@ref), implemented by every
[`AbstractFieldConverter`](@ref).

# Contract

  - It is called after [`check_compatibility`](@ref), so
    `data isa input_representation(conv)` holds and need not be checked again.
  - It returns new data that is an instance of `output_representation(conv)`, follows the
    conventions (global frame, SI units, power normalization of its port) and does not
    share mutable state with `data` unless documented.

Users call `convert_field`, not `__convert_field`.
"""
function __convert_field end
