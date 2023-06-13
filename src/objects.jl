@with_kw mutable struct SphericalCoords{T <: Real} <: Number
    t::T = 0.0
    r::T
    θ::T
    ϕ::T
end

@with_kw mutable struct Geometry{T <: Real}
    dtype::DataType
    source::SphericalCoords{T}
    observer::SphericalCoords{T}
    direction_coords::Symbol=:spherical
    getmagnification::Bool=false
    s::Integer = 2
    a::T
end
