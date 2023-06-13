"""
    setup_geometry(
        dtype::DataType=Float64;
        rsource::Real,
        θsource::Real,
        ϕsource::Real,
        robs::Real,
        θobs::Real,
        ϕobs::Real,
        a::Real,
        s::Integer=2,
        getmagnification::Bool=false,
        direction_coords::Symbol=:spherical,
    )

Setup the geometry.
"""
function setup_geometry(
    dtype::DataType=Float64;
    rsource::Real,
    θsource::Real,
    ϕsource::Real,
    robs::Real,
    θobs::Real,
    ϕobs::Real,
    a::Real,
    s::Integer=2,
    getmagnification::Bool=false,
    direction_coords::Symbol=:spherical,
)
    coords_choices = [:spherical, :shadow, :shadowpos]
    @assert direction_coords in coords_choices "`direction_coords` must be one of `$coords_choices`"

    source = SphericalCoords{dtype}(r=rsource, θ=θsource, ϕ=ϕsource)
    observer = SphericalCoords{dtype}(r=robs, θ=θobs, ϕ=ϕobs)
    return Geometry{dtype}(dtype=dtype, source=source, observer=observer, s=s, a=a,
                           direction_coords=direction_coords, getmagnification=getmagnification)
end
