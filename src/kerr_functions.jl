"""
    initial_spatial_comomentum(k::Vector{<:Real}, geometry::Geometry)

Calculate the initial covectors [p_1, p_2, p_3] for a given initial direction and geometry.
"""
function initial_spatial_comomentum(k::Vector{<:Real}, geometry::Geometry)
    @unpack r, θ = geometry.source
    a = geometry.a
    v3 = tetrad_boosting(r, geometry)
    s_θ, c_θ = sin(θ), cos(θ)


    if geometry.direction_coords == :spherical
        ψ, ρ = k
        return [sqrt((a^2*c_θ^2 + r^2)/(a^2 + r*(r - 2)))*sin(ψ)*cos(ρ), sqrt(a^2*c_θ^2 + r^2)*sin(ρ)*sin(ψ), s_θ*(a*s_θ*v3*sqrt(a^2 + r*(r - 2))*cos(ψ) + a*s_θ*sqrt(a^2 + r*(r - 2)) + (a^2 + r^2)*(v3 + cos(ψ)))/sqrt(-(v3^2 - 1)*(a^2*c_θ^2 + r^2))]
    elseif geometry.direction_coords == :shadow
        k2, k3 = k
        return [-sqrt(-(a^2*c_θ^2 + r^2)*(k2^2 + k3^2 - 1)/(a^2 + r*(r - 2))), k2*sqrt(a^2*c_θ^2 + r^2), s_θ*(a*k3*s_θ*v3*sqrt(a^2 + r*(r - 2)) + a*s_θ*sqrt(a^2 + r*(r - 2)) + (a^2 + r^2)*(k3 + v3))/sqrt(-(v3^2 - 1)*(a^2*c_θ^2 + r^2))]
    elseif geometry.direction_coords == :shadowpos
        k2, k3 = k
        return [sqrt(-(a^2*c_θ^2 + r^2)*(k2^2 + k3^2 - 1)/(a^2 + r*(r - 2))), k2*sqrt(a^2*c_θ^2 + r^2), s_θ*(a*k3*s_θ*v3*sqrt(a^2 + r*(r - 2)) + a*s_θ*sqrt(a^2 + r*(r - 2)) + (a^2 + r^2)*(k3 + v3))/sqrt(-(v3^2 - 1)*(a^2*c_θ^2 + r^2))]
    else
        return NaN
    end
end


"""
    tetrad_boosting(r::Real, geometry::Geometry)

Boosting source and observer function.
"""
function tetrad_boosting(r::Real, geometry::Geometry)
    Robs, θobs = geometry.observer.r, geometry.observer.θ
    Rsrc, θsrc = geometry.source.r, geometry.source.θ
    a = geometry.a
    return -a*exp(-(-Rsrc + r)^2)*sin(θsrc)/sqrt(Rsrc^2 - 2*Rsrc + a^2) - a*exp(-(-Robs + r)^2)*sin(θobs)/sqrt(Robs^2 - 2*Robs + a^2)
end


"""
    derivative_tetrad_boosting(r::Real, geometry::Geometry)

Derivative with respect to radius of the boosting source and observer function.
"""
function derivative_tetrad_boosting(r::Real, geometry::Geometry)
    Robs, θobs = geometry.observer.r, geometry.observer.θ
    Rsrc, θsrc = geometry.source.r, geometry.source.θ
    a = geometry.a
    return 2*a*((-Robs + r)*exp(-(-Robs + r)^2)*sin(θobs)/sqrt(Robs*(Robs - 2) + a^2) + (-Rsrc + r)*exp(-(-Rsrc + r)^2)*sin(θsrc)/sqrt(Rsrc*(Rsrc - 2) + a^2))
end


function init_values(init_direction::Vector{<:Real}, geometry::Geometry)
    @unpack t, r, θ, ϕ = geometry.source
    return [[t,r, θ,ϕ]; initial_spatial_comomentum(init_direction, geometry)]
end
