module GSHE_benchmark

# kerr_functions.jl
export initial_spatial_comomentum, tetrad_boosting, init_values
# kerr_trajectories.jl
export geodesic_odes!, gshe_odes!
# objects.jl
export SphericalCoords, Geometry
# setup.jl
export setup_geometry


import Parameters: @with_kw, @unpack

include("./objects.jl")
include("./kerr_functions.jl")
include("./kerr_trajectories.jl")
include("./setup.jl")

end
