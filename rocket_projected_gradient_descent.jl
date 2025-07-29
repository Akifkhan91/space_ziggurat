using Plots, LaTeXStrings, Test, LinearAlgebra, ReverseDiff
using CSV, DataFrames
include("../tests/run_tests.jl")
include("../tests/test_escape_velocity.jl")
include("../tests/test_circular_motion.jl")
include("../tests/test_freefall.jl")
include("../tests/test_LEO.jl")
include("../tests/test_GEO.jl")
include("../tests/test_sim_tangential_thrust_LEO.jl")
include("../tests/test_hohmann.jl")
include("../src/simulator.jl")
include("plots.jl")

# All mass units in Tonne (1 Tonne = 10^3 kg)
# All distance units in Kilo-meters (1 km = 10^3 m)
# All time units in seconds

struct units_scaling
    mass_scaling::Float64
    distance_scaling::Float64
    time_scaling::Float64
end

# Default scaling is 1 Tonnes, 1 Km, 1 seconds 
# for mass, distance and time
# Other options 1, 10^-3, 60^-1
# for Tonne, Mega-meters and minutes

# Using Tonne, Mega-meters, minutes
# t_scaling = LinRange(1e-3,1,100)
# m_scaling = LinRange(1e-3,1,100)
# d_scaling = LinRange(1e-3,1,100)

# test_pass_scaling_vector = zeros(length(m_scaling))
# for m in eachindex(m_scaling)
# scaling = units_scaling(m,1,1)

# Using Tonne, Kilo-meters, seconds
scaling = units_scaling(1,1,1)

mutable struct parameters
    T::Int64
    scaled_mu_const::Float64
    initial_mass::Float64
    final_mass::Float64
    mass_divided_by_thrust::Float64
    radius_of_the_earth::Float64
    radial_tolerance::Float64
    min_total_time::Float64
    max_total_time::Float64
    max_thrust::Float64
    ord::Int64
    coord_system::String
end

data = parameters(100, 3.983*1e5*scaling.distance_scaling^3/scaling.time_scaling^2, 
        421.3*scaling.mass_scaling, 25.6*scaling.mass_scaling, 
        0.1103*scaling.time_scaling/scaling.distance_scaling, 
        6378*scaling.distance_scaling, 0.0*scaling.distance_scaling,
        1*scaling.time_scaling, 1e5*scaling.time_scaling, 
        8.227*scaling.mass_scaling*scaling.distance_scaling/scaling.time_scaling^2,
        1, "cartesian")


# Initial positon
x_start = 6378.1*scaling.distance_scaling
y_start = 0.0*scaling.distance_scaling
z_start = 0.0*scaling.distance_scaling

# Final position
x_end = 0.0*scaling.distance_scaling
y_end = 6800.0*scaling.distance_scaling
z_end = 0.0*scaling.distance_scaling

#escape_velocity = sqrt(2 * data.scaled_mu_const/x_start)
# Orbital speed at radius r = √(GM_E/r)
v_orbital = sqrt(data.scaled_mu_const/sqrt(x_end^2+y_end^2+z_end^2))

# Initial velocity
v_x_start = 0.0*scaling.distance_scaling/scaling.time_scaling
v_y_start = 0.0*scaling.distance_scaling/scaling.time_scaling
v_z_start = 0.0*scaling.distance_scaling/scaling.time_scaling

# Final velocity of the rocket
v_x_end = -v_orbital
v_y_end = 0.0*scaling.distance_scaling/scaling.time_scaling
v_z_end = 0.0*scaling.distance_scaling/scaling.time_scaling

mutable struct InitialandBoundaryConditions
    initial_position::Vector{Float64}
    final_position::Vector{Float64}
    initial_velocity::Vector{Float64}
    final_velocity::Vector{Float64}
end

initial_and_boundary_conditions = InitialandBoundaryConditions(
    zeros(3),
    zeros(3),
    zeros(3),
    zeros(3)
    )

initial_and_boundary_conditions.initial_position .= [x_start, y_start, z_start] 
initial_and_boundary_conditions.final_position .= [x_end,y_end,z_end]
initial_and_boundary_conditions.initial_velocity .= [v_x_start, v_y_start, v_z_start]
initial_and_boundary_conditions.final_velocity .= [v_x_end, v_y_end, v_z_end]

starting_thrust_x =  data.max_thrust/sqrt(3)
starting_thrust_y =  data.max_thrust/sqrt(3)
starting_thrust_z =  data.max_thrust/sqrt(3)

starting_thrust_vector = [starting_thrust_x, starting_thrust_y, starting_thrust_z]

#run_tests(initial_and_boundary_conditions,data,scaling)
