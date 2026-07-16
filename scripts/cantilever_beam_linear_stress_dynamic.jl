using FiniteElements

include("cantilever_beam_meshes.jl")
include("plotting.jl")

## FEM Simulation
n_timesteps = 1000
Δt = 1e-8
beam_length = 1000.0 # Length of the beam
beam_height = 100.0 # Height of the beam
nx = 20 # Number of elements in x-direction
ny = 2 # Number of elements in y-direction
applied_force = 1e4 # Total force applied along the right side in -y direction
nodes, elements, forces, fixed_dofs = create_cantilever_beam_linear(beam_length,
                                                                    beam_height, nx, ny,
                                                                    applied_force)
material = Material(1e5, 0.3, 6, 7850)
u, K, f, σ = solve_fem_linear_stress_dynamic(nodes, elements, material, forces, fixed_dofs, n_timesteps, Δt)
println("-- results linear elements --")
println("number of elements: $(size(elements, 2))")
println("minimum displacement: $(minimum(u[2, :]))")
println("maximum displacement: $(maximum(u[2, :]))")
println("tip displacement (y): $(u[2, end])")  # Show tip displacement
println("Dynamic analysis completed successfully!")
