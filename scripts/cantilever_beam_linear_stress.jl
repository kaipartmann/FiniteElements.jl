using FiniteElements

include("cantilever_beam_meshes.jl")
include("plotting.jl")

## FEM Simulation
beam_length = 1000.0 # Length of the beam
beam_height = 100.0 # Height of the beam
nx = 20 # Number of elements in x-direction
ny = 2 # Number of elements in y-direction
applied_force = 1e4 # Total force applied along the right side in -y direction
nodes, elements, forces, fixed_dofs = create_cantilever_beam_linear(beam_length,
                                                                    beam_height, nx, ny,
                                                                    applied_force)
material = Material(1e5, 0.3, 6, 7850)
u, K, f, σ = solve_fem_linear_stress(nodes, elements, material, forces, fixed_dofs)
println("-- results linear elements --")
println("number of elements: $(size(elements, 2))")
println("minimum displacement: $(minimum(u[2, :]))")

## Post-Processing
fig = Figure(size=(800, 300))
ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="cantilever beam linear")
deformed_nodes = nodes .+ 2 .* u
num_nodes = size(u,2)
node_colors = zeros(Float64, num_nodes)
counts = zeros(Int, num_nodes)

for (i, element) in enumerate(eachcol(elements))
    for n in elements[:, i]
        node_colors[n] += σ[i,1]
    end
end

for n in 1:num_nodes
    if counts[n] > 0
        node_colors[n] /= counts[n]
    end
end

plot_mesh!(ax1, deformed_nodes, elements; color=node_colors)
plot_edges_linear!(ax1, deformed_nodes, elements; color=:black)
plot_nodes!(ax1, deformed_nodes; color=:black, markersize=7)
save(joinpath("img", "cantilever_beam_linear.png"), fig; px_per_unit=2)
fig
