# Main solver function
function solve_fem_linear_stress(nodes, elements, material, forces, fixed_dofs)
    # Assemble global stiffness matrix
    K = assemble_global_stiffness_linear(nodes, elements, material)

    # Force vector
    f = get_force_vector(nodes, forces)

    # Apply boundary conditions
    apply_boundary_conditions!(K, f, fixed_dofs)

    # Solve for displacements
    u = K \ f

    n_nodes = size(nodes, 2)
    σ = zeros(size(elements, 2), 3)
    D = plane_stress_stiffness(material)

    for (i,element) in enumerate(eachcol(elements))
        node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')
        strain_element = zeros(3)  # Accumulated strain

        for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
            # Get shape function derivatives in natural coordinates
            dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

            # Calculate the Jacobian matrix and its determinant
            J = SMatrix{2,4}(dN_dξ[1], dN_dη[1],
                             dN_dξ[2], dN_dη[2],
                             dN_dξ[3], dN_dη[3],
                             dN_dξ[4], dN_dη[4]) * node_positions_transposed
            invJ = inv(J)  # Inverse of the Jacobian

            # Derivatives of shape functions with respect to global coordinates (x, y)
            dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
            dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

            # Construct the B matrix (strain-displacement matrix)
            B = SMatrix{3,8}(dN_dx[1], 0.0, dN_dy[1],
                             0.0, dN_dy[1], dN_dx[1],
                             dN_dx[2], 0.0, dN_dy[2],
                             0.0, dN_dy[2], dN_dx[2],
                             dN_dx[3], 0.0, dN_dy[3],
                             0.0, dN_dy[3], dN_dx[3],
                             dN_dx[4], 0.0, dN_dy[4],
                             0.0, dN_dy[4], dN_dx[4])
            strain_element += B*u[get_dofs_linear(element)]
        end
        # Average strain over Gauss points
        strain_element /= 4
        # Calculate stress: σ = D * ε
        σ[i, :] = D * strain_element
        println("Die Spannungen für Element $i lauten", σ[i,:])
    end
    return reshape(u, 2, :), K, f, σ
end

# Main solver function for dynamic FEM analysis using Newmark-Beta method
function solve_fem_linear_stress_dynamic(nodes, elements, material, forces, fixed_dofs, n_timesteps, Δt)
    # Newmark-Beta parameters
    # β = 1/4, γ = 1/2: Unconditionally stable (average acceleration method)
    # β = 1/6, γ = 1/2: Linear acceleration method
    # β = 0,   γ = 1/2: Central difference method (conditionally stable)
    β = 1/4
    γ = 1/2

    # Assembly of global matrices
    K = assemble_global_stiffness_linear(nodes, elements, material)
    M = assemble_global_mass_linear(nodes, elements, material)

    n_dofs = 2 * size(nodes, 2)

    # Initialize displacement, velocity, and acceleration vectors
    u = zeros(n_dofs)
    v = zeros(n_dofs)
    a = zeros(n_dofs)

    # Force vector (assuming constant force for simplicity)
    f = get_force_vector(nodes, forces)

    # Apply boundary conditions to stiffness and mass matrices
    K_bc = copy(K)
    M_bc = copy(M)
    f_bc = copy(f)
    apply_boundary_conditions!(K_bc, f_bc, fixed_dofs)
    apply_boundary_conditions_mass!(M_bc, fixed_dofs)

    # Calculate initial acceleration: M*a₀ = f₀ - K*u₀
    # Since u₀ = 0 and v₀ = 0, we have: M*a₀ = f₀
    a = M_bc \ f_bc

    # Newmark-Beta constants
    a₀ = 1.0 / (β * Δt^2)
    a₁ = γ / (β * Δt)
    a₂ = 1.0 / (β * Δt)
    a₃ = 1.0 / (2.0 * β) - 1.0
    a₄ = γ / β - 1.0
    a₅ = Δt / 2.0 * (γ / β - 2.0)
    a₆ = Δt * (1.0 - γ)
    a₇ = γ * Δt

    # Effective stiffness matrix: K_eff = K + a₀*M
    K_eff = K_bc + a₀ * M_bc

    # Time stepping loop
    for t = 1:n_timesteps
        # Store old values
        a_old = copy(a)

        # Predictor phase: predict displacement and velocity
        u_pred = u + Δt * v + (Δt^2 / 2.0) * (1.0 - 2.0 * β) * a
        v_pred = v + Δt * (1.0 - γ) * a

        # Force vector at time t+1 (assuming constant forces for this implementation)
        f_t1 = f_bc

        # Effective force: f_eff = f_{t+1} + M*(a₀*u_pred + a₂*v + a₃*a)
        f_eff = f_t1 + M_bc * (a₀ * u_pred + a₂ * v + a₃ * a)

        # Solve for displacement at t+1: K_eff * u_{t+1} = f_eff
        u = K_eff \ f_eff

        # Corrector phase: calculate acceleration and velocity at t+1
        a = a₀ * (u - u_pred) - a₂ * v - a₃ * a_old
        v = v_pred + a₇ * a

        # Calculate stresses for current time step
        if t % max(1, n_timesteps ÷ 10) == 0
            σ = calculate_element_stresses(nodes, elements, material, u)
            u_nodes = reshape(u, 2, :)
            create_visualization(nodes, elements, u_nodes, σ, t)
            println("Time step $t/$n_timesteps completed")
        end
    end

    # Calculate final stress state
    σ = calculate_element_stresses(nodes, elements, material, u)

    return reshape(u, 2, :), K, f, σ
end

# Helper function to calculate element stresses
function calculate_element_stresses(nodes, elements, material, u)
    n_elements = size(elements, 2)
    σ = zeros(n_elements, 3)  # [σxx, σyy, σxy] for each element
    D = plane_stress_stiffness(material)

    for (i, element) in enumerate(eachcol(elements))
        node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')
        strain_element = zeros(3)  # Accumulated strain

        # Integration over Gauss points
        for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
            # Get shape function derivatives in natural coordinates
            dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

            # Calculate the Jacobian matrix and its determinant
            J = SMatrix{2,4}(dN_dξ[1], dN_dη[1],
                             dN_dξ[2], dN_dη[2],
                             dN_dξ[3], dN_dη[3],
                             dN_dξ[4], dN_dη[4]) * node_positions_transposed
            invJ = inv(J)  # Inverse of the Jacobian

            # Derivatives of shape functions with respect to global coordinates (x, y)
            dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
            dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

            # Construct the B matrix (strain-displacement matrix)
            B = SMatrix{3,8}(dN_dx[1], 0.0, dN_dy[1],
                             0.0, dN_dy[1], dN_dx[1],
                             dN_dx[2], 0.0, dN_dy[2],
                             0.0, dN_dy[2], dN_dx[2],
                             dN_dx[3], 0.0, dN_dy[3],
                             0.0, dN_dy[3], dN_dx[3],
                             dN_dx[4], 0.0, dN_dy[4],
                             0.0, dN_dy[4], dN_dx[4])

            # Calculate strain at this Gauss point
            strain_element += B * u[get_dofs_linear(element)]
        end

        # Average strain over Gauss points
        strain_element /= 4

        # Calculate stress: σ = D * ε
        σ[i, :] = D * strain_element
    end

    return σ
end

# Helper function to create visualization
function create_visualization(nodes, elements, u_nodes, σ, timestep)
    try
        fig = Figure(size=(800, 300))
        ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="Cantilever beam linear - Time step $timestep")
        deformed_nodes = nodes .+ 2 .* u_nodes
        num_nodes = size(u_nodes, 2)
        node_colors = zeros(Float64, num_nodes)
        counts = zeros(Int, num_nodes)

        # Average stresses to nodes
        for (i, element) in enumerate(eachcol(elements))
            for n in element
                node_colors[n] += σ[i, 1]  # Use σxx stress component
                counts[n] += 1
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
        save(joinpath("img", "cantilever_beam_linear_$timestep.png"), fig; px_per_unit=2)
    catch e
        @warn "Visualization failed at timestep $timestep: $e"
    end
end
