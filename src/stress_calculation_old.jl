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
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    σ = zeros(size(elements, 2), 3)
    for (i,element) in enumerate(eachcol(elements))
        node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')  
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
            σ[i, :] += B*u[get_dofs_linear(element)]
        end
        σ[i,:] /= 1/4
        println("Die Spannungen für Element $i lauten", σ[i,:])
    end
    return reshape(u, 2, :), K, f, σ
end

# Main solver function
function solve_fem_linear_stress_dynamic(nodes, elements, material, forces, fixed_dofs, n_timesteps, Δt)
    # Newmark-Beta-Methode
    # Spezialfälle: β = 1/4, γ = 1/2, implizit und bedingungslos stabil.
    #               β = 1/6, γ = 1/2, lineare Beschleunigungsmethode
    #               β =   0, γ = 1/2, zentrale Differenzenmethode
    # Nur bedingt stabil für γ < 1/2
    # Für γ = 1/2: second-order accurate (Fehler ~ \mathcal{O}(Δt²))
    # Für γ ≠ 1/2:  first-order accurate (Fehler ~ \mathcal{O}(Δt))
    # β ist ein Dämpfungsfaktor!
    f = []
    K = []
    σ = []
    β = 1/4
    γ = 1/2

    ndofs = size(nodes,1)*size(nodes,2)
    volume = maximum(nodes[1,:])*maximum(nodes[2,:])
    mass = material.density * volume

    u = zeros(size(nodes, 2)*2) # 2 ist die Dimension des Projekts
    v = zeros(size(nodes, 2)*2)
    a = zeros(size(nodes, 2)*2)

    a_old = a

    for t = 1:n_timesteps
        #Calculating predictors
        if β == 0.0 && γ == 0.5
            u_pred = u + Δt*v + Δt^2/2*a
            v_pred = v + Δt*(1-γ)*a
        else
            u_pred = u + Δt*v + (Δt^2)/2*(1-2*β)*a
            v_pred = v + Δt*(1-γ)*a
        end
        #Predictors finished.

        # Assemble global stiffness matrix
        K = assemble_global_stiffness_linear(nodes, elements, material)

        # Force vector
        f = get_force_vector(nodes, forces)

        # Apply boundary conditions
        apply_boundary_conditions!(K, f, fixed_dofs)

        # Bewegungsgleichung lösen: # Mü + Cu̇ + Ku = f(t), C = 0
                                    # Mü      + Ku = f(t)
                                    # Mü           = f(t) - Ku
                                    #  ü            = 1/M*(f(t) - Ku)
        #@autoinfiltrate
        a = 2*ndofs/mass*(f-K*u)

        #Calculate u_n+1 and v_n+1 from a_n+1
        if β == 0.0 && γ ≈ 0.5
            u = u_pred
            v = v_pred + Δt/2*(a_old+a)
        else
            u = u_pred + β*Δt^2*a
            v = v_pred + γ*Δt*a
        end

        n_nodes = size(nodes, 2)
        K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
        σ = zeros(size(elements, 2), 3)
        for (i,element) in enumerate(eachcol(elements)) # Elementweise Spannungsberechnung
            node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')  
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
                σ[i, :] += B*u[get_dofs_linear(element)]
            end
            σ[i,:] /= 4
            D = plane_stress_stiffness(material)
            σ[i,:] =  D * σ[i,:]
            println("Die Spannungen für Element $i lauten", σ[i,:])
        end

        u_nodes = reshape(u, 2, :)

        ## Post-Processing
        fig = Figure(size=(800, 300))
        ax1 = Axis(fig[1, 1]; aspect=DataAspect(), title="cantilever beam linear")
        deformed_nodes = nodes .+ 2 .* u_nodes
        num_nodes = size(u_nodes,2)
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
        save(joinpath("img", "cantilever_beam_linear_$t.png"), fig; px_per_unit=2)
    end

    #=n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    σ = zeros(size(elements, 2), 3)
    for (i,element) in enumerate(eachcol(elements))
        node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')  
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
            σ[i, :] += B*u[get_dofs_linear(element)]
        end
        σ[i,:] /= 4
        D = plane_stress_stiffness(material)
        σ[i,:] =  D * σ[i,:]
        println("Die Spannungen für Element $i lauten", σ[i,:])
    end=#
    return reshape(u, 2, :), K, f, σ
end