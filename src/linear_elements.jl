# 2x2 Gaussian quadrature points and weights
const GAUSS_POINTS_2x2 = ((1.0, -1/sqrt(3)), (1.0, 1/sqrt(3)))

# Shape functions and derivatives in natural coordinates (ξ, η)
function shape_functions_linear(ξ, η)
    N = SVector{4}((1 - ξ) * (1 - η) / 4,
                   (1 + ξ) * (1 - η) / 4,
                   (1 + ξ) * (1 + η) / 4,
                   (1 - ξ) * (1 + η) / 4)
    return N
end

function shape_function_derivatives_linear(ξ, η)
    dN_dξ = SVector{4}(- (1 - η) / 4,  (1 - η) / 4,  (1 + η) / 4, - (1 + η) / 4)
    dN_dη = SVector{4}(- (1 - ξ) / 4, - (1 + ξ) / 4,  (1 + ξ) / 4,  (1 - ξ) / 4)
    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadrilateral element using Gaussian quadrature
function element_stiffness_linear(nodes, element, D, thickness)
    Ke = zeros(SMatrix{8,8})  # Initialize element stiffness matrix (8x8 for a quadrilateral element)
    node_positions_transposed = SMatrix{4,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        # J = SMatrix{2,2}(sum(dN_dξ[i] * nodes[1, element[i]] for i in 1:4),
        #                  sum(dN_dξ[i] * nodes[2, element[i]] for i in 1:4),
        #                  sum(dN_dη[i] * nodes[1, element[i]] for i in 1:4),
        #                  sum(dN_dη[i] * nodes[2, element[i]] for i in 1:4))
        J = SMatrix{2,4}(dN_dξ[1], dN_dη[1],
                         dN_dξ[2], dN_dη[2],
                         dN_dξ[3], dN_dη[3],
                         dN_dξ[4], dN_dη[4]) * node_positions_transposed

        detJ = det(J)  # Determinant of the Jacobian
        invJ = inv(J)  # Inverse of the Jacobian

        # Derivatives of shape functions with respect to global coordinates (x, y)
        dN_dx = invJ[1, 1] * dN_dξ .+ invJ[1, 2] * dN_dη
        dN_dy = invJ[2, 1] * dN_dξ .+ invJ[2, 2] * dN_dη

        # Construct the B matrix (strain-displacement matrix)
        # B = zeros(MMatrix{3,8})
        # for i in 1:4
        #     B[1, 2*i-1] = dN_dx[i]
        #     B[2, 2*i]   = dN_dy[i]
        #     B[3, 2*i-1] = dN_dy[i]
        #     B[3, 2*i]   = dN_dx[i]
        # end
        B = SMatrix{3,8}(dN_dx[1], 0.0, dN_dy[1],
                         0.0, dN_dy[1], dN_dx[1],
                         dN_dx[2], 0.0, dN_dy[2],
                         0.0, dN_dy[2], dN_dx[2],
                         dN_dx[3], 0.0, dN_dy[3],
                         0.0, dN_dy[3], dN_dx[3],
                         dN_dx[4], 0.0, dN_dy[4],
                         0.0, dN_dy[4], dN_dx[4])

        # Weighting factor (Jacobian determinant and Gauss weights)
        weight = detJ * wξ * wη

        # Contribution to the element stiffness matrix
        Ke += thickness * weight * B' * D * B
    end

    return Ke
end


# Function to assemble the global stiffness matrix
function assemble_global_stiffness_linear(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness_linear(nodes, element, D, material.thickness)
        # Global node indices
        dofs = zeros(Int, 8)
        dofs[1] = get_dof(element[1], :x)
        dofs[2] = get_dof(element[1], :y)
        dofs[3] = get_dof(element[2], :x)
        dofs[4] = get_dof(element[2], :y)
        dofs[5] = get_dof(element[3], :x)
        dofs[6] = get_dof(element[3], :y)
        dofs[7] = get_dof(element[4], :x)
        dofs[8] = get_dof(element[4], :y)

        # Add element stiffness to global stiffness matrix
        for i in 1:8, j in 1:8
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
end

# Main solver function
function solve_fem_linear(nodes, elements, material, forces, fixed_dofs)
    # Assemble global stiffness matrix
    K = assemble_global_stiffness_linear(nodes, elements, material)

    # Force vector
    f = get_force_vector(nodes, forces)

    # Apply boundary conditions
    apply_boundary_conditions!(K, f, fixed_dofs)

    # Solve for displacements
    u = K \ f

    return reshape(u, 2, :), K, f
end
