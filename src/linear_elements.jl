# 2x2 Gaussian quadrature points and weights
const GAUSS_POINTS_2x2 = ((1.0, -1/sqrt(3)), (1.0, 1/sqrt(3)))

# Shape functions for linear quadrilateral elements in natural coordinates (ξ, η):
#    N1 = (1 - ξ) * (1 - η) / 4
#    N2 = (1 + ξ) * (1 - η) / 4
#    N3 = (1 + ξ) * (1 + η) / 4
#    N4 = (1 - ξ) * (1 + η) / 4
#
#    ◍---------◍
#   4|        3|
#    |         |
#    |         |
#    ◍---------◍
#    1         2
function shape_function_derivatives_linear(ξ, η)
    dN_dξ = SVector{4}(- (1 - η) / 4,  (1 - η) / 4,  (1 + η) / 4, - (1 + η) / 4)
    dN_dη = SVector{4}(- (1 - ξ) / 4, - (1 + ξ) / 4,  (1 + ξ) / 4,  (1 - ξ) / 4)
    return dN_dξ, dN_dη
end

# Function to calculate the element stiffness matrix for a quadrilateral element using
# Gaussian quadrature
function element_stiffness_linear(nodes, element, D, thickness)
    # Initialize element stiffness matrix (8x8 for a quadrilateral element)
    Ke = zeros(SMatrix{8,8})
    node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

        # Calculate the Jacobian matrix and its determinant
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

# Function to calculate the element stiffness matrix for a quadrilateral element using
# Gaussian quadrature, with u for stress calculation in iterative methods
function element_stiffness_linear(nodes, element, D, thickness, u)
    # Initialize element stiffness matrix (8x8 for a quadrilateral element)
    Ke = zeros(SMatrix{8,8})
    node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
        # Get shape function derivatives in natural coordinates
        dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

        # Calculate the Jacobian matrix and its determinant
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
        σ[i, :] += B*u[get_dofs_linear(element)]

        # Contribution to the element stiffness matrix
        Ke += thickness * weight * B' * D * B
    end
    σ[i,:] /= 1/4
    return Ke, σ
end


# Function to assemble the global stiffness matrix
function assemble_global_stiffness_linear(nodes, elements, material)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness_linear(nodes, element, D, material.thickness)

        # Global node indices
        dofs = get_dofs_linear(element)

        # Add element stiffness to global stiffness matrix
        for i in 1:8, j in 1:8
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
end

# Function to assemble the global stiffness matrix with stress calculation for iterative methods
function assemble_global_stiffness_linear(nodes, elements, material, u)
    n_nodes = size(nodes, 2)
    K = spzeros(2 * n_nodes, 2 * n_nodes)  # Global stiffness matrix
    D = plane_stress_stiffness(material)   # Material stiffness matrix

    for element in eachcol(elements)
        Ke = element_stiffness_linear(nodes, element, D, material.thickness, u)

        # Global node indices
        dofs = get_dofs_linear(element)

        # Add element stiffness to global stiffness matrix
        for i in 1:8, j in 1:8
            K[dofs[i], dofs[j]] += Ke[i, j]
        end
    end

    return K
end

# Function to calculate the element mass matrix for a quadrilateral element using
# Gaussian quadrature
function element_mass_linear(nodes, element, material)
    # Initialize element mass matrix (8x8 for a quadrilateral element)
    Me = zeros(SMatrix{8,8})
    node_positions_transposed = @views SMatrix{4,2}(nodes[:, element]')

    for (wξ, ξ) in GAUSS_POINTS_2x2, (wη, η) in GAUSS_POINTS_2x2
        # Get shape function values in natural coordinates
        N1 = (1 - ξ) * (1 - η) / 4
        N2 = (1 + ξ) * (1 - η) / 4
        N3 = (1 + ξ) * (1 + η) / 4
        N4 = (1 - ξ) * (1 + η) / 4

        N = SVector{4}(N1, N2, N3, N4)

        # Get shape function derivatives for Jacobian calculation
        dN_dξ, dN_dη = shape_function_derivatives_linear(ξ, η)

        # Calculate the Jacobian matrix and its determinant
        J = SMatrix{2,4}(dN_dξ[1], dN_dη[1],
                         dN_dξ[2], dN_dη[2],
                         dN_dξ[3], dN_dη[3],
                         dN_dξ[4], dN_dη[4]) * node_positions_transposed
        detJ = det(J)  # Determinant of the Jacobian

        # Construct the shape function matrix N for mass matrix
        N_matrix = SMatrix{2,8}(N[1], 0.0, N[2], 0.0, N[3], 0.0, N[4], 0.0,
                               0.0, N[1], 0.0, N[2], 0.0, N[3], 0.0, N[4])

        # Weighting factor (Jacobian determinant and Gauss weights)
        weight = detJ * wξ * wη

        # Contribution to the element mass matrix
        Me += material.density * material.thickness * weight * N_matrix' * N_matrix
    end

    return Me
end

# Function to assemble the global mass matrix
function assemble_global_mass_linear(nodes, elements, material)
    n_nodes = size(nodes, 2)
    M = spzeros(2 * n_nodes, 2 * n_nodes)  # Global mass matrix

    for element in eachcol(elements)
        Me = element_mass_linear(nodes, element, material)

        # Global node indices
        dofs = get_dofs_linear(element)

        # Add element mass to global mass matrix
        for i in 1:8, j in 1:8
            M[dofs[i], dofs[j]] += Me[i, j]
        end
    end

    return M
end

function get_dofs_linear(element)
    return SVector{8}(get_dof(element[1], :x), get_dof(element[1], :y),
                      get_dof(element[2], :x), get_dof(element[2], :y),
                      get_dof(element[3], :x), get_dof(element[3], :y),
                      get_dof(element[4], :x), get_dof(element[4], :y))
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
