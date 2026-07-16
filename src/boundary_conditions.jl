# Function to apply boundary conditions
function apply_boundary_conditions!(K, f, fixed_dofs::Vector{Int})
    for dof in fixed_dofs
        K[dof, :] .= 0.0
        K[:, dof] .= 0.0
        K[dof, dof] = 1.0
        f[dof] = 0.0
    end
    dropzeros!(K)
    return nothing
end

# Function to apply boundary conditions
function apply_boundary_conditions!(K, f, fixed_nodes::Dict)
    fixed_dofs = get_fixed_dofs(fixed_nodes)
    apply_boundary_conditions!(K, f, fixed_dofs)
    return nothing
end

# Function to apply boundary conditions to mass matrix
function apply_boundary_conditions_mass!(M, fixed_dofs::Vector{Int})
    for dof in fixed_dofs
        M[dof, :] .= 0.0
        M[:, dof] .= 0.0
        M[dof, dof] = 1.0  # Keep diagonal term for stability
    end
    dropzeros!(M)
    return nothing
end

# Function to apply boundary conditions to mass matrix
function apply_boundary_conditions_mass!(M, fixed_nodes::Dict)
    fixed_dofs = get_fixed_dofs(fixed_nodes)
    apply_boundary_conditions_mass!(M, fixed_dofs)
    return nothing
end

function get_fixed_dofs(fixed_nodes::Dict{Int,Symbol})
    fixed_dofs = Vector{Int}()
    for (node, dimsymb) in fixed_nodes
        if dimsymb === :x
            push!(fixed_dofs, get_dof(node, :x))
        elseif dimsymb === :y
            push!(fixed_dofs, get_dof(node, :y))
        elseif dimsymb === :xy
            push!(fixed_dofs, get_dof(node, :x), get_dof(node, :y))
        else
            msg = "unknown dimension specifier :$(dimsymb)!\n"
            throw(ArgumentError(msg))
        end
    end
    return fixed_dofs
end

function get_force_vector(nodes, forces)
    f = zeros(2 * size(nodes, 2))
    for (node, force) in forces
        dofx, dofy = get_dof(node, :x), get_dof(node, :y)
        f[dofx] = force[1]
        f[dofy] = force[2]
    end
    return f
end
