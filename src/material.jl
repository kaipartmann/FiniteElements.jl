# Define the material properties (Plane stress)
struct Material
    E::Float64  # Young's modulus
    ν::Float64  # Poisson's ratio
    thickness::Float64  # Thickness of the structure
    density::Float64 # Density of Material
end

# Constructor with default density for backward compatibility
Material(E, ν, thickness) = Material(E, ν, thickness, 7850.0)

function plane_stress_stiffness(material::Material)
    (; E, ν) = material
    factor = E / (1.0 - ν^2)
    D = factor * SMatrix{3,3}(1.0, ν, 0.0, ν, 1.0, 0.0, 0.0, 0.0, (1.0 - ν) / 2.0)
    return D
end
