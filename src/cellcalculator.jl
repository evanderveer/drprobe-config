include("index_generator.jl")

struct CellParameters
    a::Real
    b::Real
    c::Real
    α::Real
    β::Real
    γ::Real
end

function orthonormal_basis_matrix(
    cell_parameters::CellParameters
)
    (; a, b, c, α, β, γ) = cell_parameters
    cx = c*cosd(β)
    cy = c*(cosd(α)-cosd(β)*cosd(γ))/sind(γ)
    cz = sqrt(c^2 - cx^2 - cy^2)
    
    #Return the CoB matrix
    [a b*cosd(γ) cx; 
     0 b*sind(γ) cy; 
     0 0         cz]
end

function change_to_orthonormal_basis(
    vector::Vector{<:Real}, 
    cell_parameters::CellParameters
    )
    CoBmatrix = orthonormal_basis_matrix(cell_parameters)
    CoBmatrix * vector
end

function change_from_orthonormal_basis(
    vector::Vector{<:Real}, 
    cell_parameters::CellParameters
    )
    CoBmatrix = orthonormal_basis_matrix(cell_parameters)
    inv(CoBmatrix) * vector
end

(yield_index, reset_index) = index_generator()

function find_lower_index_zone(
    zone_axis::Vector{<:Real}, 
    cell_parameters::CellParameters, 
    tolerance::Real, 
    max_iterations::Int
    )
    
    zone_axis_orthonormal_basis = change_to_orthonormal_basis(zone_axis, cell_parameters)
    reset_index()

    for _ in 1:max_iterations
        trial_zone_axis = collect(yield_index())

        if trial_zone_axis == zone_axis
            return zone_axis
        end

        trial_orthonormal = change_to_orthonormal_basis(trial_zone_axis, cell_parameters)

        if angle(trial_orthonormal, zone_axis_orthonormal_basis) < tolerance
            return trial_zone_axis
        end
    end
    
    zone_axis
end

function find_orthogonal_axis(
    zone_axis, 
    cell_parameters, 
    tolerance, 
    find_lower_index,
    max_iterations
    )
    proj_vector_ortho = change_to_orthonormal_basis(zone_axis, cell_parameters)
    reset_index()
    for _ in 1:max_iterations
        index = yield_index()
        vec = change_to_orthonormal_basis([index...], cell_parameters)
        angle = acosd(clamp(dot(proj_vector_ortho, vec) / (norm(proj_vector_ortho) * norm(vec)), -1, 1))
        if abs(angle .- 90) < tolerance
            return [index...]
        end
    end
    println("No lower index zone within tolerance")
end

function find_third_vector(
    vector_1,
    vector_2, 
    cell_parameters, 
    tolerance, 
    find_lower_index, 
    max_iterations
    )

    vec_1_ortho = change_to_orthonormal_basis(vector_1, cell_parameters)
    vec_2_ortho = change_to_orthonormal_basis(vector_2, cell_parameters)

    third_vector = change_from_orthonormal_basis(cross(vec_1_ortho, vec_2_ortho), cell_parameters)
    if sqrt(sum(third_vector .^ 2)) > 30
        println("High index lattice vector, trying to find lower index vector")
    end
    third_vector = round.(Int32, find_lower_index_zone(third_vector, cell_parameters, tolerance, max_iterations))
    println("New lattice vector: $third_vector")
    third_vector
end

function find_orthogonal_cell(
    zone_axis::Vector{<:Int},
    cell_parameters::CellParameters; 
    tolerance::Real = 1, 
    find_lower_index::Bool = true, 
    max_iterations::Int = 1E7
    )
    if find_lower_index && sqrt(sum(zone_axis .^ 2)) > 30
        println("High index zone axis, trying to find lower index zone")
        zone_axis = find_lower_index_zone(zone_axis, cell_parameters, tolerance, max_iterations)
        println("New zone axis: $zone_axis")
    end

    orthogonal_vector = find_orthogonal_axis(zone_axis, cell_parameters, tolerance, find_lower_index, max_iterations)
    third_vector = find_third_vector(zone_axis, orthogonal_vector, cell_parameters, tolerance, find_lower_index, max_iterations)

    println("\nNew orthogonal cell found. a = $(Tuple(orthogonal_vector)), b = $(Tuple(third_vector)), c (ZA) = $(Tuple(zone_axis))")

    angles = round.(90 .- find_vector_angles([orthogonal_vector, third_vector, zone_axis], cell_parameters), digits=4)
    println("Orthogonal cell angle error: α = $(angles[1])°, β = $(angles[2])°, γ = $(angles[3])°")

    #Return the change of basis matrix
    [orthogonal_vector third_vector zone_axis]
end

function find_vector_angles(
    vectors::Vector{<:Vector{<:Real}}, 
    cell_parameters::CellParameters
    )
    orthonormal_vectors = change_to_orthonormal_basis.(vectors, cell_parameters)

    α = angle(orthonormal_vectors[2], orthonormal_vectors[3])
    β = angle(orthonormal_vectors[1], orthonormal_vectors[3])
    γ = angle(orthonormal_vectors[1], orthonormal_vectors[2])

    [α, β, γ]
end

angle(a, b) = acosd(clamp(a ⋅ b / (norm(a) * norm(b)), -1, 1))

function load_cell(
    filename::String
)
    if splitext(filename)[2] != ".cel"
        throw(ArgumentError("filename must be a .cel file"))
    end
    f = open(filename)
    readline(f) #Skip header line
    cell_parameters = [c for c in split(readline(f), " ") if c != ""]
    data = readdlm(f)
    close(f)
    cell_parameters = parse.(Float64, cell_parameters[2:end])
    return (CellParameters(cell_parameters...), data)
end



abs2