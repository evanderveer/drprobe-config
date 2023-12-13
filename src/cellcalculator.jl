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
    
    #Transform the zone axis into orthonormal basis
    zone_axis_orthonormal_basis = change_to_orthonormal_basis(zone_axis, cell_parameters)
    reset_index()

    for _ in 1:max_iterations
        #Try out potential parallel vectors of increasing index
        trial_zone_axis = collect(yield_index())

        if trial_zone_axis == zone_axis
            return zone_axis
        end

        trial_orthonormal = change_to_orthonormal_basis(trial_zone_axis, cell_parameters)

        if angle(trial_orthonormal, zone_axis_orthonormal_basis) < tolerance
            return trial_zone_axis
        end
    end
    
    #If no lower index vector, return the input rounded to nearest integer
    round.(Int, zone_axis)
end

function find_orthogonal_axis(
    zone_axis::Vector{<:Real}, 
    cell_parameters::CellParameters, 
    tolerance::Real, 
    max_iterations::Int
    )

    #Transform the zone axis into orthonormal basis
    zone_axis_orthonormal = change_to_orthonormal_basis(zone_axis, cell_parameters)
    reset_index()

    for _ in 1:max_iterations
        #Try out potential orthogonal vectors of increasing index
        trial_zone_axis = collect(yield_index())
        trial_orthonormal = change_to_orthonormal_basis(trial_zone_axis, cell_parameters)
        if abs(angle(zone_axis_orthonormal, trial_orthonormal) .- 90) < tolerance
            return trial_zone_axis
        end
    end

end

function find_orthogonal_axis(
    vector_1::Vector{<:Real}, 
    vector_2::Vector{<:Real}, 
    cell_parameters::CellParameters, 
    )
    vector_1_orthogonal = change_to_orthonormal_basis(vector_1, cell_parameters)
    vector_2_orthogonal = change_to_orthonormal_basis(vector_2, cell_parameters)

    change_from_orthonormal_basis(vector_1_orthogonal × vector_2_orthogonal, cell_parameters)
end

function find_orthogonal_cell(
    zone_axis::Vector{<:Int},
    cell_parameters::CellParameters; 
    tolerance::Real = 1, 
    max_iterations::Int = 10_000_000
    )

    orthogonal_vector = find_orthogonal_axis(
        zone_axis, 
        cell_parameters, 
        tolerance, 
        max_iterations)

    if isnothing(orthogonal_vector)
        error("Could not find orthogonal axis. Increase tolerance or maximum number of iterations.")
    end

    third_vector = find_orthogonal_axis(
        zone_axis, 
        orthogonal_vector, 
        cell_parameters)

    zone_axis = find_lower_index_zone(zone_axis, cell_parameters, tolerance, max_iterations)
    third_vector = find_lower_index_zone(third_vector, cell_parameters, tolerance, max_iterations)

    print_results(orthogonal_vector, third_vector, zone_axis, cell_parameters)

    #Return the change of basis matrix
    [orthogonal_vector third_vector zone_axis]
end

function print_results(a, b, c, cell_parameters)
    println("\nNew orthogonal cell found. a = $(Tuple(a)), b = $(Tuple(b)), c (ZA) = $(Tuple(c))")

    angles = round.(90 .- find_vector_angles([a, b, c], cell_parameters), digits=4)
    println("Orthogonal cell angle deviation from 90°: α = $(angles[1])°, β = $(angles[2])°, γ = $(angles[3])°")
end

function find_vector_angles(
    vectors::Vector{<:Vector{<:Real}}, 
    cell_parameters::CellParameters
    )
    orthonormal_vectors = change_to_orthonormal_basis.(vectors, Ref(cell_parameters))

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