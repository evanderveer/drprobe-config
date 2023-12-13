include("index_generator.jl")

"""
    struct CellParameters
    Represents the parameters of a crystallographic cell.

    Fields:
    - `a::Real`: Length of the first cell vector.
    - `b::Real`: Length of the second cell vector.
    - `c::Real`: Length of the third cell vector.
    - `α::Real`: Angle between the second and third cell vectors (in degrees).
    - `β::Real`: Angle between the first and third cell vectors (in degrees).
    - `γ::Real`: Angle between the first and second cell vectors (in degrees).
"""
struct CellParameters
    a::Real
    b::Real
    c::Real
    α::Real
    β::Real
    γ::Real
end

"""
    orthonormal_basis_matrix(cell_parameters::CellParameters)
    Compute the change of basis matrix to an orthonormal basis for the given cell parameters.

    Returns a 3x3 matrix representing the orthonormal basis.
"""
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

"""
    change_to_orthonormal_basis(vector::Vector{<:Real}, cell_parameters::CellParameters)
    Transform a vector from the lattice vector basis to the orthonormal basis.

    Returns the vector in the orthonormal basis.
"""
function change_to_orthonormal_basis(
    vector::Vector{<:Real}, 
    cell_parameters::CellParameters
    )
    CoBmatrix = orthonormal_basis_matrix(cell_parameters)
    CoBmatrix * vector
end

"""
    change_from_orthonormal_basis(vector::Vector{<:Real}, cell_parameters::CellParameters)
    Transform a vector from the orthonormal basis back to the lattice vector basis.

    Returns the vector in the lattice vector basis.
"""
function change_from_orthonormal_basis(
    vector::Vector{<:Real}, 
    cell_parameters::CellParameters
    )
    CoBmatrix = orthonormal_basis_matrix(cell_parameters)
    inv(CoBmatrix) * vector
end

(yield_index, reset_index) = index_generator()

"""
    find_lower_index_zone(
        zone_axis::Vector{<:Real},
        cell_parameters::CellParameters,
        tolerance::Real,
        max_iterations::Int
    )
    Find a lower index zone axis given the initial zone axis and cell parameters.

    Returns the lower index zone axis within the specified tolerance and maximum iterations.
"""
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

"""
    find_orthogonal_axis(
        zone_axis::Vector{<:Real},
        cell_parameters::CellParameters,
        tolerance::Real,
        max_iterations::Int
    )
    Find an orthogonal axis to the given zone axis within the specified tolerance and maximum iterations.

    Returns the orthogonal axis.
"""
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

"""
    find_orthogonal_axis(
        vector_1::Vector{<:Real},
        vector_2::Vector{<:Real},
        cell_parameters::CellParameters
    )
    Find an orthogonal axis given two vectors and cell parameters.

    Returns the orthogonal axis.
"""
function find_orthogonal_axis(
    vector_1::Vector{<:Real}, 
    vector_2::Vector{<:Real}, 
    cell_parameters::CellParameters, 
    )
    vector_1_orthogonal = change_to_orthonormal_basis(vector_1, cell_parameters)
    vector_2_orthogonal = change_to_orthonormal_basis(vector_2, cell_parameters)

    change_from_orthonormal_basis(vector_1_orthogonal × vector_2_orthogonal, cell_parameters)
end

"""
    find_orthogonal_cell(
        zone_axis::Vector{<:Int},
        cell_parameters::CellParameters;
        tolerance::Real = 1,
        max_iterations::Int = 10_000_000
    )
    Tries to find an orthogonal unit cell given a zone axis and cell parameters.

    Returns the change of basis matrix from the given unit cell to the new orthogonal cell.
"""
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

"""
    print_results(a, b, c, cell_parameters)
    Print information about the new orthogonal cell.

    Displays the lengths and angles of the orthogonal cell vectors.
"""
function print_results(a, b, c, cell_parameters)
    println("\nNew orthogonal cell found. a = $(Tuple(a)), b = $(Tuple(b)), c (ZA) = $(Tuple(c))")

    angles = round.(90 .- find_vector_angles([a, b, c], cell_parameters), digits=4)
    println("Orthogonal cell angle deviation from 90°: α = $(angles[1])°, β = $(angles[2])°, γ = $(angles[3])°")
end

"""
    find_vector_angles(vectors::Vector{<:Vector{<:Real}}, cell_parameters::CellParameters)
    Calculate the angles between vectors in the orthonormal basis.

    Returns a vector containing the angles α, β, and γ (in degrees).
"""
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

"""
    angle(a::Vector{<:Real}, b::Vector{<:Real})
    Calculate the angle in degrees between two vectors a and b.
"""
function angle(
    a::Vector{<:Real}, 
    b::Vector{<:Real}
    )
    if length(a) != length(b)
        throw(ArgumentError("Vectors must have the same length"))
    end

    dot_product = a ⋅ b
    magnitude_product = norm(a) * norm(b)
    
    if magnitude_product ≈ 0.0
        throw(ArgumentError("Cannot calculate angle with zero-length vector"))
    end

    dot_clamped = clamp(dot_product / magnitude_product, -1, 1)
    return acosd(dot_clamped)
end

"""
    load_cell(filename::String)
    Load cell parameters and atom positions in the unit cell from a .cel file.

    Returns a tuple containing the CellParameters and the atom positions.
"""
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