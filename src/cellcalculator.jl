include("index_generator.jl")

const ELEMENTS = begin
    lowercase.(readdlm("src/element_list.txt"))
end

const CELL_ID_STRING = "# Orthogonalized unit cell\n"

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

struct CellException <: Exception
    message::String
end

"""
    orthonormal_basis_matrix(cell_parameters::CellParameters)
    Compute the change of basis matrix to an orthonormal basis for the given cell parameters.
    Taken from *Essentials of Crystallography, Vol. 1* by McKie & McKie, p. 156. 
    Returns a 3x3 matrix representing the orthonormal basis.
"""
function orthonormal_basis_matrix(
    cell_parameters::CellParameters
)
    (; a, b, c, α, β, γ) = cell_parameters
    γstar = acosd( ( cosd(α)*cosd(β) - cosd(γ) )/( sind(α)*sind(β) ) )

    #Return the CoB matrix
    [ a*sind(β)*sind(γstar) 0         0; 
     -a*sind(β)*cosd(γstar) b*sind(α) 0; 
      a*cosd(β)             b*cosd(α) c]
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
    maximum_iterations::Int = 10_000_000,
    maximum_index::Int = 100
    )

    orthogonal_vector = find_orthogonal_axis(
        zone_axis, 
        cell_parameters, 
        tolerance, 
        maximum_iterations)

    if isnothing(orthogonal_vector)
        error("Could not find orthogonal axis. Increase tolerance or maximum number of iterations.")
    elseif maximum(orthogonal_vector) > maximum_index
        throw(CellException("orthogonal axis index exceeds maximum, increase tolerance"))
    end

    third_vector = find_orthogonal_axis(
        zone_axis, 
        orthogonal_vector, 
        cell_parameters)

    if maximum(third_vector) > maximum_index
        throw(CellException("orthogonal axis index exceeds maximum, increase tolerance"))
    end

    zone_axis = find_lower_index_zone(zone_axis, cell_parameters, tolerance, maximum_iterations)
    third_vector = find_lower_index_zone(third_vector, cell_parameters, tolerance, maximum_iterations)

    print_results(orthogonal_vector, third_vector, zone_axis, cell_parameters)



    #Return the change of basis matrix
    transpose([orthogonal_vector third_vector zone_axis])
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

    println("The new orthogonal unit cell is $(det([a b c])) times larger than the old unit cell.")
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
        throw(ArgumentError("Cannot calculate angle with zero-length vector: a = $a, b = $b"))
    end

    dot_clamped = clamp(dot_product / magnitude_product, -1, 1)
    return acosd(dot_clamped)
end

"""
    transform_atom_positions(
        atom_positions::AbstractMatrix{<:Real}, 
        cell_parameters::CellParameters)
    Transform the atomic positions in the unit cell defined by `cell_parameters` to
    equivalent positions in the new unit cell given a change of basis matrix between
    the two. 
    Taken from *Essentials of Crystallography, Vol. 1* by McKie & McKie, p. 164. 
"""
function transform_atom_positions(
    position_basis::AbstractVecOrMat, 
    CoBmatrix::AbstractMatrix{<:Real}
    )
    transformation_matrix = transpose(inv(CoBmatrix))
    positions = position_basis[2:4, :]

    transformed_positions = hcat([transformation_matrix * col for col in eachcol(positions)]...)
    [
        permutedims(position_basis[1,:]); 
        transformed_positions; 
        position_basis[5:end, :]
    ]
end

function transform_basis(
    position_basis::AbstractVecOrMat,
    CoBmatrix::AbstractMatrix{<:Real}
    )
    extended_positions = extend_atom_positions(position_basis, CoBmatrix)
    transformed_positions = transform_atom_positions(extended_positions, CoBmatrix)
    filter_positions(transformed_positions)
end

function extend_atom_positions(
    position_basis::AbstractVecOrMat,
    CoBmatrix::AbstractMatrix{<:Real}
    )
    for (axis, range) in zip(([1,0,0], [0,1,0], [0,0,1]), 
                            expansion_ranges(CoBmatrix))
        position_basis = expand_along_axis(position_basis, range, axis)
    end
    position_basis
end

function expand_along_axis(
    position_basis::AbstractVecOrMat,
    range::UnitRange,
    axis::AbstractVector{<:Int}
    )
    new_basis = deepcopy(position_basis)
    for translation_magnitude in range
        translation_vector = translation_magnitude * axis
        added_positions = [permutedims(position_basis[1, :]); 
                           position_basis[2:4, :] .+ translation_vector; 
                           position_basis[5:end,:]]
        new_basis = [new_basis added_positions]
    end
    new_basis
end

function expansion_ranges(
    CoBmatrix::AbstractMatrix{<:Real}
    )
    ranges = Vector{UnitRange}()

    for col in eachcol(CoBmatrix)
        min_to_fill = floor(Int64, minimum([col..., sum(col), 0]))
        max_to_fill = ceil(Int64, maximum([col..., sum(col)]))

        #Adding a few extra uc's seems to be necessary for some reason
        push!(ranges, min_to_fill-5:max_to_fill+5)
    end
    ranges
end

function filter_positions(
    positions::AbstractVecOrMat
    )
    positions[2:4, :] .= round.(positions[2:4, :], digits=5)

    positions_inside_unit_cell = positions[:, [all(0 .<= position[2:4] .<= 1) 
                                               for position in eachcol(positions)]]
    positions_without_duplicates = []
    for position in eachcol(positions_inside_unit_cell)
        if !(position ∈ positions_without_duplicates)
            push!(positions_without_duplicates, position)
        end
    end
    hcat(positions_without_duplicates...)
end

function new_cell_parameters(
    cell_parameters::CellParameters,
    CoBmatrix::AbstractMatrix{<:Real}
    )
    (; a, b, c, α, β, γ) = cell_parameters

    #Basically McKie, p. 163
    lattice_parameters = round.(sqrt.(CoBmatrix .^2 * [a, b, c] .^ 2)', digits=5)
    #Now pretend that all angles are 90°
    CellParameters(lattice_parameters..., 90, 90, 90)
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
    data = data_to_matrix(data)
    return (CellParameters(cell_parameters...), data)
end

function data_to_matrix(data)
    out_data = []
    for row in eachrow(data)
        if !(lowercase(row[1]) ∈ ELEMENTS)
            continue
        end
        push!(out_data, row)
    end
    hcat(out_data...)
end

function save_cell(
    filename::String,
    cell_parameters::CellParameters,
    data::AbstractMatrix
    )
    f = open(filename, "w")
    write(f, CELL_ID_STRING)
    writedlm(f,
             permutedims([0, 
                          cell_parameters.a, 
                          cell_parameters.b, 
                          cell_parameters.c, 
                          cell_parameters.α, 
                          cell_parameters.β, 
                          cell_parameters.γ]),
             ' ')
    writedlm(f, permutedims(data))
    close(f)
end

function make_block(
    cell_parameters::CellParameters,
    basis,
    block_size
)
    (; a, b, c, α, β, γ) = cell_parameters
    CoBmatrix = [
                block_size[1]/a 0 0;
                0 block_size[2]/b 0;
                0 0 block_size[3]/c
                ]
    new_cps = new_cell_parameters(cell_parameters, CoBmatrix)
    println(new_cps)
    new_basis = DrProbeConfig.transform_basis(basis, CoBmatrix)
    (new_cps, new_basis)
end