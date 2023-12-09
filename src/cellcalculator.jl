macro index_list(max_index)
    indices = vcat(
                collect(
                    Iterators.product(
                        -max_index:max_index, 
                        -max_index:max_index, 
                        -max_index:max_index)
                        )
                    ...)
    sort(indices, by=x->sum(x.^2))[2:end]
end

const INDEX_LIST = @index_list 50

function change_basis(vector, a, b, c, α, β, γ)
    cx = c*cosd(β)
    cy = c*(cosd(α)-cosd(β)*cosd(γ))/sind(γ)
    cz = sqrt(c^2 - cx^2 - cy^2)

    A = [a b*cosd(γ) cx; 
         0 b*sind(γ) cy; 
         0 0         cz] * vector
end

function find_lower_index_zone(zone_axis, cell_parameters, tolerance)
    proj_vector_ortho = change_basis(zone_axis, cell_parameters...)
    for index in INDEX_LIST
        vec = change_basis([index...], cell_parameters...)
        angle = acosd(clamp(dot(proj_vector_ortho, vec)/(norm(proj_vector_ortho)*norm(vec)), -1, 1))
        if angle < tolerance
            return [index...]
        end
    end
    println("No lower index zone within tolerance")
end

function find_orthogonal_axis(zone_axis, cell_parameters, tolerance)
    proj_vector_ortho = change_basis(zone_axis, cell_parameters...)

    for index in INDEX_LIST
        vec = change_basis([index...], cell_parameters...)
        angle = acosd(clamp(dot(proj_vector_ortho, vec)/(norm(proj_vector_ortho)*norm(vec)), -1, 1))
        if abs(angle .- 90) < tolerance
            return [index...]
        end
    end
    println("No lower index zone within tolerance")
end

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
    return (parse.(Float64, cell_parameters[2:end]), data)
end


