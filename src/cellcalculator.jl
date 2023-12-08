function change_basis(vector, a, b, c, α, β, γ)
    cx = c*cosd(β)
    cy = c*(cosd(α)-cosd(β)*cosd(γ))/sind(γ)
    cz = sqrt(c^2 - cx^2 - cy^2)

    A = [a b*cosd(γ) cx; 
         0 b*sind(γ) cy; 
         0 0         cz] * vector
end

function find_orthogonal_axis(projection_vector, cell_parameters, max_index)
    proj_vector_ortho = change_basis(projection_vector, cell_parameters...)
    angles = []
    indices = []

    for h in -max_index:max_index
    for k in -max_index:max_index
    for l in -max_index:max_index
        if h==k==l==0; continue; end
        vec = change_basis([h,k,l], cell_parameters...)
        angle = acosd(clamp(dot(proj_vector_ortho, vec)/(norm(proj_vector_ortho)*norm(vec)), -1, 1))
        push!(angles, angle)
        push!(indices, (h,k,l))
    end
    end
    end

    sort_perm = sortperm(abs.(angles .- 90))

    indices[sort_perm][1]
end