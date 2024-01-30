function make_gaussian_kernel(
    config::Dict
    )

    resolution = (config["scan-frame"]["resolution"]["y"], config["scan-frame"]["resolution"]["x"])
    size_nm = config["scan-frame"]["size"]["y"], config["scan-frame"]["size"]["x"]
    source_radius = config["source-radius"]

    sample_rates = size_nm ./ resolution
    kernel_cutoff = 3 * source_radius

    kernel_dimensions = ceil.(kernel_cutoff ./ sample_rates) 

    kernel = zeros(Float32, (round.(Int, 2 .* kernel_dimensions .+ 1))...)
    kernel_offset = Origin(-kernel_dimensions[1], -kernel_dimensions[2])(kernel)
    

    kernprm = -log(2)/source_radius^2

    for index in CartesianIndices(kernel_offset)
        sum_of_squares_index = sum((sample_rates .* Tuple(index)) .^ 2)

        if sum_of_squares_index > kernel_cutoff^2
            continue
        end

        kernel_offset[index] = exp(sum_of_squares_index * kernprm)
    end

    kernel_offset ./ sum(kernel_offset)
end

function apply_spatial_coherence(
    config::Dict,
    image::AbstractMatrix{<:Real}
    )
    type = config["convolutions"]["source-size"]["profile"]
    if type != "gaussian"
        error("only gaussian source profile is currently supported")
    end

    kernel = make_gaussian_kernel(config)

    imfilter(image, kernel, "circular")
end


