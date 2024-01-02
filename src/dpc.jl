function construct_iDPC(
    images::Tuple{T, T, T, T},
    size_nm::Tuple{<:Real, <:Real}
    ) where T<:AbstractMatrix{<:Real}

    KX = similar(images[1])
    KY = similar(images[1])
    Ksquared = similar(images[1])

    fill_reciprocal!(KY, KX, Ksquared, size_nm)

    DPCx = images[1] .- images[3]
    DPCy = images[2] .- images[4]

    fDPCx = fftshift(fft(DPCx))
    fDPCy = fftshift(fft(DPCy))

    fftiDPC = (KX .* fDPCx .+ KY .* fDPCy)  ./ (2 * π * 1im * Ksquared)

    real(ifft(ifftshift(fftiDPC)))
end

function reciprocalGrid(metadata, size_px, size_nm)
    # Central pixel coordinate (x, y)
    center = size_px ./2 

    # Reciprocal space grid and magnitude [1/nm]
    kx = (1 / size_nm[2]) * ((1:metadata.nx) .- center[1])
    ky = (1 / size_nm[2]) * ((1:metadata.ny) .- center[2])

    KX, KY = meshgrid(kx, ky)

    return KX, KY, center
end

function fill_reciprocal!(
    KY,
    KX,
    Ksquared,
    size_nm
    )

    reciprocal_lattice_vectors = 1 ./ size_nm

    for index in CartesianIndices(KY)
        index_tuple = Tuple(index)
        reciprocal_position = index_tuple .* reciprocal_lattice_vectors
        KY[index] = reciprocal_position[1]
        KX[index] = reciprocal_position[2]
        Ksquared[index] = sum(reciprocal_position .^ 2)
    end
end

function make_square(image)
    im_size = size(image)
    if im_size[1] < im_size[2]
        new_im = [image; image][1:im_size[2], 1:im_size[2]]
    else
        new_im = [image image][1:im_size[1], 1:im_size[1]]
    end
    new_im
end