using DrProbeConfig
zone_axis = [1,2,-5]
cell_parameters, basis = DrProbeConfig.load_cell("test/BTO.cel")

cob = DrProbeConfig.find_orthogonal_cell(zone_axis, cell_parameters, tolerance = 0.3)

new_parameters = DrProbeConfig.new_cell_parameters(cell_parameters, cob)

new_basis = DrProbeConfig.transform_basis(basis, cob)

(block_parameters, block_basis) = DrProbeConfig.make_block(new_parameters, new_basis, [1,1,10])

DrProbeConfig.save_cell("test/BTOblock.cel", block_parameters, block_basis)
