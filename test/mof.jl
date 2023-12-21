using DrProbeConfig
zone_axis = [-44910,52138,-36471]
cell_parameters, basis = DrProbeConfig.load_cell("test/mof/mof.cel")

cob = DrProbeConfig.find_orthogonal_cell(zone_axis, cell_parameters, tolerance = 0.3)

new_parameters = DrProbeConfig.new_cell_parameters(cell_parameters, cob)

new_basis = DrProbeConfig.transform_basis(basis, cob)

(block_parameters, block_basis) = DrProbeConfig.make_block(new_parameters, new_basis)

DrProbeConfig.save_cell("test/mof/mof.cel", new_parameters, new_basis)
