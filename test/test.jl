using DrProbeConfig
zone_axis = [1,1,0]
cell_parameters, basis = DrProbeConfig.load_cell("test/BTO.cel")
cob = DrProbeConfig.find_orthogonal_cell(zone_axis, cell_parameters, tolerance = 1)

new_parameters = DrProbeConfig.new_cell_parameters(cell_parameters, cob)
new_basis = DrProbeConfig.transform_basis(basis, cob)
DrProbeConfig.save_cell("test/testBTO.cel", new_parameters, new_basis)
