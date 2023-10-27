using Gridap.CellData
using Infiltrator

function _get_hat_function_cellfield(i, basis_data, model)
  Ω = Triangulation(model)
  cell_to_ith_node(c) = c[i]
  A = lazy_map(cell_to_ith_node, basis_data)
  GenericCellField(A, Ω, ReferenceDomain())
end

function _get_hat_functions_on_cells(model)
  # Always order 1
  reffe = ReferenceFE(lagrangian, Float64, 1)
  V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  fe_basis = get_fe_basis(V0)
  bd = Gridap.CellData.get_data(fe_basis)
  #@show typeof(bd)
  bd
end

function _build_patch_RHS_vectors(𝐀ₕ, f, model, weight, dvRT, dvp, Qₕ)
  hat_fns_on_cells = _get_hat_functions_on_cells(model)
  RHS_RT_form(ψ) = ∫(ψ * (weight ⋅ 𝐀ₕ ⋅ dvRT))*Qₕ
  RHS_L²_form(ψ) = ∫((f * ψ + 𝐀ₕ ⋅ ∇(ψ))*dvp)*Qₕ
  cur_num_cells = num_cells(model)
  # Hardcoded for triangles
  nodes_per_cell = 3
  # TODO: find eltype of the vectors
  cell_RHS_RTs = Matrix{Vector{Float64}}(undef, cur_num_cells, nodes_per_cell)
  cell_RHS_L²s = Matrix{Vector{Float64}}(undef, cur_num_cells, nodes_per_cell)
  Ω = Triangulation(model)
  for i = 1:nodes_per_cell
    ψᵢ = _get_hat_function_cellfield(i, hat_fns_on_cells, model)
    #writevtk(Ω, "psi$(i)", cellfields=["ψᵢ" => ψᵢ])
    cell_RHS_RT = RHS_RT_form(ψᵢ)
    cell_RHS_L² = RHS_L²_form(ψᵢ)
    #@show typeof(parallel_smart_collect(cell_RHS_RT))
    cell_RHS_RTs[:, i] = lazy_map(vec, parallel_smart_collect(cell_RHS_RT))
    cell_RHS_L²s[:, i] = lazy_map(vec, parallel_smart_collect(cell_RHS_L²))
  end
  cell_RHS_RTs, cell_RHS_L²s
end

function _build_cellwise_matrices(duRT, weight, dvRT, dvp, Qₕ)
  cell_mass_mats = ∫(weight ⋅ duRT ⋅ dvRT) * Qₕ
  cell_mixed_mats = ∫((∇ ⋅ duRT) * dvp) * Qₕ
  @exfiltrate
  cell_mass_mats = parallel_smart_collect(cell_mass_mats)
  cell_mixed_mats = parallel_smart_collect(cell_mixed_mats)
  cell_mass_mats, cell_mixed_mats
end

function _build_lagange_row(dvp, Qₕ)
  #cur_num_cells = num_cells(model)
  #cell_Λ_vecs = Matrix{Vector{Float64}}(undef, cur_num_cells, 1)
  cell_Λ_vecs = ∫(1 * dvp) * Qₕ
  #cell_Λ_vecs = cell_Λ_cells 
  cell_Λ_vecs = map(vec, parallel_smart_collect(cell_Λ_vecs))
end

function build_all_cellwise_objects(𝐀ₕ, f, weight, spaces, model, RT_order, measure)
  Tₕ = Triangulation(model)
  if isnothing(measure)
    measure = Measure(Tₕ, 2 * RT_order + 2)
  end
  Qₕ = CellData.get_cell_quadrature(measure)
  #Qₕ = CellQuadrature(Tₕ, quad_order)
  dvp = get_trial_fe_basis(spaces.L²_space)
  duRT = get_fe_basis(spaces.RT_space)
  dvRT = get_trial_fe_basis(spaces.RT_space)
  args = (weight, dvRT, dvp, Qₕ)
  cell_mass_mats, cell_mixed_mats = _build_cellwise_matrices(duRT, args...)
  cell_RHS_RTs, cell_RHS_L²s = _build_patch_RHS_vectors(𝐀ₕ, f, model, args...)
  cell_Λ_vecs = _build_lagange_row(dvp, Qₕ)
  (; cell_mass_mats, cell_Λ_vecs, cell_mixed_mats, cell_RHS_RTs, cell_RHS_L²s)
end
