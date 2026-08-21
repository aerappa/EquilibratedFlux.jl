using Gridap.FESpaces
using Gridap.ReferenceFEs

#=
For a Neumann-type patch centered at vertex a, the equilibrated flux
reconstruction σₐ built on the patch must have its RT normal-trace dofs on
a's own incident ΓN edges *prescribed* (essential condition) to the
ψₐ-weighted moments of the Neumann datum, i.e. dof_i(ψₐ σN) for each RT dof
functional dof_i associated with that edge. Summing the contributions from
an edge's two endpoint patches exactly reconstructs the RT interpolant of
σN on that edge, since ψₐ + ψᵦ ≡ 1 there (partition of unity).

We compute this via Gridap's own Hdiv Dirichlet-value machinery
(`TrialFESpace`), which correctly handles the Piola map / edge orientation
that would otherwise need to be reproduced by hand. Rather than calling it
once per Neumann vertex, we reuse the same "local hat function" trick
already used for the volume RHS in CellwiseAssembler.jl: for i = 1:3 (local
vertex position within a cell), compute the lift using the field that
equals, on each cell, that cell's i-th local P1 basis function. Restricted
to any given patch's own cells, and picking out the correct local position
per cell via `node_to_offsets`, this reproduces exactly ψₐ.
=#
function compute_neumann_lift(patches, model, spaces, RT_order, neumann_tags, neumann_data)
  if isempty(neumann_tags)
    return patches
  end
  reffeRT = ReferenceFE(raviart_thomas, Float64, RT_order)
  V_N = FESpace(model, reffeRT; conformity = :Hdiv, dirichlet_tags = neumann_tags)
  Ω = Triangulation(model)
  gN = CellField(neumann_data, Ω)

  hat_fns_on_cells = _get_hat_functions_on_cells(model)
  nodes_per_cell = 3
  dirichlet_values = Vector{Vector{Float64}}(undef, nodes_per_cell)
  for i = 1:nodes_per_cell
    ψᵢ = _get_hat_function_cellfield(i, hat_fns_on_cells, model)
    U_i = TrialFESpace(V_N, ψᵢ * gN)
    dirichlet_values[i] = U_i.dirichlet_values
  end

  topo = get_grid_topology(model)
  cell_to_edge = Geometry.get_faces(topo, 2, 1)
  cell_to_edge_cache = array_cache(cell_to_edge)
  concrete_reffeRT = ReferenceFEs.ReferenceFE(TRI, raviart_thomas, Float64, RT_order)
  local_edge_to_dofs = get_face_own_dofs(concrete_reffeRT)[4:6]

  main_cell_dofs = get_cell_dof_ids(spaces.RT_space)
  vn_cell_dofs = get_cell_dof_ids(V_N)
  main_cache = array_cache(main_cell_dofs)
  vn_cache = array_cache(vn_cell_dofs)

  new_patches = Vector{Patch}(undef, length(patches))
  for (idx, patch) in enumerate(patches)
    if !(patch isa NeumannPatch)
      new_patches[idx] = patch
      continue
    end
    dof_ids = Int[]
    dof_values = Float64[]
    for (k, cellid) in enumerate(patch.data.patch_cell_ids)
      offset = patch.data.node_to_offsets[k]
      cell_edges = getindex!(cell_to_edge_cache, cell_to_edge, cellid)
      target_local_edges = findall(e -> e in patch.edge_ids, cell_edges)
      isempty(target_local_edges) && continue
      vn_ids = getindex!(vn_cache, vn_cell_dofs, cellid)
      main_ids = getindex!(main_cache, main_cell_dofs, cellid)
      dvals = dirichlet_values[offset]
      for local_e in target_local_edges
        for ℓ in local_edge_to_dofs[local_e]
          @assert vn_ids[ℓ] < 0 "expected Neumann-tagged edge dof to be Dirichlet in V_N"
          push!(dof_ids, main_ids[ℓ])
          push!(dof_values, dvals[-vn_ids[ℓ]])
        end
      end
    end
    new_patches[idx] = NeumannPatch(patch.data, patch.edge_ids, dof_ids, dof_values)
  end
  new_patches
end
