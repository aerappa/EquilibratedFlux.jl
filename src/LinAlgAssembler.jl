using SparseArrays
using Gridap.FESpaces
using Gridap.Arrays
using FastLapackInterface

struct LinAlgObjs end

instantiate_M(RT_dofs_per_cell, max_patch_cells) =
  zeros(RT_dofs_per_cell * max_patch_cells, RT_dofs_per_cell * max_patch_cells)

instantiate_B(RT_dofs_per_cell, L²_dofs_per_cell, max_patch_cells) =
  zeros(RT_dofs_per_cell * max_patch_cells, L²_dofs_per_cell * max_patch_cells)

instantiate_A(M, B) = [M B; transpose(B) zeros(size(B)[2], size(B)[2])]

instantiate_RHS_RT(RT_dofs_per_cell, max_patch_cells) =
  zeros(max_patch_cells * RT_dofs_per_cell)

instantiate_RHS_L²(L²_dofs_per_cell, max_patch_cells) =
  zeros(max_patch_cells * L²_dofs_per_cell)

instantiate_RHS(RHS_RT, RHS_L²) = [RHS_RT; RHS_L²]

instantiate_Λ(RHS_L²) = zeros(length(RHS_L²))

function get_dofs_per_cell(k, d)
  RT_dofs_per_cell = (k + d + 1) * binomial(k + d - 1, k)
  L²_dofs_per_cell = binomial(k + d, k)
  (RT_dofs_per_cell, L²_dofs_per_cell)
end

function test_dofs_per_cell(k, d, spaces)
  RT_cell_dof_ids = get_cell_dof_ids(spaces.RT_space)
  L²_cell_dof_ids = get_cell_dof_ids(spaces.L²_space)
  num_dofs_per_cell_RTs = lazy_map(length, RT_cell_dof_ids)
  num_dofs_per_cell_L²s = lazy_map(length, L²_cell_dof_ids)
  RT_dofs_per_cell_formula = (k + d + 1) * binomial(k + d - 1, k)
  L²_dofs_per_cell_formula = binomial(k + d, k)
  cache_RT = array_cache(num_dofs_per_cell_RTs)
  cache_L² = array_cache(num_dofs_per_cell_L²s)
  for i = 1:length(num_dofs_per_cell_RTs)
    cell_dofs_RT = getindex!(cache_RT, num_dofs_per_cell_RTs, i)
    @assert cell_dofs_RT == RT_dofs_per_cell_formula
  end
  for i = 1:length(num_dofs_per_cell_L²s)
    cell_dofs_L² = getindex!(cache_L², num_dofs_per_cell_L²s, i)
    @assert cell_dofs_L² == L²_dofs_per_cell_formula
  end
end

function instantiate_linalg(RT_order, dim, metadata)
  (RT_dofs_per_cell, L²_dofs_per_cell) = get_dofs_per_cell(RT_order, dim)
  max_patch_cells = metadata.max_num_patch_cells
  M = instantiate_M(RT_dofs_per_cell, max_patch_cells)
  B = instantiate_B(RT_dofs_per_cell, L²_dofs_per_cell, max_patch_cells)
  A = instantiate_A(M, B)
  # Pre-allocate the pivot vector for the matrix A
  ws = LUWs(A)
  RHS_RT = instantiate_RHS_RT(RT_dofs_per_cell, max_patch_cells)
  RHS_L² = instantiate_RHS_L²(L²_dofs_per_cell, max_patch_cells)
  RHS = instantiate_RHS(RHS_RT, RHS_L²)
  Λ = similar(RHS_L²)
  σ_loc = similar(RHS)
  (; M, B, A, ws, Λ, RHS_RT, RHS_L², RHS, σ_loc)
end
