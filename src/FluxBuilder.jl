using LinearAlgebra
using Gridap
using Gridap.Adaptivity
using ChunkSplitters
#using TimerOutputs
#
#macro timeit(timer, title, body)
#          esc(body)
#end

#const to = TimerOutput()
#
function build_equilibrated_flux(𝐀ₕ, f, model::AdaptedDiscreteModel, RT_order; 
    measure = nothing, weight = 1.0, return_contributions=false)
  build_equilibrated_flux(𝐀ₕ, f, model.model, RT_order, measure = measure, weight = weight, return_contributions = return_contributions)
end


"""
    build_equilibrated_flux(𝐀ₕ, f, model, RT_order; measure = nothing, weight= 1.0)

TODO: relevant docstring
"""
function build_equilibrated_flux(𝐀ₕ, f, model, RT_order; measure = nothing, weight= 1.0, return_contributions=false)
  topo = get_grid_topology(model)
  @assert all(p->p==TRI,get_polytopes(topo))
  patches, metadata = create_patches(model, RT_order)
  spaces = build_global_spaces(model, RT_order)
  cell_objects = build_all_cellwise_objects(𝐀ₕ, f, weight, spaces, model, RT_order, measure)
  linalgs = [instantiate_linalg(RT_order, 2, metadata) for i = 1:Threads.nthreads()]
  dms = build_DOFManagers(spaces)
  σ_gl = zero(spaces.RT_space)
  diri_patches::Vector{DirichletPatch{Int32}} =
    filter(patch -> patch isa DirichletPatch, patches)
  int_patches::Vector{InteriorPatch{Int32}} =
    filter(patch -> patch isa InteriorPatch, patches)
  η_loc = return_contributions ? zeros(length(patches)) : nothing
  build_equilibrated_flux(patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms; η=η_loc)
  if return_contributions
    return η_loc, σ_gl # η_loc contains the squared norms 0.5 * (σ_loc - Ah)^2 for each patch
  else
    return σ_gl
  end
end

function matrix_scatter!(patch_mat, cell_mats, dm_col, dm_row, patch_data)
  fill!(patch_mat, 0)
  for cellid in patch_data.patch_cell_ids
    update_cell_local_dofs!(dm_row, cellid)
    update_cell_local_dofs!(dm_col, cellid)
    cell_mat = cell_mats[cellid]
    for i in axes(cell_mat, 1)
      for j in axes(cell_mat, 2)
        patch_mat[dm_row.cell_dofs_loc[i], dm_col.cell_dofs_loc[j]] += cell_mat[i, j]
      end
    end
  end
end

function single_vector_scatter!(patch_vec, cell_vecs, dm, patch_data)
  fill!(patch_vec, 0)
  for cellid in patch_data.patch_cell_ids
    update_cell_local_dofs!(dm, cellid)
    cell_vec = cell_vecs[cellid]
    for i in axes(cell_vec, 1)
      patch_vec[dm.cell_dofs_loc[i]] += cell_vec[i]
    end
  end
end

function vector_scatter!(patch_vec, cell_vecs, dm, patch_data)
  fill!(patch_vec, 0)
  node_to_offsets = patch_data.node_to_offsets
  for (i, cellid) in enumerate(patch_data.patch_cell_ids)
    cell_vec_all_nodes = @view cell_vecs[cellid, :]
    offset = node_to_offsets[i]
    update_cell_local_dofs!(dm, cellid)
    cell_vec = cell_vec_all_nodes[offset]
    for i in axes(cell_vec, 1)
      patch_vec[dm.cell_dofs_loc[i]] += cell_vec[i]
    end
  end
end

function setup_patch_system!(A, RHS, linalg, dm_RT, dm_L²)
  A .= 0
  RHS .= 0
  free_patch_dofs_RT = dm_RT.free_patch_dofs_loc
  free_patch_dofs_L² = dm_L².free_patch_dofs_loc
  ndofs_RT = length(free_patch_dofs_RT)
  ndofs_L² = length(free_patch_dofs_L²)
  ndofs = ndofs_RT + ndofs_L²
  RHS[1:ndofs_RT] = @view linalg.RHS_RT[free_patch_dofs_RT]
  RHS[ndofs_RT+1:ndofs] = @view linalg.RHS_L²[free_patch_dofs_L²]
  A[1:ndofs_RT, 1:ndofs_RT] = @view linalg.M[free_patch_dofs_RT, free_patch_dofs_RT]
  A[1:ndofs_RT, ndofs_RT+1:ndofs] = @view linalg.B[free_patch_dofs_RT, free_patch_dofs_L²]
  A[ndofs_RT+1:ndofs, 1:ndofs_RT] =
    transpose(@view linalg.B[free_patch_dofs_RT, free_patch_dofs_L²])
end

function patch_hat_Ah2_contribution(patch_data, cell_Ah2_diags)
  s = 0.0
  node_to_offsets = patch_data.node_to_offsets
  for (i, cellid) in enumerate(patch_data.patch_cell_ids)
    s += cell_Ah2_diags[cellid, node_to_offsets[i]]
  end
  s
end

function count_free_dofs(dm::DOFManager)
  free_patch_dofs = dm.free_patch_dofs_loc
  length(free_patch_dofs)
end

function add_lagrange!(A, dm_RT, dm_L², Λ)
  ndofs_RT = count_free_dofs(dm_RT)
  ndofs_L² = count_free_dofs(dm_L²)
  ndofs = ndofs_RT + ndofs_L²
  A[ndofs+1, ndofs_RT+1:ndofs] = @view Λ[1:ndofs_L²]
  A[ndofs_RT+1:ndofs, ndofs+1] = transpose(@view Λ[1:ndofs_L²])
end

function scatter_to_global_σ!(σ_gl, dm_RT, σ_patch, n_free_dofs_RT)
  for i = 1:n_free_dofs_RT
    σ_gl[dm_RT.patch_dofs_gl[i]] += σ_patch[i]
  end
end

function solve_patch!(linalg, n_free_dofs)
  # Extract free dofs slices of the monolithic objects
  # Solve the patch problem in place
  A_free_dofs = @view linalg.A[1:n_free_dofs, 1:n_free_dofs]
  RHS_free_dofs = @view linalg.RHS[1:n_free_dofs]
  σ_free_dofs = @view linalg.σ_loc[1:n_free_dofs]
  ldiv!(σ_free_dofs, LU(LAPACK.getrf!(linalg.ws, A_free_dofs)...), RHS_free_dofs)
end

function build_DOFManagers(spaces)
  dm_RTs = [DOFManager(spaces.RT_space) for i = 1:Threads.nthreads()]
  dm_L²s = [DOFManager(spaces.L²_space) for i = 1:Threads.nthreads()]
  (dm_RTs, dm_L²s)
end

#The convention will be that
#  global := dof ordering for all the dofs
#  local  := dof ordering for the patch
#  an underscore _gl at the end of a variable indicates global dof index
#  and similarly for _loc always patch local
#
function build_equilibrated_flux(
  patches::AbstractVector{<:Patch},
  σ_gl,
  linalgs,
  cell_objects,
  RT_order,
  (dm_RTs, dm_L²s);
  nchunks=Threads.nthreads(), # TODO: Consider different chunk sizes?
  η=nothing
)
  #println("Loop on patches")
  co = cell_objects
  BLAS_nthreads = BLAS.get_num_threads()
  BLAS.set_num_threads(1)
  σ_gls = [zeros(size(σ_gl)) for i = 1:nchunks]
  #for patch in patches # Serial
  Threads.@threads for (patchid_range, ichunk) in chunks(1:length(patches), nchunks)
    for patchid in patchid_range
      #tid = Threads.threadid()
      linalg = linalgs[ichunk]
      dm_L² = dm_L²s[ichunk]
      dm_RT = dm_RTs[ichunk]
      σ_gl_chunk = σ_gls[ichunk]
      patch = patches[patchid]
      # Change the local numbering for the current patch
      update_patch_dofs!(dm_RT, patch.data)
      update_patch_dofs!(dm_L², patch.data)
      ## Scatter the cell based matrices in cell_objects to the reused matrices in
      ## linalg
      matrix_scatter!(linalg.M, co.cell_mass_mats, dm_RT, dm_RT, patch.data)
      matrix_scatter!(linalg.B, co.cell_mixed_mats, dm_L², dm_RT, patch.data)
      ## Idem for vectors
      vector_scatter!(linalg.RHS_RT, co.cell_RHS_RTs, dm_RT, patch.data)
      vector_scatter!(linalg.RHS_L², co.cell_RHS_L²s, dm_L², patch.data)
      single_vector_scatter!(linalg.Λ, co.cell_Λ_vecs, dm_L², patch.data)
      ## Now that scatter to local system is complete, remove fixed dofs
      remove_homogeneous_neumann_dofs!(dm_RT, patch.data, RT_order)
      ## Count the free dofs for this patch once the BCs are imposed
      n_free_dofs_RT = count_free_dofs(dm_RT)
      n_free_dofs_L² = count_free_dofs(dm_L²)
      n_free_dofs = n_free_dofs_RT + n_free_dofs_L²
      ## Use the sub-matrices and vectors generated from the scatters to build
      ## the monolithic objects
      setup_patch_system!(linalg.A, linalg.RHS, linalg, dm_RT, dm_L²)
      ## Handle the pure Neumann case
      if patch isa InteriorPatch
        add_lagrange!(linalg.A, dm_RT, dm_L², linalg.Λ)
        n_free_dofs += 1
      end
      solve_patch!(linalg, n_free_dofs)
      ## Scatter to the global FE object's free_values
      scatter_to_global_σ!(σ_gl_chunk, dm_RT, linalg.σ_loc, n_free_dofs_RT)

      if !isnothing(η)
        ## compute 0.5 * norm(σ_loc  - Ah)^2 contribution
        free_patch_dofs_RT = dm_RT.free_patch_dofs_loc
        σ_RT = linalg.σ_loc[1:n_free_dofs_RT]
        M_RT = @view linalg.M[free_patch_dofs_RT, free_patch_dofs_RT]
        RHS_RT = @view linalg.RHS_RT[free_patch_dofs_RT]
        η[patchid] = 0.5 * dot(M_RT* σ_RT, σ_RT)
        @assert η[patchid] >= -1e-12
        η[patchid] -= dot(σ_RT, RHS_RT)
        η[patchid] += 0.5 * patch_hat_Ah2_contribution(patch.data, co.cell_Ah2_diags)
        @assert η[patchid] >= -1e-12
      end
    end
  end
  σ_gl .+= sum(σ_gls)
  BLAS.set_num_threads(BLAS_nthreads)
  η, σ_gl
end