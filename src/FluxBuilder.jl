using LinearAlgebra
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
    measure = nothing, weight = 1.0, neumann_tags = String[], neumann_data = nothing)
  build_equilibrated_flux(𝐀ₕ, f, model.model, RT_order, measure = measure, weight = weight,
    neumann_tags = neumann_tags, neumann_data = neumann_data)
end


"""
    build_equilibrated_flux(𝐀ₕ, f, model, RT_order; measure = nothing, weight= 1.0,
                             neumann_tags = String[], neumann_data = nothing)

TODO: relevant docstring

`neumann_tags` marks the part ΓN of the boundary (as model face labeling
tags) where an inhomogeneous Neumann condition `-∇u⋅n = σN` holds; the rest
of the boundary is assumed to be (possibly inhomogeneous) Dirichlet. When
given, `neumann_data` must be a function `x -> VectorValue` giving the
physical flux `σN(x)⋅n(x)` (e.g. `-∇u_exact` for a manufactured solution),
matching the `𝐀ₕ = -∇(uh)` convention used elsewhere in this package.
Vertices lying on the interface between the Dirichlet and Neumann parts of
the boundary are not yet supported.
"""
function build_equilibrated_flux(𝐀ₕ, f, model, RT_order; measure = nothing, weight= 1.0,
    neumann_tags = String[], neumann_data = nothing)
  topo = get_grid_topology(model)
  @assert all(p->p==TRI,get_polytopes(topo))
  @assert isempty(neumann_tags) || !isnothing(neumann_data) "neumann_data must be given when neumann_tags is non-empty"
  patches, metadata = create_patches(model, RT_order; neumann_tags = neumann_tags)
  spaces = build_global_spaces(model, RT_order)
  patches = compute_neumann_lift(patches, model, spaces, RT_order, neumann_tags, neumann_data)
  cell_objects = build_all_cellwise_objects(𝐀ₕ, f, weight, spaces, model, RT_order, measure)
  linalgs = [instantiate_linalg(RT_order, 2, metadata) for i = 1:Threads.nthreads()]
  dms = build_DOFManagers(spaces)
  σ_gl = zero(spaces.RT_space)
  diri_patches::Vector{DirichletPatch{Int32}} =
    filter(patch -> patch isa DirichletPatch, patches)
  int_patches::Vector{InteriorPatch{Int32}} =
    filter(patch -> patch isa InteriorPatch, patches)
  neu_patches::Vector{NeumannPatch{Int32}} =
    filter(patch -> patch isa NeumannPatch, patches)
  build_equilibrated_flux(diri_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms)
  build_equilibrated_flux(int_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms)
  build_equilibrated_flux(neu_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms)
  σ_gl
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
  nchunks=Threads.nthreads() # TODO: Consider different chunk sizes?
)
  #println("Loop on patches")
  co = cell_objects
  BLAS_nthreads = BLAS.get_num_threads()
  BLAS.set_num_threads(1)
  σ_gls = [zeros(size(σ_gl)) for i = 1:nchunks]
  #for patch in patches # Serial
  # Threads.@threads for (patchid_range, ichunk) in chunks(1:length(patches), nchunks)
  Threads.@threads for (ichunk, patchid_range) in enumerate(index_chunks(1:length(patches); n=nchunks))
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
      ## Neumann patches: their own directly-incident ΓN edges carry a
      ## prescribed (essential), generally nonzero, RT normal trace. Locate
      ## their local positions *before* any dof is removed from dm_RT (both
      ## this and the homogeneous removal below shrink dm_RT.patch_dofs_gl,
      ## which would otherwise invalidate position-based lookups), lift
      ## their contribution into the RHS, then remove them afterwards.
      if patch isa NeumannPatch
        neu_local = find_local_dof_positions(dm_RT, patch.dof_ids)
        g_fixed = patch.dof_values
        linalg.RHS_RT .-= linalg.M[:, neu_local] * g_fixed
        linalg.RHS_L² .-= transpose(linalg.B[neu_local, :]) * g_fixed
      end
      ## Now that scatter to local system is complete, remove fixed dofs
      remove_homogeneous_neumann_dofs!(dm_RT, patch.data, RT_order)
      if patch isa NeumannPatch
        remove_prescribed_dofs!(dm_RT, patch.dof_ids, neu_local)
      end
      ## Count the free dofs for this patch once the BCs are imposed
      n_free_dofs_RT = count_free_dofs(dm_RT)
      n_free_dofs_L² = count_free_dofs(dm_L²)
      n_free_dofs = n_free_dofs_RT + n_free_dofs_L²
      ## Use the sub-matrices and vectors generated from the scatters to build
      ## the monolithic objects
      setup_patch_system!(linalg.A, linalg.RHS, linalg, dm_RT, dm_L²)
      ## Handle the pure Neumann case (fully-interior or fully-ΓN patches)
      if patch isa InteriorPatch || patch isa NeumannPatch
        add_lagrange!(linalg.A, dm_RT, dm_L², linalg.Λ)
        n_free_dofs += 1
      end
      solve_patch!(linalg, n_free_dofs)
      ## Scatter to the global FE object's free_values
      scatter_to_global_σ!(σ_gl_chunk, dm_RT, linalg.σ_loc, n_free_dofs_RT)
      ## Directly assign the prescribed (Neumann-lifted) dof contributions;
      ## each shared ΓN edge accumulates the ψₐ-weighted share from both of
      ## its endpoint patches, recovering the full σN moment on that edge.
      if patch isa NeumannPatch
        for (gdof, val) in zip(patch.dof_ids, patch.dof_values)
          σ_gl_chunk[gdof] += val
        end
      end
    end
  end
  σ_gl .+= sum(σ_gls)
  #patch::DirichletPatch = patches[1]
  #@show to
  BLAS.set_num_threads(BLAS_nthreads)
  σ_gl
end
