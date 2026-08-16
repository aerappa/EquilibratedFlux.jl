#using TimerOutputs

struct DOFManager{T, M <: Matrix{T}, V <: Vector{T}} # TODO: get rid of these types
  # All the cell dofs stored as a matrix
  all_cell_dofs_gl::M
  # The current patch dofs in the global ennumeration
  patch_dofs_gl::V
  # The current free dofs in patch local ennumeration for slicing
  # into the patch local objects
  free_patch_dofs_loc::V
  # The current cell's dofs in the patch local ennumeration
  cell_dofs_loc::V
end

function DOFManager(space)
  all_cell_dofs_gl = make_matrix_from_lazy(get_cell_dof_ids(space))
  patch_dofs_RT_gl = sizehint!(eltype(all_cell_dofs_gl)[], 256)
  push!(patch_dofs_RT_gl, 0)
  cell_dofs_loc = sizehint!(eltype(all_cell_dofs_gl)[], 256)
  push!(cell_dofs_loc, 0)
  free_patch_dofs_loc = sizehint!(eltype(all_cell_dofs_gl)[], 256)
  push!(free_patch_dofs_loc, 0)
  DOFManager(all_cell_dofs_gl, patch_dofs_RT_gl, free_patch_dofs_loc, cell_dofs_loc)
end

DOFManager(dm::DOFManager) = DOFManager(
  copy(dm.all_cell_dofs_gl),
  copy(dm.patch_dofs_gl),
  copy(dm.free_patch_dofs_loc),
  copy(dm.cell_dofs_loc),
)

function update_cell_local_dofs!(dm::DOFManager, cellid)
  cur_cell_dofs_gl = @view dm.all_cell_dofs_gl[cellid, :]
  empty!(dm.cell_dofs_loc)
  for id in cur_cell_dofs_gl
    new_id = findfirst(n -> n == id, dm.patch_dofs_gl)
    new_id isa Nothing && error("Cannot update cell local dofs!")
    push!(dm.cell_dofs_loc, new_id)
  end
end

function make_matrix_from_lazy(lazy)
  rows = length(lazy)
  cols = length(lazy[1])
  matrix = zeros(eltype(eltype(lazy)), rows, cols)
  lazy_cache = array_cache(lazy)
  for i = 1:rows
    matrix[i, :] = getindex!(lazy_cache, lazy, i)
  end
  matrix
end

#function _get_edge_dofs(edge_ids, dofs_per_edge)
#  edge_dofs = eltype(edge_ids)[]
#  for edge_id in edge_ids
#    # Enumerate all dofs on the given edge
#    for i = 0:(dofs_per_edge-1)
#      edge_dof = (edge_id * dofs_per_edge) - i
#      push!(edge_dofs, edge_dof)
#    end
#  end
#  edge_dofs
#end


function remove_homogeneous_neumann_dofs!(dm, patch_data, RT_order)
  dofs_per_edge = RT_order + 1
  # TODO: for the moment, there are two loops, first to remove the local
  # edge dofs, and then to remove the global edge dofs. The problem is that
  # the global edge dofs cannot be modied in the first loop because their
  # positions are being checked to get the local_edge_dofs.
  patch_dofs_gl = dm.patch_dofs_gl
  free_patch_dofs_loc = dm.free_patch_dofs_loc
  bdry_edge_ids = patch_data.bdry_edge_ids
  for edge_id in bdry_edge_ids
  #for edge_id in patch_data.bdry_edge_ids
    # Enumerate all dofs on the given edge
    for i = 0:(dofs_per_edge-1)
      bdry_edge_dof = (edge_id * dofs_per_edge) - i
      local_edge_dof = findfirst(n -> n == bdry_edge_dof, patch_dofs_gl)
      local_edge_dof isa Nothing && error("local_edge_dof cannot be computed!")
      local_edge_dof_idx = findfirst(n -> n == local_edge_dof, free_patch_dofs_loc)
      local_edge_dof_idx isa Nothing && error("local_edge_dof cannot be computed!")
      deleteat!(free_patch_dofs_loc,  local_edge_dof_idx)
      #@timeit to "filter" filter!(n -> n ≠ local_edge_dof, free_patch_dofs_loc)
    end
  end
  for edge_id in bdry_edge_ids
  #for edge_id in patch_data.bdry_edge_ids
    # Enumerate all dofs on the given edge
    for i = 0:(dofs_per_edge-1)
      bdry_edge_dof = (edge_id * dofs_per_edge) - i
      edge_dof_idx = findfirst(n -> n == bdry_edge_dof, patch_dofs_gl)
      edge_dof_idx isa Nothing && error("local_edge_dof cannot be computed!")
      deleteat!(patch_dofs_gl, edge_dof_idx)
      #filter!(n -> n ≠ edge_dof, patch_dofs_gl)
    end
  end
  #@show length(dm.patch_dofs_gl) 
  #@show length(dm.free_patch_dofs_loc)
  #@assert length(dm.patch_dofs_gl) == length(dm.free_patch_dofs_loc)
end

#=
Finds, for each of the given GLOBAL RT dof ids, its position in the
*current* (assumed not-yet-shrunk-by-any-other-removal) dm.patch_dofs_gl.
This is a pure lookup: it is the caller's responsibility to call this
*before* any dof removal (homogeneous or prescribed) touches patch_dofs_gl,
since it relies on patch_dofs_gl still reflecting the original local
numbering used by e.g. linalg.M/linalg.B/linalg.RHS_RT. The returned
indices remain valid references into those fixed-size arrays regardless of
what happens to patch_dofs_gl afterwards.
=#
function find_local_dof_positions(dm::DOFManager, dof_ids)
  patch_dofs_gl = dm.patch_dofs_gl
  local_indices = Vector{Int}(undef, length(dof_ids))
  for k in eachindex(dof_ids)
    loc = findfirst(n -> n == dof_ids[k], patch_dofs_gl)
    loc isa Nothing && error("prescribed dof cannot be located in the current patch!")
    local_indices[k] = loc
  end
  local_indices
end

#=
Removes the given GLOBAL RT dof ids from the patch's free system, exactly
as remove_homogeneous_neumann_dofs! does for the (always-zero) internal
patch boundary. Unlike that function, the dofs being removed here are
generally *nonzero* (prescribed from Neumann data), so their original local
positions (`local_indices`, from `find_local_dof_positions` called on a
still-pristine patch_dofs_gl, i.e. *before* this function or
remove_homogeneous_neumann_dofs! ran) must be used by the caller beforehand
to lift the prescribed values into the patch RHS. This function only needs
to be safe to call *after* other removals have already shrunk
patch_dofs_gl: it locates targets purely by value, never by position, so it
does not depend on patch_dofs_gl being pristine.
=#
function remove_prescribed_dofs!(dm::DOFManager, dof_ids, local_indices)
  patch_dofs_gl = dm.patch_dofs_gl
  free_patch_dofs_loc = dm.free_patch_dofs_loc
  for loc in local_indices
    idx = findfirst(n -> n == loc, free_patch_dofs_loc)
    idx isa Nothing && error("prescribed dof already removed from the free set!")
    deleteat!(free_patch_dofs_loc, idx)
  end
  for gdof in dof_ids
    idx = findfirst(n -> n == gdof, patch_dofs_gl)
    idx isa Nothing && error("prescribed dof cannot be located while shrinking patch_dofs_gl!")
    deleteat!(patch_dofs_gl, idx)
  end
  nothing
end

function update_global_patch_dofs!(dm::DOFManager, patch_data)
  empty!(dm.patch_dofs_gl)
  patch_cell_ids = patch_data.patch_cell_ids
  #RT_dofs_all = unique(foldl(union, patch_cell_dof_RT))
  for i ∈ patch_cell_ids
    cell_dofs = @view dm.all_cell_dofs_gl[i, :]
    for j ∈ cell_dofs
      if j ∉ dm.patch_dofs_gl
        push!(dm.patch_dofs_gl, j)
      end
    end
  end
end

function update_free_patch_dofs!(dm, patch_data)
  empty!(dm.free_patch_dofs_loc)
  # TODO: possibly for RT only handle removing edge DOFs
  for i = 1:length(dm.patch_dofs_gl)
    push!(dm.free_patch_dofs_loc, i)
  end
end

function update_patch_dofs!(dm::DOFManager, patch_data)
  update_global_patch_dofs!(dm, patch_data)
  update_free_patch_dofs!(dm, patch_data)
  nothing
end
