using Gridap.Geometry

abstract type Patch{T} end

struct PatchData{T}
  node_to_offsets::Vector{T}
  patch_cell_ids::Vector{T}
  bdry_edge_ids::Vector{T}
  all_edge_ids::Vector{T}
end

struct DirichletPatch{T} <: Patch{T}
  data::PatchData{T}
end

struct InteriorPatch{T} <: Patch{T}
  data::PatchData{T}
end

# TODO: Implement
struct NeumannPatch{T} <: Patch{T}
  data::PatchData{T}
  bdry_data::Vector{Float64}
end

function _is_interior_object(i, d, labels)
  entity_idx_i = labels.d_to_dface_to_entity[d][i]
  labels.tag_to_name[entity_idx_i] == "interior"
end

_is_boundary_node_labels(i, labels) = !_is_interior_object(i, 1, labels)

_is_boundary_edge_labels(e, labels) = !_is_interior_object(e, 2, labels)

_is_boundary_edge_patch_i_e2n(e, i, edge_to_node) = i ∉ edge_to_node[e]

#=
This function constructs the map from a given node to the indices
in the cell_to_node array which it appears. For example with
cell_to_node = [[1, 2, 4], [2, 4, 5], [2, 3, 5], [3, 5, 6],
                [4, 5, 7], [5, 7, 8], [5, 6, 8], [6, 8, 9]]
The first few entries of node_to_offsets are
node_to_offsets = [[1], [2, 1, 1], [2, 1], [3, 2, 1],..
The idea is that eventually we will want to slice into a cell_RHS arrays
with certain cells e.g.
patch_cell_RHS_RTs = cell_RHS_RTs[patch_cells]
and then next extract the current patch via
patch_RT_cell_RHS_RTs = patch_cell_RHS_RTs[node_to_offsets[patchid]]
=#
function _get_nodes_to_offsets(model)
  topo = get_grid_topology(model)
  cell_to_node = Geometry.get_faces(topo, 2, 0)
  node_to_cell = Geometry.get_faces(topo, 0, 2)
  cell_to_node_cache = array_cache(cell_to_node)
  node_to_cell_cache = array_cache(node_to_cell)
  num_nodes = length(node_to_cell)
  node_to_offsets = [[] for i = 1:num_nodes]
  for nodeid = 1:num_nodes
    cellids = getindex!(node_to_cell_cache, node_to_cell, nodeid)
    for cellid in cellids
      nodeids = getindex!(cell_to_node_cache, cell_to_node, cellid)
      node_offset = findfirst(i -> i == nodeid, nodeids)
      push!(node_to_offsets[nodeid], node_offset)
    end
  end
  node_to_offsets
end

function _get_node_to_offsets(
  nodeid,
  cell_to_node,
  node_to_cell,
  cell_to_node_cache,
  node_to_cell_cache,
)
  node_to_offsets = typeof(nodeid)[]
  cellids = getindex!(node_to_cell_cache, node_to_cell, nodeid)
  for cellid in cellids
    nodeids = getindex!(cell_to_node_cache, cell_to_node, cellid)
    node_offset = findfirst(i -> i == nodeid, nodeids)
    push!(node_to_offsets, node_offset)
  end
  node_to_offsets
end

function _get_patch_edge_ids(patch_cell_ids, cell_to_edge, cell_to_edge_cache)
  patch_edge_ids = getindex!(cell_to_edge_cache, cell_to_edge, patch_cell_ids)
  unique!(patch_edge_ids.data)
end

function _get_boundary_edges(patch_edge_ids, patch_cell_ids, edge_to_cell, edge_to_cell_cache, is_boundary_patch)
  bdry_edge_ids = eltype(patch_edge_ids)[]
  for e in patch_edge_ids
    # internal boundary edges have an adjacent cell outside the patch
    cells = getindex!(edge_to_cell_cache,edge_to_cell,e)
    is_internal = !all([cell in patch_cell_ids for cell in cells])
    # for internal vertices, global boundary edges are also constrained
    if is_internal || (!is_boundary_patch && length(cells)==1)
      push!(bdry_edge_ids, e)
    end
  end
  bdry_edge_ids
end

function test_edge_dof_consistency(model, RT_order, RT_space)
  topo = get_grid_topology(model)
  dofs_per_edge = RT_order + 1
  node_to_cell = Geometry.get_faces(topo, 0, 2)
  cell_to_edge = Geometry.get_faces(topo, 2, 1)
  node_to_cell_cache = array_cache(node_to_cell)
  cell_to_edge_cache = array_cache(cell_to_edge)
  cell_dofs_RT = get_cell_dof_ids(RT_space)
  num_nodes = length(node_to_cell)
  for i = 1:num_nodes
    patch_cell_ids = getindex!(node_to_cell_cache, node_to_cell, i)
    patch_edge_ids = _get_patch_edge_ids(patch_cell_ids, cell_to_edge, cell_to_edge_cache)
    patch_edge_dofs = _get_edge_dofs(patch_edge_ids, dofs_per_edge)
    sort!(patch_edge_dofs)
    RT_dofs = unique!(cell_dofs_RT[patch_cell_ids].data)
    sort!(RT_dofs)
    nedofs = length(patch_edge_dofs)
    @assert RT_dofs[1:nedofs] == patch_edge_dofs
  end
end

function create_patches(model, RT_order)
  labels = get_face_labeling(model)
  #is_boundary_edge(e) = _is_boundary_edge_labels(e, labels)
  is_boundary_node(i) = _is_boundary_node_labels(i, labels)
  dofs_per_edge = RT_order + 1
  #dofs_per_cell = RT_order * (RT_order + 1)
  topo = get_grid_topology(model)
  node_to_cell = Geometry.get_faces(topo, 0, 2)
  cell_to_node = Geometry.get_faces(topo, 2, 0)
  cell_to_edge = Geometry.get_faces(topo, 2, 1)
  edge_to_cell = Geometry.get_faces(topo, 1, 2)
  #nodes_to_offsets = _get_nodes_to_offsets(model)
  # Type of the table, Int32 it seems
  T = eltype(cell_to_edge[1])
  patches = Patch[]
  # Need to declare caches first to prevent reinstatiation
  node_to_cell_cache = array_cache(node_to_cell)
  cell_to_edge_cache = array_cache(cell_to_edge)
  cell_to_node_cache = array_cache(cell_to_node)
  edge_to_cell_cache = array_cache(edge_to_cell)
  num_nodes = length(node_to_cell)
  max_num_patch_cells = 0
  for i::T = 1:num_nodes
    is_boundary = is_boundary_node(i)
    patch_cell_ids = getindex!(node_to_cell_cache, node_to_cell, i)
    all_edge_ids = _get_patch_edge_ids(patch_cell_ids, cell_to_edge, cell_to_edge_cache)
    max_num_patch_cells = max(max_num_patch_cells, length(patch_cell_ids))
    bdry_edge_ids = _get_boundary_edges(all_edge_ids, patch_cell_ids, edge_to_cell, edge_to_cell_cache, is_boundary)
    #all_dofs = _get_edge_dofs(patch_edge_ids, dofs_per_edge)
    node_to_offsets = _get_node_to_offsets(
      i,
      cell_to_node,
      node_to_cell,
      cell_to_node_cache,
      node_to_cell_cache,
    )
    # Copy needed because otherwise modifying the underlying array
    data = PatchData(node_to_offsets, copy(patch_cell_ids), bdry_edge_ids, all_edge_ids)
    if is_boundary
      patch = DirichletPatch(data)
      #filter!(!is_boundary_edge, bdry_edge_ids)
    else
      patch = InteriorPatch(data)
    end
    push!(patches, patch)
  end
  #metadata = (; nodes_to_offsets, max_num_patch_cells)
  # Can probably delete this is max_num_patch_cells can be known a priori
  metadata = (; max_num_patch_cells)
  patches, metadata
end
