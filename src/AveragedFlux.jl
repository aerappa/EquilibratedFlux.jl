using Gridap.Geometry
using Gridap
using FillArrays

average(xs) = sum(xs) / length(xs)

midpoint(xs) = average(xs)

function get_positive_normals(trian)
    n = get_normal_vector(trian)
    edge_endpoints = get_cell_points(trian)
    n.plus(edge_endpoints)
    # plus == ⁺ is the positively oriented normal: outward on boundary and from
    # low to high GID in the interior
   normals = evaluate(n.plus, edge_endpoints)
   lazy_map(midpoint, normals) |> parallel_smart_collect
end

function get_cell_field_at_midpoints(𝓣ₕ, 𝐀ₕ)
  reffe_mid = VectorValue(0.5, 0.5)
  mid_cellpoint = CellPoint(Fill(reffe_mid, num_cells(𝓣ₕ)), 𝓣ₕ, ReferenceDomain())
  evaluate(𝐀ₕ, mid_cellpoint) |> parallel_smart_collect
end

function get_edge_lengths(edge_to_nodes, node_coords)
	num_edges = length(edge_to_nodes)
  edge_to_nodes_cache = array_cache(edge_to_nodes)
  edge_lengths = zeros(num_edges)
  for edge_ind in 1:num_edges 
    node = getindex!(edge_to_nodes_cache, edge_to_nodes, edge_ind)
    edge_lengths[edge_ind] = norm(node_coords[node[1]] - node_coords[node[2]])
  end
  edge_lengths
end

function scatter_to_σ!(σ_arr, edge_inds, edge_to_cells, edge_lengths, normals, 𝐀ₕ_at_midpoints)
  normals_cache = array_cache(normals)
  edge_to_cells_cache = array_cache(edge_to_cells)
  #𝐀ₕ_average = zero(𝐀ₕ_at_midpoints[1])
  for (i, edge_ind) in enumerate(edge_inds)
    # Cells that share this vertex
    cell_inds = getindex!(edge_to_cells_cache, edge_to_cells, edge_ind)
    #n_F = n_all[edge_ind]
    normal = getindex!(normals_cache, normals, i)
    #𝐀ₕ_average += getindex!(𝐀ₕ_at_midpoints_cache, 𝐀ₕ_at_midpoints, 1)
    #for cell_ind in cell_inds
    #  𝐀ₕ_average += getindex!(𝐀ₕ_at_midpoints_cache, 𝐀ₕ_at_midpoints, cell_ind)
    #end
    #𝐀ₕ_average /= length(𝐀ₕ_average)
    𝐀ₕ_average  = average(𝐀ₕ_at_midpoints[cell_inds])
    σ_arr[edge_ind] = edge_lengths[edge_ind]*(𝐀ₕ_average ⋅ normal)  # -|F|{∇uh} ⋅ n_F
  end
end


"""
    build_averaged_flux(𝐀ₕ, model)

TODO: relevant docstring
"""
function build_averaged_flux(𝐀ₕ, model)
  𝓣ₕ = Triangulation(model)
	𝓢ₕ = SkeletonTriangulation(model)
	𝓑ₕ = BoundaryTriangulation(model)
	grid = get_grid(𝓣ₕ)
	topo = Geometry.GridTopology(grid)
	edge_to_cells = Geometry.get_faces(topo, 1, 2)
	edge_to_nodes = Geometry.get_faces(topo, 1, 0)
	num_edges = length(edge_to_nodes)
  node_coords = Geometry.get_node_coordinates(model)
  edge_lengths = get_edge_lengths(edge_to_nodes, node_coords)
	n_skel = get_positive_normals(𝓢ₕ)
	n_bdry = get_positive_normals(𝓑ₕ)
  bdry_edge_inds = findall(map(cells -> length(cells)==1, edge_to_cells))
	skel_edge_inds = setdiff(1:length(edge_to_cells), bdry_edge_inds)
	reffe_RT₀ = ReferenceFE(raviart_thomas, Float64, 0)
	RT₀ = FESpace(model, reffe_RT₀)
  𝓣ₕ = Triangulation(model)
  𝐀ₕ_at_midpoints  = get_cell_field_at_midpoints(𝓣ₕ, 𝐀ₕ)
  σ_arr = zeros(num_edges)
  # Bdry
  scatter_to_σ!(σ_arr, bdry_edge_inds, edge_to_cells, edge_lengths, n_bdry, 𝐀ₕ_at_midpoints)
  # Skeleton
  scatter_to_σ!(σ_arr, skel_edge_inds, edge_to_cells, edge_lengths, n_skel, 𝐀ₕ_at_midpoints)
 	FEFunction(RT₀, σ_arr)
end
