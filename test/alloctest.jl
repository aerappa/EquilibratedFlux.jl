using EquilibratedFlux
using Gridap
using Gridap.Geometry
const EF = EquilibratedFlux

#=
Zero-allocation regression test for the performance-critical per-patch
solve loop (matrix_scatter!/vector_scatter!/solve_patch!/etc., invoked via
the internal, 6-argument-vector `build_equilibrated_flux` method), run in
serial (nchunks=1, bypassing Threads.@threads entirely). This mirrors the
exact setup the top-level, user-facing
`build_equilibrated_flux(𝐀ₕ, f, model, RT_order)` performs internally,
including the type annotations on the filtered patch vectors (without
which the loop variable is abstractly typed and every call boxes,
producing spurious allocations that also disguise real ones).

Two meshes are used since a single mesh can only be pure-Dirichlet or
pure-Neumann (mixed boundaries are not yet supported, see NeumannLift.jl):
one exercises DirichletPatch + InteriorPatch, the other InteriorPatch +
NeumannPatch.
=#
function _check_zero_alloc(patch_vec, σ_gl, linalgs, cell_objects, RT_order, dms, σ_gls)
  isempty(patch_vec) && return
  # warm up (JIT compile) before measuring
  build_equilibrated_flux(patch_vec, σ_gl, linalgs, cell_objects, RT_order, dms, σ_gls; nchunks = 1)
  fill!(σ_gl, 0)
  allocated = @allocated build_equilibrated_flux(
    patch_vec, σ_gl, linalgs, cell_objects, RT_order, dms, σ_gls; nchunks = 1)
  @test allocated == 0
end

let
  u_exact(x) = sin(2 * pi * x[1]) * sin(pi * x[2])
  f(x) = 5 * pi^2 * u_exact(x)

  n = 10
  order = 2
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition) |> simplexify
  degree = 2 * order + 2
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  U = TrialFESpace(V0, u_exact)
  a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
  b(v) = ∫(v * f) * dΩ
  op = AffineFEOperator(a, b, U, V0)
  uh = solve(op)
  Ah = -∇(uh)

  RT_order = order
  patches, metadata = EF.create_patches(model, RT_order)
  spaces = EF.build_global_spaces(model, RT_order)
  cell_objects = EF.build_all_cellwise_objects(Ah, f, 1.0, spaces, model, RT_order, nothing)
  linalgs = [EF.instantiate_linalg(RT_order, 2, metadata) for i = 1:1]
  dms = EF.build_DOFManagers(spaces)
  σ_gl = zero(spaces.RT_space)
  σ_gls = [zeros(size(σ_gl.free_values)) for i = 1:1]

  diri_patches::Vector{EF.DirichletPatch{Int32}} =
    filter(patch -> patch isa EF.DirichletPatch, patches)
  int_patches::Vector{EF.InteriorPatch{Int32}} =
    filter(patch -> patch isa EF.InteriorPatch, patches)
  @test !isempty(diri_patches)
  @test !isempty(int_patches)

  _check_zero_alloc(diri_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms, σ_gls)
  _check_zero_alloc(int_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms, σ_gls)
end

let
  gradu(x) = VectorValue(3 * x[1]^2, 3 * x[2]^2)
  f(x) = -(2.0 - 4.0)
  Aex(x) = -gradu(x)

  n = 10
  order = 2
  domain = (0, 1, 0, 1)
  partition = (n, n)
  model = CartesianDiscreteModel(domain, partition) |> simplexify
  labels = get_face_labeling(model)
  add_tag_from_tags!(labels, "neumann", ["boundary"])
  degree = 2 * order + 2
  Ω = Triangulation(model)
  dΩ = Measure(Ω, degree)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V0 = TestFESpace(model, reffe; conformity = :H1)
  V0cf = FESpaceWithConstantFixed(V0, true)
  U = TrialFESpace(V0cf)
  ΓN = BoundaryTriangulation(model, tags = ["neumann"])
  dΓN = Measure(ΓN, degree)
  nΓN = get_normal_vector(ΓN)
  σN = -(CellField(gradu, ΓN)) ⋅ nΓN
  a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dΩ
  b(v) = ∫(v * f) * dΩ - ∫(v * σN) * dΓN
  op = AffineFEOperator(a, b, U, V0cf)
  uh = solve(op)
  Ah = -∇(uh)

  RT_order = order
  patches, metadata = EF.create_patches(model, RT_order; neumann_tags = ["neumann"])
  spaces = EF.build_global_spaces(model, RT_order)
  patches = EF.compute_neumann_lift(patches, model, spaces, RT_order, ["neumann"], Aex)
  cell_objects = EF.build_all_cellwise_objects(Ah, f, 1.0, spaces, model, RT_order, nothing)
  linalgs = [EF.instantiate_linalg(RT_order, 2, metadata) for i = 1:1]
  dms = EF.build_DOFManagers(spaces)
  σ_gl = zero(spaces.RT_space)
  σ_gls = [zeros(size(σ_gl.free_values)) for i = 1:1]

  int_patches::Vector{EF.InteriorPatch{Int32}} =
    filter(patch -> patch isa EF.InteriorPatch, patches)
  neu_patches::Vector{EF.NeumannPatch{Int32}} =
    filter(patch -> patch isa EF.NeumannPatch, patches)
  @test !isempty(int_patches)
  @test !isempty(neu_patches)

  _check_zero_alloc(int_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms, σ_gls)
  _check_zero_alloc(neu_patches, σ_gl.free_values, linalgs, cell_objects, RT_order, dms, σ_gls)
end
