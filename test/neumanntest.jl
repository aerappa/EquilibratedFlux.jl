using EquilibratedFlux
using Gridap
using Gridap.FESpaces

let
  # Pure Neumann problem (ΓN = ∂Ω) with a manufactured cubic solution. Its
  # gradient (hence σN) and its Laplacian (hence f) are low-degree
  # polynomials, exactly representable at the chosen polynomial order, so
  # the equilibrium property and the Neumann-datum match can be
  # checked to near machine precision, exactly as for the pure-Dirichlet
  # case in dirichlettest.jl. u itself is degree 3, not exactly
  # representable at order 2, so uh retains a genuine (mesh-dependent)
  # discretization error and the effectivity index is a meaningful,
  # well-conditioned ratio.
  u(x) = x[1]^3 + x[2]^3
  gradu(x) = VectorValue(3 * x[1]^2, 3 * x[2]^2)
  f(x) = -(6 * x[1] + 6 * x[2]) # = -Δu

  order = 2
  for n = 5:5:20
    domain = (0, 1, 0, 1)
    partition = (n, n)
    model = CartesianDiscreteModel(domain, partition) |> simplexify
    labels = get_face_labeling(model)
    add_tag_from_tags!(labels, "neumann", ["boundary"])

    degree = 2 * order + 2
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)
    reffe = ReferenceFE(lagrangian, Float64, order)

    # Pure Neumann primal problem is only unique up to an additive constant;
    # pin one dof to remove the null space (the H1-seminorm error and the
    # flux reconstruction are both insensitive to this additive constant).
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

    H1err = sqrt(sum(∫(∇(u - uh) ⋅ ∇(u - uh)) * dΩ))

    neumann_data(x) = -gradu(x)
    σ_eq = build_equilibrated_flux(-∇(uh), f, model, order;
      neumann_tags = ["neumann"], neumann_data = neumann_data)

    η_eq = sqrt(sum(∫((σ_eq + ∇(uh)) ⋅ (σ_eq + ∇(uh))) * dΩ))
    @test η_eq > H1err
    eff = η_eq / H1err
    @test isapprox(eff, 1.0, atol = 5e-2)

    # equilibrium property: ∇⋅σ_eq should exactly reproduce f
    reffeL² = ReferenceFE(lagrangian, Float64, order)
    VL² = FESpace(model, reffeL²; conformity = :L2)
    mL²(u, v) = ∫(u * v) * dΩ
    bL²(v) = ∫(v * f) * dΩ
    f_proj = solve(AffineFEOperator(mL², bL², VL², VL²))
    div_misfit = sqrt(sum(∫((∇ ⋅ σ_eq - f_proj) ⋅ (∇ ⋅ σ_eq - f_proj)) * dΩ))
    @test isapprox(div_misfit, 0.0, atol = 1e-8)

    # Neumann datum: for a polynomial σN, σ_eq should reproduce it exactly
    # up to the RT_order polynomial moments on each ΓN face.
    σ_eq_dot_n = σ_eq ⋅ nΓN
    neumann_mismatch = sqrt(sum(∫((σ_eq_dot_n - σN) * (σ_eq_dot_n - σN)) * dΓN))
    @test isapprox(neumann_mismatch, 0.0, atol = 1e-8)
  end
end
