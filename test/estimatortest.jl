using EquilibratedFlux
using Gridap

let
  # Define some helper functions
  L2_inner_product(f, g, dx) = ∫(f ⋅ g) * dx

  L2_norm_squared(f, dx) = L2_inner_product(f, f, dx)

  function L2_norm_squared(f, model, order)
    degree = 2 * order + 2
    Ω = Triangulation(model)
    dx = Measure(Ω, degree)
    L2_norm_squared(f, dx)
  end

  function L²_projection_f(model, reffe, f, dx)
    V = TestFESpace(model, reffe; conformity = :L2)
    m(u, v) = ∫(u * v)*dx
    b(v) = ∫(v * f) * dx
    op_proj = AffineFEOperator(m, b, V, V)
    solve(op_proj)
  end

  dim = 2

  for order = 1:5
    for n = 5:5:40
      #n = 10 # Number of elements in x and y for square mesh
      domain = (0,1,0,1)
      partition = (n, n)
      model = CartesianDiscreteModel(domain, partition)
      # Change to triangles
      model = simplexify(model)
      𝓣ₕ = Triangulation(model)
      u(x) = sin(2*pi*x[1]) * sin(pi*x[2])
      f(x) = 5 * pi^2 * u(x) # = -Δu

      # Polynomial order
      degree = 2 * order + 2
      dx = Measure(𝓣ₕ, degree)
      reffe = ReferenceFE(lagrangian, Float64, order)
      V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
      U = TrialFESpace(V0, u)
      a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dx
      b(v) = ∫(v * f) * dx
      op = AffineFEOperator(a, b, U, V0)
      uh = solve(op)
      H1err = L2_norm_squared(∇(u - uh), dx) |> sum |> sqrt
      
      # perform different flux equilibrations
      η_loc, σ_eq = build_equilibrated_flux(-∇(uh), f, model, order; return_contributions=true)
      σ_ave = build_averaged_flux(-∇(uh), model)
      η_vertex = sqrt(sum(η_loc) * (dim+1))
      η_ave = L2_norm_squared(σ_ave + ∇(uh), dx) |> sum |> sqrt
      η_eq = L2_norm_squared(σ_eq + ∇(uh), dx) |> sum |> sqrt

      @test η_vertex > H1err
      @test η_eq > H1err
      eff = η_eq / H1err
      @test isapprox(eff, 1.0, atol=1e-2)
      nothing
    end
  end
end
