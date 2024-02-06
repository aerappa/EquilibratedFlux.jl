using EquilibratedFlux 
using LinearAlgebra 
using ArgParse
using Gridap

L2_norm(f, dx) = L2_norm_int(f, dx) |> sum |> sqrt

L2_norm_int(f, dx) = L2_inner_product(f, f, dx)

L2_inner_product(f, g, dx) = ∫(f ⋅ g) * dx


function solve_poisson(model, u, f, order)
  degree = 2 * order + 2
  Ω = Triangulation(model)
  dx = Measure(Ω, degree)
  reffe = ReferenceFE(lagrangian, Float64, order)
  V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
  dofs = num_free_dofs(V0)
  U = TrialFESpace(V0, u)
  a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dx
  b(v) = ∫(v * f) * dx
  A = assemble_matrix(a, U, V0)
  luA = lu(A)
  rhs = assemble_vector(b, V0)
  uh_dofs = luA \ rhs
  uh = FEFunction(V0, uh_dofs)
  uh, dofs
end

function compute_estimator(σ, 𝐀ₕ, model, degree)
  order = 2 * degree + 2
  Ω = Triangulation(model)
  dx = Measure(Ω, order)
  L2_norm(σ - 𝐀ₕ, dx)
end

function compute_error(u, uh, model, degree)
  order = 2 * degree + 2
  Ω = Triangulation(model)
  dx = Measure(Ω, order)
  L2_norm(∇(u - uh), dx)
end

function run_timing_test(ns, soln_degree, RT_degree, foldername)
  times_poisson = Float64[]
  times_cell_loop = Float64[]
  times_bdry_loops = Float64[]
  times_dm = Float64[]
  times_int_loops = Float64[]
  times_flux_total = Float64[]
  times_flux_ave = Float64[]
  all_dofs = Float64[]
  #JLD2.save(joinpath(foldername, "results.jld2"), Base.@locals)
  #return
  let
    # Solving -Δu = f
    #u(x) = 30*x[1] * (x[1] - 1) * x[2] * (x[2] - 1)
    #∇u(x) = 30*VectorValue((2 * x[1] - 1) * x[2] * (x[2] - 1), (2 * x[2] - 1) * x[1] * (x[1] - 1))
    #f(x) = 30*(-2 * (x[1] * x[1] + x[2] * x[2]) + 2 * (x[1] + x[2]))
    u(x) = sin(2*pi*x[1])*sin(2*pi*x[2])
    ∇u(x) = 2*pi*VectorValue(cos(2*pi*x[1])*sin(2*pi*x[2]),sin(2*pi*x[1])*cos(2*pi*x[2]))
    f(x) = 8*pi^2*u(x)
    for n in ns
      model = CartesianDiscreteModel((0,1,0,1), (n,n)) |> simplexify
      println("Assemble and solve")
      GC.gc()
      @show Threads.nthreads()
      BLAS.set_num_threads(Threads.nthreads())
      time_poisson = @elapsed begin
        println("time poisson")
        @time uh, dofs = solve_poisson(model, u, f, soln_degree)
      end
      @show time_poisson
      @show dofs
      push!(all_dofs, dofs)
      push!(times_poisson, time_poisson)
      ∇uh = ∇(uh)
      println("Flux assemble")
      GC.gc()
      time_flux = @elapsed begin
        #σ, time_cell_loop, time_bdry_loop, time_int_loop, time_dm = build_flux(∇uh, f, model, RT_degree, weight =1.0)
        σ= build_equilibrated_flux(-∇uh, f, model, RT_degree)
        η = compute_estimator(σ, -∇uh, model, soln_degree)
        e = compute_error(u, uh, model, soln_degree)
        @show η / e
      end
      @show time_flux
      #push!(times_dm, time_dm)
      #push!(times_flux_total, time_flux)
      #push!(times_cell_loop, time_cell_loop)
      #push!(times_bdry_loops, time_bdry_loop)
      #push!(times_int_loops, time_int_loop)
      #time_ave = @elapsed σ_ave = build_averaged_flux(model, ∇uh)
      #push!(times_flux_ave, time_ave)
      #@time (H1err, est) = error_estimate(model, u, uh, ∇uh, σ, f, soln_degree)
      println("-----------------------------")
      #(H1err_ave, est_ave) = error_estimate(model, u, uh, ∇uh, σ_ave, f, soln_degree)
      #println("-----------------------------")
    end
  end
  num_threads = Threads.nthreads()
  #JLD2.save(joinpath(foldername, "results.jld2"), Base.@locals)
  return Base.@locals
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output-dir"
            help = "Directory to save to"
            default = "."
        "--degree", "-p"
            help = "Polynomial degree"
            arg_type = Int
            default = 1
        "--dofs", "-d"
        help = "Total desired max dofs (approximate)"
            arg_type = Int
            default = 5e6
    end
    return parse_args(s)
end

# For getting the n to get m dofs for a given polynomial degree k
function quadratic(a, b, c)
	discr = b^2 - 4*a*c
	discr >= 0 ?   ( (-b + sqrt(discr))/(2a), (-b - sqrt(discr))/(2a) ) : error("Only complex roots")
end

function get_n_for_m_dofs_degree_k(m, k)
  result = quadratic((2*binomial(k-1, 2) + 3*(k-1) + 1), -2*k, 1 - m)
  if result isa Number
    return Int64(ceil(result))
  else
    return Int64(ceil(result[1]))
  end
end

function main_CLI()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
    k = parsed_args["degree"]
    desired_max_dofs = parsed_args["dofs"]
    dofs = [Int64(floor(desired_max_dofs / 4^k)) for k = 0:7]
    ns  = [get_n_for_m_dofs_degree_k(dof,k) for dof in dofs]
    @show ns
    run_timing_test(ns, k, k, parsed_args["output-dir"])
end

function dummy_main()
    k = 3
    desired_max_dofs = 100000
    contraction_factor = 4
    min_dof = 0
    dofs = [Int64(floor(desired_max_dofs / contraction_factor^i)) for i = 0:min_dof]
    ns  = [get_n_for_m_dofs_degree_k(dof,k) for dof in dofs]
    @show ns
    run_timing_test(ns, k, k, "test")
end

#main_CLI()
dummy_main()
