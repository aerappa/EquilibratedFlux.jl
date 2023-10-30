#=
function solve_poisson(model, u, f, order)
degree = 2 * order + 2
Ω = Triangulation(model)
dx = Measure(Ω, degree)
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe; conformity = :H1, dirichlet_tags = "boundary")
U = TrialFESpace(V0, u)
a(u, v) = ∫(∇(v) ⊙ ∇(u)) * dx
b(v) = ∫(v * f) * dx
op = AffineFEOperator(a, b, U, V0)
uh = solve(op)
writevtk(Ω, "Lshaped_soln", cellfields = ["uh" => uh, "∇uh" => ∇(uh)])
uh
end

function compute_estimator(σ, 𝐀ₕ, model, order)
degree = 2 * order + 1
Ω = Triangulation(model)
dx = Measure(Ω, degree)
#η = ∫((σ + 𝐀ₕ) ⋅ (σ + 𝐀ₕ)) * dx
η = L2_norm_int(σ + 𝐀ₕ, dx)
writevtk(Ω, "Estimator", cellfields = ["η" => σ + 𝐀ₕ])
get_array(η)
end

function error_estimate(model, u, uh, 𝐀ₕ, σ, f, order)
degree = 2 * order
Ω = Triangulation(model)
dx = Measure(Ω, degree)
#H1err = sum(∫(∇(u - uh) ⋅ ∇(u - uh)) * dx)
H1err = L2_norm(∇(u - uh), dx)
div_check = L2_norm(∇ ⋅ σ - f, dx)
est = L2_norm(σ + 𝐀ₕ, dx)
@show div_check
@show H1err
@show est
@show est / H1err
(H1err, est)
end
=#

