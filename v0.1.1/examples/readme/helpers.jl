using Gridap
using LaTeXStrings
using GridapMakie, CairoMakie

# Define some helper functions
L2_inner_product(f, g, dx) = ∫(f ⋅ g) * dx

L2_norm_squared(f, dx) = L2_inner_product(f, f, dx)

function L2_norm_squared(f, model, order)
  degree = 2 * order + 2
  trian = Triangulation(model)
  dx = Measure(trian, degree)
  L2_norm_squared(f, dx)
end

function L²_projection(model, reffe, f, dx)
  V = TestFESpace(model, reffe; conformity = :L2)
  m(u, v) = ∫(u * v)*dx
  b(v) = ∫(v * f) * dx
  op_proj = AffineFEOperator(m, b, V, V)
  solve(op_proj)
end

function plot_error_and_estimator(trian, ηave_arr, ηeq_arr, H1err_arr)
    max_val = maximum([ηave_arr..., ηeq_arr..., H1err_arr...])
    min_val = minimum([ηave_arr..., ηeq_arr..., H1err_arr...])
    ηave_vis = CellField(ηave_arr, trian)
    ηeq_vis = CellField(ηeq_arr, trian)
    H1err_vis = CellField(H1err_arr, trian)
    fig = Figure(resolution = (700, 600))
    ga = fig[1, 1] = GridLayout()
    axerr = Axis(ga[1, 1], xlabel = L"x", ylabel = L"y", title = L"$H_0^1$ seminorm error")
    axeq  = Axis(ga[2, 1], xlabel = L"x", ylabel = L"y", title = L"$$Equilibrated flux estimator")
    axave = Axis(ga[2, 2], xlabel = L"x", ylabel = L"y", title = L"$$Averaged flux estimator")
    plot_error = plot!(axerr, trian, H1err_vis, colorrange=(min_val,max_val), colormap=:viridis)
    plot_eq    = plot!(axeq,  trian, ηeq_vis,    colorrange=(min_val,max_val), colormap=:viridis)
    plot_aver  = plot!(axave, trian, ηave_vis,    colorrange=(min_val,max_val), colormap=:viridis)
    Colorbar(fig[1,2], limits=(min_val, max_val), colormap=:viridis)
    fig
end
#save("comparison.png", fig) #src
#
function plot_divergence_misfit(trian, eq_div_vis, ave_div_vis)
  fig = Figure(resolution = (1600, 800))
  ga = fig[1, 1] = GridLayout()
  axdiveq = Axis(ga[1, 1], xlabel = L"x", ylabel = L"y", title = L"$$ Divergence misfit equilibrated flux")
  axdivave= Axis(ga[1, 3], xlabel = L"x", ylabel = L"y", title = L"$$ Divergence misfit averaged flux")
  plot_div_eq = plot!(axdiveq, trian, eq_div_vis,    colormap=:viridis)
  Colorbar(ga[1,2], plot_div_eq)
  plot_div_ave = plot!(axdivave, trian, ave_div_vis, colormap=:viridis)
  Colorbar(ga[1,4], plot_div_ave)
  fig
end
