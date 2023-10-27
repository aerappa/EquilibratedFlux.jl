using Gridap

function build_triangle_model(n)
  domain = (0, 1, 0, 1)
  partition = (n, n)
  CartesianDiscreteModel(domain, partition) |> simplexify
end

function build_global_spaces(model, k)
  reffeRT = ReferenceFE(raviart_thomas, Float64, k)
  reffeP = ReferenceFE(lagrangian, Float64, k)
  RT_space = FESpace(model, reffeRT, conformity = :HDiv)
  L²_space = FESpace(model, reffeP, conformity = :L2)
  (; RT_space, L²_space)
end
