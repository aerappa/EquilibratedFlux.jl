using GLMakie
using Gridap
using GridapMakie
using Gridap.CellData
using Gridap.Adaptivity

let
  domain = (0,1,0,1)
  fields = CellField[]
  for i = 1:5
    n = 2^i
    # This works
    #n = 8
    parition = (n, n)
    model = CartesianDiscreteModel(domain, parition) |> simplexify
    #model = model.model
    Ω = Triangulation(model)
    field = CellField(rand(num_cells(model)), Ω)
    push!(fields, field)
  end

  idx = Observable(1)
  field_plot = lift(idx) do idx
    fields[idx]
  end

  fig, ax, plt = plot(field_plot)
  framerate = 5
  idxs = 1:length(fields)
  record(fig, "animation.gif", idxs; framerate=framerate, compression=0) do this_idx
    idx[] = this_idx
  end
end
