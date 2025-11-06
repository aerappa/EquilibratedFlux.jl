using ChunkSplitters

function parallel_smart_collect(lazy_array; nchunks=Threads.nthreads())
  len = length(lazy_array)
  ntids = Threads.nthreads()
  non_lazy = Array{eltype(lazy_array)}(undef, len)
  caches = [array_cache(lazy_array) for _ = 1:ntids]
  #Threads.@threads  for i = 1:n
  # Threads.@threads for (id_range, ichunk) in chunks(1:n; n=nchunks)
  Threads.@threads for (ichunk, id_range) in enumerate(index_chunks(1:len; n=nchunks))
    for i in id_range
      cache = caches[ichunk]
      non_lazy[i] = copy(getindex!(cache, lazy_array, i))
    end
  end
  non_lazy
end

function smart_collect(lazy_array)
  n = length(lazy_array)
  non_lazy = Array{eltype(lazy_array)}(undef, n)
  cache = array_cache(lazy_array)
  for i = 1:n
    non_lazy[i] = copy(getindex!(cache, lazy_array, i))
    #getindex!(cache, lazy_array, i)
  end
  non_lazy
end
