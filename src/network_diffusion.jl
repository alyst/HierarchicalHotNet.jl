"""
    stepmatrix(g::Union{AbstractSimpleWeightedGraph, AbstractMatrix};
               inedge_weights::Union{AbstractVector, Nothing} = nothing,
               normalize_weights::Bool = true) -> AbstractMatrix

Construct the step probabilities matrix for the random walk on graph `g`.

Returns matrix ``S``, where ``s_{i,j}`` is the probability to travel along
the edge ``(v_j, v_i) ∈ g``.

* `inedge_weights::AbstractVector`: (optional) an array of factors
   for each vertex to scale the weights of incoming edges
* `normalize_weights=true`: normalize the weights of the columns in
   the final output, so that it sums to 1
"""
function stepmatrix(adjmtx::AbstractMatrix{<:Number};
                    inedge_weights::Union{AbstractVector, Nothing} = nothing,
                    normalize_weights::Bool=true)
    check_square(adjmtx, "Adjacency matrix")
    adj_mtx = !isnothing(inedge_weights) ?
        Diagonal(inedge_weights) * adjmtx :
        copy(adjmtx)
    if normalize_weights
        node_scales = dropdims(sum(adj_mtx, dims=1), dims=1)
        @inbounds for (i, w) in enumerate(node_scales)
            node_scales[i] = ifelse(w == 0, 1.0, inv(w))
        end
        rmul!(adj_mtx, Diagonal(node_scales))
    end
    return adj_mtx
end

stepmatrix(g::AbstractSimpleWeightedGraph; kwargs...) =
    stepmatrix(Graphs.adjacency_matrix(g, dir=:in); kwargs...)

"""
    random_walk_matrix(g::AbstractSimpleWeightedGraph,
                       restart_probability::Number=0.1;
                       stepmtx_kwargs...) -> Matrix{Float64}

Calculate the stationary distribution of vertex transition probabilities
for a random walk with restart on graph `g`.

Returns the matrix ``W``, where ``w_{i, j}`` is the ``v_j → v_i`` transition probability
along the edges of the graph `g`.
"""
random_walk_matrix(g::AbstractSimpleWeightedGraph,
                   restart_probability::Number=0.1;
                   stepmtx_kwargs...) =
    random_walk_matrix(Matrix(stepmatrix(g; stepmtx_kwargs...)),
                       restart_probability)

function random_walk_matrix(adjmtx::AbstractMatrix,
                            restart_probability::Number)
    (0 <= restart_probability <= 1) ||
        throw(DomainError(restart_probability, "restart_probability should be in [0, 1] range"))
    check_square(adjmtx, "Adjacency matrix")
    return restart_probability * inv(convert(Matrix,
            Diagonal(I, size(adjmtx, 1)) -
            (1-restart_probability) * adjmtx))
end

"""
    similarity_matrix(g::AbstractSimpleWeightedGraph,
                      node_weights::AbstractVector;
                      restart_probability::Union{Number, Nothing} = nothing,
                      stepmtx_kwargs...) -> Matrix{Float64}

Calculate the matrix of node transition probabilities for the random walk with restart
taking into account the node restarting probabilities `node_weights`.

Returns the matrix ``Ŵ``, where ``ŵ_{i, j}`` is the probability of ``v_j → v_i`` transition
considering all possible starting vertices of graph `g`.

### Parameters
  * `node_weights::AbstractVector`: vector of node restart probabilities
  * `restart_probability`: random walk restart probability
  * [`stepmatrix()`](@ref HierarchicalHotNet.stepmatrix) keyword arguments

### See also
[`HierarchicalHotNet.random_walk_matrix`](@ref)
"""
similarity_matrix(g::Union{AbstractMatrix, AbstractSimpleWeightedGraph},
    node_weights::AbstractVector;
    restart_probability::Union{Number, Nothing} = nothing,
    stepmtx_kwargs...) =
        random_walk_matrix(g, restart_probability; stepmtx_kwargs...) *
        Diagonal(node_weights)

"""
    stabilized_stepmatrix(g::AbstractSimpleWeightedGraph,
                          node_weights::AbstractVector;
                          restart_probability::Union{Number, Nothing} = nothing,
                          stepmtx_kwargs...) -> Matrix{Float64}

Calculate the matrix for the node transition probabilities at each iteration
of the random walk with restart taking into account the probabilities of visiting the nodes.

Returns matrix ``Ŝ``. If ``w`` is the node starting probabilities (`node_weights`),
``S`` is the step matrix, and ``W`` is the corresponding transition probabilities matrix
for a random walk with restart, then
```math
    Ŝ = (1-r) S × \\mathrm{diag}(W × w) + r \\mathrm{diag}(w).
```
### Parameters
  * `node_weights::AbstractVector`: vector of node restart probabilities
  * `restart_probability`: random walk restart probability
  * `node_weights_diffused`: precalculated node visiting weights for a random walk with restart
    (in case the random matrix is already calculated before)
  * [`stepmatrix()`](@ref HierarchicalHotNet.stepmatrix) keyword arguments
"""
stabilized_stepmatrix(
            g::Union{AbstractMatrix, AbstractSimpleWeightedGraph},
            node_weights::AbstractVector;
            restart_probability::Union{Number, Nothing} = nothing,
            node_weights_diffused::Union{Nothing, AbstractVector} = nothing,
            stepmtx_kwargs...) =
    stabilized_stepmatrix(stepmatrix(g; stepmtx_kwargs...), node_weights, restart_probability,
                          node_weights_diffused = node_weights_diffused)

function stabilized_stepmatrix(
            stepmtx::AbstractMatrix,
            node_weights::AbstractVector,
            restart_probability::Number,
            node_weights_diffused::Union{Nothing, AbstractVector} = nothing)
    diffused_node_weights = node_weights_diffused !== nothing ? node_weights_diffused :
        random_walk_matrix(stepmtx, restart_probability) * node_weights
    diffsum = sum(diffused_node_weights)
    k = diffsum != 0 ? sum(node_weights)/diffsum : 0.0 # k to maintain the weights sum after diffusion
    return stepmtx * (k * (1 - restart_probability) * Diagonal(diffused_node_weights)) +
           convert(typeof(stepmtx), restart_probability * Diagonal(node_weights))
end

"""
    neighborhood_weights(adjmtx::AbstractMatrix, g::AbstractGraph) -> Vector{Float64}

Given the `adjmtx` adjacency matrix for the graph ``g'`` and the graph `g`, calculate
the sum of the weights of the ``g'`` outgoing edges that are also present in `g`.
Returns the vector with the weights sum for each vertex.

The graphs have to have the same number of nodes.
"""
function neighborhood_weights(adjmtx::AbstractMatrix, g::AbstractGraph)
    check_square(adjmtx, "Adjacency matrix")
    n = size(adjmtx, 1)
    (n == nv(g)) || throw(DimensionMismatch("Adjacency matrix columns ($n) do not match the number of graph vertices ($(nv(g)))"))
    return [sum(j -> j != i ? weights[j] : 0.0, neighborhood(g, i, 1))
            for (i, weights) in enumerate(eachcol(adjmtx))]
end
