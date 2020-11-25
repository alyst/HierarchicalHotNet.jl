"""
    stepmatrix(g::Union{AbstractSimpleWeightedGraph, AbstractMatrix};
               inedge_weights::Union{AbstractVector, Nothing} = nothing,
               normalize_weights::Bool = true) -> AbstractMatrix

Construct the step matrix for the random walk.

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
    stepmatrix(LightGraphs.weights(g); kwargs...)

"""
    random_walk_matrix(g::AbstractSimpleWeightedGraph,
                       restart_probability::Number=0.1;
                       stepmtx_kwargs...)

Find the matrix that transforms given initial node probabilities into
stationary distribution of visiting probabilities of a random walk with restart.
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
                      stepmtx_kwargs...)

Calculate the matrix for the random walk with restart.

### Parameters
* `node_weights::AbstractVector`: vector of node restart probabilities
* `restart_probability`: random walk restart probability
* [`HierarchicalHotNet.stepmatrix`](@ref) keyword arguments
"""
similarity_matrix(g::Union{AbstractMatrix, AbstractSimpleWeightedGraph},
    node_weights::AbstractVector;
    restart_probability::Union{Number, Nothing} = nothing,
    stepmtx_kwargs...) =
        random_walk_matrix(g, restart_probability; stepmtx_kwargs...) *
        Diagonal(node_weights)

function neighborhood_weights(adjmtx::AbstractMatrix, g::AbstractGraph)
    check_square(adjmtx, "Adjacency matrix")
    n = size(adjmtx, 1)
    (n == nv(g)) || throw(DimensionMismatch("Adjacency matrix columns ($n) do not match the number of graph vertices ($(nv(g)))"))
    return [sum(j -> j != i ? weights[j] : 0.0, neighborhood(g, i, 1))
            for (i, weights) in enumerate(eachcol(adjmtx))]
end
