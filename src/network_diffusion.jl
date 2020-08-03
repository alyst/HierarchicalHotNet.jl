"""
Construct the step matrix for the random walk.
"""
function stepmatrix(adjmtx::AbstractMatrix{<:Number};
                     normalize_weights::Bool=true)
    check_square(adjmtx, "Adjacency matrix")
    adj_mtx = copy(adjmtx)
    if normalize_weights
        node_degs = sum(adj_mtx, dims=1)
        node_degs[node_degs .== 0] .= 1
        adj_mtx ./= node_degs
    end
    return adj_mtx
end

stepmatrix(g::AbstractSimpleWeightedGraph; kwargs...) =
    stepmatrix(LightGraphs.weights(g); kwargs...)

"""
Find the matrix that transforms given initial node probabilities into
stationary distribution of visiting probabilities of a random walk with restart.
"""
random_walk_matrix(g::AbstractSimpleWeightedGraph,
                   restart_probability::Number=0.1;
                   normalize_weights::Bool=true) =
    random_walk_matrix(Matrix(stepmatrix(g, normalize_weights=normalize_weights)),
                       restart_probability)

function random_walk_matrix(adjmtx::AbstractMatrix,
                            restart_probability::Number = 0.1)
    check_square(adjmtx, "Adjacency matrix")
    return restart_probability * inv(
            diagm(0 => fill(1.0, size(adjmtx, 1))) -
            (1-restart_probability) * adjmtx)
end

function similarity_matrix(g::AbstractSimpleWeightedGraph,
                           node_weights::AbstractVector;
                           restart_probability::Number = 0.1)
    return random_walk_matrix(g, restart_probability) *
           Diagonal(node_weights)
end

function neighborhood_weights(adjmtx::AbstractMatrix, g::AbstractGraph)
    check_square(adjmtx, "Adjacency matrix")
    n = size(adjmtx, 1)
    (n == nv(g)) || throw(DimensionMismatch("Adjacency matrix columns ($n) do not match the number of graph vertices ($(nv(g)))"))
    return [sum(j -> j != i ? weights[j] : 0.0, neighborhood(g, i, 1))
            for (i, weights) in enumerate(eachcol(adjmtx))]
end
