"""
Construct the walk matrix for the random walk.
"""
function walk_matrix(g::AbstractSimpleWeightedGraph;
                     normalize_weights::Bool=true)
    adj_mtx = copy(weights(g))
    if normalize_weights
        node_degs = sum(adj_mtx, dims=1)
        node_degs[node_degs .== 0] .= 1
        adj_mtx ./= node_degs
    end
    return adj_mtx
end

"""
Find the matrix that transforms given initial node probabilities into
stationary distribution of visiting probabilities of a random walk with restart.
"""
random_walk_matrix(g::AbstractSimpleWeightedGraph;
                   restart_probability::Number=0.1,
                   normalize_weights::Bool=true) =
    random_walk_matrix(Matrix(walk_matrix(g, normalize_weights=normalize_weights)),
                       restart_probability=restart_probability)

function random_walk_matrix(adjmtx::AbstractMatrix;
                            restart_probability::Number = 0.1)
    size(adjmtx, 1) == size(adjmtx, 2) || throw(DimensionMismatch("Square matrix required"))
    return restart_probability * inv(
            diagm(0 => fill(1.0, size(adjmtx, 1))) -
            (1-restart_probability) * adjmtx)
end

function similarity_matrix(g::AbstractSimpleWeightedGraph,
                           node_weights::AbstractVector;
                           restart_probability::Number = 0.1)
    return random_walk_matrix(g, restart_probability=restart_probability) *
           Diagonal(node_weights)
end

function neighborhood_weights(adjmtx::AbstractMatrix, g::AbstractGraph)
    size(adjmtx, 1) == size(adjmtx, 2) || throw(DimensionMismatch("Square matrix expected, $(size(adjmtx)) given"))
    n = size(adjmtx, 1)
    n == nv(g) || throw(DimensionMismatch("Adjacency matrix cols ($n) don't match the graph vertices $(nv(g))"))
    return [sum(j -> j != i ? weights[j] : 0.0, neighborhood(g, i, 1))
            for (i, weights) in enumerate(eachcol(adjmtx))]
end
