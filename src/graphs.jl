function import_digraph(edges_df::DataFrame;
                        src_col::Symbol, dest_col::Symbol, weight_col::Symbol)
    edges_df = copy(edges_df, copycols=false)
    nodes = sort!(unique(vcat(edges_df[!, src_col], edges_df[!, dest_col])))
    edges_df[!, :__src_node__] = searchsortedfirst.(Ref(nodes), edges_df[!, src_col])
    edges_df[!, :__dest_node__] = searchsortedfirst.(Ref(nodes), edges_df[!, dest_col])
    res = SimpleWeightedDiGraph{Int, Float64}(length(nodes))
    for r in eachrow(edges_df)
        add_edge!(res, r.__src_node__, r.__dest_node__, r[weight_col])
    end
    return res, nodes
end

function hhotnet_example_graph()
    g = SimpleWeightedGraph(25)
    gene_scores_a = Vector{Float64}()
    gene_scores_b = Vector{Float64}()

    # Add a clique
    clique_nodes = 1:8
    for i in clique_nodes
        push!(gene_scores_a, 5)
        push!(gene_scores_b, 0.1)
        # pos[letters[i]] = (math.cos(2*math.pi*i/8)-2, math.sin(2*math.pi*i/8))
        for j in (i+1):last(clique_nodes)
            add_edge!(g, i, j)
        end
    end

    # Add an approximate clique
    approx_clique_nodes = 9:16
    nlinks = 2
    for i in first(approx_clique_nodes):(last(approx_clique_nodes)+nlinks)
        isrc = mod(i, approx_clique_nodes)
        # pos[letters[i]] = (math.cos(2*math.pi*i/8)+2, math.sin(2*math.pi*i/8))
        push!(gene_scores_a, 10)
        for j in (i+1):(i+nlinks)
            add_edge!(g, isrc, mod(j, approx_clique_nodes))
        end
    end
    append!(gene_scores_b, [0.7, 0.4, 0.5, 0.6, 0.6, 0.6, 0.5, 0.4])

    # Add a cycle
    cycle_nodes = 17:24
    for i in cycle_nodes
        push!(gene_scores_a, 3)
        add_edge!(g, i, mod(i+1, cycle_nodes))
        # pos[letters[i]] = (math.cos(2*math.pi*i/8), math.sin(2*math.pi*i/8)-2)
    end
    append!(gene_scores_b, [0.1, 0.2, 0.1, 0.2, 0.1, 0.2, 0.1, 0.2])

    # Add a point connecting the clique, approximate clique, and cycle, add an isolated point.
    hub_node = 25
    push!(gene_scores_a, 5)
    push!(gene_scores_b, 1.0)
    add_edge!(g, first(clique_nodes), hub_node)
    add_edge!(g, first(approx_clique_nodes), hub_node)
    add_edge!(g, first(cycle_nodes), hub_node)
    #pos['y'] = (0.0, 0.0)

    return g
end

"""
Generate the graph from the example in Figure 1 of Tarjan (1983).
"""
function tarjan1983_example_graph()
    g = SimpleWeightedDiGraph(7)
    for (u, v, w) in [('a', 'b', 10), ('b', 'a', 12), ('b', 'c', 30), ('d', 'c',  6), ('d', 'e', 16),
                      ('e', 'd', 13), ('e', 'f',  8), ('f', 'a', 26), ('a', 'g', 15), ('g', 'b', 35),
                      ('c', 'g', 45), ('g', 'c', 22), ('d', 'g', 14), ('g', 'e', 50), ('f', 'g', 20)]
        add_edge!(g, u-'a'+1, v-'a'+1, w)
    end
    return g
end
