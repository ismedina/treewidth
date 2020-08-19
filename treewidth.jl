using LightGraphs
using StatsBase: sample

function ⊂(A,B)
    for a in A
        (a ∉ B) ? (return false) : nothing
    end
    return true
end

"""
    random_graph(Nn, Ne)

Random graph with Nn nodes and at Ne edges.
"""
function random_graph(Nn, Ne)
    Ne ≤ Nn*(Nn-1)÷2 || error("Number of edges must be smaller or equal than N(N-1)/2, with N the number of vertices")
    G = Graph(Nn)
    possible_edges = Tuple{Int, Int}[]
    for i in 1:Nn, j in 1:i-1
        push!(possible_edges, (i,j))
    end
    edges = sample(possible_edges, Ne, replace = false)
    for e in edges
        add_edge!(G, e[1], e[2])
    end
    return G
end

"""
    line_graph(G)

Line graph of G.
"""
function line_graph(G::Graph)
    LG = Graph()
    nodeinfo = []

    # add nodes
    for i in 1:nv(G)
        neigh = neighbors(G, i)
        for j in neigh
            # add the node if not present
            edge = (i < j)  ? (i, j) : (j, i)
            if edge ∉ nodeinfo
                add_vertex!(LG)
                push!(nodeinfo, edge)
            end
        end
    end

    # add edges
    for (i, e1) in enumerate(nodeinfo)
        for (j, e2) in enumerate(nodeinfo[i+1:end])
            common = e1 ∩ e2
            if length(common) > 0
                LightGraphs.add_edge!(LG, i, j + i)
            end
        end
    end
    return LG, nodeinfo
end

"""
    lacking_for_clique_neigh(G::Graph, i)

Finds the edges that are lacking to `G` for the neighborhood of vertex `i`
to be a clique. Returns this number of edges and the edges themselves.
"""
function lacking_for_clique_neigh(G::Graph, i::Int)
    neigh = neighbors(G,i)
    lacking = Tuple{Int,Int}[]
    for j in eachindex(neigh)
        for k in 1:j-1
            if ! has_edge(G, neigh[j], neigh[k])
                push!(lacking, (neigh[j], neigh[k]))
            end
        end
    end
    return length(lacking), lacking
end

function n_lacking_for_clique_neigh(G::Graph, i::Int)
    neigh = neighbors(G,i)
    #lacking = Tuple{Int,Int}[]
    n_lacking = 0
    for j in eachindex(neigh)
        for k in 1:j-1
            if ! has_edge(G, neigh[j], neigh[k])
                n_lacking += 1
            end
        end
    end
    return n_lacking
end

"""
    rem_vertex_fill!(G, i, lacking, ordering, vertex_label)

Fills `G` with the lacking edges for the neighborhood of `i` to be a clique,
removes vertex `i` of `G`, pushes it to `ordering` and updates `vertex_label`.
LightGraphs changes the order of the vertices when one is removed like this:

removed                  moved
   |                       |
  [v1    v2    v3    v4    v5]    -->    [v5    v2    v3    v4]
"""
function rem_vertex_fill!(G::Graph, i::Int, lacking::Vector{Tuple{Int,Int}},
                          ordering::Vector{Int}, vertex_label)
    for e in lacking
        LightGraphs.add_edge!(G, e[1], e[2])
    end
    push!(ordering, vertex_label[i])
    rem_vertex!(G, i)
    v = pop!(vertex_label)
    if i ≤ nv(G)
        vertex_label[i] = v
    end
end

function rem_vertex_no_fill!(G::Graph, i::Int,
                          ordering::Vector{Int}, vertex_label)
    push!(ordering, vertex_label[i])
    rem_vertex!(G, i)
    v = pop!(vertex_label)
    if i ≤ nv(G)
        vertex_label[i] = v
    end
end


"""
    min_fill_ordering(G)

Find an ordering of the vertices of `G` using the min-fill heuristic
(cfr. [Bodlaender, _Discovering Treewidth_](http://webdoc.sub.gwdg.de/ebook/serien/ah/UU-CS/2005-018.pdf))
"""
function min_fill_ordering(G::Graph)
    H = copy(G)
    ordering = Int[]
    #no_lacking = Tuple{Int,Int}[] # to pass when there is no edge lacking
    vertex_label = collect(1:nv(H))
    while nv(H) > 0
        success = false
        for i in nv(H):-1:1
            if degree(H, i) == 0
                rem_vertex_no_fill!(H, i, ordering, vertex_label)
                success = true
            end
        end
        for i in nv(H):-1:1
            if degree(H, i) == 1
                rem_vertex_no_fill!(H, i, ordering, vertex_label)
                success = true
            end
        end

        if ! success
            degrees = degree(H)
            J = sortperm(degrees)

            found_clique = false
            v = 0
            best_n_lacking = Inf
            #best_lacking = Tuple{Int,Int}[]
            for i in 1:nv(H)
                j = J[i]
                #n_lacking, lacking = lacking_for_clique_neigh(H, j)
                n_lacking = n_lacking_for_clique_neigh(H, j)
                if n_lacking == 0
                    rem_vertex_no_fill!(H, j, ordering, vertex_label)
                    found_clique = true
                    break
                elseif n_lacking < best_n_lacking
                    v = j
                    best_n_lacking = n_lacking
                    #best_lacking = lacking
                end
            end
            # if running until here, remove v
            if ! found_clique
                best_n_lacking, best_lacking = lacking_for_clique_neigh(H, v)
                rem_vertex_fill!(H, v, best_lacking, ordering, vertex_label)
                #print("h")
            end
        end
    end
    ordering
end

"""
    triangulation(G::Graph, ordering)

Chordal completion of `G` following order `ordering`.
For each vertex add edges between its higher numbered neighbors.
"""
function triangulation(G::Graph, ordering)
    H = copy(G)
    for (i, v) in enumerate(ordering)
        neigh = neighbors(H, v)
        high_neigh = neigh ∩ ordering[i+1:end] # get higher numbered neighbors
        for (j, i1) in enumerate(high_neigh)
            for i2 in high_neigh[1:j-1]
                LightGraphs.add_edge!(H, i1, i2)
            end
        end
    end
    return H
end

"""
    tree_decomposition(G)

Finds a tree decomposition of G using the min-fill heuristic. Returns the
width of the decomposition, the tree and the bags.
"""
function tree_decomposition(G::Graph)
    ordering = min_fill_ordering(G)
    H = triangulation(G, ordering)
    up_neighs = [neighbors(H, ordering[i]) ∩ ordering[i+1:end] for i in 1:nv(H)]
    # first bag is the first node within a maximum clique following the ordering
    clique_neigh_idx = findfirst(length.(up_neighs) .== nv(H)-1:-1:0)
    first_bag = up_neighs[clique_neigh_idx] ∪ ordering[clique_neigh_idx]
    decomp = Graph(1)
    bags = [first_bag]
    tw = length(first_bag) - 1

    for i in clique_neigh_idx-1:-1:1
        neigh = up_neighs[i]

        # find a bag all neighbors are in
        old_bag_idx = 0
        for (j,bag) in enumerate(bags)
            if neigh ⊂ bag
                old_bag_idx = j
                break
            end
        end
        (old_bag_idx == 0) && (old_bag_idx = 1) # no old_bag was found: just connect to the first_bag

        # create new node
        add_vertex!(decomp)
        new_bag = [neigh; ordering[i]]
        push!(bags, new_bag)

        # update treewidth
        tw = max(tw, length(new_bag)-1)

        # add edge to decomposition
        LightGraphs.add_edge!(decomp, old_bag_idx, nv(decomp))

     end
     return tw, decomp, bags
end

"""
    is_tree_decomposition(G, tree, bags)

Checks if `(tree, bags)` forms a tree decomposition of `G`.
"""
function is_tree_decomposition(G, tree, bags)
    # 1. check union property
    if sort(∪(bags...)) != collect(1:nv(G))
        @warn ("Union of bags is not equal to union of vertices")
        return false
    end

    # 2. check edge property
    for edge in edges(G)
        edge_found = false
        e = (edge.src, edge.dst)
        for B in bags
            edge_found = edge_found | (e ⊂ B)
        end
        if ! edge_found
            @warn("Edge $e not found in any bag")
            return false
        end
    end

    # 3. check subtree property
    subgraphs = [Int[] for i in 1:nv(G)] # subgraph[i] are the bags that contain `i`
    for (i,b) in enumerate(bags)
        for j in b
            push!(subgraphs[j], i)
        end
    end

    for (v, s) in enumerate(subgraphs)
        if length(s) > 0
            subtree, _ = induced_subgraph(tree, s)
            if ! is_connected(subtree)
                @warn("Subgraph for vertex $v not connected")
                return false
            end
        end
    end
    return true
end
