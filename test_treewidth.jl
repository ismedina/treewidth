using Test
using TestSetExtensions
using LightGraphs
include("treewidth.jl")

# for extensive check of the random tests increase the following parameter
samples_per_test = 1

@testset ExtendedTestSet "subset function" begin
    A = [1,2,3]
    B = [1,2,3,4]

    @test A ⊂ B
    @test !(B ⊂ A)

    # works for all iterable objects
    @test (1,2) ⊂ A
    @test A ⊂ (1,2,3,4)
    @test 1 ⊂ A
    @test 1:3 ⊂ A
    @test "ba" ⊂ "abc"

    # doesn't work for non-iterable objects
    @test_throws MethodError (complete_graph(3) ⊂ complete_graph(3))
end

@testset ExtendedTestSet "random_graph" begin
    Nn = 10
    Ne = 20
    G = random_graph(Nn, Ne)
    @test nv(G) == Nn
    @test ne(G) == Ne

    Ne = Nn*(Nn-1)÷2 + 1
    @test_throws ErrorException("Number of edges must be smaller or equal than N(N-1)/2, with N the number of vertices") random_graph(Nn, Ne)
end

@testset ExtendedTestSet "line_graph" begin
    G = random_regular_graph(4,2) # a square
    LG, _ = line_graph(G)
    @test (nv(LG) == 4) & prod(degree(LG) .== 2) # in this case G and LG are isomorphic
end

# Tree decomposition subroutines
@testset ExtendedTestSet "lacking_for_clique_neigh" begin
    G = complete_graph(5)
    @test (0,[]) == lacking_for_clique_neigh(G,1)
    rem_edge!(G,2,3)
    @test (1,[(2,3)]) == lacking_for_clique_neigh(G,1)
end

@testset ExtendedTestSet "rem_vertex_fill!" begin
    ## `rem_vertex_fill!`: eliminates a vertexs in the graph and fills
    # their neighborhood
    G = complete_graph(5)
    rem_edge!(G,2,3)
    ordering = Int[]
    vertex_label = [1, 2, 3, 4, 5]
    rem_vertex_fill!(G, 1, [(2,3)], ordering, vertex_label)

    @test G == complete_graph(4)
    @test ordering == [1]
    @test vertex_label == [5, 2, 3, 4]
end

@testset ExtendedTestSet "triangulation" begin
    ## `triangulation`: from a graph and an ordering get a chordal completion
    # (i. e. for every n-cycle (n>3) there is a 3-cycle that traverses
    # only three of the original nodes)
    # since enumerating the cycles on a graph grows faster than
    # exponentially in the number of edges, the parameters must be kept small
    for i in 1:samples_per_test
        Nn = rand(4:8)
        Ne = rand(Nn:(Nn*(Nn-1))÷2)
        G = random_graph(Nn, Ne)
        ordering = min_fill_ordering(G)
        H = triangulation(G, ordering)

        dg = DiGraph(H) # cycle enumerating only works on digraphs
        cycles = simplecycles(dg)
        length_3_cycles = cycles[length.(cycles) .== 3]
        success = true
        for cycle in cycles[length.(cycles) .> 3]
            success = success & any(map(x -> x ⊂ cycle, length_3_cycles))
        end
        @test success
    end
end

@testset ExtendedTestSet "tree decomposition" begin
    ## tree decomposition of a complete graph
    # A complete graph has itself as its only chordal completion.
    # Since the min_fill_ordering works by constructing a chordal completion,
    # and for complete graphs there is only one, the approximated treewidth is exact.
    for n in [10, 25, 50]
        G = complete_graph(n)
        tw, _ = tree_decomposition(G)
        @test tw == n-1
    end

    # TODO: find more examples in which the heuristic is supposed to yield the correct result
end

@testset ExtendedTestSet "is_tree_decomposition" begin
    # We consider the following graph and tree decomposition
    #    5
    #  /   \            tree                   [2, 3, 4]
    # 4 ——— 3            -->                    /     \
    # |     |       decomposition       [2, 4, 1]     [3, 4, 5]
    # 1 ——— 2


    # `ìs_tree_decomposition`: validates tree decompositions
    G = Graph(5)
    add_edge!(G, 1, 2)
    add_edge!(G, 1, 4)
    add_edge!(G, 2, 3)
    add_edge!(G, 4, 3)
    add_edge!(G, 5, 3)
    add_edge!(G, 4, 5)

    tree = Graph(3)
    add_edge!(tree, 1, 2)
    add_edge!(tree, 1, 3)

    bags = [[2, 3, 4],
            [2, 4, 1],
            [3, 4, 5]]

    @test is_tree_decomposition(G, tree, bags)

    # removing node 1 of the second bag leads to an invalid
    # tree decomposition (not all vertices present)
    bags = [[2, 3, 4],
            [2, 4],
            [3, 4, 5]]

    @test ! (@test_logs (:warn, "Union of bags is not equal to union of vertices") is_tree_decomposition(G, tree, bags))


    # deleting node 2 of the first bag leads to an invalid
    # tree decomposition (edge (2, 3) not contained in any bag)
    bags = [[3, 4],
            [2, 4, 1],
            [3, 4, 5]]
    e = (2,3)
    @test ! (@test_logs (:warn, "Edge $e not found in any bag") is_tree_decomposition(G, tree, bags))

    # removing node 4 of the first bag leads to an invalid
    # tree decomposition (subgraph for vertex 4 not connected)
    bags = [[2, 3],
            [2, 4, 1],
            [3, 4, 5]]

    @test ! (@test_logs (:warn, "Subgraph for vertex 4 not connected") is_tree_decomposition(G, tree, bags))


    # validate tree decomposition algorithm on random graphs
    for i in 1:samples_per_test
        Nn = rand(20:50)
        Ne = rand(3*Nn:Nn*(Nn-1)÷2)
        G = random_graph(Nn, Ne)
        tw, tree, bags = tree_decomposition(G)
        @test is_tree_decomposition(G, tree, bags)
    end
end
