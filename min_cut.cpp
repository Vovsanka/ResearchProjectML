#include "min_cut.hpp"


int MinCut::solveMinCut(int vertices, std::vector<std::tuple<int,int,int>> edges, int s, int t) {
    Digraph g(vertices); 
    auto capacity = get(edge_capacity, g);
    auto rev = get(edge_reverse, g);
    auto residual_capacity = get(edge_residual_capacity, g);

    // Add edges with capacities
    for (auto& e : edges) {
        auto [u, v, c] = e;
        if (c < 0) throw std::runtime_error("MinCutProblem: negative edges are not allowed!");
        auto e1 = add_edge(u, v, g).first;
        auto e1r = add_edge(v, u, g).first;
        auto e2 = add_edge(v, u, g).first; 
        auto e2r = add_edge(u, v, g).first;
        capacity[e1] = capacity[e2] = c; 
        capacity[e1r] = capacity[e2r] = 0;
        rev[e1] = e1r;
        rev[e1r] = e1;
        rev[e2] = e2r;
        rev[e2r] = e2;
    };

    // Compute max flow
    int minCut = push_relabel_max_flow(g, s, t);

    std::vector<bool> visited(vertices, false);

    std::function<void(int)> dfs = [&](int v) -> void {
        if (visited[v]) return;
        visited[v] = true;
        for (auto [ei, e_end] = out_edges(v, g); ei != e_end; ++ei) {
            if (residual_capacity[*ei] > 0) {
                dfs(target(*ei, g));
            }
        }
    };

    // Run DFS on the residual network to determine the min cut subset
    dfs(s);

    return minCut;
}

int MinCut::solveGlobalMinCut(std::vector<std::tuple<int,int,int>> edges) {
    Bigraph g;
    for (auto& e : edges) {
        auto [u, v, c] = e;
        if (c < 0) throw std::runtime_error("MinCutProblem: negative edges are not allowed!");
        add_edge(u, v, c, g);
    };

    WeightMap weight_map = get(edge_weight, g);

    std::vector<bool> parity(num_vertices(g));
    auto parity_map = make_iterator_property_map(parity.begin(), get(vertex_index, g));

    int minCut = stoer_wagner_min_cut(g, weight_map, boost::parity_map(parity_map));

    return minCut;
}