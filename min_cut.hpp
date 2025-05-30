#ifndef MIN_CUT_HPP
#define MIN_CUT_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>
#include <boost/property_map/property_map.hpp>


namespace MinCut {

    using namespace boost;

    typedef adjacency_list<vecS, vecS, directedS,
        property<vertex_name_t, std::string>,
        property<edge_capacity_t, int,
            property<edge_residual_capacity_t, int,
                property<edge_reverse_t, adjacency_list<>::edge_descriptor>>>> Digraph;

    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, int>> Bigraph;
    typedef property_map<Bigraph, edge_weight_t>::type WeightMap;

    int solveMinCut(int vertices, std::vector<std::tuple<int,int,int>> edges, int s, int t); // handled as bigraph!

    int solveGlobalMinCut(std::vector<std::tuple<int,int,int>> edges);
}

#endif