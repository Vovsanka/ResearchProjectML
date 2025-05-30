#include "clustering_problem.hpp"

template<typename S>
ClusteringProblem<S>::ClusteringProblem(
    const std::vector<S> &samples,
    const std::function<int(Utuple<3,S>)> &tripleCostCB,
    const std::function<int(Utuple<2,S>)> &pairCostCB
) : samples(samples) {
    int sampleCount = samples.size();
    // init relevant pairs
    relevantPairs.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            Upair indexPair({i, j});
            Utuple<2,S> samplePair({samples[i], samples[j]});
            int c = pairCostCB(samplePair);
            if (c) {
                pairCosts[indexPair] = c;
                relevantPairs[i].push_back(j);
                relevantPairs[j].push_back(i);
            }
        }
    }
    // init relevant triples
    relevantTriples.resize(sampleCount, {});
    for (int i = 0; i < sampleCount; i++) {
        for (int j = i + 1; j < sampleCount; j++) {
            for (int k = j + 1; k < sampleCount; k++) {
                Utriple indexTriple({i, j, k});
                Utuple<3,S> sampleTriple({samples[i], samples[j], samples[k]});
                int c = tripleCostCB(sampleTriple);
                if (c) {
                    tripleCosts[indexTriple] = c;
                    relevantTriples[i].push_back(std::make_pair(j, k)); // j < k
                    relevantTriples[j].push_back(std::make_pair(i, k)); // i < k
                    relevantTriples[k].push_back(std::make_pair(i, j)); // j < k
                }
            }
        }
    }
}
