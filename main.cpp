#include <iostream>

#include "clustering_problem.hpp"



//     ClusteringProblem<S> createIndependentCutSubproblem(std::vector<int> subsamples) {
//         int subSampleCount = subsamples.size();
        
//         // subsamples as a set
//         std::vector<bool> subsetR(sampleCount, false);
//         for (auto i : subsamples) {
//             subsetR[i] = true;
//         }

//         // index of samples in the subproblem
//         std::vector<int> indexOf(sampleCount);
//         for (int ind = 0; ind < subSampleCount; ind++) {
//             indexOf[subsamples[ind]] = ind;
//         }

//         // filter relevant pairs for th subproblem
//         std::vector<std::vector<int>> subRelevantPairs(subSampleCount);
//         std::map<UnorderedPair<>, int> subPairCosts;
//         for (int i : subsamples) {
//             for (int j : relevantPairs[i]) {
//                 if (subsetR[j]) {
//                     subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
//                     subPairCosts[UnorderedPair<>(indexOf[i], indexOf[j])] = pairCosts[UnorderedPair<>(i, j)];
//                 }
//             }
//         }
//         // filter relevant triples for the subproblem
//         std::vector<std::vector<std::pair<int, int>>> subRelevantTriples(subSampleCount);
//         std::map<UnorderedTriple<>, int> subTripleCosts;
//         for (int i : subsamples) {
//             for (auto [j, k] : relevantTriples[i]) {
//                 if (subsetR[j] && subsetR[k]) {
//                     subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
//                     subTripleCosts[UnorderedTriple(indexOf[i], indexOf[j], indexOf[k])] = tripleCosts[UnorderedTriple<>(i, j, k)];
//                 }
//             }
//         }
//         // create the subproblem
//         return ClusteringProblem(
//             subSampleCount,
//             subRelevantTriples,
//             subTripleCosts,
//             subRelevantPairs,
//             subPairCosts
//         );
//     }

//     bool applyIndependentSubproblemCut() {
//         std::vector<std::vector<int>> partition;
//         std::vector<bool> processed(sampleCount, false);
//         for (int i = 0; i < sampleCount; i++) {
//             if (processed[i]) continue;
//             // start the BFS from i
//             std::vector<int> chosen = {i};
//             std::queue <int> q;
//             processed[i] = true;
//             q.push(i);
//             while (!q.empty()) {
//                 int current = q.front();
//                 q.pop();
//                 for (auto j : relevantPairs[current]) {
//                     UnorderedPair<> indexPair(current, j);
//                     int c = pairCosts[indexPair];
//                     if (c >= 0) continue; 
//                     if (!processed[j]) {
//                         processed[j] = true;
//                         chosen.push_back(j);
//                         q.push(j);
//                     }
//                 }
//                 for (auto [j, k] : relevantTriples[current]) {
//                     UnorderedTriple<> indexTriple(current, j, k);
//                     int c = tripleCosts[indexTriple];
//                     if (c >= 0) continue; 
//                     if (!processed[j]) {
//                         processed[j] = true;
//                         chosen.push_back(j);
//                         q.push(j);
//                     }
//                     if (!processed[k]) {
//                         processed[k] = true;
//                         chosen.push_back(k);
//                         q.push(k);
//                     }
//                 }
//             }
//             partition.push_back(chosen);
//         }

//         if (partition.size() == 1) return false; // no cuts => no smaller subproblems
//         std::cout << "Applying the independent subproblem cut (proposition 3.1)" << std::endl;

//         // fix the labels
//         for (int subsetI = 0; subsetI < partition.size(); subsetI++) {
//             for (int subsetJ = subsetI + 1; subsetJ < partition.size(); subsetJ++) {
//                 for (auto i : partition[subsetI]) {
//                     for (auto j : partition[subsetJ]) {
//                         labelFixed[i][j] = labelFixed[j][i] = true;
//                         labelValue[i][j] = labelValue[j][i] = false;
//                     }
//                 }
//             }
//         }

//         // create, solve the subproblems, accumulate the results
//         int clusterOffset = 0;
//         for (auto backIndexing : partition) {
//             // create and solve the independent subproblems
//             auto subproblem = createIndependentCutSubproblem(backIndexing);
//             subproblem.solve();
//             // accumulate the results
//             resultingCost += subproblem.getResultingCost();
//             // accumulate the labels
//             auto subLabelFixed = subproblem.getIndexLabelFixed();
//             auto subLabelValue = subproblem.getIndexLabelValue();
//             for (int i = 0; i < subLabelFixed.size(); i++) {
//                 for (int j = i + 1; j < subLabelFixed.size(); j++) {
//                     int originalI = backIndexing[i], originalJ = backIndexing[j];
//                     labelFixed[originalI][originalJ] = labelFixed[originalJ][originalI] = subLabelFixed[i][j];
//                     labelValue[originalI][originalJ] = labelValue[originalJ][originalI] = subLabelValue[i][j];
//                 }
//             } 
//         }
        
//         return true;
//     }

//     void createSolveAccumulateJoinSubproblem(std::vector<bool> subsetR) {
//         std::vector<int> joinSamples;
//         for (int i = 0; i < sampleCount; i++) {
//             if (subsetR[i]) joinSamples.push_back(i);
//         }

//         // aplying proposition 5.1
//         int subSampleCount = sampleCount - joinSamples.size() + 1;
//         int indexOfJoint = subSampleCount - 1;

//         // index of samples in the subproblem
//         std::vector<std::vector<int>> backIndexing(subSampleCount);
//         std::vector<int> indexOf(sampleCount);
//         for (int i = 0, ind = 0; i < sampleCount; i++) {
//             if (subsetR[i]) {
//                 indexOf[i] = indexOfJoint;
//                 backIndexing[indexOfJoint].push_back(i);
//             } else {
//                 indexOf[i] = ind;
//                 backIndexing[ind].push_back(i);
//                 ind++;
//             }
//         }

//         // filter relevant pairs for the subproblem which have no relations to the samples being joint
//         std::vector<std::vector<int>> subRelevantPairs(subSampleCount);
//         std::map<UnorderedPair<>, int> subPairCosts;
//         for (int i = 0; i < sampleCount; i++) {
//             if (subsetR[i]) continue; // skip the relations for the joint set
//             for (int j : relevantPairs[i]) {
//                 if (subsetR[j]) continue; // skip the relations for the joint set
//                 subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
//                 subPairCosts[UnorderedPair<>(indexOf[i], indexOf[j])] = pairCosts[UnorderedPair<>(i, j)];
//             }
//         }
        
//         // filter relevant triples for the subproblem which have no relations to the samples being joint
//         std::vector<std::vector<std::pair<int, int>>> subRelevantTriples(subSampleCount);
//         std::map<UnorderedTriple<>, int> subTripleCosts;
//         for (int i = 0; i < sampleCount; i++) {
//             if (subsetR[i]) continue; // skip the relations for the joint set
//             for (auto [j, k] : relevantTriples[i]) {
//                 if (subsetR[j] || subsetR[k]) continue; // skip the relations for the joint set
//                 subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
//                 subTripleCosts[UnorderedTriple(indexOf[i], indexOf[j], indexOf[k])] = tripleCosts[UnorderedTriple<>(i, j, k)];
//             }
//         }

//         // compute the costs to the joint subset as well as the inner joining cost
//         for (int i : joinSamples) {
//             for (int j : relevantPairs[i]) {
//                 if (subsetR[j]) {
//                     if (i < j) resultingCost += pairCosts[UnorderedPair<>(i, j)]; // consider one direction (i, j) and skip (j, i)
//                 } else {
//                     UnorderedPair<> indexPair(indexOf[i], indexOf[j]);
//                     if (!subPairCosts[indexPair]) {
//                         subPairCosts[indexPair] = 0;
//                         subRelevantPairs[indexOf[i]].push_back(indexOf[j]);
//                         subRelevantPairs[indexOf[j]].push_back(indexOf[i]);
//                     }
//                     subPairCosts[indexPair] += pairCosts[UnorderedPair<>(i, j)];
//                 }
//             }
//         }
//         for (int i : joinSamples) {
//             for (auto [j, k] : relevantTriples[i]) {
//                 bool innerJ = subsetR[j], innerK = subsetR[k];
//                 if (innerJ && innerK) {
//                     // (i < j < k) consider only (i, j, k) and skip the others
//                     if (i < j) resultingCost += tripleCosts[UnorderedTriple<>(i, j, k)];
//                 } else if (innerJ || innerK) {
//                     int outer, inner;
//                     if (innerJ) {
//                         inner = j;
//                         outer = k;
//                     }
//                     if (innerK) {
//                         inner = k;
//                         outer = j;
//                     }
//                     if (i > inner) continue; // consider one direction (i, inner) and skip (inner, i)
//                     UnorderedPair<> indexPair(indexOf[i], indexOf[outer]);
//                     if (!subPairCosts[indexPair]) {
//                         subPairCosts[indexPair] = 0;
//                         subRelevantPairs[indexOf[i]].push_back(indexOf[outer]);
//                         subRelevantPairs[indexOf[outer]].push_back(indexOf[i]);
//                     }
//                     subPairCosts[indexPair] += tripleCosts[UnorderedTriple<>(i, j, k)];
//                 } else { // j and k are outer
//                     UnorderedTriple<> indexTriple(indexOf[i], indexOf[j], indexOf[k]);
//                     if (!subTripleCosts[indexTriple]) {
//                         subTripleCosts[indexTriple] = 0;
//                         subRelevantTriples[indexOf[i]].push_back(std::make_pair(indexOf[j], indexOf[k]));
//                         subRelevantTriples[indexOf[j]].push_back(std::make_pair(indexOf[k], indexOf[i])); // indexOf[i] > indexOf[k] by definition above
//                         subRelevantTriples[indexOf[k]].push_back(std::make_pair(indexOf[j], indexOf[i])); // indexOf[i] > indexOff[j] by definition above
//                     }
//                     subTripleCosts[indexTriple] += tripleCosts[UnorderedTriple<>(i, j, k)];
//                 }
//             } 
//         }

//         // create and solve the subproblem
//         auto subproblem = ClusteringProblem(
//             subSampleCount,
//             subRelevantTriples,
//             subTripleCosts,
//             subRelevantPairs,
//             subPairCosts
//         );
//         subproblem.solve();

//         // accumulate the results (resulting cost for the join has been accumulated)
//         resultingCost += subproblem.getResultingCost();
//         // fix the labels for join
//         for (int indI = 0; indI < joinSamples.size(); indI++) {
//             for (int indJ = indI + 1; indJ < joinSamples.size(); indJ++) {
//                 int i = joinSamples[indI], j = joinSamples[indJ];
//                 labelFixed[i][j] = labelFixed[j][i] = true;
//                 labelValue[i][j] = labelValue[j][i] = true;
//             }
//         }
//         // accumulate the labels
//         auto subLabelFixed = subproblem.getIndexLabelFixed();
//         auto subLabelValue = subproblem.getIndexLabelValue();
//         for (int i = 0; i < subLabelFixed.size(); i++) {
//             for (int j = i + 1; j < subLabelFixed.size(); j++) {
//                 for (int originalI : backIndexing[i]) {
//                     for (int originalJ : backIndexing[j]) {
//                         labelFixed[originalI][originalJ] = labelFixed[originalJ][originalI] = subLabelFixed[i][j];
//                         labelValue[originalI][originalJ] = labelValue[originalJ][originalI] = subLabelValue[i][j];
//                     }
//                 }   
//             }
//         } 
//     }

//     bool checkSubsetJoinForIndexSubset(std::vector<bool> &indexSubset) {   
//         // compute rhs
//         int lhsLowerBound = 0; // avoid MinCut computation for lhs>rhs (in particular the edge cases with lhs=0, rhs<0 because of 3.1 applied before)
//         int singleR = 0, doubleR = 0;
//         for (int i = 0; i < sampleCount; i++) {
//             if (!indexSubset[i]) continue; // i is in R
//             for (auto j : relevantPairs[i]) {
//                 int c = pairCosts[UnorderedPair<>(i, j)];
//                 if (indexSubset[j]) {
//                     if (i < j) lhsLowerBound += c;
//                 } else if (c < 0) {
//                     singleR += c;
//                 }
//             }
//             for (auto [j, k] : relevantTriples[i]) {
//                 int c = tripleCosts[UnorderedTriple<>(i, j, k)];
//                 if (indexSubset[j] && indexSubset[k]) {
//                     if (i < j) lhsLowerBound += c;
//                     continue; // j or k must be not in R
//                 }
//                 if (c > 0) continue; // omit positive costs
//                 if (!indexSubset[k] && !indexSubset[j]) {
//                     singleR += c;
//                 } else {
//                     doubleR += c; // since two elements of the triple are in R, the cost will be added twice
//                 }
//             }
//         }
//         int rhs = singleR + doubleR/2;
//         if (lhsLowerBound > rhs) return false;

//         int lhs = -solveMinCutForIndexSubset(indexSubset, true, false, true);
//         return (lhs <= rhs);
//     }

//     bool applySubsetJoin() {
//         // subset join criterion 3.11
//         // heuristically construct and check the candidate sets R for possible joining
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j = i + 1; j < sampleCount; j++) {
//                 if (pairCosts[UnorderedPair<>(i, j)] > 0) continue;

//                 std::vector<bool> subsetR(sampleCount, false); // R doesn't have to be a connected component
//                 subsetR[i] = subsetR[j] = true;
//                 std::vector<bool> joinR;
//                 if (checkSubsetJoinForIndexSubset(subsetR)) joinR = subsetR;

//                 std::set<int> candidates;
//                 for (int k = 0; k < sampleCount; k++) {
//                     if (!subsetR[k]) {
//                         candidates.insert(k);
//                     }
//                 }

//                 std::function<int(int)> computeOffset = [&](int i) {
//                     // positive offset if not mergeable
//                     if (subsetR[i]) return 1;
//                     int offset = 0;
//                     for (auto j : relevantPairs[i]) {
//                         if (!subsetR[j]) continue;
//                         int c = pairCosts[UnorderedPair<>(i, j)];
//                         if (c > 0) {
//                             return 1;
//                         } else {
//                             offset += c;
//                         }
//                     }
//                     for (auto [j, k] : relevantTriples[i]) {
//                         if (!subsetR[j] || !subsetR[k]) continue; // 2 of 3 triple elements are already in R
//                         int c = tripleCosts[UnorderedTriple<>(i, j, k)];
//                         if (c > 0) {
//                             return 1;
//                         } else {
//                             offset += c;
//                         }
//                     }
//                     return offset;
//                 };

//                 while(!candidates.empty()) {
//                     std::vector<int> badCandidates;
//                     int bestCandidate = -1, bestOffset = 0;
//                     for (auto k : candidates) {
//                         int offset = computeOffset(k);
//                         if (offset > 0) {
//                             badCandidates.push_back(k);
//                         }
//                         if (offset <= bestOffset) {
//                             bestOffset = offset;
//                             bestCandidate = k;
//                         }
//                     }
//                     for (auto k : badCandidates) {
//                         candidates.erase(k);
//                     }
//                     if (bestCandidate != -1) {
//                         subsetR[bestCandidate] = true;
//                         candidates.erase(bestCandidate);
//                         if (checkSubsetJoinForIndexSubset(subsetR)) joinR = subsetR;
//                     }
//                 }
//                 if (!joinR.empty()) {
//                     std::cout << "Applying the subset join (proposition 3.11)" << std::endl;
//                     createSolveAccumulateJoinSubproblem(joinR);
//                     return true;
//                 }
//             }
//         }
//         return false;
//     }
    
//     bool applyPairJoin() {
//         // proposition 3.4
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j = i + 1; j < sampleCount; j++) {
//                 // compute lhs
//                 int lhs = 0;
//                 UnorderedPair<> indexPair(i, j);
//                 int c = pairCosts[indexPair];
//                 if (c < 0) lhs += -2*c;
//                 for (auto [k1, k2] : relevantTriples[i]) {
//                     if (k1 != j && k2 != j) continue; // k1 or k2 is j
//                     int c = tripleCosts[UnorderedTriple<> (i, k1, k2)];
//                     if (c < 0) lhs += -c;
//                 }
//                 // compute rhs (assume that the underlying graph is connected after applying the proposition 3.1)
//                 int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, i, {j});

//                 if (lhs >= rhs) {
//                     std::vector<bool> subsetR(sampleCount, false);
//                     subsetR[i] = subsetR[j] = true;
//                     std::cout << "Applying the pair join (proposition 3.4)" << std::endl;
//                     // std::cout << i << " " << j << std::endl;
//                     // std::cout << lhs << " vs " << rhs << std::endl;
//                     createSolveAccumulateJoinSubproblem(subsetR);
//                     return true;
//                 }
//             }
//         }
//         return false;
//     }

//     bool applyComplexPairJoin() {
//         // proposition 3.6
//         for (int i = 0; i < sampleCount; i++) {
//             for (int k = i + 1; k < sampleCount; k++) {
//                 for (int j = 0; j < sampleCount; j++) {
//                     if (j == i || j == k) continue;
//                     // compute lhs1, lhs2, lhs3
//                     int lhs1 = 0, lhs2 = 0, lhs3 = 0;
//                     UnorderedTriple<> indexTriple(i, j, k);
//                     {
//                         int c = tripleCosts[indexTriple];
//                         lhs3 += c;
//                         if (c < 0) {
//                             lhs1 += -c;
//                             lhs2 += -c; 
//                         }
//                     }
//                     UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
//                     {
//                         int c = pairCosts[indexPairIJ];
//                         lhs3 += c;
//                         if (c < 0) lhs1 += -2*c; 
//                     }
//                     {
//                         int c = pairCosts[indexPairJK];
//                         lhs3 += c;
//                         if (c < 0) lhs2 += -2*c; 
//                     }
//                     {
//                         int c = pairCosts[indexPairIK];
//                         lhs3 += c;
//                         if (c < 0) {
//                             lhs1 += -2*c;
//                             lhs2 += -2*c;
//                         }
//                     }
//                     for (auto [p, q] : relevantTriples[i]) {
//                         if (p == j || p == k || q == j || q == k) {
//                             int c = tripleCosts[UnorderedTriple<>(i, p, q)];
//                             if (c < 0) lhs1 += -c;
//                         }
//                     }
//                     for (auto [p, q] : relevantTriples[k]) {
//                         if (p == j || p == i || q == j || q == i) {
//                             int c = tripleCosts[UnorderedTriple<>(k, p, q)];
//                             if (c < 0) lhs2 += -c;
//                         }
//                     }
//                     // compute rhs3
//                     int rhs3 = 0;
//                     for (auto [p, q] : relevantTriples[i]) {
//                         if (p == j || p == k || q == j || q == k) continue;
//                         int c = tripleCosts[UnorderedTriple<>(i, p, q)];
//                         rhs3 -= abs(c);
//                     }
//                     for (auto [p, q] : relevantTriples[j]) {
//                         if (p == i || p == k || q == i || q == k) continue;
//                         int c = tripleCosts[UnorderedTriple<>(j, p, q)];
//                         rhs3 -= abs(c);
//                     }
//                     for (auto [p, q] : relevantTriples[k]) {
//                         if (p == i || p == j || q == i || q == j) continue;
//                         int c = tripleCosts[UnorderedTriple<>(k, p, q)];
//                         rhs3 -= abs(c);
//                     }
//                     for (auto p : relevantPairs[i]) {
//                         if (p == j || p == k) continue;
//                         int c = pairCosts[UnorderedPair(i, p)];
//                         rhs3 -= abs(c);
//                     }
//                     for (auto p : relevantPairs[j]) {
//                         if (p == i || p == k) continue;
//                         int c = pairCosts[UnorderedPair(j, p)];
//                         rhs3 -= abs(c);
//                     }
//                     for (auto p : relevantPairs[k]) {
//                         if (p == i || p == j) continue;
//                         int c = pairCosts[UnorderedPair(k, p)];
//                         rhs3 -= abs(c);
//                     }
//                     // check the 3d condition
//                     if (!(lhs3 <= rhs3)) continue; 
//                     // check the 1st and the 2d condition
//                     int rhs1 = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, i, {j, k});
//                     int rhs2 = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, true, false, k, {i, j});
//                     if (lhs1 >= rhs1 && lhs2 >= rhs2) {
//                         std::vector<bool> subsetR(sampleCount, false);
//                         subsetR[i] = subsetR[k] = true; // join i and k
//                         std::cout << "Applying the complex pair join (proposition 3.6)" << std::endl;
//                         // std::cout << i << " " << j << " " << k << std::endl;
//                         // std::cout << "(" << lhs1 << "," << lhs2 << "," << lhs3 << ") vs ";
//                         // std:: cout << "(" << rhs1 << "," << rhs2 << "," << rhs3 << ")" << std::endl;
//                         createSolveAccumulateJoinSubproblem(subsetR);
//                         return true;
//                     }
//                 }
//             }
//         }
//         return false;
//     }

//     bool applyExplicitPairJoin() {
//         // proposition 3.8
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j = i + 1; j < sampleCount; j++) {
//                 int lhs = 0;
//                 UnorderedPair indexPair(i, j);
//                 if (pairCosts[indexPair]) lhs += pairCosts[indexPair];
//                 // compute rhs
//                 int rhs = 0;
//                 for (int k : relevantPairs[i]) {
//                     if (k == j) continue;
//                     int c = pairCosts[UnorderedPair<>(i, k)];
//                     if (c < 0) rhs += c;
//                 }
//                 for (int k : relevantPairs[j]) {
//                     if (k == i) continue;
//                     int c = pairCosts[UnorderedPair<>(j, k)];
//                     if (c < 0) rhs += c;
//                 }
//                 for (auto [k1, k2] : relevantTriples[i]) {
//                     int c = tripleCosts[UnorderedTriple<>(i, k1, k2)];
//                     if (c < 0) rhs += c;
//                 }
//                 for (auto [k1, k2] : relevantTriples[j]) {
//                     if (k1 == i || k2 == i) continue; // the triples (i, j, *) have already been considered in the previos for-loop
//                     int c = tripleCosts[UnorderedTriple<>(j, k1, k2)];
//                     if (c < 0) rhs += c;
//                 }
//                 // check the condition
//                 if (lhs <= rhs) {
//                     std::vector<bool> subsetR(sampleCount, false);
//                     subsetR[i] = subsetR[j] = true;
//                     std::cout << "Applying the explicit pair join (proposition 3.8)" << std::endl;
//                     std::cout << i << " " << j << std::endl;
//                     std::cout << lhs << " vs " << rhs << std::endl;
//                     createSolveAccumulateJoinSubproblem(subsetR);
//                     return true;
//                 }
//             }
//         }
//         return false;
//     }

//     bool applyExplicitPairJoinViaTriple() {
//         // proposition 3.9
//         for (int i = 0; i < sampleCount; i++) {
//             for (int k = i + 1; k < sampleCount; k++) {
//                 for (int j = 0; j < sampleCount; j++) {
//                     if (j == i || j == k) continue;
//                     // compute costs for the triple
//                     UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
//                     UnorderedTriple<> indexTriple(i, j, k);
//                     int cIJ = pairCosts[indexPairIJ];
//                     int cIK = pairCosts[indexPairIK];
//                     int cJK = pairCosts[indexPairJK];
//                     int cIJK = tripleCosts[indexTriple];
//                     // compute rhs
//                     int singleR = 0, doubleR = 0;
//                     std::vector<int> elementsR = {i, j, k};
//                     std::vector<bool> subsetR(sampleCount, false);
//                     subsetR[i] = subsetR[j] = subsetR[k] = true;
//                     for (int i1 : elementsR) {
//                         for (int j1 : relevantPairs[i1]) {
//                             if (subsetR[j1]) continue;
//                             int c = pairCosts[UnorderedPair<>(i1, j1)];
//                             if (c < 0) singleR += c;
//                         }
//                         for (auto [j1, k1] : relevantTriples[i1]) {
//                             if (subsetR[j1] && subsetR[k1]) continue;
//                             int c = tripleCosts[UnorderedTriple<>(i1, j1, k1)];
//                             if (c > 0) continue;
//                             if (subsetR[j1] xor subsetR[k1]) {
//                                 doubleR += c;
//                             } else { // only i1 is in R
//                                 singleR += c;
//                             }
//                         }
//                     }
//                     int rhs = singleR + doubleR/2;
//                     // check the conditions
//                     if (
//                         cIJ + cIK <= 0 &&
//                         cIJ + cJK <= 0 &&
//                         cIK + cJK <= 0 &&
//                         cIJ + cIK + cJK <= 0 &&
//                         2*cIJ + 2*cIK + 2*cJK + cIJK <= 0 &&
//                         cIJ + cIK + cIJK <= rhs &&
//                         cJK + cIK + cIJK <= rhs
//                     ) {
//                         std::vector<bool> subsetR(sampleCount, false);
//                         subsetR[i] = subsetR[k] = true; // join i and k
//                         std::cout << "Applying the complex pair join (proposition 3.9)" << std::endl;
//                         std::cout << i << " " << j << " " << k << std::endl;
//                         std::cout << rhs << std::endl;
//                         createSolveAccumulateJoinSubproblem(subsetR);
//                         return true;
//                     }
//                 }
//             }
//         }
//         return false;
//     }

//     bool applyTripleJoin() {
//         // proposition 3.5
//         // compute lhs base value
//         int lhsBase = 0;
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j : relevantPairs[i]) {
//                 if (i > j) continue;
//                 int c = pairCosts[UnorderedPair<>(i, j)];
//                 if (c > 0) lhsBase -= c;
//             }
//             for (auto [j, k] : relevantTriples[i]) {
//                 if (i > j) continue;
//                 int c = tripleCosts[UnorderedTriple<>(i, j, k)];
//                 if (c > 0) lhsBase -= c;
//             }
//         }
//         // iterate over all unordered triples (i, j, k)
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j = i + 1; j < sampleCount; j++) {
//                 for (int k = j + 1; k < sampleCount; k++) {
//                     // compute lhs
//                     int lhs = lhsBase;
//                     UnorderedPair<> indexPairIJ(i, j), indexPairIK(i, k), indexPairJK(j, k);
//                     UnorderedTriple<> indexTriple(i, j, k);
//                     int cIJ = pairCosts[indexPairIJ];
//                     if (cIJ < 0) lhs += -2*cIJ;
//                     int cIK = pairCosts[indexPairIK];
//                     if (cIK < 0) lhs += -2*cIK;
//                     int cJK = pairCosts[indexPairJK];
//                     if (cJK < 0) lhs += -cJK;
//                     int cIJK = tripleCosts[indexTriple];
//                     if (cIJK < 0) lhs += -2*cIJK;
//                     lhs += std::min(std::min(0, cIJ), std::min(cIK, cJK)); // min for 4 cases
//                     // compute rhs
//                     int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, false, false, i, {j, k});
//                     if (lhs >= rhs) {
//                         std::vector<bool> subsetR(sampleCount, false);
//                         subsetR[i] = subsetR[j] = subsetR[k] = true; // join ijk
//                         std::cout << "Applying the triple join (proposition 3.5)" << std::endl;
//                         std::cout << i << " " << j << " " << k << std::endl;
//                         std::cout << lhs << " vs " << rhs << std::endl;
//                         createSolveAccumulateJoinSubproblem(subsetR);
//                         return true;
//                     }
//                 }
//             }
//         }
//         return false;
//     }

//     void applyPairCuts() {
//         // proposition 3.2
//         for (int i = 0; i < sampleCount; i++) {
//             for (int j = i + 1; j < sampleCount; j++) {
//                 // check if already cut
//                 if (labelFixed[i][j]) {
//                     if (labelValue[i][j]) throw std::runtime_error("Error in the program detected: cannot cut the joint elements! ");
//                     else continue;
//                 }
//                 // check the condition
//                 int lhs = 0;
//                 UnorderedPair<> indexPair(i, j);
//                 {
//                     int c = pairCosts[indexPair];
//                     if (c > 0) lhs += c;
//                 }
//                 int rhs = solveMinCutForIndexSubset(std::vector<bool>(sampleCount, true), true, false, false, i, {j});
//                 if (lhs >= rhs) {
//                     std::cout << "Applying the pair cut (proposition 3.2)" << std::endl;
//                     std::cout << i << " " << j << std::endl;
//                     std::cout << lhs << " vs " << rhs << std::endl;
//                     // assume no joint labels in this problems, because it would be solved by a join-subproblem!
//                     labelFixed[i][j] = labelFixed[j][i] = true;
//                     labelValue[i][j] = labelValue[j][i] = false;
//                 }
//             }
//         }
//     }

//     void applyTripleCuts() {
//         // TODO
//     }



int cost(Utuple<3,char> t) {
    // if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
    //     return -5;
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
        return -22; // -7
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
        return 10;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e')
        return 10;
    return 0;
}

// int cost(UnorderedTriple<char> t) {
//     // example: 3.1 and 3.11 are sufficient
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -1;
//     if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -15;
//     if (t[0] == 'd' && t[1] == 'e' && t[2] == 'h') return 50;
//     if (t[0] == 'e' && t[1] == 'f' && t[2] == 'h') return -50;
//     if (t[0] == 'f' && t[1] == 'g' && t[2] == 'i') return -30;
//     if (t[0] == 'd' && t[1] == 'h' && t[2] == 'k') return -2;
//     if (t[0] == 'i' && t[1] == 'k' && t[2] == 'l') return -4;
//     if (t[0] == 'j' && t[1] == 'l' && t[2] == 'm') return -10;
//     return 0;
// }


// int cost(UnorderedTriple<char> t) {
//     // // pyramid example below (3.1 + 3.11 are not sufficient) (3.4 is sufficient)
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e') return -75;
//     // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10) (3.5 is sufficient for 100)
//     if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 10; // 10
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
//     if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
//     return 0;
// }

int main() {
    std::vector<char> samples = {'a', 'b', 'c', 'd', 'e'};
    ClusteringProblem<char> problem(samples, cost);
    problem.solve();
    problem.printResults();
    
    return 0;
}