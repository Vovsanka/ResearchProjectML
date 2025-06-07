#include <iostream>

#include "clustering_problem.hpp"





//     bool applyExplicitPairJoinViaTriple() {
//         // 3.9
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
//                         std::cout << "Applying the complex pair join (3.9)" << std::endl;
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
//         // 3.5
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
//                         std::cout << "Applying the triple join (3.5)" << std::endl;
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
//         // 3.2
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
//                     std::cout << "Applying the pair cut (3.2)" << std::endl;
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



// int cost(Utuple<3,char> t) {
//     // if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c')
//     //     return -5;
//     if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd')
//         return -22; // -7
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd')
//         return 10;
//     if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e')
//         return 10;
//     return 0;
// }

// int cost(Utuple<3,char> t) {
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


int cost(Utuple<3,char> t) {
    // // pyramid example below (3.1 + 3.11 are not sufficient) (3.4 is sufficient)
    // if (t[0] == 'a' && t[1] == 'b' && t[2] == 'e') return -75;
    // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10) (3.5 is sufficient for 100)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 100; // 10 or 100
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'c') return -50;
    if (t[0] == 'a' && t[1] == 'b' && t[2] == 'd') return -50;
    if (t[0] == 'a' && t[1] == 'c' && t[2] == 'd') return -50;
    return 0;
}

int main() {
    // std::vector<char> samples = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm'};
    std::vector<char> samples = {'a', 'b', 'c', 'd', 'e'};
    ClusteringProblem<char> problem(samples, cost);
    problem.solve();
    return 0;
}