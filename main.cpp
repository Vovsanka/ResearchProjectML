#include <iostream>

#include "clustering_problem.hpp"



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
    // pyramid example below (3.1 + 3.11 + 3.4 are not sufficient) (3.6 is sufficient for 10) (everything is insufficient for 100)
    if (t[0] == 'b' && t[1] == 'c' && t[2] == 'd') return 1000; // 10 or 100
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