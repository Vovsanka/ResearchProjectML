#include <iostream>

template <typename S = int> // domain S of samples
class CubicSetPartitionProblem {
    double (*cost)(S, S, S);  // cost function for triples SxSxS
public: 
    explicit CubicSetPartitionProblem(double (*costCB)(S, S, S)) {
        this->cost = costCB;
    }

    double getCost(S s1, S s2, S s3) {
        return cost(s1, s2, s3);
    }
};

double cost(char s1, char s2, char s3) {
    if (s1 == 'a' && s2 == 'b' && s3 == 'c')
        return -5;
    if (s1 == 'b' && s2 == 'c' && s3 == 'd')
        return -7;
    if (s1 == 'a' && s2 == 'b' && s3 == 'd')
        return 10;
    return 0;
}

int main() {
    CubicSetPartitionProblem<char> problem(cost);
    std::cout << problem.getCost('a', 'b', 'c') << std::endl;
    return 0;
}