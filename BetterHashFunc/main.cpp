#include "types.cpp"
#include <iostream>

int main() {
    ElasticHashTable<int, int> eht(100, 0.1, 99999999);
    std::vector<std::pair<int, int>> elts = {{1, 10}, {2, 20}, {3, 30}, {4, 40}, {5, 50}, {6, 60}, {7, 70}, {8, 80}, {9, 90}, {10, 100}};

    eht.insertBatch(elts);

    eht.printTable();
    return 0;
};

