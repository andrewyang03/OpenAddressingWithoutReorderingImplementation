#include <stdio.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <cmath>
#include <bitset>
#include <functional>

template <typename K, typename V>
class ElasticHashTable {
    private:
        int c; // C is a large positive constant
        int n; // Number of Arrays
        double delta; // Defined delta
        std::vector<std::vector<std::pair<K, V>>>arrays; // Vector for each array instance
        
        int phi(int i, int j) {
            std::bitset<32> bin_i(i);
            std::bitset<32> bin_j(j);

            std::string merged = "";

            for (int k = 31; k >=0 ; --k) {
                // Interweaves the values of binary j digits with a preceding 1
                merged += "1";
                merged += bin_j[k] ? "1": "0";

            }
            merged += "0";

            for (int k = 31; k >= 0; --k) {
                // Concatenates the values of binary i digits
                merged += bin_i[k] ? "1": "0";
            }
            
            // Return and recast as an int
            return static_cast<int>(std::bitset<64>(merged).to_ulong());
        }


        int hashFunction(const K& key, int i, int j) {
            std::hash<K> hashFunc;
            size_t hashVal = hashFunc(key);
            int returnVal = static_cast<int>((hashVal ^ phi(i, j)) % arrays[i].size());
            return std::abs(returnVal);
        }

        int f(double epsilon) {
            return static_cast<int>(c * 
                std::min(
                    std::pow(std::log(std::pow(epsilon, -1)), 2), 
                    std::log(std::pow(delta, -1))
                )
            );
        }
    
    public:
        ElasticHashTable(int size, double d, int c): n(size), delta(d), c(c) {
            // Break A into size n disjoint Arrays satisfying A_n+1 = abs(abs(A_n) / (2 + - 1))
            int numArrays = std::ceil(log2(size));
            arrays.resize(numArrays);

            int currrentSize = n;

            for (int i = 0; i < numArrays; ++i) {
                arrays[i].resize(size / std::pow(2, i));

                // Initialize each slot with a sentinal value
                for (auto& slot: arrays[i]) {
                    slot = std::make_pair(K(-1), V());
                }

                // Initialize the first array with 75% occupancy
                if (i == 0) {
                    int fillCt = static_cast<int>(0.75 * arrays[i].size());
                    for (int j = 0; j < fillCt; ++j) {
                        arrays[i][j] = std::make_pair(K(j), V(j));
                    }
                }

                // Adjust for the next one
                currrentSize = currrentSize / 2 + (std::rand() % 3 - 1); //adjust by +- 1
            }
        }

        void insertBatch(const std::vector<std::pair<K, V>>& elements) {
            int maxInsertions = n - std::floor(delta * n);
            int insertCt = 0;

            for (int batch = 0; batch < arrays.size(); ++batch) {
                int i = batch;

                int batchSize = calcBatchSize(i);

                for (int j = 0; j < batchSize && insertCt < elements.size() && insertCt < maxInsertions; ++j) {
                    const auto& element = elements[insertCt];
                    //Insert the key value pair at index i
                    insertElement(element.first, element.second, i);
                    std::cout << "Inserting key: " << element.first << ", value: " << element.second << ", into array " << i << std::endl;
                    insertCt ++;
                }
            }
        }

        void printTable() {
            // Print the contents of the hash table for verification
            std::cout << "Hash Table after insertion" <<std::endl << std::flush;
            for (size_t i = 0; i < arrays.size(); ++ i) {
                std::cout <<"Array A_" << i + 1 << ":" << std::endl << std::flush;
                for (const auto& slot : arrays[i]) {
                    if (slot.first != K()) {
                        std::cout << " Key: " << slot.first << ", Value: " << slot.second << std::endl << std::flush;
                    }
                }
            }
        }
    
    private:
        int calcBatchSize(int i) {
            if (i == 0) {
                // is ceil(0.75 * |A_1|)
                int batchSize = static_cast<int>(std::ceil(0.75 * arrays[0].size()));
                return std::max(1, batchSize);
            } else {
                // |A_i| - floor(delta |A_i| / 2) - ceiling(0.75 |A_i|) + ceil(0.75 |A_i+1|)
                int A_i = arrays[i].size();
                int A_iPlusOne = arrays[i + 1].size();
                int batchSize =  A_i - std::floor(delta * A_i / 2) - std::ceil(0.75 * A_i) + std::ceil(0.75 * A_iPlusOne);
                return std::max(1, batchSize);
            }
        }

        int countOccupiedSlots(const std::vector<std::pair<K, V>>& array) {
            int count = 0;
            // Iterate and check if each slot is occupied by a key
            for (const auto& slot : array) {
                if (slot.first != K(-1)) {
                    count++;
                }
            }
            return count;
        }


        void insertElement(const K& key, const V& value, int i) {
            // Determines the fraction of free spaces for A_i and A_i+1

            // epsilon 1 = (abs(A_i) - expected batch size) / abs(A_i)
            double epsilon1 = 1.0 - (static_cast<double>(countOccupiedSlots(arrays[i])) / arrays[i].size());

            // epsilon 2 = (abs(A_i+1) - expected batch size) / abs(A_i + 1)
            double epsilon2 = 1.0 - (static_cast<double>(countOccupiedSlots(arrays[i + 1])) / arrays[i + 1].size());

            bool inserted = false;


            if (epsilon1 > delta / 2 && epsilon2 > 0.25) {
                // Goes into first free slot of the positions hash_i(x), ...

                for (int j = 1; j <= f(epsilon1); ++j) {
                    int slot = hashFunction(key, i, j) % arrays[i].size();
                    
                    // Free slot
                    if (arrays[i][slot].first == K(-1)) {
                        arrays[i][slot] = std::make_pair(key, value);
                        inserted = true;
                        break;
                    }
                }
                
                if (!inserted) {
                    // Otherwise goes into first freeslot of hash_i+1(x), ...
                    for (int j = 1; j <= f(epsilon2); ++j) {
                        int slot = hashFunction(key, i, j);

                        // Free slot
                        if (arrays[i + 1][slot].first == K(-1)) {
                            arrays[i + 1][slot] = std::make_pair(key, value);
                            inserted = true;
                            break;
                        }
                    }
                }

            } else if (epsilon1 <= delta / 2) {
                // Goes into first freeslot of hash_i+1(x), ...
                for (int j = 1; j <= f(epsilon1); ++j) {
                    int slot = hashFunction(key, i, j);

                    // Free slot
                    if (arrays[i + 1][slot].first == K(-1)) {
                        arrays[i + 1][slot] = std::make_pair(key, value);
                        inserted = true;
                        break;
                    }
                }
            } else if (epsilon2 <= 0.25) {
                // Goes into first free slot of the positions hash_i(x), ...
                for (int j = 1; j <= f(epsilon2); ++j) {
                    int slot = hashFunction(key, i, j);

                    // Free slot
                    if (arrays[i][slot].first == K(-1)) {
                        arrays[i][slot] = std::make_pair(key, value);
                        inserted = true;
                        break;
                    }
                }
            }
        }
};