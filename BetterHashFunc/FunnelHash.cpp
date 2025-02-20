#include <stdio.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <cmath>
#include <bitset>
#include <functional>

template <typename K, typename V>

class FunnelHashTable {
    private:
    double delta;
    int alpha; 
    int beta;
    std::vector<std::vector<std::vector<std::pair<K, V>>>> arrays;

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
    // Addume delta <= 1/8
    FunnelHashTable(double d, int n): delta(d) {
        // a = ceil(4 log d^-1 + 10)
        alpha = std::ceil(4 * std::log2(1 / delta) + 10);
        beta = std::ceil(2 * std::log2(1 / delta));

        // Create a + 1 arrays
        arrays.resize(alpha + 1)

        // Determine special array size
        int specialArraySize = std::max(std::floor(3 * delta * n / 4), std::ceil(delta * n / 2));
        arrays[alpha].resize(specialArraySize);
        
        // Determine the size of A' such that it is divisible by b
        int aPrimeSize = n - specialArraySize

        while (aPrimeSize % beta != 0) {
            aPrimeSize --;
        }

        // Readjust special array size
        int specialArraySize = n - aPrimeSize;
        arrays[alpha].resize(specialArraySize);

        // Initialize each slot with a sentinal value
        for (auto& slot: arrays[alpha]) {
            slot = std::make_pair(K(-1), V());
        }
        
        // Split A' into alpha arrays geometrically decreasing with rate 3/4
        int curSize = aPrimeSize;
        for (int i = 0; i < alpha; ++i) {
            arrays[i].resize(curSize);
             // Divide each A_i into arrays of A_i_j, each of size beta
            for (int j = 0; j < arrays[i].size(); ++j) {
                arrays[i][j].resize(beta);

                // Initialize each slot with a sentinal value
                for (auto& slot: arrays[i][j]) {
                    slot = std::make_pair(K(-1), V());
                }
            }
            curSize = static_cast<int>(curSize * 3 / 4);
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
        for (int j = 0; j < arrays[i].size(); ++j) {
            // Hash the key to determine the subarray index
            int subarrayIndex = hashFunction(key, i , j);

            // Iterate through A_i_j for an empty slot and insert into first empty slot
            for (k = 0; k < beta; ++k) {
                if (arrays[i][subarrayIndex][k] == K(-1)) {
                    arrays[i][subarrayIndex] = std::make_pair(key, value);
                    // return on success
                    return;
                }
            }
        }
        // Upon failure, try special array
        for (int k = 0; k < arrays[alpha].size(); ++k) {
            if (arrays[i][alpha][k] == K(-1)) {
                arrays[i][subarrayIndex] = std::make_pair(key, value);
                return;
            }
        }

        // If all slots are full, the insertion fails
        std::cerr << "Insertion failed: No available slots for key " << key << std::endl;
    }
}
