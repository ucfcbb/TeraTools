
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

const char BASES[] = {'A', 'C', 'G', 'T'};
const char COMPLEMENT[] = {'T', 'G', 'C', 'A'};
const int BASES_COUNT = 4;

int main(int argc, char* argv[]) {
    // Seed handling
    unsigned int seed;
    if (argc > 1) {
        seed = std::stoul(argv[1]);
    } else {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 rng(seed);

    // Random number of strings (between 1 and 20)
    std::uniform_int_distribution<int> num_strings_dist(100, 200);
    int num_strings = num_strings_dist(rng);

    // For each string, random length (between 10 and 100), random content
    std::uniform_int_distribution<int> len_dist(1000, 10000);
    std::uniform_int_distribution<int> base_dist(0, BASES_COUNT - 1);

    // std::cout << num_strings << std::endl;
    for (int i = 0; i < num_strings; ++i) {
        int len = len_dist(rng);
        std::string s, comp;
        s.reserve(len), comp.reserve(len);
        for (int j = 0; j < len; ++j) {
            int idx = base_dist(rng);
            s += BASES[idx];
            comp += COMPLEMENT[idx];
        }
        std::reverse(comp.begin(), comp.end());
        std::cout << s << std::endl;
        std::cout << comp << std::endl;
    }

    std::cout << "SEED: " << seed << std::endl;
    std::cout << "NUM_STRINGS: " << num_strings << std::endl;
    return 0;
}
