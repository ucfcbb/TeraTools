
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

const char BASES[] = {'A', 'C', 'G', 'T', 'N'};
const int BASES_COUNT = 5;

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
    std::cout << "SEED: " << seed << std::endl;

    // Random number of strings (between 1 and 20)
    std::uniform_int_distribution<int> num_strings_dist(1, 20);
    int num_strings = num_strings_dist(rng);
    std::cout << "NUM_STRINGS: " << num_strings << std::endl;

    // For each string, random length (between 10 and 100), random content
    std::uniform_int_distribution<int> len_dist(10, 100);
    std::uniform_int_distribution<int> base_dist(0, BASES_COUNT - 1);

    for (int i = 0; i < num_strings; ++i) {
        int len = len_dist(rng);
        std::string s;
        s.reserve(len);
        for (int j = 0; j < len; ++j) {
            s += BASES[base_dist(rng)];
        }
        std::cout << s << std::endl;
    }

    return 0;
}
