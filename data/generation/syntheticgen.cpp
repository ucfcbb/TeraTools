
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
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <num_strings_low> <num_strings_high> <len_low> <len_high> [seed]"
                  << std::endl;
        return 1;
    }
    unsigned int seed;
    int num_strings_low = std::stoi(argv[1]);
    int num_strings_high = std::stoi(argv[2]);
    int len_low = std::stoi(argv[3]);
    int len_high = std::stoi(argv[4]);
    if (argc > 5) {
        seed = std::stoul(argv[5]);
    } else {
        std::random_device rd;
        seed = rd();
    }
    std::mt19937 rng(seed);

    std::uniform_int_distribution<int> num_strings_dist(num_strings_low, num_strings_high);
    int num_strings = num_strings_dist(rng);
    std::uniform_int_distribution<int> len_dist(len_low, len_high);
    std::uniform_int_distribution<int> base_dist(0, BASES_COUNT - 1);

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

    std::cout << num_strings_low << " " << num_strings_high << " " << len_low << " "
              << len_high << " " << seed << std::endl;
    return 0;
}
