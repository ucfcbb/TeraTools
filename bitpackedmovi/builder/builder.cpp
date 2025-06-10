#include <iostream>
#include <sdsl/vectors.hpp>

struct uRange {
    uint64_t min, max;
};

std::ostream& operator<<(std::ostream& os, uRange range) {
    os << '[' << range.min << ',' << range.max << ']';
    return os;
}

int main() {
    /*
     * recieves from cin: 
     * # runs
     * # bits needed for alphabet
     * run length encoded runs, each run is first the character, then the length
    */
    uint64_t runs, alphbits, maxalph, lenbits = 8, maxLen = (1 << lenbits) - 1;
    std::cin >> runs >> alphbits;
    if (alphbits > 64) {
        std::cerr << "Number of bits per alphabet symbol must be at most 64, recieved " << alphbits << std::endl;
        return 1;
    }
    if (runs == 0) {
        std::cerr << "There must be at least one run, recieved 0 for the number of runs parameter" << std::endl;
        return 1;
    }
    maxalph = ((alphbits == 64)? 0 : (1 << alphbits)) - 1;
    std::cout << "Number of runs: " << runs 
        << "\nNumber of bits per symbol in rlbwt: " << alphbits 
        << "\nNumber of bits per run for encoding length: " << lenbits
        << "\nCurrent maximum possible length of a run: " << maxLen 
        << std::endl;

    uint64_t alph = 0, len = 0, numRun = 0;
    sdsl::int_vector<> rlbwt = sdsl::int_vector<>(runs, 0, alphbits);
    sdsl::int_vector<> runlens = sdsl::int_vector<>(runs, 0, lenbits);
    
    uRange alphRange, lenRange;

    if (std::cin >> alph >> len) {
        alphRange = {alph,alph};
        lenRange = {len,len};
        ++numRun;
    }
    else {
        std::cerr << "Failed to read first run's character and length" << std::endl;
        return 1;
    }

    while (std::cin >> alph >> len) { 
        if (alph > maxalph) {
            std::cerr << "Recieved symbol larger than maximum alphabet symbol. Recieved " << alph 
                << ". Maximum possible: " << maxalph 
                << ". According to " << alphbits << " bits per alphabet symbol" << std::endl;
            return 1;
        }
        if (len > maxLen) {
            std::cout << "Expanding lenbits from " << lenbits << " to " << lenbits*2 << std::endl;
            lenbits *= 2;
            maxLen = ((lenbits == 64)? 0 : (1 << lenbits)) - 1;
            sdsl::util::expand_width(runlens, lenbits);
        }

        alphRange.min = std::min(alphRange.min, alph);
        alphRange.max = std::max(alphRange.max, alph);
        lenRange.min = std::min(lenRange.min, len);
        lenRange.max = std::max(lenRange.max, len);

        rlbwt[numRun] = alph;
        runlens[numRun] = len;
        ++numRun;
    }
    if (numRun != runs) {
        std::cerr << "Expected " << runs << " runs, recieved " << numRun << std::endl;
        return 1;
    }

    std::cout << "Alphabet range: " << alphRange << "\nRun Lengths range: " << lenRange << std::endl;
    return 0;
}
