#include <iostream>
#include"fm-index.h"
#include <sdsl/vectors.hpp>

struct uRange {
    uint64_t min, max;
};

std::ostream& operator<<(std::ostream& os, uRange range) {
    os << '[' << range.min << ',' << range.max << ']';
    return os;
}

int main(int argc, char *argv[]) {
    //const char characters[7] = "$ACGTN";
    if (argc != 2 && argc != 3) {
        std::cerr << "Usage: builder [-mmap] <input.fmd>\n\nIf -mmap is passed, the file is memory mapped before reading, otherwise traditional file io is used\n";
        std::cerr << argc-1 << " arguments passed instead of 1 or 2" << std::endl;
        return 1;
    }
    if (argc == 3 && strcmp(argv[1], "-mmap") !=0) {
        std::cerr << "If two arguments are passed, the first must be \"-mmap\"\n";
        std::cerr << "Currently passed: \"" << argv[1] << "\"\n";
        return 1;
    }
    char* inputfmd = argv[argc-1];
    int use_mmap = (argc == 3);

    rb3_fmi_t fmi;
    rb3_fmi_restore(&fmi, inputfmd, use_mmap);
    if (fmi.e == 0 && fmi.r == 0) {
        std::cerr << "ERROR: failed to load fmd from index file " << inputfmd << std::endl;
        return 1;
    }

    uint64_t runs = 0, alphbits, maxalph, lenbits = 1, maxLen = (1 << lenbits) - 1, totalLen = 0;
    uint64_t alph = 0, len = 0;
    uRange alphRange, lenRange;

    if (!fmi.e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return 1;
    }

    alphbits = (uint64_t)fmi.e->abits;
    maxalph = ((alphbits == 64)? 0 : (1 << alphbits)) - 1;
    std::cout << "Number of bits per symbol in rlbwt: " << alphbits << std::endl;

    int64_t l;
    int c = 0;
    rlditr_t itr;
    rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?

    sdsl::int_vector<> rlbwt = sdsl::int_vector<>(1, 0, alphbits);
    sdsl::int_vector<> runlens = sdsl::int_vector<>(1, 0, lenbits);

    auto setAlphLenCheckBounds = [&alph, &maxalph, &alphbits, &len, &maxLen, &lenbits, &runlens, &c, &l, &runs] () {
        alph = (uint64_t)c;
        len = (uint64_t)l;
        if (alph > maxalph) {
            std::cerr << "Recieved symbol larger than maximum alphabet symbol. Recieved " << alph 
                << ". Maximum possible: " << maxalph 
                << ". According to " << alphbits << " bits per alphabet symbol" << std::endl;
            exit(1);
        }
        while (len > maxLen) {
            std::cout << "Recieved symbol larger than current maximum length. Recieved " << len
                << ". Current run index: " << runs << std::endl;;
            std::cout << "Expanding lenbits from " << lenbits << " to " << lenbits*2 << std::endl;
            lenbits *= 2;
            maxLen = ((lenbits == 64)? 0 : (1 << lenbits)) - 1;
            sdsl::util::expand_width(runlens, lenbits);
        }
    };

    if ((l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
        //std::cout << c << ' ' << l << std::endl;
        setAlphLenCheckBounds();
        alphRange = {alph,alph};
        lenRange = {len,len};
        totalLen += len;
        ++runs;
    }
    else {
        std::cerr << "Failed to read first run's character and length" << std::endl;
        return 1;
    }

    while ((l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
        //std::cout << c << ' ' << l << std::endl;
        if (runs == rlbwt.size()) {
            std::cout << "Resizing sdsl::int_vectors from " << runs << " elements (runs) to "
                << runs*2 << " elements (runs).\n";
            rlbwt.resize(runs*2);
            runlens.resize(runs*2);
        }

        setAlphLenCheckBounds();

        alphRange.min = std::min(alphRange.min, alph);
        alphRange.max = std::max(alphRange.max, alph);
        lenRange.min = std::min(lenRange.min, len);
        lenRange.max = std::max(lenRange.max, len);
        totalLen += len;

        rlbwt[runs] = alph;
        runlens[runs] = len;

        ++runs;
    }
    //std::cerr << "Final BWT:\n" << s << std::endl;


    std::cout << "Number of runs: " << runs 
        << "\nNumber of bits per symbol in rlbwt: " << alphbits 
        << "\nNumber of bits per run for encoding length: " << lenbits
        << "\nCurrent maximum possible length of a run: " << maxLen 
        << std::endl;


    std::cout << "Alphabet range: " << alphRange << "\nRun Lengths range: " << lenRange << std::endl;
    std::cout << "Total BWT length: " << totalLen << std::endl;

    std::cout << "Shrinking rlbwt to size" << std::endl;
    std::cout << "Previous element width in bits: " << (int)rlbwt.width() << std::endl;
    sdsl::util::bit_compress(rlbwt);
    std::cout << "New element width in bits: " << (int)rlbwt.width() << std::endl;
    std::cout << "Previous length in number of elements: " << rlbwt.size() << std::endl;
    std::cout << "Previous capacity by number of bits: " << rlbwt.capacity() << std::endl;
    rlbwt.resize(runs);
    std::cout << "New length in number of elements: " << rlbwt.size() << std::endl;
    std::cout << "New capacity by number of bits: " << rlbwt.capacity() << std::endl;
    std::cout << std::endl;
    std::cout << "Shrinking runlens to size" << std::endl;
    std::cout << "Previous element width in bits: " << (int)runlens.width() << std::endl;
    sdsl::util::bit_compress(runlens);
    std::cout << "New element width in bits: " << (int)runlens.width() << std::endl;
    std::cout << "Previous length in number of elements: " << runlens.size() << std::endl;
    std::cout << "Previous capacity by number of bits: " << runlens.capacity() << std::endl;
    runlens.resize(runs);
    std::cout << "New length in number of elements: " << runlens.size() << std::endl;
    std::cout << "New capacity by number of bits: " << runlens.capacity() << std::endl;

    std::cout << "Final rlbwt size in bytes: " << sdsl::size_in_bytes(rlbwt) << std::endl;
    std::cout << "Final runlens size in bytes: " << sdsl::size_in_bytes(runlens) << std::endl;
    return 0;
}
