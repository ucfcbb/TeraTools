#include <iostream>
#include"fm-index.h"
#include <sdsl/vectors.hpp>

#include<chrono>
#include<stack>
#include<string>

class Timer {
    std::stack<std::pair<std::chrono::time_point<std::chrono::high_resolution_clock>, std::string>> processes; 
    const char scopeChar;
    const std::string prefix;
    std::ostream& os;
    public:

    void start(std::string processName) {
        os << prefix << std::string(1+processes.size(), scopeChar) << "Starting '" << processName << "'\n";
        processes.emplace(std::chrono::high_resolution_clock::now(), processName);
    }
    void stop() {
        auto end = std::chrono::high_resolution_clock::now();
        auto beginNamePair = processes.top();
        processes.pop();
        os << prefix << std::string(1+processes.size(), scopeChar) << "Ending '" << beginNamePair.second 
            << "'. It took " << std::chrono::duration<double>(end - beginNamePair.first).count() << " seconds\n";
    }
    void stopAllProcesses() {
        while (processes.size())
            stop();
    }

    Timer() = delete;
    Timer(const char c, const std::string& pref, std::ostream& s): scopeChar(c), prefix(pref), os(s) {}
    ~Timer() {
        stopAllProcesses();
    }
}Timer('|', "Timer:", std::cout);

struct uRange {
    uint64_t min, max;
};

std::ostream& operator<<(std::ostream& os, uRange range) {
    os << '[' << range.min << ',' << range.max << ']';
    return os;
}

int main(int argc, char *argv[]) {
    Timer.start("Reading RLBWT from FMD");
    sdsl::int_vector<> rlbwt, runlens;
    {
        Timer.start("Reading Arguments");
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

        Timer.stop(); //Reading Arguments
        Timer.start((use_mmap)? "Loading fmd with mmap" : "Loading fmd");

        rb3_fmi_t fmi;
        rb3_fmi_restore(&fmi, inputfmd, use_mmap);
        if (fmi.e == 0 && fmi.r == 0) {
            std::cerr << "ERROR: failed to load fmd from index file " << inputfmd << std::endl;
            return 1;
        }

        Timer.stop(); //(use_mmap)? "Loading fmd with mmap" : "Loading fmd"

        if (!fmi.e) {
            std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
            return 1;
        }

        Timer.start("Reading fmd for parameters");
        uint64_t runs = 0, totalLen = 0, alphbits, lenbits;
        uRange alphRange, lenRange;

        rlditr_t itr1;
        rld_itr_init(fmi.e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
        int64_t l;
        int c = 0;


        if ((l = rld_dec(fmi.e, &itr1, &c, 0)) > 0) {
            alphRange = {(uint64_t)c,(uint64_t)c};
            lenRange = {(uint64_t)l,(uint64_t)l};
            totalLen += (uint64_t)l;
            ++runs;
        }
        else {
            std::cerr << "Failed to read first run's character and length" << std::endl;
            return 1;
        }

        while ((l = rld_dec(fmi.e, &itr1, &c, 0)) > 0) {
            alphRange.min = std::min(alphRange.min, (uint64_t)c);
            alphRange.max = std::max(alphRange.max, (uint64_t)c);
            lenRange.min = std::min(lenRange.min, (uint64_t)l);
            lenRange.max = std::max(lenRange.max, (uint64_t)l);

            totalLen += (uint64_t)l;
            ++runs;
        }

        alphbits = sdsl::bits::hi(alphRange.max) + 1;
        lenbits = sdsl::bits::hi(lenRange.max) + 1;
        if (alphbits != (uint64_t)fmi.e->abits) 
            std::cout << "WARNING: computed bits per symbol not equal to bits used in fmd. Computed: " 
                << alphbits << ", ropebwt3: " << (uint64_t)fmi.e->abits << std::endl;

        std::cout << "Number of runs: " << runs 
            << "\nNumber of bits per symbol in rlbwt: " << alphbits 
            << "\nNumber of bits per run for encoding length: " << lenbits
            << std::endl;


        std::cout << "Alphabet range: " << alphRange << "\nRun Lengths range: " << lenRange << std::endl;
        std::cout << "Total BWT length: " << totalLen << std::endl;

        Timer.stop(); //Reading fmd for parameters
        Timer.start("Reading fmd into sdsl");

        uint64_t alph = 0, len = 0, newTotalLen = 0, newRuns = 0;


        rlditr_t itr;
        rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?

        rlbwt = sdsl::int_vector<>(runs, 0, alphbits);
        runlens = sdsl::int_vector<>(runs, 0, lenbits);

        while ((l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
            alph = (uint64_t)c;
            len = (uint64_t)l;
            if (alph < alphRange.min || alph > alphRange.max) {
                std::cerr << "ERROR: Run symbol outside of previously found range, symbol: " << alph << ", Range: " << alphRange << std::endl;
                return 1;
            }
            if (len < lenRange.min || len > lenRange.max) {
                std::cerr << "ERROR: Run length outside of previously found range, length: " << len << ", Range: " << lenRange << std::endl;
                return 1;
            }
            newTotalLen += len;

            rlbwt[newRuns] = alph;
            runlens[newRuns] = len;

            ++newRuns;
        }
        Timer.stop(); //Reading fmd into sdsl

        if (newRuns != runs) {
            std::cerr << "ERROR: Number of runs found is different in first and second reads. First: " << runs << ", second: " << newRuns << std::endl;
            return 1;
        }
        if (newTotalLen != totalLen) {
            std::cerr << "ERROR: Total bwt length found is different in first and second reads. First: " << totalLen << ", second: " << newTotalLen << std::endl;
            return 1;
        }


        /*
        Timer.start("Shrinking sdsl");
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
        Timer.stop(); //Shrinking sdsl
        */

        Timer.start("Computing sdsl size");
        std::cout << "Final rlbwt size in bytes: " << sdsl::size_in_bytes(rlbwt) << std::endl;
        std::cout << "Final runlens size in bytes: " << sdsl::size_in_bytes(runlens) << std::endl;
        Timer.stop(); //Computing sdsl size
        rb3_fmi_free(&fmi);
    }
    Timer.stop(); //Reading RLBWT from FMD
    return 0;
}
