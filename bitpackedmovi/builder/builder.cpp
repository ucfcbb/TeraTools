#include<iostream>
#include"fm-index.h"
#include<sdsl/vectors.hpp>
#include<vector>

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

//returns whether an RLBWT represented by a chars int vec and a lens int vec is equivalent to an fmi
bool areEqual(sdsl::int_vector<> chars, sdsl::int_vector<> lens, const rb3_fmi_t& fmi) {
    if (chars.size() != lens.size()) {
        std::cerr << "areEqual called for a symbol array and length array of different sizes" << std::endl;
        exit(1);
    }

    rlditr_t itr;
    rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
    uint64_t currentRun = 0;
    int64_t l = 0;
    int c = 0;
    while (currentRun < chars.size() && (l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
        if ((uint64_t)l != lens[currentRun] || (uint64_t)c != chars[currentRun])
            return false;
        ++currentRun;
    }
    return currentRun == chars.size();
}

struct IntervalPoint {
    /*
    //represents a position in a range [0,n-1] that is composed of x intervals
    //[i_0,i_1-1],[i_1,i_2-1],[i_2,i_3-1],...,[i_{x-1},n-1]
    //a position p in [0,n-1] in this range is represented by 
    //position, interval, offset s.t.
    // - position = p
    // - interval = j s.t. i_j <= p and i_{j+1} > p
    // - offset   = k s.t. i_j + k = p (therefore, k in [0,i_{j+1}-i_j-1]
    */
    uint64_t position, interval, offset;
};
/*
bool operator==(const IntervalPoint& lhs, const IntervalPoint& rhs) {
    return lhs.position == rhs.position && lhs.interval == rhs.interval && lhs.offset == rhs.offset;
}
*/

bool operator!=(const IntervalPoint& lhs, const IntervalPoint& rhs) {
    return lhs.position != rhs.position || lhs.interval != rhs.interval || lhs.offset != rhs.offset;
}

//assumptions:
//no runs of length 0
//the same runlens vector is passed for every call and unmodified
//distance provided added to current position doesn't result in an out-of-bounds IntervalPoint 
//  (except for position = n, the first position after the range
void AdvanceIntervalPoint(IntervalPoint& intPoint, uint64_t distance, const sdsl::int_vector<>& runlens) {
    uint64_t remaining = intPoint.offset + distance;
    while (remaining && remaining >= runlens[intPoint.interval])
        remaining -= runlens[intPoint.interval++];
    intPoint.offset = remaining;
    intPoint.position += distance;
}

//NOTE: interval points returned by mapLF don't have valid position fields, they are set to -1
//assumptions: inputs are valid and correspond to each other
IntervalPoint mapLF(const IntervalPoint& intPoint, 
        const sdsl::int_vector<>& runlens, 
        const sdsl::int_vector<> & toRun, 
        const sdsl::int_vector<> & toOffset) {
    IntervalPoint res;
    res.position = (uint64_t)-1;
    res.interval = toRun[intPoint.interval];
    res.offset = toOffset[intPoint.interval] + intPoint.offset;
    while (runlens[res.interval] <= res.offset)
        res.offset -= runlens[res.interval++];
    return res;
}

void printStructures(
        const sdsl::int_vector<>& rlbwt,
        const sdsl::int_vector<>& runlens,
        const sdsl::int_vector<>& toRun,
        const sdsl::int_vector<>& toOffset) {
    std::cout << "Printing run-length data structures. Format:\n"
        << "\tsymbol\tlength\trun that head of this run maps to with LF\toffset within mapped to run of the LF mapping of head of this run\n";
    uint64_t runs = rlbwt.size();
    if (runs != runlens.size() || runs != toRun.size() || runs != toOffset.size()) {
        std::cerr << "ERROR: length of passed run-length compressed data structures for rlbwt, runlens, and LF don't match.\n"
            << "ERROR:\trlbwt length: " << rlbwt.size() << '\n'
            << "ERROR:\trunlens length: " << runlens.size() << '\n'
            << "ERROR:\ttoRun length: " << toRun.size() << '\n'
            << "ERROR:\ttoOffset length: " << toOffset.size() << '\n';
        exit(1);
    }
    for (uint64_t i = 0; i < runs; ++i) {
        std::cout << '\t' << rlbwt[i]
            << '\t' << runlens[i]
            << '\t' << toRun[i]
            << '\t' << toOffset[i] << '\n';
    }
}

int main(int argc, char *argv[]) {
    Timer.start("builder");
    Timer.start("Program Initialization");
    int use_mmap;
    rb3_fmi_t fmi;
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
        use_mmap = (argc == 3);

        Timer.stop(); //Reading Arguments
        Timer.start((use_mmap)? "Loading fmd with mmap" : "Loading fmd");

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
    }
    Timer.stop(); //Program Initialization

    Timer.start("Constructing RLBWT from FMD");
    sdsl::int_vector<> rlbwt, runlens;
    uint64_t runs = 0, totalLen = 0, alphbits, lenbits;
    uRange alphRange, lenRange;
    std::vector<uint64_t> alphCounts, alphRuns;
    {
        Timer.start("Reading fmd for parameters");

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

        if (alphRange.max == (uint64_t)-1) {
            std::cerr << "Maximum alphabet symbol is 2^64 - 1. This program assumes this is not the case." << std::endl;
            return 1;
        }

        alphCounts.resize(alphRange.max+1);
        alphRuns.resize(alphRange.max+1);
        for (uint64_t i = 0; i <= alphRange.max; ++i)
            alphCounts[i] = alphRuns[i] = 0;

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
            alphCounts[alph] += len;
            ++alphRuns[alph];

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

        std::cout << "Number of runs and number of total occurrences in text for each character is printed below. Format:\n\tsymbol\truns\toccurrences\n";
        for (uint64_t i = 0; i <= alphRange.max; ++i) 
            std::cout << '\t' << i << '\t' << alphRuns[i] << '\t' << alphCounts[i] << '\n';

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

        Timer.start("Verifying correctness of constructed RLBWT");

        std::cout << "Built RLBWT and input FMD are " 
            << ((areEqual(rlbwt, runlens, fmi))? "equal" : "not equal!") << std::endl;
        Timer.stop(); //Verifying correctness of constructed RLBWT

    }
    Timer.stop(); //Constructing RLBWT from FMD

    Timer.start("Constructing LF from RLBWT");
    sdsl::int_vector<> toRun(runs, 0, sdsl::bits::hi(runs) + 1);
    sdsl::int_vector<> toOffset(runs, 0, runlens.width());
    std::vector<IntervalPoint> alphStarts(alphRange.max+2);
    {
        Timer.start("Computing alphStarts");
        {
            alphStarts[0] = {0, 0, 0};
            uint64_t thisAlph = 0, thisAlphLeft = alphCounts[0];
            uint64_t thisRunStart = 0;

            //at beginning of one iteration of this loop,
            //we are at the beginning of a run in the bwt with thisAlphLeft
            //characters of thisAlph to traverse
            for (uint64_t i = 0; i < runs; ++i) {
                uint64_t l = runlens[i];
                //at the beginning of one iteration of this loop, we are at a position
                //within a run (possibly the first) such that the last position in the 
                //F array of thisAlph occurs in this run (run i).
                //!!!!!!!!!!!!!
                //thisAlphLeft is the number of thisAlph characters left to traverse + the
                //number of non thisAlph characters before the thisAlph characters in the 
                //current run
                //!!!!!!!!!!!!!
                while (thisAlphLeft < l) {
                    ++thisAlph;
                    alphStarts[thisAlph] = { thisRunStart + thisAlphLeft, i, thisAlphLeft };
                    thisAlphLeft += alphCounts[thisAlph];
                }

                thisAlphLeft -= l;
                thisRunStart += l;
            }

            if (thisAlph != alphRange.max) {
                std::cerr << "ERROR: When computing alphStarts, didn't end up on largest alphabet value! Largest value: "
                    << alphRange.max << ". Ended up on: " << thisAlph << std::endl;
                return 1;
            }
            if (thisAlphLeft) {
                std::cerr << "ERROR: When computing alphStarts, there are a nonzero amount of "
                    << alphRange.max << " characters left to look for, but the end of the rlbwt has been reached." << std::endl;
                return 1;
            }
            if (thisRunStart != totalLen) {
                std::cerr << "Final position is not total length of bwt in alphStart computation! Final position: " 
                    << thisRunStart << ". Total length: " << totalLen;
                return 1;
            }

            //equivalent to alphStarts[thisAlph] = { thisRunStart + thisAlphLeft, i, thisAlphLeft };
            //since thisAlphLeft = 0
            alphStarts[++thisAlph] = { thisRunStart, runs, 0 };
        }
        Timer.stop(); //Computing alphStarts

        std::cout << "Computed alphStarts printed below."
            << " Note: ENDBWT is not a symbol, but is stored as a sentinel in the alphStarts array at the endFormat:\n"
            << "\tsymbol\tposition\trun\toffset\n";
        for (uint64_t i = 0; i <= alphRange.max; ++i)
            std::cout << '\t' << i << '\t' << alphStarts[i].position
                << '\t' << alphStarts[i].interval
                << '\t' << alphStarts[i].offset << '\n';
        std::cout << '\t' << "ENDBWT" << '\t' << alphStarts[alphRange.max+1].position
            << '\t' << alphStarts[alphRange.max+1].interval
            << '\t' << alphStarts[alphRange.max+1].offset << '\n';

        Timer.start("Computing run and offsets for LF from run starts");
        {
            std::vector<IntervalPoint> currentAlphLFs(alphRange.max+1);
            for (uint64_t i = 0; i <= alphRange.max; ++i) 
                currentAlphLFs[i] = alphStarts[i];

            for (uint64_t i = 0; i < runs; ++i) {
                uint64_t l = runlens[i], c = rlbwt[i];

                toRun[i] = currentAlphLFs[c].interval;
                toOffset[i] = currentAlphLFs[c].offset;

                AdvanceIntervalPoint(currentAlphLFs[c], l, runlens);
            }

            //verify all currentAlphLFs ended up at the start of the next alph
            bool allGood = true;
            for (uint64_t i = 0; i <= alphRange.max; ++i) {
                if (currentAlphLFs[i] != alphStarts[i+1]) {
                    allGood = false;
                    std::cerr << "Symbol " << i << " LF didn't end up at the start of "
                        << " symbol " << i + 1 << " in F array.\nIt ended up at: "
                        << "{ position: " << currentAlphLFs[i].position
                        << ", interval: " << currentAlphLFs[i].interval
                        << ", offset: " << currentAlphLFs[i].offset << " }.\n"
                        << "The next symbol starts at " 
                        << "{ position: " << alphStarts[i+1].position
                        << ", interval: " << alphStarts[i+1].interval
                        << ", offset: " << alphStarts[i+1].offset << " }."
                        << std::endl;
                }
            }

            if (allGood)
                std::cout << "The LFs of all symbols ended up at the start of the next symbol.\n";
            else {
                std::cerr << "ERROR: The LFs of some symbols didn't end up at the start of the next symbol. (See above)." << std::endl;
                return 1;
            }
        }
        Timer.stop(); //Computing run and offsets for LF from run starts
        
        //printStructures(rlbwt, runlens, toRun, toOffset);

        Timer.start("Verifying computed LF by sum of sequence lengths");
        uint64_t lfsDone = 0;
        IntervalPoint start{ (uint64_t)-1, 0, 0}, current{(uint64_t)-1, 0, 0};
        //In theory, this is a multi-dollar BWT. I.E. for strings, s_0, s_1, s_2,..., s_k,
        //the text is s_0,$_0,s_1,$_1,s_2,$_2,...,s_k,$_k concatenated in that order,
        //where $_i < $_j for i < j and $_i < a for all non $_* characters in the text
        //HOWEVER: In practice, the multi-dollars incur a big cost by increasing the alphabet size
        //therefore they are usually all treated as $. This costs us something: we are no longer able
        //to perform LF on the BWT locations with $. If desired, the value i in BWT[x] = $_i must be recovered.
        //then LF[x] = i. We forgo this ability in implementation because we don't care about it and costs index size

        //However this does make verifying the LF function more difficult. (Typically, we could just LF n times
        //and verify that we return to the origin, then LF is a permutation with one cycle and likely correct.)
        //Here instead, we indpentently LF through each string contained in the text, sum up their lengths, 
        //and claim that if the sum of the lengths is accurate then the LF is likely correct. We could
        //also keep a bool array and make sure every position of the BWT is traversed exactly once,
        //but this would be costly in terms of memory for large compressed BWTs (150 GB for human472)
        const uint64_t outputtingInterval = 3e7;
        Timer.start("Timing LFs 0 to " + std::to_string(std::min(totalLen, outputtingInterval - 1)));
        //count the length of sequence seq
        for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
            if (seq) {
                ++start.offset;
                if (start.offset == runlens[start.interval]) {
                    start.offset = 0;
                    ++start.interval;
                }
            }
            current = start;
            do {
                if ((lfsDone+1)%outputtingInterval == 0) {
                    Timer.stop();
                    Timer.start("Timing LFs " + std::to_string(lfsDone+1) + " to " + std::to_string(std::min(totalLen, lfsDone+outputtingInterval)));
                }
                ++lfsDone;

                /*
                   std::cout << "current: { position: " << current.position
                   << ", interval: " << current.interval
                   << ", offset: " << current.offset << " }.\n";
                 */

                current = mapLF(current, runlens, toRun, toOffset);
            } while (lfsDone <= totalLen && rlbwt[current.interval] != 0);
            ++lfsDone; //counting LF to endmarker, which isn't performed in actuality (simulated by incrementing start.offset...)
        }
        Timer.stop();
        /*
        std::cout << "current: { position: " << current.position
            << ", interval: " << current.interval
            << ", offset: " << current.offset << " }.\n";
        */
        
        if (lfsDone != totalLen) {
            std::cerr << "ERROR: Sum of sequence length not equal to the length of the BWT! Sequence length is computed by LF.\n"
                << "ERROR: A total of " << lfsDone << " LFs were computed. This is not equal to the length of the BWT, which is "
                << totalLen << ". NOTE: if lfsDone = length + 1, LF (very likely) was terminated early and didn't compute "
                << "the full sequence lengths. This is because verification is automatically terminated after length + 1"
                << " LFs because the LF function is then known to be incorrect." << std::endl;
            return 1;
        }
        std::cout << "LF is likely correct, the sum of sequence lengths computed by it is " << lfsDone << ". " 
            << lfsDone << " LFs were computed in order to verify this.\n";
        Timer.stop(); //Verifying computed LF by sum of sequence lengths
    }
    Timer.stop(); //Constructing LF from RLBWT

    rb3_fmi_free(&fmi);

    Timer.stop(); //builder
    return 0;
}
