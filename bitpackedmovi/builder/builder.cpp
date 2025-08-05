#include<iostream>
#include"fm-index.h"
#include<sdsl/vectors.hpp>
#include<vector>

#include<chrono>
#include<stack>
#include<string>
#include<iomanip>

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
        uint64_t thisRunLength = lens[currentRun], thisRunChar = chars[currentRun];
        //reconcatenating endmarker runs for comparison
        while (thisRunChar == 0 && currentRun + 1 < chars.size() && chars[currentRun + 1] == 0) {
            ++currentRun;
            thisRunLength += lens[currentRun];
        }

        if ((uint64_t)l != thisRunLength || (uint64_t)c != thisRunChar)
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

struct InvertibleMoveStructure {
    //suffix array samples at
    sdsl::int_vector<> SeqAt;
    sdsl::int_vector<> PosAt;

    //characters between this sample and the next (including this sample)
    //unnecessary, but makes coding cleaner
    sdsl::int_vector<> IntLen;

    sdsl::int_vector<> AboveToInterval;
    sdsl::int_vector<> AboveToOffset;
    sdsl::int_vector<> BelowToInterval;
    sdsl::int_vector<> BelowToOffset;

    sdsl::int_vector<> AboveLCP;

    IntervalPoint mapPhi(const IntervalPoint& intPoint) {
        IntervalPoint res;
        res.position = (uint64_t) -1;
        res.interval = AboveToInterval[intPoint.interval];
        res.offset = AboveToOffset[intPoint.interval] + intPoint.offset;
        while (IntLen[res.interval] <= res.offset)
            res.offset -= IntLen[res.interval++];
        return res;
    }

    IntervalPoint mapInvPhi(const IntervalPoint& intPoint) {
        IntervalPoint res;
        res.position = (uint64_t) -1;
        res.interval = BelowToInterval[intPoint.interval];
        res.offset = BelowToOffset[intPoint.interval] + intPoint.offset;
        while(IntLen[res.interval] <= res.offset)
            res.offset -= IntLen[res.interval++];
        return res;
    }
};

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
        uint64_t field_width,
        const sdsl::int_vector<>& rlbwt,
        const sdsl::int_vector<>& runlens,
        const sdsl::int_vector<>& toRun,
        const sdsl::int_vector<>& toOffset,
        const sdsl::int_vector<>& SATopRunInt,
        const sdsl::int_vector<>& SABotRunInt,
        const InvertibleMoveStructure& PhiInvPhi,
        const sdsl::int_vector<>& runSampledAt,
        const std::vector<IntervalPoint>& revEquivLF) {
    std::cout << "Printing run-length SA order data structures. Format:\n"
        << R"("runInd"	"symbol"	"length"	"run that head of this run maps to with LF"	"offset within mapped to run of the LF mapping of head of this run"	"interval that suffix of top run maps to in phiinvphi"	"interval that suffix of the bottom of the run maps to in phiinvphi")"
        << '\n';
    uint64_t runs = rlbwt.size();
    if (runs != runlens.size() || runs != toRun.size() || runs != toOffset.size() || runs != SATopRunInt.size() || runs != SABotRunInt.size()) {
        std::cerr << "ERROR: length of passed run-length compressed data structures for rlbwt, runlens, LF, SATopRunInt, and SABotRunInt don't match.\n"
            << "ERROR:\trlbwt length: " << rlbwt.size() << '\n'
            << "ERROR:\trunlens length: " << runlens.size() << '\n'
            << "ERROR:\ttoRun length: " << toRun.size() << '\n'
            << "ERROR:\ttoOffset length: " << toOffset.size() << '\n'
            << "ERROR:\tSATopRunInt length: " << SATopRunInt.size() << '\n'
            << "ERROR:\tSABotRunInt length: " << SABotRunInt.size() << '\n';
        exit(1);
    }
    for (uint64_t i = 0; i < runs; ++i) {
        std::cout << std::setw(field_width) << i 
            << std::setw(field_width) << rlbwt[i]
            << std::setw(field_width) << runlens[i]
            << std::setw(field_width) << toRun[i]
            << std::setw(field_width) << toOffset[i]
            << std::setw(field_width) << SATopRunInt[i]
            << std::setw(field_width) << SABotRunInt[i] << '\n';
    }

    std::cout 
        << "------------------------------------------------------------------------------------------\n\n\n"
        << "------------------------------------------------------------------------------------------\n"
        << "Printing run-length Text order data structures. Format:\n"
        << R"("interval Ind"	"SeqAt"	"PosAt"	"IntLen"	"runSampledAt"	"AboveToInterval"	"AboveToOffset"	"BelowToInterval"	"BelowToOffset"	"RevEquivLF Interval"	"RevEquivLF Offset"	"AboveLCP")" << '\n';
    if (2*runs < PhiInvPhi.SeqAt.size() || PhiInvPhi.SeqAt.size() != PhiInvPhi.PosAt.size() || PhiInvPhi.PosAt.size() != PhiInvPhi.IntLen.size() || PhiInvPhi.IntLen.size() != runSampledAt.size()
            || PhiInvPhi.AboveToInterval.size() != runSampledAt.size()
            || PhiInvPhi.AboveToInterval.size() != PhiInvPhi.AboveToOffset.size()
            || PhiInvPhi.AboveToOffset.size() != PhiInvPhi.BelowToInterval.size()
            || PhiInvPhi.BelowToInterval.size() != PhiInvPhi.BelowToOffset.size()
            || PhiInvPhi.BelowToOffset.size() != revEquivLF.size()) {
        std::cerr << "ERROR: length of passed run-length compressed PhiInvPhi data structures not consistent for SeqAt PosAt and IntLen. They should have equal lengths and length <= 2*runs\n"
            << "ERROR:\truns: " << runs << '\n'
            << "ERROR:\tSeqAt length: " << PhiInvPhi.SeqAt.size() << '\n'
            << "ERROR:\tPosAt length: " << PhiInvPhi.PosAt.size() << '\n'
            << "ERROR:\tIntLen length: " << PhiInvPhi.IntLen.size() << '\n'
            << "ERROR:\trunSampledAt length: " << runSampledAt.size() << '\n'
            << "ERROR:\tAboveToInterval length: " << PhiInvPhi.AboveToInterval.size() << '\n'
            << "ERROR:\tAboveToOffset length: " << PhiInvPhi.AboveToOffset.size() << '\n'
            << "ERROR:\tBelowToInterval length: " << PhiInvPhi.BelowToInterval.size() << '\n'
            << "ERROR:\tBelowToOffset length: " << PhiInvPhi.BelowToOffset.size() << '\n'
            << "ERROR:\trevEquivLF length: " << revEquivLF.size() << '\n'
            << "ERROR:\tAboveLCP length: " << PhiInvPhi.AboveLCP.size() << '\n';
        exit(1);
    }
    for (uint64_t i = 0; i < PhiInvPhi.SeqAt.size(); ++i) {
        std::cout << std::setw(field_width) << i 
            << std::setw(field_width) << PhiInvPhi.SeqAt[i]
            << std::setw(field_width) << PhiInvPhi.PosAt[i]
            << std::setw(field_width) << PhiInvPhi.IntLen[i]
            << std::setw(field_width) << runSampledAt[i]
            << std::setw(field_width) << PhiInvPhi.AboveToInterval[i]
            << std::setw(field_width) << PhiInvPhi.AboveToOffset[i]
            << std::setw(field_width) << PhiInvPhi.BelowToInterval[i]
            << std::setw(field_width) << PhiInvPhi.BelowToOffset[i]
            << std::setw(field_width) << revEquivLF[i].interval
            << std::setw(field_width) << revEquivLF[i].offset 
            << std::setw(field_width) << PhiInvPhi.AboveLCP[i] << '\n';
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

        if (!rb3_fmi_is_symmetric(&fmi)) {
            std::cerr << "ERROR: fmd is not 'symmetric' (does this mean bidirectional?) fmd must be bidirectional" << std::endl;
            return 1;
        }
    }
    Timer.stop(); //Program Initialization

    Timer.start("Constructing RLBWT from FMD");
    sdsl::int_vector<> rlbwt, runlens;

    //'repaired' values. These are the values of our constructed index.
    //They may differ from ropebwt3's values because we split endmarkers into separate runs
    uint64_t runs = 0, totalLen = 0, alphbits, lenbits;
    uRange alphRange, lenRange;
    std::vector<uint64_t> alphCounts, alphRuns;
    {
        Timer.start("Reading fmd for parameters");
        {
            //original ropebwt3 values
            uint64_t RB3_runs = 0, RB3_lenbits;
            uRange RB3_lenRange;



            rlditr_t itr1;
            rld_itr_init(fmi.e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;


            if ((l = rld_dec(fmi.e, &itr1, &c, 0)) > 0) {
                alphRange = {(uint64_t)c,(uint64_t)c};

                RB3_lenRange = {(uint64_t)l,(uint64_t)l};
                lenRange = (c == 0)? uRange{(uint64_t)1, (uint64_t)1} : RB3_lenRange;

                runs += (c == 0)? (uint64_t)l : 1;
                ++RB3_runs;

                totalLen += (uint64_t)l;
            }
            else {
                std::cerr << "Failed to read first run's character and length" << std::endl;
                return 1;
            }

            while ((l = rld_dec(fmi.e, &itr1, &c, 0)) > 0) {
                alphRange.min = std::min(alphRange.min, (uint64_t)c);
                alphRange.max = std::max(alphRange.max, (uint64_t)c);

                RB3_lenRange.min = std::min(RB3_lenRange.min, (uint64_t)l);
                RB3_lenRange.max = std::max(RB3_lenRange.max, (uint64_t)l);
                lenRange.min = std::min(lenRange.min, (uint64_t)((c == 0)? 1 : l));
                lenRange.max = std::max(lenRange.max, (uint64_t)((c == 0)? 1 : l));

                ++RB3_runs;
                runs += (c == 0)? (uint64_t)l : 1;

                totalLen += (uint64_t)l;
            }

            if (alphRange.max == (uint64_t)-1) {
                std::cerr << "Maximum alphabet symbol is 2^64 - 1. "
                    << "This program assumes this is not the case (it can only handle alphabet <= (2^64) - 2." << std::endl;
                return 1;
            }

            std::cout << "INFO: The parameters for our constructed BWT (i.e. #runs, max length, etc.) may be "
                << "different from those of the input (ropebwt3).\nINFO: This is because in our constructed BWT, "
                << "each endmarker is contained in its own run.\n";

            alphbits = sdsl::bits::hi(alphRange.max) + 1;
            RB3_lenbits = sdsl::bits::hi(RB3_lenRange.max) + 1;
            lenbits = sdsl::bits::hi(lenRange.max) + 1;
            if (alphbits != (uint64_t)fmi.e->abits) 
                std::cout << "WARNING: computed bits per symbol not equal to bits used in fmd. Computed: " 
                    << alphbits << ", ropebwt3: " << (uint64_t)fmi.e->abits << std::endl;

            std::cout << "Input number of runs (i.e. before splitting endmarker runs): " << RB3_runs 
                << "\nThis index number of runs (i.e. after splitting endmarker runs): " << runs 
                << "\nNumber of bits per symbol in rlbwt: " << alphbits 
                << "\nInput number of bits per run for encoding length (i.e. before splitting endmarker runs): " << RB3_lenbits
                << "\nThis index number of bits per run for encoding length (i.e. after splitting endmarker runs): " << lenbits
                << std::endl;


            std::cout << "Alphabet range: " << alphRange 
                << "\nInput run lengths range (i.e. before splitting endmarker runs): " << RB3_lenRange
                << "\nThis index run lengths range (i.e. after splitting endmarker runs): " << lenRange << std::endl;
            std::cout << "Total BWT length: " << totalLen << std::endl;
        }
        Timer.stop(); //Reading fmd for parameters
        Timer.start("Reading fmd into sdsl");
        {
            alphCounts.resize(alphRange.max+1, 0);
            alphRuns.resize(alphRange.max+1, 0);
            std::vector<uint64_t> RB3_alphRuns(alphRange.max+1, 0);

            uint64_t alph = 0, len = 0, newTotalLen = 0, newRuns = 0;


            rlditr_t itr;
            rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?

            rlbwt = sdsl::int_vector<>(runs, 0, alphbits);
            runlens = sdsl::int_vector<>(runs, 0, lenbits);

            int64_t l;
            int c = 0;
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

                ++RB3_alphRuns[alph];
                alphRuns[alph] += (alph == 0)? len : 1;

                if (alph == 0) {
                    while (len) {
                        rlbwt[newRuns] = alph;
                        runlens[newRuns] = 1;
                        ++newRuns;
                        --len;
                    }
                }
                else {
                    rlbwt[newRuns] = alph;
                    runlens[newRuns] = len;
                    ++newRuns;
                }
            }

            if (newRuns != runs) {
                std::cerr << "ERROR: Number of runs found is different in first and second reads. First: " << runs << ", second: " << newRuns << std::endl;
                return 1;
            }
            if (newTotalLen != totalLen) {
                std::cerr << "ERROR: Total bwt length found is different in first and second reads. First: " << totalLen << ", second: " << newTotalLen << std::endl;
                return 1;
            }

            std::cout << "Number of runs and number of total occurrences in text for each character is printed below. Format:\n"
                << "\tsymbol\truns in this index\truns in input\toccurrences\n";
            for (uint64_t i = 0; i <= alphRange.max; ++i) 
                std::cout << '\t' << i << '\t' << alphRuns[i] << '\t' << RB3_alphRuns[i] << '\t' << alphCounts[i] << '\n';
        }
        Timer.stop(); //Reading fmd into sdsl

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

        Timer.start("Computing RLBWT size");
        std::cout << "Final rlbwt size in bytes: " << sdsl::size_in_bytes(rlbwt) << std::endl;
        std::cout << "Final runlens size in bytes: " << sdsl::size_in_bytes(runlens) << std::endl;
        Timer.stop(); //Computing sdsl size

        Timer.start("Verifying correctness of constructed RLBWT");

        std::cout << "Built RLBWT and input FMD are " 
            << ((areEqual(rlbwt, runlens, fmi))? "equal" : "not equal!") << std::endl;
        Timer.stop(); //Verifying correctness of constructed RLBWT

    }
    Timer.stop(); //Constructing RLBWT from FMD

    rb3_fmi_free(&fmi);

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
        
        Timer.start("Computing LF size");
        std::cout << "Final toRun size in bytes: " << sdsl::size_in_bytes(toRun) << std::endl;
        std::cout << "Final toOffset size in bytes: " << sdsl::size_in_bytes(toOffset) << std::endl;
        Timer.stop(); //Computing LF size

        /*
        Timer.start("Verifying computed LF by sum of sequence lengths");
        {
            uint64_t sumSeqLengths = 0;
            //This below comment is outdated. It is kept for informational purposes and posterity.
            //We now keep each dollar in a separate run. We still verify by sequence
            //because it is easily parallelizable and complete LF traversal is by far the slowest part of construction.
            //(Naturally, since it is the only O(n) instead of O(r) part. It's not even O(n) because I haven't
            //split the intervals of the move data structure yet.)
            //-----------------------------------------------------------------------------------------------------------
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
            //compute starts
            //-----------------------------------------------------------------------------------------------------------
            std::vector<IntervalPoint> starts(alphCounts[0]);
            IntervalPoint start{ (uint64_t)-1, 0, 0};
            starts[0] = start;
            for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                ++start.offset;
                if (start.offset == runlens[start.interval]) {
                    start.offset = 0;
                    ++start.interval;
                }
                starts[seq] = start;
            }

            //count the length of sequence seq
            #pragma omp parallel for schedule(guided)
            for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                IntervalPoint current{starts[seq]};
                uint64_t seqLen = 1; //counting LF to endmarker, which isn't performed in actuality (simulated by incrementing start.offset...)
                while (rlbwt[current.interval] != 0) {
                    ++seqLen;
                    current = mapLF(current, runlens, toRun, toOffset);
                }
                #pragma omp critical 
                {
                    sumSeqLengths += seqLen;
                    stringStarts[seq] = current;
                }
            }

            if (sumSeqLengths != totalLen) {
                std::cerr << "ERROR: Sum of sequence length not equal to the length of the BWT! Sequence length is computed by LF.\n"
                    << "ERROR: A total of " << sumSeqLengths << " LFs were computed. This is not equal to the length of the BWT, which is "
                    << totalLen << "." << std::endl;
                return 1;
            }
            std::cout << "LF is likely correct, the sum of sequence lengths computed by it is " << sumSeqLengths << ". " 
                << sumSeqLengths - alphCounts[0] << " LFs were computed in order to verify this.\n";
        }
        Timer.stop(); //Verifying computed LF by sum of sequence lengths
        */
    }
    Timer.stop(); //Constructing LF from RLBWT

    Timer.start("SA sampling");
    std::vector<IntervalPoint> stringStarts(alphCounts[0]);

    //suffixes are 0-indexed
    //sequences are 0-indexed

    //suffix array samples at the top of runs
    //the suffix at the top of run i is the SATopRunInt[i] interval in the text order
    sdsl::int_vector<> SATopRunInt;
    //SABotRunInt is not needed after the data structure is built but is needed to build the data strcutre (for sampling)
    sdsl::int_vector<> SABotRunInt;

    sdsl::int_vector<> runSampledAt;

    //Samples in Text order
    InvertibleMoveStructure PhiInvPhi;

    uint64_t maxSeqLen = 0;
    std::vector<uint64_t> seqLens(alphCounts[0]); //seq lengths, counts endmarker so the empty string has length 1
    //number of times each string is at the top (and bottom respectively) of a run, 
    //counts endmarker, so value for the empty string would be 1
    std::vector<uint64_t> seqNumsTopRun(alphCounts[0]), seqNumsBotRun(alphCounts[0]), seqNumsTopOrBotRun(alphCounts[0]); 
    {
        uint64_t sumSeqLengths = 0;
        uint64_t maxIntLen = 0;
        Timer.start("Auxiliary info computation (seqLens, seqNumsTopRun, seqNumsBotRun, seqNumsTopOrBotRun)");
        {
            std::vector<IntervalPoint> starts(alphCounts[0]);
            IntervalPoint start{ (uint64_t)-1, 0, 0};
            starts[0] = start;
            for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                ++start.offset;
                if (start.offset == runlens[start.interval]) {
                    start.offset = 0;
                    ++start.interval;
                }
                starts[seq] = start;
            }

            #pragma omp parallel for schedule(guided)
            for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                IntervalPoint current{starts[seq]};
                uint64_t seqLen = 0, seqNumTopRun = 0, seqNumBotRun = 0, seqNumTopOrBotRun = 0, seqNumOneRun = 0;
                uint64_t currentTopOrBotIntervalLen = 0; //number of characters in the current interval on the sequence since the last time the sequence was at the top or bottom of a run
                uint64_t maxTopOrBotIntervalLen = 0;
                while (rlbwt[current.interval] != 0) {
                    //add suffix at current in SA
                    ++seqLen;
                    ++currentTopOrBotIntervalLen;
                    seqNumOneRun += runlens[current.interval] == 1;
                    seqNumTopRun += (current.offset == 0);
                    seqNumBotRun += (current.offset == runlens[current.interval] - 1);
                    if ((current.offset == 0) || (current.offset == runlens[current.interval] - 1)) {
                        seqNumTopOrBotRun++;
                        maxTopOrBotIntervalLen = std::max(maxTopOrBotIntervalLen, currentTopOrBotIntervalLen);
                        currentTopOrBotIntervalLen = 0;
                    }
                    
                    current = mapLF(current, runlens, toRun, toOffset);
                } 
                //add suffix 0 of seq
                ++seqLen;
                ++currentTopOrBotIntervalLen;
                seqNumOneRun += runlens[current.interval] == 1;
                seqNumTopRun += (current.offset == 0);
                seqNumBotRun += (current.offset == runlens[current.interval] - 1);
                if ((current.offset == 0) || (current.offset == runlens[current.interval] - 1)) {
                    seqNumTopOrBotRun++;
                    maxTopOrBotIntervalLen = std::max(maxTopOrBotIntervalLen, currentTopOrBotIntervalLen);
                    currentTopOrBotIntervalLen = 0;
                }
                else {
                    std::cerr << "ERROR: suffix 0 of seq " << seq << " not at top or bottom of run (should be both, since in a run of length 1!\n";
                    exit(1);
                }

                //save results
                #pragma omp critical
                {
                    sumSeqLengths += seqLen;
                    stringStarts[seq] = current;

                    seqLens[seq] = seqLen;
                    seqNumsTopRun[seq] = seqNumTopRun;
                    seqNumsBotRun[seq] = seqNumBotRun;
                    seqNumsTopOrBotRun[seq] = seqNumTopOrBotRun;
                    
                    maxIntLen = std::max(maxIntLen, maxTopOrBotIntervalLen);

                    if (seqNumTopOrBotRun != seqNumTopRun + seqNumBotRun - seqNumOneRun) {
                        std::cerr << "Number of times at bottom, top, boundary, and number runs of length one not consistent for seq " 
                            << seq << "!\nTop: " << seqNumTopRun
                            << "\nBot: " << seqNumBotRun 
                            << "\nTop or Bot (Boundary): " << seqNumTopOrBotRun
                            << "\nRuns of length one: " << seqNumOneRun << "\n";
                        exit(1);
                    }
                }
            }
        }
        Timer.stop(); //Auxiliary info computation (seqLens, seqNumsTopRun, seqNumsBotRun)

        if (sumSeqLengths != totalLen) {
            std::cerr << "ERROR: Sum of sequence length not equal to the length of the BWT! Sequence length is computed by LF.\n"
                << "ERROR: A total of " << sumSeqLengths << " LFs were computed. This is not equal to the length of the BWT, which is "
                << totalLen << "." << std::endl;
            return 1;
        }
        std::cout << "LF is likely correct, the sum of sequence lengths computed by it is " << sumSeqLengths << ". " 
            << sumSeqLengths - alphCounts[0] << " LFs were computed in order to verify this.\n";

        Timer.start("Prefix summing auxiliary info");
        maxSeqLen = seqLens[0];
        for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
            maxSeqLen = std::max(maxSeqLen, seqLens[seq]);

            seqNumsTopRun[seq] += seqNumsTopRun[seq-1];
            seqNumsBotRun[seq] += seqNumsBotRun[seq-1];
            seqNumsTopOrBotRun[seq] += seqNumsTopOrBotRun[seq-1];
        }
        Timer.stop(); //Prefix summing auxiliary info

        std::cout << "Maximum sequence length: " << maxSeqLen << '\n';
        std::cout << "Maximum interval length: " << maxIntLen << '\n';

        if (seqNumsTopRun.back() != runs) {
            std::cerr << "ERROR:seqNumsTopRun.back() != runs, " << seqNumsTopRun.back() << " != " << runs << '\n';
            exit(1);
        }
        if (seqNumsBotRun.back() != runs) {
            std::cerr << "ERROR: seqNumsBotRun.back() != runs, " << seqNumsBotRun.back() << " != " << runs << '\n';
            exit(1);
        }
        std::cout << "Total run tops: " << seqNumsTopRun.back() << '\n'
            << "Total run bots: " << seqNumsBotRun.back() << '\n'
            << "Total run top and bots: " << seqNumsTopOrBotRun.back() << '\n'
            << "Therefore, total runs of length 1 [(tops + bots) - (tops and bots)]: " << seqNumsTopRun.back() + seqNumsBotRun.back() - seqNumsTopOrBotRun.back() << '\n';

        SATopRunInt = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
        SABotRunInt = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);

        PhiInvPhi.SeqAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(alphCounts[0] - 1) + 1);
        PhiInvPhi.PosAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxSeqLen - 1) + 1);

        PhiInvPhi.IntLen = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen) + 1);

        PhiInvPhi.AboveToInterval = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
        PhiInvPhi.AboveToOffset = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen - 1) + 1);
        PhiInvPhi.BelowToInterval = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
        PhiInvPhi.BelowToOffset = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen - 1) + 1);


        runSampledAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(runs - 1) + 1);

        Timer.start("Sampling");
        {
            std::vector<IntervalPoint> starts(alphCounts[0]);
            IntervalPoint start{ (uint64_t)-1, 0, 0};
            starts[0] = start;
            for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                ++start.offset;
                if (start.offset == runlens[start.interval]) {
                    start.offset = 0;
                    ++start.interval;
                }
                starts[seq] = start;
            }

            Timer.start("Sampling the SA order");
            #pragma omp parallel for schedule(guided)
            for(uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                IntervalPoint current{starts[seq]};
                uint64_t pos = seqLens[seq] - 1;
                uint64_t prevPos = seqLens[seq];
                uint64_t posTopRun = seqNumsTopRun[seq] - 1;
                uint64_t posBotRun = seqNumsBotRun[seq] - 1;
                uint64_t posRun = seqNumsTopOrBotRun[seq] - 1;
                uint64_t finalTopRun = (seq == 0)? 0 : seqNumsTopRun[seq-1];
                uint64_t finalBotRun = (seq == 0)? 0 : seqNumsBotRun[seq-1];
                uint64_t finalRun = (seq == 0)? 0 : seqNumsTopOrBotRun[seq-1];
                while (rlbwt[current.interval] != 0) {
                    if (current.offset == 0 || current.offset == runlens[current.interval] - 1) {
                        if (current.offset == 0) {
                            #pragma omp critical
                            {
                                if (SATopRunInt[current.interval] != 0) {
                                    std::cerr << "ERROR: this run's top sample has already been set!\n";
                                }
                                SATopRunInt[current.interval] = posRun;
                                if (posTopRun == finalTopRun) {
                                    std::cerr << "ERROR: posTopRun to reach last position before endmarker found!\n";
                                }
                            }
                            --posTopRun;
                        }
                        if (current.offset == runlens[current.interval] - 1) {
                            #pragma omp critical
                            {
                                if (SABotRunInt[current.interval] != 0) {
                                    std::cerr << "ERROR: this run's bot sample has already been set!\n";
                                }
                                SABotRunInt[current.interval] = posRun;
                                if (posBotRun == finalBotRun) {
                                    std::cerr << "ERROR: posBotRun to reach last position before endmarker found!\n";
                                }
                            }
                            --posBotRun;
                        }

                        #pragma omp critical
                        {
                            if (PhiInvPhi.SeqAt[posRun] != 0 || PhiInvPhi.PosAt[posRun] != 0) {
                                std::cerr << "ERROR: this interval's phi sample has already been set!\n";
                            }
                            PhiInvPhi.SeqAt[posRun] = seq;
                            PhiInvPhi.PosAt[posRun] = pos;
                            PhiInvPhi.IntLen[posRun] = prevPos - pos;
                            runSampledAt[posRun] = current.interval;
                            if (posRun == finalRun) {
                                std::cerr << "ERROR: posRun to reach last position before endmarker found!\n";
                            }
                        }
                        prevPos = pos;
                        --posRun;

                    }
                    if (pos == 0) {
                        std::cerr << "ERROR: pos in sequence to reach -1 before endmarker found!\n";
                    }
                    --pos;
                    current = mapLF(current,runlens, toRun, toOffset);
                }
                if (pos != 0) {
                    std::cerr << "ERROR: full sequence not traversed!\n";
                }
                if (posTopRun != finalTopRun) {
                    std::cerr << "ERROR: didn't traverse all top samples of this sequence!\n";
                }
                if (posBotRun != finalBotRun) {
                    std::cerr << "ERROR: didn't traverse all bottom samples of this sequence!\n";
                }
                if (posRun != finalRun) {
                    std::cerr << "ERROR: didn't traverse all interval samples of this sequence!\n";
                }
                if (current.offset == 0) {
                    #pragma omp critical
                    {
                        SATopRunInt[current.interval] = posRun;
                    }
                    --posTopRun;

                    PhiInvPhi.SeqAt[posRun] = seq;
                    PhiInvPhi.PosAt[posRun] = pos;
                    PhiInvPhi.IntLen[posRun] = prevPos - pos;
                    runSampledAt[posRun] = current.interval;
                }
                else {
                    std::cerr << "ERROR: offset != 0 in last run, endmarker run!\n";
                }
                if (current.offset == runlens[current.interval] - 1) {
                    #pragma omp critical
                    {
                        SABotRunInt[current.interval] = posRun;
                    }
                    --posBotRun;
                }
                else {
                    std::cerr << "ERROR: offset not at end of run in last run, endmarker run!\n";
                }
            }
            Timer.stop(); //Sampling in SA order
            
            Timer.start("Sampling the Text order");
            //#pragma omp parallel for schedule(guided)
            for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                //go from beginning of sequence to end, by interval in phi
                uint64_t currentIntervalIndex = (seq == 0)? 0 :  seqNumsTopOrBotRun[seq-1];
                uint64_t seqTraversed = 0;

                //do first interval
                if (runlens[runSampledAt[currentIntervalIndex]] != 1) {
                    std::cerr << "ERROR: first interval of sequences starts at a run of length not equal to 1"
                       << " (should be 1 since bwt value should be endmarker of previous sequence)!\n";
                }
                if (PhiInvPhi.SeqAt[currentIntervalIndex] != seq)
                    std::cerr << "ERROR: seq at interval doesn't match seq in Phi!\n";
                if (SATopRunInt[runSampledAt[currentIntervalIndex]] != SABotRunInt[runSampledAt[currentIntervalIndex]])
                    std::cerr << "ERROR: top and bottom run sequence samples don't match in endmarker run (should be length 1)!\n";

                PhiInvPhi.AboveToInterval[currentIntervalIndex] = SABotRunInt[(runSampledAt[currentIntervalIndex] == 0)? SABotRunInt.size()-1 : runSampledAt[currentIntervalIndex]-1];
                PhiInvPhi.AboveToOffset[currentIntervalIndex] = 0;
                PhiInvPhi.BelowToInterval[currentIntervalIndex] = SATopRunInt[(runSampledAt[currentIntervalIndex] == SATopRunInt.size() - 1)? 0 : runSampledAt[currentIntervalIndex]+1];
                PhiInvPhi.BelowToOffset[currentIntervalIndex] = 0;

                //std::cout << "PhiInvPhi.AboveToInterval[currentIntervalIndex]: " << PhiInvPhi.AboveToInterval[currentIntervalIndex] << std::endl;
                //std::cout << "PhiInvPhi.BelowToInterval[currentIntervalIndex]: " << PhiInvPhi.BelowToInterval[currentIntervalIndex] << std::endl;

                seqTraversed += PhiInvPhi.IntLen[currentIntervalIndex];
                ++currentIntervalIndex;

                //std::cout << "Here" << std::endl;
                //std::cout << "seqTraversed: " << seqTraversed << ". currentIntervalIndex: " << currentIntervalIndex << std::endl;

                while (seqTraversed < seqLens[seq]) {
                    //std::cout << "In here" << std::endl;
                    //add next interval
                    uint64_t runIndex = runSampledAt[currentIntervalIndex];
                    uint64_t topRunInt = SATopRunInt[runIndex];
                    uint64_t botRunInt = SABotRunInt[runIndex];
                    if (!((seq == PhiInvPhi.SeqAt[topRunInt] && seqTraversed == PhiInvPhi.PosAt[topRunInt]) ||
                                (seq == PhiInvPhi.SeqAt[botRunInt] && seqTraversed == PhiInvPhi.PosAt[botRunInt])))
                        std::cerr << "ERROR: Beginning of run interval in sequence is not equal to the sample at"
                           << " the beginning or the end of the corresponding interval!\n";
                    //std::cout << "AFter this?" << std::endl;
                    //computing above sample
                    if (seq == PhiInvPhi.SeqAt[topRunInt] && seqTraversed == PhiInvPhi.PosAt[topRunInt]) {
                        uint64_t runAboveIndex = (runIndex == 0)? runs - 1 : runIndex - 1;
                        #pragma omp critical
                        {
                            PhiInvPhi.AboveToInterval[currentIntervalIndex] = SABotRunInt[runAboveIndex];
                            PhiInvPhi.AboveToOffset[currentIntervalIndex] = 0;
                        }
                    }
                    else {
                        //use previous above sample
                        uint64_t prevInt = PhiInvPhi.AboveToInterval[currentIntervalIndex - 1];
                        uint64_t prevOffset = PhiInvPhi.AboveToOffset[currentIntervalIndex - 1];
                        prevOffset += PhiInvPhi.IntLen[currentIntervalIndex - 1];
                        while (prevOffset >= PhiInvPhi.IntLen[prevInt]) {
                            prevOffset -= PhiInvPhi.IntLen[prevInt];
                            ++prevInt;
                        }
                        #pragma omp critical
                        {
                            PhiInvPhi.AboveToInterval[currentIntervalIndex] = prevInt;
                            PhiInvPhi.AboveToOffset[currentIntervalIndex] = prevOffset;
                        }
                    }

                    //std::cout << "How about this?" << std::endl;

                    //computing below sample
                    if (seq == PhiInvPhi.SeqAt[botRunInt] && seqTraversed == PhiInvPhi.PosAt[botRunInt]) {
                        //std::cout << "In if" << std::endl;
                        uint64_t runBelowIndex = (runIndex == runs - 1)? 0 : runIndex + 1;
                        #pragma omp critical
                        {
                            PhiInvPhi.BelowToInterval[currentIntervalIndex] = SATopRunInt[runBelowIndex];
                            PhiInvPhi.BelowToOffset[currentIntervalIndex] = 0;
                        }
                    }
                    else {
                        //std::cout << "In else" << std::endl;
                        //use previous below sample
                        uint64_t prevInt = PhiInvPhi.BelowToInterval[currentIntervalIndex - 1];
                        uint64_t prevOffset = PhiInvPhi.BelowToOffset[currentIntervalIndex - 1];
                        prevOffset += PhiInvPhi.IntLen[currentIntervalIndex - 1];
                        //std::cout << "Got here" << std::endl;
                        while (prevOffset >= PhiInvPhi.IntLen[prevInt]) {
                            //std::cout << "prevInt: " << prevInt << ". prevOffset: " << prevOffset << std::endl;
                            prevOffset -= PhiInvPhi.IntLen[prevInt];
                            ++prevInt;
                        }
                        //std::cout << "prevInt: " << prevInt << ". prevOffset: " << prevOffset << std::endl;
                        //std::cout << "passed while loop" << std::endl;
                        #pragma omp critical
                        {
                            PhiInvPhi.BelowToInterval[currentIntervalIndex] = prevInt;
                            PhiInvPhi.BelowToOffset[currentIntervalIndex] = prevOffset;
                        }
                    }

                    //std::cout << "How aboutttt this?" << std::endl;

                    seqTraversed += PhiInvPhi.IntLen[currentIntervalIndex];
                    ++currentIntervalIndex;
                }

                //std::cout << "Here after" << std::endl;

                #pragma omp critical
                if (seqTraversed != seqLens[seq]) {
                    std::cerr << "ERROR: Traversed a sequence some length not equal to it's actual length!\n";
                    std::cerr << "ERROR: Traversed seq " << seq << " " << seqTraversed << " characters. Actual length: " << seqLens[seq] << '\n';
                }
                #pragma omp critical
                if (currentIntervalIndex != seqNumsTopOrBotRun[seq]) {
                    std::cerr << "ERROR: Traversed sequence but ended up on some interval in PhiInvPhi other than the first interval of the next sequence!\n";
                    std::cerr << "ERROR: Ended up on interval " << currentIntervalIndex << ". Should have ended up on " << seqNumsTopOrBotRun[seq] << '\n';
                }
            }
            Timer.stop(); //Sampling in Textorder
        }
        Timer.stop(); //Sampling
    }
    Timer.stop(); //SA sampling
    
    Timer.start("Verifying Phi");
    {
        bool pass = true;
        #pragma omp parallel for schedule(guided)
        for (uint64_t run = 0; run < runs; ++run) {
            //check if bottom of run phis to bottom of next run
            IntervalPoint start{uint64_t(-1), SABotRunInt[run], 0},
                          end{uint64_t(-1), SABotRunInt[(run == 0)? SABotRunInt.size() - 1: (run - 1)], 0};
            for(uint64_t ops = 0; ops < runlens[run]; ++ops)
                start = PhiInvPhi.mapPhi(start);
            if (start != end) {
                #pragma omp critical
                {
                    pass = false;
                    std::cerr << "ERROR: suffix did not end up at the end of previous run in verification of phi when phiing runlen times from the bottom of a run!\n";
                }
            }
        }
        if (pass)
            std::cout << "Phi move data structure is likely correct, it is a permutation with one cycle\n";
    }
    Timer.stop(); //Verifying Phi
    
    Timer.start("Verifying InvPhi");
    {
        bool pass = true;
        #pragma omp parallel for schedule(guided)
        for (uint64_t run = 0; run < runs; ++run) {
            //check if bottom of run phis to bottom of next run
            IntervalPoint start{uint64_t(-1), SATopRunInt[run], 0},
                          end{uint64_t(-1), SATopRunInt[(run == SABotRunInt.size()-1)? 0 : (run + 1)], 0};
            for(uint64_t ops = 0; ops < runlens[run]; ++ops)
                start = PhiInvPhi.mapInvPhi(start);
            if (start != end) {
                #pragma omp critical
                {
                    pass = false;
                    std::cerr << "ERROR: suffix did not end up at the end of previous run in verification of phi when phiing runlen times from the bottom of a run!\n";
                }
            }
        }
        if (pass)
            std::cout << "InvPhi move data structure is likely correct, it is a permutation with one cycle\n";
    }
    Timer.stop(); //Verifying InvPhi
            
    Timer.start("LCP computation");
    {
        std::vector<IntervalPoint> starts(alphCounts[0]);
        std::vector<IntervalPoint> revEquivLF(seqNumsTopOrBotRun.back());
        IntervalPoint start{ (uint64_t)-1, 0, 0};
        starts[0] = start;
        for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
            ++start.offset;
            if (start.offset == runlens[start.interval]) {
                start.offset = 0;
                ++start.interval;
            }
            starts[seq] = start;
        }

        Timer.start("Reverse sampling");
        #pragma omp parallel for schedule(guided)
        for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
            //traverse seq from end to beginning, storing samples of LF position of seq in the intervals of rev(seq) in the PhiInvPHi data structure
            uint64_t revSeq = (seq%2 == 0)? seq + 1 : seq - 1;
            IntervalPoint current{starts[seq]};
            //uint64_t posSeq = seqLens[seq] - 1;
            uint64_t revSeqIntervalIndex = (revSeq == 0)? 0 : seqNumsTopOrBotRun[revSeq - 1];
            uint64_t revSeqIntervalOffset = 0;
            uint64_t finalRevSeqIntervalIndex = seqNumsTopOrBotRun[revSeq];

            //check if seq and revSeq lengths are equal
            #pragma omp critical
            if (seqLens[seq] != seqLens[revSeq]) {
                std::cerr << "ERROR: Length of sequence (" << seq << ") and its reverse (" << revSeq << ") are not equal! Respectively: " 
                    << seqLens[seq] << " and " << seqLens[revSeq] << "!\n";
            }

            while (revSeqIntervalIndex != finalRevSeqIntervalIndex) {
                if (revSeqIntervalOffset == 0) {
                    #pragma omp critical
                    {
                        if (revEquivLF[revSeqIntervalIndex].interval != 0 || revEquivLF[revSeqIntervalIndex].offset != 0)
                            std::cerr << "ERROR: setting already set revEquivLF in reverse sampling!\n";
                        revEquivLF[revSeqIntervalIndex] = current;
                    }
                }

                current = mapLF(current, runlens, toRun, toOffset);
                ++revSeqIntervalOffset;
                if (revSeqIntervalOffset == PhiInvPhi.IntLen[revSeqIntervalIndex]) {
                    revSeqIntervalOffset = 0;
                    ++revSeqIntervalIndex;
                }
            }
        }
        Timer.stop(); //Reverse sampling

        Timer.start("LCP Sampling");
        {
            PhiInvPhi.AboveLCP.resize(seqNumsTopOrBotRun.back());

            Timer.start("locating alphStarts in PhiInvPhi");
            uint64_t sampledAlphStarts = 0;
            std::vector<IntervalPoint> phiAlphStarts(alphStarts.size()-1);
            for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
                //std::cout << "Hre" << std::endl;
                phiAlphStarts[al] = {uint64_t(-1), SATopRunInt[alphStarts[al].interval], 0};
                for (uint64_t i = 0; i < alphStarts[al].offset; ++i)
                    phiAlphStarts[al] = PhiInvPhi.mapInvPhi(phiAlphStarts[al]);
                sampledAlphStarts += phiAlphStarts[al].offset == 0;
            }
            Timer.stop(); //locating alphStarts in PhiInvPhi
            //std::cout << "Hredone" << std::endl;

            std::cout << "Of " << phiAlphStarts.size() << ", " << sampledAlphStarts << " alphStarts are at the beginning or end of a run in the BWT\n";
            std::cout << "Locations of alphStarts in PhiInvPhi intervals:\n\tsymbol\tinterval\toffset\n";
            for (uint64_t al = 0; al < phiAlphStarts.size(); ++al)
                std::cout << '\t' << al
                    << '\t' << phiAlphStarts[al].interval
                    << '\t' << phiAlphStarts[al].offset << '\n';

            
            std::vector<bool> alphStartFound(phiAlphStarts.size());
            for (uint64_t al = 0; al < alphStartFound.size(); ++al)
                alphStartFound[al] = phiAlphStarts[al].offset != 0;
            //compute only for those with AboveToOffset = 0
            #pragma omp parallel for schedule(guided)
            for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                uint64_t startInterval = (seq == 0)? 0 : seqNumsTopOrBotRun[seq - 1];
                for (uint64_t currentInterval = startInterval; currentInterval < seqNumsTopOrBotRun[seq]; ++currentInterval) {
                    if (PhiInvPhi.AboveToOffset[currentInterval] != 0) {
                        #pragma omp critical
                        {
                            PhiInvPhi.AboveLCP[currentInterval] = PhiInvPhi.AboveLCP[currentInterval -1] - PhiInvPhi.IntLen[currentInterval - 1];
                            if (PhiInvPhi.AboveLCP[currentInterval] < PhiInvPhi.IntLen[currentInterval]){
                                bool found = false;
                                if (PhiInvPhi.AboveLCP[currentInterval] == 0 && PhiInvPhi.IntLen[currentInterval] == 1) {
                                    //possibly at the beginning of alphabet in F column, LCP of 0 is correct
                                    for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
                                        if (phiAlphStarts[al].interval != currentInterval || phiAlphStarts[al].offset != 0)
                                            continue;
                                        if (alphStartFound[al])
                                            std::cerr << "ERROR: already found this alplhStart in phiinv data structure!\n";
                                        alphStartFound[al] = true;
                                        found = true;
                                    }
                                }
                                if (currentInterval == seqNumsTopOrBotRun[seq] - 1) {
                                    //LCP should be interval length - 1
                                    if (PhiInvPhi.AboveLCP[currentInterval] != PhiInvPhi.IntLen[currentInterval] - 1) {
                                        std::cerr << "ERROR: Last interval of sequence has weird matching length, maybe including endmarker!\n";
                                    }
                                    else 
                                        found = true;
                                }
                                if (!found) 
                                    std::cerr << "ERROR: forwarded LCP computed smaller than interval length and is not one of alphStarts!\n"
                                        << "ERROR: forwarded: " << PhiInvPhi.AboveLCP[currentInterval] << ". Interval length: " << PhiInvPhi.IntLen[currentInterval] << '\n';
                            }
                        }
                        continue;
                    }

                    uint64_t matchingLength = 0;
                    IntervalPoint revSeq{revEquivLF[currentInterval]};
                    IntervalPoint revSeqAbove{revEquivLF[PhiInvPhi.AboveToInterval[currentInterval]]};

                    while (rlbwt[revSeq.interval] != 0 && rlbwt[revSeqAbove.interval] != 0 && rlbwt[revSeq.interval] == rlbwt[revSeqAbove.interval]) {
                        revSeq = mapLF(revSeq, runlens, toRun, toOffset);
                        revSeqAbove = mapLF(revSeqAbove, runlens, toRun, toOffset);
                        ++matchingLength;
                    }

                    #pragma omp critical
                    {
                        PhiInvPhi.AboveLCP[currentInterval] = matchingLength;
                        if (matchingLength < PhiInvPhi.IntLen[currentInterval]) {
                            bool found = false;
                            if (matchingLength == 0 && PhiInvPhi.IntLen[currentInterval] == 1) {
                                //possibly at the beginning of alphabet in F column, LCP of 0 is correct
                                for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
                                    if (phiAlphStarts[al].interval != currentInterval || phiAlphStarts[al].offset != 0)
                                        continue;
                                    if (alphStartFound[al])
                                        std::cerr << "ERROR: already found this alphStart in the PhiInvPhi data structure!\n";
                                    alphStartFound[al] = true;
                                    found = true;
                                }
                            }
                            if (currentInterval == seqNumsTopOrBotRun[seq] - 1) {
                                //LCP should be interval length - 1
                                if (matchingLength !=  PhiInvPhi.IntLen[currentInterval] - 1) {
                                    std::cerr << "ERROR: Last interval of sequence has weird matching that might include the endmarker!\n";
                                }
                                else 
                                    found = true;
                            }
                            if (!found)
                                std::cerr << "ERROR: Computed LCP smaller than interval length and is not one of alphStarts!\n"
                                    << "ERROR: Computed: " << matchingLength <<". Interval length: " << PhiInvPhi.IntLen[currentInterval] << '\n';
                        }
                    }
                }
            }

            bool re = true;
            for (auto t : alphStartFound)
                re = re && t;
            if (re)
                std::cout << "Found all " << sampledAlphStarts << " sampled alplh starts.\n";
            else 
                std::cerr << "ERROR: Didn't find all alph starts!\n";

            Timer.stop(); //LCP Sampling
            /*
            Timer.start("Printing Structures");
            printStructures(10, rlbwt, runlens, toRun, toOffset, SATopRunInt, SABotRunInt, PhiInvPhi, runSampledAt, revEquivLF);
            Timer.stop(); //Printing Structures
            */
        }
    }
    Timer.stop(); //LCP computation

    Timer.start("Shrinking LCP");
    std::cout << "Shrinking AboveLCP to size\n";
    std::cout << "Previous element width in bits: " << (int)PhiInvPhi.AboveLCP.width() << '\n';
    sdsl::util::bit_compress(PhiInvPhi.AboveLCP);
    std::cout << "New element width in bits: " << (int)PhiInvPhi.AboveLCP.width() << '\n';
    Timer.stop(); //Shrinking LCP
    
    Timer.start("RLBWT Repair");
    {
        Timer.start("Detecting endmarkers in runs in RLBWT");
        //set of runIDs that are runs of multiple endmarkers in the RLBWT
        for (const auto& intPoint : stringStarts) {
            if (intPoint.offset) {
                std::cerr << "ERROR: There is a run of endmarkers in the RLBWT with length > 1!\n";
                return 1;
            }
        }
        std::cout << "Every endmarker in the RLBWT is in a run of length 1.\n";
        Timer.stop(); //Detecting endmarkers in runs in RLBWT

        Timer.start("Correcting LFs of endmarkers");
        toRun[stringStarts[0].interval] = alphCounts[0] - 1;
        toOffset[stringStarts[0].interval] = 0;
        for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
            toRun[stringStarts[seq].interval] = seq - 1;
            toOffset[stringStarts[seq].interval] = 0;
        }
        Timer.stop(); //Correcting LFs of endmarkers

        /*
        //This code is correct, but is commented out because it takes very long for BWTs with large n (O(n) non-constant LF calls)
        Timer.start("Testing all LFs by permutation with one cycle");
        IntervalPoint start{ (uint64_t)-1, 0, 0}, current{ (uint64_t)-1, 0, 0};
        uint64_t lfsDone = 0;
        uint64_t outputInterval = 1e8;

        Timer.start("Timing LFs 0 to " + std::to_string(std::min(totalLen, outputInterval - 1)));
        do {
            current = mapLF(current, runlens, toRun, toOffset);
            ++lfsDone;
            if (lfsDone % outputInterval == outputInterval - 1) {
                Timer.stop();
                Timer.start("Timing LFs " + std::to_string(lfsDone + 1) + " to " + std::to_string(std::min(totalLen, lfsDone + outputInterval)));
            }
        } while (lfsDone <= totalLen && current != start);
        Timer.stop();

        if (lfsDone != totalLen) {
            std::cerr << "ERROR: LF Invalid after endmarker repair. " << lfsDone 
                << " LFs done in an attempt to traverse one cycle of the BWT. Length is actually " << totalLen
                << ". NOTE: if sumSeqLengths = length + 1, LF (very likely) was terminated early and didn't compute "
                << "the full sequence lengths. This is because verification is automatically terminated after length + 1"
                << " LFs because the LF function is then known to be incorrect." << std::endl;
            return 1;
        }

        std::cout << "LF is a permutation with one cycle, therefore it is very likely correct.\n";
        Timer.stop(); //Testing all LFs by permutation with one cycle
        */
    }
    Timer.stop(); //RLBWT Repair

    Timer.stop(); //builder
    return 0;
}
