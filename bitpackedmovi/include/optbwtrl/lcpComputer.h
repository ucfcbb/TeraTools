#include"util.h"
#include"moveStructure.h"
#include"fm-index.h"
#include<queue>

enum FFMethod { LINEAR, EXPONENTIAL };

struct packedTripleVector {
    typedef uint64_t size_type;
    uint8_t a,b,c,width;
    sdsl::bit_vector data;

    packedTripleVector() = default;

    packedTripleVector(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t elements): a(a), b(b), c(c), width(a+b+c), data(width*elements,0) {
        if (a > 64 || b > 64 || c > 64) {
            std::cerr << "ELEMENTS MUST BE <= 64 BITS!" << std::endl;
            exit(1);
        }
    }

    //packedTripleVector(packedTripleVector&& m): a(m.a), b(m.b), c(m.c), width(m.width), data(std::move(m.data)) {}

    template<unsigned el>
    inline std::pair<uint8_t, uint8_t> offW() const ;

    template<unsigned el>
    uint64_t get(const uint64_t ind) const {
        static_assert(el < 3, "The element accessed in a packed triple must be between 0 and 2 inclusive");
        auto [off, w] = offW<el>();
        return data.get_int(ind*width + off, w);
    }

    template<unsigned el>
    void set(const uint64_t ind, const uint64_t val) {
        static_assert(el < 3, "The element accessed in a packed triple must be between 0 and 2 inclusive");
        auto [off, w] = offW<el>();
        //std::cout << "off " << static_cast<unsigned>(off) << " w " << static_cast<unsigned>(w) << " ind " << ind << " ind*width " << ind*width
            //<< " ind*width + off " << ind*width + off << " val " << val << std::endl;
        data.set_int(ind*width + off, val, w);
        //std::cout << "leaving set" << std::endl;
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(a, out, child, "a");
        bytes += sdsl::serialize(b, out, child, "b");
        bytes += sdsl::serialize(c, out, child, "c");
        bytes += sdsl::serialize(width, out, child, "width");
        bytes += sdsl::serialize(data, out, child, "data");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(a, in);
        sdsl::load(b, in);
        sdsl::load(c, in);
        sdsl::load(width, in);
        sdsl::load(data, in);
    }

    uint64_t size() const {
        return data.size()/width;
    }
};
            
template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<0>()  const {
    return {0, a};
}

template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<1>()  const {
    return {a, b};
}

template<>
inline std::pair<uint8_t, uint8_t> packedTripleVector::offW<2>()  const {
    return {a + b, c};
}

//copied from bench.cpp, clean up structure of code later
struct MoveStructureTable {
    typedef uint64_t size_type;
    //D_index, D_offset, intlens
    packedTripleVector data;

    MoveStructureTable() = default;

    MoveStructureTable(MoveStructure& mv): data(mv.D_index.width(), mv.D_offset.width(), mv.intLens->width(), mv.D_index.size()) {
        //std::cout << "calc: " << ((mv.D_index.width()+mv.D_offset.width() + mv.intLens->width())*mv.D_index.size()+7)/8 << std::endl;
        //std::cout << "Input index size: " << sdsl::size_in_bytes(mv) << std::endl;
        //std::cout << "This index size: " << sdsl::size_in_bytes(*this) << std::endl;
        //std::cout << "In constructor" << std::endl;
        //std::cout << data.data.size() << std::endl;
        for (uint64_t i = 0; i < mv.D_index.size(); ++i) {
            //std::cout << "i " << i << std::endl;
            data.set<0>(i, mv.D_index[i]);
            data.set<1>(i, mv.D_offset[i]);
            data.set<2>(i, (*mv.intLens)[i]);
        }
        //std::cout << "Out constructor" << std::endl;
    }

    //MoveStructureTable(MoveStructureTable&& mv): data(std::move(mv.data)) {}

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

        bool operator!=(const IntervalPoint& rhs) const {
            return position != rhs.position || interval != rhs.interval || offset != rhs.offset;
        }
    };

    IntervalPoint map(const IntervalPoint& intPoint) const {
        IntervalPoint res;
        res.position = static_cast<uint64_t>(-1);
        res.interval = data.get<0>(intPoint.interval);
        res.offset = data.get<1>(intPoint.interval) + intPoint.offset;
        while (data.get<2>(res.interval) <= res.offset)
            res.offset -= data.get<2>(res.interval++);
        return res;
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(data, out, child, "data");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(data,in);
    }

    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{static_cast<uint64_t>(-1), 0, 0};

        do {
            curr = map(curr);
            ++totalOps;
        } while ((curr.interval || curr.offset) && totalOps < N + 1);

        return totalOps == N;
    }

    //returns whether the move structure is a permutation of length N
    //uses numIntervals log numIntervals auxiliary bits
    bool permutationLengthN(size_type N) const {
        size_type runs = data.size();
        sdsl::int_vector<> nextInt(runs, runs, sdsl::bits::hi(runs) + 1);
        size_type totalOps = 0;

        #pragma omp parallel for schedule(dynamic, 1024)
        for (uint64_t i = 0; i < runs; ++i) {
            IntervalPoint curr{static_cast<uint64_t>(-1), i, 0};
            uint64_t ops = 0;
            do {
                curr = map(curr);
                ++ops;
            } while (curr.offset);

            #pragma omp critical 
            {
                nextInt[i] = curr.interval;
                totalOps += ops;
            }
        }

        if (totalOps != N)
            return false;

        uint64_t traversed = 1, curr = 0;
        while (nextInt[curr] && nextInt[curr] != runs && traversed < runs) {
            curr = nextInt[curr];
            ++traversed;
        }

        return traversed == runs && nextInt[curr] == 0;
    }
};

//copied from bench.cpp, clean up structure of code later
struct MoveStructureStartTable {
    typedef uint64_t size_type;
    //D_index, D_offset, startPos
    packedTripleVector data;
    
    MoveStructureStartTable() = default;

    MoveStructureStartTable(const MoveStructure& mv) {
        uint64_t largestStartPos = 0;
        for (uint64_t i = 1; i < mv.D_index.size(); ++i)
            largestStartPos += (*mv.intLens)[i - 1];

        data = packedTripleVector(mv.D_index.width(), mv.D_offset.width(), sdsl::bits::hi(largestStartPos) + 1, mv.D_index.size() + 1);
        //std::cout << "In constructor" << std::endl;
        //std::cout << data.data.size() << std::endl;
        largestStartPos = 0;
        for (uint64_t i = 0; i < mv.D_index.size(); ++i) {
            //std::cout << "i " << i << std::endl;
            data.set<0>(i, mv.D_index[i]);
            data.set<1>(i, mv.D_offset[i]);
            data.set<2>(i, largestStartPos);
            largestStartPos += (*mv.intLens)[i];
        }
        data.set<2>(mv.D_index.size(), largestStartPos);
        //std::cout << "Out constructor" << std::endl;
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

        bool operator!=(const IntervalPoint& rhs) const {
            return position != rhs.position || interval != rhs.interval || offset != rhs.offset;
        }
    };
    
    template <FFMethod M = EXPONENTIAL>
    IntervalPoint map(const IntervalPoint& intPoint) const;

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(data, out, child, "data");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(data, in);
    }

    template<FFMethod M>
    //returns whether the move structure is a permutation of length N with exactly one cycle
    //uses numIntervals log numIntervals auxiliary bits
    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{0, 0, 0};

        do {
            curr = map<M>(curr);
            ++totalOps;
        } while (curr.position && totalOps < N + 1);

        return totalOps == N;

    }
    
    template<FFMethod M>
    //returns whether the move structure is a permutation of length N with exactly one cycle
    //uses numIntervals log numIntervals auxiliary bits
    bool permutationLengthN(size_type N) const {
        size_type runs = data.size() - 1;
        sdsl::int_vector<> nextInt(runs, runs, sdsl::bits::hi(runs) + 1);
        size_type totalOps = 0;

        #pragma omp parallel for schedule(dynamic, 1024)
        for (uint64_t i = 0; i < runs; ++i) {
            IntervalPoint curr{static_cast<uint64_t>(-1), i, 0};
            uint64_t ops = 0;
            do {
                curr = map<M>(curr);
                ++ops;
            } while (curr.offset);

            #pragma omp critical 
            {
                nextInt[i] = curr.interval;
                totalOps += ops;
            }
        }

        if (totalOps != N)
            return false;

        uint64_t traversed = 1, curr = 0;
        while (nextInt[curr] && nextInt[curr] != runs && traversed < runs) {
            curr = nextInt[curr];
            ++traversed;
        }

        return traversed == runs && nextInt[curr] == 0;
    }
};

template<> 
MoveStructureStartTable::IntervalPoint MoveStructureStartTable::map<LINEAR>(const IntervalPoint& intPoint) const {
    //std::cout << "enter map" << std::endl;
    IntervalPoint res;
    res.interval = data.get<0>(intPoint.interval);
    res.offset = data.get<1>(intPoint.interval) + intPoint.offset;
    res.position = data.get<2>(res.interval) + res.offset;
    if (data.get<2>(res.interval + 1) > res.position)
        return res;
    //std::cout << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
    //res.interval should never be > data.get<2>.size()
    while (data.get<2>(res.interval + 1) <= res.position)
        ++res.interval;
    res.offset = res.position - data.get<2>(res.interval);
    //std::cout << "exit map" << std::endl;
    return res;
}

template<> 
MoveStructureStartTable::IntervalPoint MoveStructureStartTable::map<EXPONENTIAL>(const IntervalPoint& intPoint) const {
    //std::cout << "enter map" << std::endl;
    IntervalPoint res;
    res.interval = data.get<0>(intPoint.interval);
    res.offset = data.get<1>(intPoint.interval) + intPoint.offset;
    res.position = data.get<2>(res.interval) + res.offset;
    if (data.get<2>(res.interval + 1) > res.position)
        return res;
    uint64_t dist = 2;
    //std::cout << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
    //res.interval should never be > data.size()
    while (dist < data.size() - 1 - res.interval && data.get<2>(res.interval + dist) <= res.position)
        dist *= 2;
    res.interval += dist/2;
    
    dist = dist/2 - 1;
    using myint = uint64_t;

    myint lo = res.interval;
    myint hi = std::min(lo + dist, data.size() - 2);
    myint ans = -1;
    while(lo <= hi) {
        myint mid = (lo + hi) >> 1;
        myint here = data.get<2>(mid);
        if (here <= res.position){
            lo = (ans = mid) + 1;
        } else {
            hi = mid - 1;
        }
    }
    
    res.interval = ans;
    res.offset = res.position - data.get<2>(res.interval);
    //std::cout << "exit map" << std::endl;
    return res;
}

class LCPComputer {
    uint64_t totalLen;

    sdsl::int_vector<> F;
    
    //Flens is now part of Psi
    //sdsl::int_vector<> Flens;
    MoveStructureTable Psi;

    sdsl::int_vector<> intAtTop;

    //sdsl::int_vector<> PhiIntLen;
    MoveStructureStartTable PhiNEWONEHE;

    sdsl::int_vector<> PLCPsamples;

    void ConstructPsi(const rb3_fmi_t* rb3, sdsl::int_vector<> & F, MoveStructureTable &Psi, uint64_t & numSequences) {
        Timer.start("Constructing Psi from FMD");

        uint64_t runs = 0, alphbits, lenbits;
        totalLen = 0;
        uRange lenRange;
        std::vector<uint64_t> alphRuns;

        uRange RB3_lenRange;
        Timer.start("Reading fmd for parameters");
        {
            //original ropebwt3 values
            uint64_t RB3_runs = 0, RB3_lenbits;

            uRange alphRange;



            rlditr_t itr1;
            rld_itr_init(rb3->e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;
            //int prevc;


            if ((l = rld_dec(rb3->e, &itr1, &c, 0)) > 0) {
                alphRange = {static_cast<uint64_t>(c),static_cast<uint64_t>(c)};

                RB3_lenRange = {static_cast<uint64_t>(l),static_cast<uint64_t>(l)};
                lenRange = (c == 0)? uRange{static_cast<uint64_t>(1), static_cast<uint64_t>(1)} : RB3_lenRange;

                ++RB3_runs;
                runs += (c == 0)? static_cast<uint64_t>(l) : 1;
                if (static_cast<uint64_t>(c) >= alphRuns.size())
                    alphRuns.resize(c+1);
                alphRuns[c] += (c == 0)? static_cast<uint64_t>(l) : 1;

                totalLen += static_cast<uint64_t>(l);
            }
            else {
                std::cerr << "Failed to read first run's character and length" << std::endl;
                exit(1);
            }
            //prevc = c;

            while ((l = rld_dec(rb3->e, &itr1, &c, 0)) > 0) {
                alphRange.min = std::min(alphRange.min, static_cast<uint64_t>(c));
                alphRange.max = std::max(alphRange.max, static_cast<uint64_t>(c));

                RB3_lenRange.min = std::min(RB3_lenRange.min, static_cast<uint64_t>(l));
                RB3_lenRange.max = std::max(RB3_lenRange.max, static_cast<uint64_t>(l));
                lenRange.min = std::min(lenRange.min, static_cast<uint64_t>((c == 0)? 1 : l));
                lenRange.max = std::max(lenRange.max, static_cast<uint64_t>((c == 0)? 1 : l));

                ++RB3_runs;
                runs += (c == 0)? static_cast<uint64_t>(l) : 1;
                if (static_cast<uint64_t>(c) >= alphRuns.size())
                    alphRuns.resize(c+1);
                alphRuns[c] += (c == 0)? static_cast<uint64_t>(l) : 1;

                totalLen += static_cast<uint64_t>(l);

                //if (c == prevc) {
                //    std::cerr << "ERROR: Two runs of the same character follow each other from ropebwt3!" << std::endl;
                //    exit(1);
                //}
                //prevc = c;
            }

            numSequences = alphRuns[0];

            if (alphRange.max == static_cast<uint64_t>(-1)) {
                std::cerr << "Maximum alphabet symbol is 2^64 - 1. "
                    << "This program assumes this is not the case (it can only handle alphabet <= (2^64) - 2." << std::endl;
                exit(1);
            }

            std::cout << "INFO: The parameters for our constructed BWT (i.e. #runs, max length, etc.) may be "
                << "different from those of the input (ropebwt3).\nINFO: This is because in our constructed BWT, "
                << "each endmarker is contained in its own run.\n";

            alphbits = sdsl::bits::hi(alphRange.max) + 1;
            RB3_lenbits = sdsl::bits::hi(RB3_lenRange.max) + 1;
            lenbits = sdsl::bits::hi(lenRange.max) + 1;
            if (alphbits != static_cast<uint64_t>(rb3->e->abits)) 
                std::cout << "WARNING: computed bits per symbol not equal to bits used in fmd. Computed: " 
                    << alphbits << ", ropebwt3: " << static_cast<uint64_t>(rb3->e->abits) << std::endl;

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

        MoveStructure tempPsi;
        Timer.start("Reading fmd to create Flens");
        {
            tempPsi.intLens = new sdsl::int_vector<>(runs, 0, lenbits);

            std::vector<uint64_t> alphFRunStarts(alphRuns.size());
            for (uint64_t i = 1; i < alphRuns.size(); ++i) {
                alphFRunStarts[i] = alphFRunStarts[i-1] + alphRuns[i-1];
            }

            rlditr_t itr1;
            rld_itr_init(rb3->e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;

            while ((l = rld_dec(rb3->e, &itr1, &c, 0)) > 0) {
                if (c == 0) {
                    for (uint64_t i = 0; i<static_cast<uint64_t>(l); ++i)
                        (*tempPsi.intLens)[alphFRunStarts[c]++] = 1;
                }
                else {
                    (*tempPsi.intLens)[alphFRunStarts[c]++] = l;
                }
            }

            //error checking:
            for (uint64_t i = 0, prevStart = 0; i < alphRuns.size(); ++i) {
                if (alphFRunStarts[i] - prevStart != alphRuns[i]) {
                    std::cerr << "alphFRunStarts[" << i << "] did not end up at correct run!" << std::endl;
                    exit(1);
                }
                prevStart = alphFRunStarts[i];
            }
        }
        Timer.stop(); //Reading fmd to create Flens

        //NOTE: Psi of endmarker runs will be incorrect, will fix in a later step
        Timer.start("Reading fmd to create D_index and D_offset");
        {
            tempPsi.D_index = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(runs - 1) + 1);
            tempPsi.D_offset = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(lenRange.max - 1) + 1);

            std::vector<uint64_t> alphFRunStarts(alphRuns.size());
            for (uint64_t i = 1; i < alphRuns.size(); ++i) {
                alphFRunStarts[i] = alphFRunStarts[i-1] + alphRuns[i-1];
            }

            uint64_t currentRun = 0, currentOffset = 0;

            rlditr_t itr1;
            rld_itr_init(rb3->e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;

            while ((l = rld_dec(rb3->e, &itr1, &c, 0)) > 0) {
                if (c == 0) {
                    for (uint64_t i = 0; i<static_cast<uint64_t>(l); ++i) {
                        tempPsi.D_index[alphFRunStarts[c]] = currentRun;
                        tempPsi.D_offset[alphFRunStarts[c]] = currentOffset;

                        ++alphFRunStarts[c];

                        ++currentOffset;
                        currentOffset %= (*tempPsi.intLens)[currentRun];
                        currentRun += (currentOffset == 0);
                    }
                }
                else {
                    tempPsi.D_index[alphFRunStarts[c]] = currentRun;
                    tempPsi.D_offset[alphFRunStarts[c]] = currentOffset;

                    ++alphFRunStarts[c];

                    currentOffset += l;
                    while (currentOffset && currentOffset >= (*tempPsi.intLens)[currentRun])
                        currentOffset -= (*tempPsi.intLens)[currentRun++];
                }
            }

            //error checking:
            for (uint64_t i = 0, prevStart = 0; i < alphRuns.size(); ++i) {
                if (alphFRunStarts[i] - prevStart != alphRuns[i]) {
                    std::cerr << "alphFRunStarts[" << i << "] did not end up at correct run!" << std::endl;
                    exit(1);
                }
                prevStart = alphFRunStarts[i];
            }
        }
        Timer.stop(); //Reading fmd to create D_index and D_offset

        //obviously, we could just replace F with a bit vector of length r with sigma set bits
        //and use rank on the bitvector
        //for the same time complexity but r bits instead of r log sigma
        //I'm not sure how much slower that is. I haven't tried it yet.
        Timer.start("Constructing F");
        {
            F = sdsl::int_vector<>(runs, 7, alphbits);
            uint64_t curr = 0;
            for (uint64_t alph = 0; alph < alphRuns.size(); ++alph)
                for (uint64_t i = 0; i < alphRuns[alph]; ++i)
                    F[curr++] = alph;
        }
        Timer.stop(); //Constructing F

        Psi = std::move(MoveStructureTable(tempPsi));
        delete tempPsi.intLens;
        Timer.stop(); //Constructing Psi from FMD
    }

    //O(n) time, not needed for LCP computation
    //LCP computation is O(n) anyways and we can repair while we do it
    //although we don't use psi after so may as well not repair it
    //this function fixes the psi mappings for the endmarker runs (of F)
    //It also computes numTopRuns and seqLens
    //numTopRuns and seqLens are vectors of length equal to the number of sequences
    //the i+1-th value of numTopRuns is the number of times suffixes of sequence i
    //are at the (SA position corresponding to the) top of a run in the BWT. 
    //(includes the termination symbol at the end of sequence i, $_i)
    //the i+1-th value of seqLens is the length of seq i,
    //(includes the termination symbol) 
    //Actually, numTopRuns and seqLens are the prefix sums of the above definitions
    //
    //
    //for every suffix x at the top of a run in the BWT, 
    //there is an input interval of psi, j, where 
    //suffix x-1 is at the top of the input interval
    //intAtTop stores, for every input interval of psi
    //with suffix x-1 at the top of the input interval,
    //how many suffixes < x are at the top of a run in the BWT
    void ComputeAuxAndRepairPsi(uint64_t& maxPhiIntLen, std::vector<uint64_t> & numTopRuns, std::vector<uint64_t> & seqLens, sdsl::int_vector<> & intAtTop, 
            const sdsl::int_vector<>& F, MoveStructureTable& Psi, sdsl::int_vector<>& PhiIntLen, const uint64_t numSequences) {
        Timer.start("Computing numTopRuns, seqLens, and repairing Psi of endmarkers in F");

        numTopRuns.resize(numSequences + 1);
        seqLens = std::vector<uint64_t>(numSequences + 1);

        intAtTop = sdsl::int_vector<>(F.size(), -1, sdsl::bits::hi(F.size() - 1) + 1);

        //computing seqLens, numTopRuns, correctSeqPsis, and maxPhiIntLen
        Timer.start("Parallel seq traversal");
        {
            std::vector<uint64_t> maxPhiIntLenPerSeq(numSequences);
            std::vector<MoveStructureTable::IntervalPoint> correctSeqPsis(numSequences);
            #pragma omp parallel for schedule(dynamic, 1)
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                MoveStructureTable::IntervalPoint start = {static_cast<uint64_t>(-1), seq, 0}, curr;
                start = Psi.map(start);
                curr = start;

                //suffix 0 of seq seq is at the top of a run in the bwt since its bwt value is a termination string
                uint64_t numTop = 1;
                uint64_t len = 1;
                uint64_t maxPhiIntLenThisSeq = 0;
                uint64_t currPhiIntLen = 1;

                while (curr.interval >= numSequences) {
                    //suffix x is at the top of a run in the bwt iff suffix x - 1 is at the start
                    //of an input interval of psi
                    numTop += (curr.offset == 0);
                    if (curr.offset == 0) {
                        maxPhiIntLenThisSeq = std::max(maxPhiIntLenThisSeq, currPhiIntLen);
                        currPhiIntLen = 0;
                    }
                    ++len;
                    ++currPhiIntLen;
                    curr = Psi.map(curr);
                }

                //for last interval in sequence, if condition is unnecessary, guaranteed to be true.
                if (curr.offset == 0) {
                    maxPhiIntLenThisSeq = std::max(maxPhiIntLenThisSeq, currPhiIntLen);
                    currPhiIntLen = 0;
                }

                if (curr.offset) {
                    std::cerr << "ERROR: Run of endmarkers in F of length more than 1!" << std::endl;
                    exit(1);
                }

                uint64_t seqStartingAtStart = curr.interval;

                //should not need omp critical
                correctSeqPsis[(seqStartingAtStart)? seqStartingAtStart - 1 : numSequences - 1] = start;
                seqLens[seqStartingAtStart + 1] = len;
                numTopRuns[seqStartingAtStart + 1] = numTop;
                maxPhiIntLenPerSeq[seqStartingAtStart] = maxPhiIntLenThisSeq;
            }

            maxPhiIntLen = 0;
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                Psi.data.set<0>(seq, correctSeqPsis[seq].interval);
                Psi.data.set<1>(seq, correctSeqPsis[seq].offset);
                maxPhiIntLen = std::max(maxPhiIntLen, maxPhiIntLenPerSeq[seq]);
            }
        }
        Timer.stop(); //Parallel seq traversal

        //from here on, numTopRuns and seqLens are the exclusive prefix sums of their previous definition
        Timer.start("Prefix summing auxiliary data");
        for (uint64_t i = 1; i < seqLens.size(); ++i) {
            seqLens[i] += seqLens[i-1];
            numTopRuns[i] += numTopRuns[i-1];
        }

        if (numTopRuns.back() != F.size()) {
            std::cerr << "ERROR: Number of runs in numTopRuns doesn't sum to total number of runs!" << std::endl;
            exit(1);
        }
        if (seqLens.back() != totalLen) {
            std::cerr << "ERROR: Lengths in seqLens doesn't sum to total length!" << std::endl;
            exit(1);
        }
        Timer.stop(); //Prefix summing auxiliary data

        std::cout << "INFO: maximum interval length for phi data structure: " << maxPhiIntLen << std::endl;

        PhiIntLen = sdsl::int_vector<>(F.size(), 0, sdsl::bits::hi(maxPhiIntLen) + 1);
        //computing intAtTop and PhiIntLen
        Timer.start("Second Parallel seq traversal");
        #pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t seq = 0; seq < numSequences; ++seq) {
            uint64_t prevSeq = (seq)? seq - 1 : numSequences - 1;
            MoveStructureTable::IntervalPoint curr = {static_cast<uint64_t>(-1), prevSeq, 0};
            curr = Psi.map(curr);

            uint64_t currentInt = numTopRuns[seq];
            uint64_t currIntLen = 1;
            do {
                if (curr.offset == 0) {
                    #pragma omp critical
                    {
                        PhiIntLen[currentInt] = currIntLen;
                        currIntLen = 0;
                        intAtTop[curr.interval] = currentInt++;
                    }
                }
                curr = Psi.map(curr);
                ++currIntLen;
                //std::cout << "curr.interval " << curr.interval << " curr.offset " << curr.offset << std::endl;
            } while (curr.interval >= numSequences);

            if (curr.offset == 0) {
                #pragma omp critical
                {
                    PhiIntLen[currentInt] = currIntLen;
                    currIntLen = 0;
                    intAtTop[curr.interval] = currentInt++;
                }
            }

            if (currentInt != numTopRuns[seq + 1]) {
                std::cerr << "ERROR: Didn't reach beginning of next sequence in numTopRuns!" << std::endl;
                //std::cout << currentInt << " " << numTopRuns[seq + 1] << std::endl;
                //for (uint64_t i = 0; i < F.size(); ++i) {
                    //std::cout << intAtTop[i] << std::endl;
                //}
                exit(1);
            }
        }
        Timer.stop(); //Second Parallel seq traversal
        Timer.stop(); //Computing numTopRuns, seqLens, and repairing Psi of endmarkers in F
    }

    void ConstructPhiAndSamples(const sdsl::int_vector<>& F, const MoveStructureTable& Psi, sdsl::int_vector<>& PhiIntLen, const uint64_t FlensBits,
            const std::vector<uint64_t>& numTopRuns, const std::vector<uint64_t>& seqLens, const sdsl::int_vector<>& intAtTop, const uint64_t numSequences, const uint64_t maxPhiIntLen,
            uint64_t & sampleInterval, sdsl::int_vector<> &Psi_Index_Samples, sdsl::int_vector<> &Psi_Offset_Samples) {
        Timer.start("Construct Phi and Samples");
        //MoveStructure tempPhi;
        //tempPhi.intLens = &PhiIntLen;
        //computing, for each seq i, 
        // 1. whenever suffix j of i is at the bottom of a run, the input interval and offset of suffix j of i in Phi (D_index and D_offset of the top of the run below it)
        // 2. Psi_Index_Samples and Psi_Offset_Samples: ISA samples (in psi move data structure), for every suffix that is a multiple of n/r (of the original text)
        uint64_t numRuns = F.size();
        sampleInterval = totalLen/numRuns;
        uint64_t numSamples = (totalLen%sampleInterval != 0) + (totalLen/sampleInterval);
        //std::cout << "sampleInterval " << sampleInterval << std::endl;
        //std::cout << "numSamples " << numSamples << std::endl;
        //tempPhi.D_index = sdsl::int_vector<>(numRuns, 0, sdsl::bits::hi(numRuns - 1) + 1);
        //tempPhi.D_offset = sdsl::int_vector<>(numRuns, 0, sdsl::bits::hi(maxPhiIntLen - 1) + 1);
        {
            uint64_t sum = 0;
            for (uint64_t i = 0; i < PhiIntLen.size(); ++i)
                sum += PhiIntLen[i];
            PhiNEWONEHE.data = packedTripleVector(sdsl::bits::hi(numRuns - 1) + 1, sdsl::bits::hi(maxPhiIntLen - 1) + 1, sdsl::bits::hi(sum) + 1, PhiIntLen.size() + 1);
            sum = 0;
            for (uint64_t i = 0; i < PhiIntLen.size(); ++i){
                PhiNEWONEHE.data.set<2>(i, sum);
                sum += PhiIntLen[i];
            }
            PhiNEWONEHE.data.set<2>(PhiIntLen.size(), sum);
            PhiIntLen = sdsl::int_vector<>();
        }
        Psi_Index_Samples = sdsl::int_vector<>(numSamples, 0, sdsl::bits::hi(numRuns - 1) + 1);
        Psi_Offset_Samples = sdsl::int_vector<>(numSamples, 0, FlensBits);
        std::cout << "FlensBits: " << FlensBits << std::endl;
        std::cout << "Psi_Offset_Samples.width(): " << static_cast<uint64_t>(Psi_Offset_Samples.width()) << std::endl;
        #pragma omp parallel for schedule (dynamic, 1)
        for (uint64_t seq = 0; seq < numSequences; ++seq) {
            //curr is the interval point in the psi move data structure of suffix suff
            //if curr is at the top of a psi interval, then suff+1 is at the top of an rlbwt interval
            //if curr is at the bottom of a psi interval, then suff+1 is at the bottom of an rlbwt interval
            MoveStructureTable::IntervalPoint curr = {static_cast<uint64_t>(-1), (seq)? seq - 1 : numSequences - 1, 0};
            curr = Psi.map(curr);
            uint64_t suff = seqLens[seq];
            //phiPoint is the interval point in the move structure of suff
            MoveStructureStartTable::IntervalPoint phiPoint = {PhiNEWONEHE.data.get<2>(numTopRuns[seq]), numTopRuns[seq], 0};
            MoveStructureStartTable::IntervalPoint phiPointAtPhiOutputIntervalStart = phiPoint;

            while (suff < seqLens[seq+1]) {
                //2.
                //store psi samples if needed
                if (suff % sampleInterval == 0) {
                    #pragma omp critical
                    {
                        Psi_Index_Samples[suff/sampleInterval] = curr.interval;
                        Psi_Offset_Samples[suff/sampleInterval] = curr.offset;
                    }
                }

                //update phiPoint to suff + 1
                ++phiPoint.offset;
                ++phiPoint.position;
                //phiPoint.offset == (*Phi.intLens)[phiPoint.interval] iff curr.offset == 0)
                assert((phiPoint.position == PhiNEWONEHE.data.get<2>(phiPoint.interval + 1)) ==
                        (curr.offset == 0));
                //suff + 1 starts an input interval iff suff ends an input interval
                if (curr.offset == 0) {
                    if (phiPoint.position != PhiNEWONEHE.data.get<2>(phiPoint.interval + 1)) {
                        std::cerr << "Offset != length at the end of phi interval!" << std::endl;
                        exit(1);
                    }
                    phiPoint.offset = 0;
                    ++phiPoint.interval;
                }

                //if suff is the end of an output interval
                if (curr.offset == Psi.data.get<2>(curr.interval) - 1) {
                    uint64_t runBelow = (curr.interval+1)%numRuns;
                    uint64_t phiInterval = intAtTop[runBelow];
                    #pragma omp critical 
                    {
                        PhiNEWONEHE.data.set<0>(phiInterval, phiPointAtPhiOutputIntervalStart.interval);
                        PhiNEWONEHE.data.set<1>(phiInterval, phiPointAtPhiOutputIntervalStart.offset);
                    }
                    phiPointAtPhiOutputIntervalStart = phiPoint;
                }
                
                //now, curr is the position of suff + 1
                curr = Psi.map(curr);
                //update suff
                ++suff;
                /*
                std::cout << "suff " << suff << std::endl;
                std::cout << "curr: " << curr.interval << ' ' << curr.offset << std::endl;
                std::cout << "intLen: " << (*Psi.intLens)[curr.interval] << std::endl;
                //1. 
                //when curr at the bottom of a psi input interval, suff is the start of a phi output interval
                //then suff-1 is the end of a phi output interval, then invPhi(suff-1) is the end of a phi input interval
                //so, update the start of the previous input interval to phiPointAtPhiOutputIntervalStart
                if (curr.offset == (*Psi.intLens)[curr.interval] - 1)
                    std::cout << "starts output interval" << std::endl;
                if (curr.offset == (*Psi.intLens)[curr.interval] - 1) {
                    std::cout << "in if" << std::endl;
                    uint64_t runBelow = (curr.interval+1)%numRuns;
                    uint64_t phiInterval = intAtTop[runBelow];
                    if (phiInterval == 25) {
                        std::cout << "25: " << seq << " " << suff << std::endl;
                    }
                    #pragma omp critical 
                    {
                        Phi.D_index[phiInterval] = phiPointAtPhiOutputIntervalStart.interval;
                        Phi.D_offset[phiInterval] = phiPointAtPhiOutputIntervalStart.offset;
                    }
                    phiPointAtPhiOutputIntervalStart = phiPoint;
                }

                bool suffAtTop = (curr.offset == 0);
                ++suff;
                ++phiPoint.offset;

                if (suff == 30) {
                    std::cout << suffAtTop << std::endl;
                }
                //2.
                //start new interval in phi if suff is at the top of a run
                if (suffAtTop) {
                    #pragma omp critical
                    {
                        ++phiPoint.interval;
                        phiPoint.offset = 0;
                    }
                }
                */
            }
            /*
            std::cout << "done while" << std::endl;
            //condition in the following if statement is redundant, should never be false.
            //1. 
            //when at the bottom of a psi input interval, suff is the start of a phi output interval
            //so, update the start of the previous input interval to phiPointAtPhiOutputIntervalStart
            if (curr.offset == (*Psi.intLens)[curr.interval] - 1) {
                uint64_t runBelow = (curr.interval+1)%numRuns;
                uint64_t phiInterval = intAtTop[runBelow];
                #pragma omp critical 
                {
                    Phi.D_index[phiInterval] = phiPointAtPhiOutputIntervalStart.interval;
                    Phi.D_offset[phiInterval] = phiPointAtPhiOutputIntervalStart.offset;
                }
                //phiPointAtPhiOutputIntervalStart = phiPoint;
            }
            else {
                std::cerr << "ERROR: didn't set phi values of last interval of seq!" << std::endl;
                exit(1);
            }
            */

            if (suff != seqLens[seq+1]) {
                std::cerr << "ERROR: suff didn't end up at the beginning of the next sequence!" << std::endl;
                exit(1);
            }
            if (phiPoint.interval != numTopRuns[seq+1]) {
                std::cerr << "ERROR: phiPoint didn't end up at the beginning of the next sequence!" << std::endl;
                exit(1);
            }
        }
        Timer.stop(); //Construct Phi and Samples
    }

    void ComputePLCPSamples(const sdsl::int_vector<>& intAtEnd, const uint64_t numSequences, const sdsl::int_vector<>& F, 
            const MoveStructureTable& Psi, const std::vector<uint64_t>& numTopRuns, const std::vector<uint64_t>& seqLens,
            const uint64_t sampleInterval, const sdsl::int_vector<>& Psi_Index_Samples, const sdsl::int_vector<>& Psi_Offset_Samples) {
        Timer.start("LCP Computation");
        PLCPsamples = sdsl::int_vector<>(F.size(), 0, 1);

        /*
        sdsl::int_vector<> suffStartingPhiInt(F.size(), 0, sdsl::bits::hi(seqLens.back() - 1) + 1);
        for (uint64_t i = 1; i < F.size(); ++i)
            suffStartingPhiInt[i] = suffStartingPhiInt[i-1] + PhiNEWONEHE.data.get<2>(i-1);
         */

        #pragma omp parallel for schedule (dynamic, 1)
        for (uint64_t seq = 0; seq < numSequences; ++seq) {
            uint64_t suffMatchEnd = seqLens[seq], currIntStart = seqLens[seq];
            MoveStructureTable::IntervalPoint suffMatchEndIntPoint = Psi.map({static_cast<uint64_t>(-1), ((seq)? seq - 1 : numSequences - 1), 0});
            for (uint64_t phiInt = numTopRuns[seq]; phiInt < numTopRuns[seq+1]; ++phiInt) { 
                //get suffix above phiInt
                MoveStructureStartTable::IntervalPoint suffMatchingToPhiIntPoint = PhiNEWONEHE.map({static_cast<uint64_t>(-1), phiInt, 0});
                uint64_t suffMatchingTo = PhiNEWONEHE.data.get<2>(suffMatchingToPhiIntPoint.interval) + suffMatchingToPhiIntPoint.offset;
                uint64_t psiIntAbove = intAtEnd[(phiInt)? phiInt - 1 : F.size() - 1];
                MoveStructureTable::IntervalPoint coordAbove = Psi.map({static_cast<uint64_t>(-1), psiIntAbove, 0});
                if (coordAbove.offset == 0) {
                    coordAbove.interval = (coordAbove.interval)? coordAbove.interval - 1 : F.size() - 1;
                    coordAbove.offset = Psi.data.get<2>(coordAbove.interval);
                }
                --coordAbove.offset;
                if (suffMatchEnd - currIntStart) {
                    uint64_t destinationSuff = suffMatchingTo + (suffMatchEnd - currIntStart);
                    uint64_t sampleNum = destinationSuff/sampleInterval;
                    uint64_t closestSampledSuff = sampleNum * sampleInterval;
                    //by default, use above run
                    uint64_t dist = suffMatchEnd - currIntStart;
                    //otherwise, if sample is closer
                    if (dist > destinationSuff - closestSampledSuff) {
                        coordAbove = {static_cast<uint64_t>(-1), Psi_Index_Samples[sampleNum], Psi_Offset_Samples[sampleNum]};
                        dist = destinationSuff - closestSampledSuff;
                    }

                    for (uint64_t i = 0; i < dist; ++i) 
                        coordAbove = Psi.map(coordAbove);
                }
                //compute lcp of interval phiInt
                //the LCP value at the start of this interval is at least the interval length - 1
                //(it is only the interval length - 1 when the suffix at the start of the interval starts with c
                //and is the lexicographically smallest suffix that starts with c.
                //happens by default for all endmarkers since multidollar bwt

                //In our application, LCP can be less than the interval length when Ns are present in the interval
                //since we don't allow matches to include Ns. If this wasn't the case, the currLCP could be initialized to at least (*Phi.intLens)[phiInt] - 1

                while (F[suffMatchEndIntPoint.interval] != 0 && F[suffMatchEndIntPoint.interval] != 5 &&
                        F[suffMatchEndIntPoint.interval] == F[coordAbove.interval]) {
                    suffMatchEndIntPoint = Psi.map(suffMatchEndIntPoint);
                    coordAbove = Psi.map(coordAbove);
                    ++suffMatchEnd;
                }

                #pragma omp critical
                {
                    sdsl::util::expand_width(PLCPsamples, sdsl::bits::hi(suffMatchEnd - currIntStart) + 1);
                    PLCPsamples[phiInt] = suffMatchEnd - currIntStart;
                }


                currIntStart += PhiNEWONEHE.data.get<2>(phiInt);
                if (suffMatchEnd < currIntStart) {
                    suffMatchEnd = currIntStart;
                    suffMatchEndIntPoint = Psi.map({static_cast<uint64_t>(-1), intAtEnd[phiInt], 0});
                }
            }
        }

        sdsl::util::bit_compress(PLCPsamples);
        std::cout << "PLCP width: " << static_cast<uint64_t>(PLCPsamples.width()) << std::endl;
        Timer.stop(); //LCP Computation
    }

    void ComputeMinLCPRun(const sdsl::int_vector<>& intAtTop, const sdsl::int_vector<>& F, const MoveStructureTable& Psi, std::ofstream& out) {
        //O(r log sigma) time, can be skipped if RLBWT is maintained or re-read
        Timer.start("Computing Min LCP per Run");
        struct MappedPositionRunPair {
            uint64_t psiInputInt;
            MoveStructureTable::IntervalPoint runStart;

            MappedPositionRunPair(uint64_t interval, const MoveStructureTable& psi): psiInputInt(interval) {
                runStart = psi.map({static_cast<uint64_t>(-1), psiInputInt, 0});
            }
        };

        struct ComparePositionPair {
            bool operator()(const MappedPositionRunPair& lhs, const MappedPositionRunPair& rhs) const {
                return lhs.runStart.interval > rhs.runStart.interval || (lhs.runStart.interval == rhs.runStart.interval && lhs.runStart.offset > rhs.runStart.offset);
            }
        };


        std::priority_queue<MappedPositionRunPair,std::vector<MappedPositionRunPair>,ComparePositionPair> nextAlph;
        std::vector<uint64_t> endmarkerOrder;
        {
            uint64_t seq = 0;
            while (F[seq] == 0)
                nextAlph.emplace(seq++,Psi);
            endmarkerOrder.reserve(nextAlph.size());
            while(nextAlph.size()) {
                endmarkerOrder.push_back(nextAlph.top().psiInputInt);
                nextAlph.pop();
            }
        }

        //get first of each character
        nextAlph.emplace(endmarkerOrder[0], Psi);
        for (uint64_t i = endmarkerOrder.size(); i < F.size(); ++i)
            if (F[i] != F[i-1])
                nextAlph.emplace(i, Psi);

        MappedPositionRunPair firstRun = nextAlph.top();

        uint64_t runs = 0, endmarkerPosition = 0;
        while (nextAlph.size()) {
            MappedPositionRunPair t = nextAlph.top();
            nextAlph.pop();
            if (F[t.psiInputInt] != 0 && t.psiInputInt != F.size() - 1 && F[t.psiInputInt] == F[t.psiInputInt+1])
                nextAlph.emplace(t.psiInputInt+1, Psi);
            if (F[t.psiInputInt] == 0 && endmarkerPosition != endmarkerOrder.size() - 1)
                nextAlph.emplace(endmarkerOrder[++endmarkerPosition], Psi);
            ++runs;

            uint64_t runLen = Psi.data.get<2>(t.psiInputInt);
            //out << F[t.psiInputInt] << ' ' << runLen << ' ';
            t = (nextAlph.size())? nextAlph.top() : firstRun;
            uint64_t phiInt = (intAtTop[t.psiInputInt]+1) % F.size();
            MoveStructureStartTable::IntervalPoint pPoint = {static_cast<uint64_t>(-1), phiInt, 0};
            uint64_t minL = static_cast<uint64_t>(-1), minLloc = static_cast<uint64_t>(-1);
            for (uint64_t i = 0; i < runLen; ++i) {
                pPoint = PhiNEWONEHE.map(pPoint);
                uint64_t l = PLCPsamples[pPoint.interval] - pPoint.offset;
                //use <= for first index
                if (l < minL) {
                    minL = l;
                    minLloc = runLen - 1 - i;
                }
            }
            out << "( " << minLloc << ", " << minL << ")\n";
        }

        if (runs != F.size()) {
            std::cerr << "ERROR: runs found by recovering RLBWT not equal to runs in F!\n";
            exit(1);
        }
        Timer.stop(); //Computing Min LCP per Run
    }

    public:
    typedef uint64_t size_type;

    //input: a run length encoding of a multidollar BWT where all dollars are represented by 0
    //All characters between (and including) 0 and max_char are assumed to have more than 0 occurrences
    //in the text. max_char is the maximum character in the text
    LCPComputer(rb3_fmi_t* rb3, const std::string& safeTempName, std::ofstream& lcpOut) {
        std::ofstream tempOutFile(safeTempName);
        if (!tempOutFile.is_open()) {
            std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' failed to open for writing!" << std::endl;
            exit(1);
        }

        sdsl::memory_monitor::start();

        //Psi.intLens = &Flens;
        uint64_t numSequences;
        {
            auto event = sdsl::memory_monitor::event("Construct Psi");
            ConstructPsi(rb3, F, Psi, numSequences);

            rb3_fmi_free(rb3);
        }

        /*
        for (uint64_t i = 0; i < Psi.D_index.size(); ++i) {
            std::cout 
                << F[i] << '\t'
                << (*Psi.intLens)[i] << '\t'
                << Psi.D_index[i] << '\t'
                << Psi.D_offset[i] << '\n';
        }
        {
            MoveStructure::IntervalPoint end{static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            curr = end;
            char con[] = "$ACGTN";
            uint64_t count = 0;
            do {
                curr = Psi.map(curr);
                std::cout << con[F[curr.interval]];
                ++count;
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }
        */

        //numTopRuns and seqLens are vectors of length equal to the number of sequences+1
        //the i+1-th value of numTopRuns is the number of times suffixes of sequence i
        //are at the (SA position corresponding to the) top of a run in the BWT. 
        //(includes the termination symbol at the end of sequence i, $_i)
        //the i+1-th value of seqLens is the length of seq i,
        //(includes the termination symbol) 
        std::vector<uint64_t> numTopRuns, seqLens;
        //for every suffix x at the top of a run in the BWT, 
        //there is an input interval of psi, j, where 
        //suffix x-1 is at the top of the input interval
        //intAtTop stores, for every input interval of psi
        //with suffix x-1 at the top of the input interval,
        //how many suffixes < x - 1 are at the top of a run in the BWT
        uint64_t maxPhiIntLen;
        sdsl::int_vector<> PhiIntLen;
        {
            auto event = sdsl::memory_monitor::event("Compute Auxiliary Data and Repair Endmarker Psis");
            ComputeAuxAndRepairPsi(maxPhiIntLen, numTopRuns, seqLens, intAtTop, F, Psi, PhiIntLen, numSequences);
        }

        {
            auto event = sdsl::memory_monitor::event("Computing and storing intAtEnd");
            Timer.start("Computing and storing intAtEnd");
            sdsl::int_vector<> intAtEnd(F.size(), 0, sdsl::bits::hi(F.size() - 1) + 1);
            //intAtEnd[j] is the input interval of Psi where the suffix at the end of input interval j of Phi occurs
            //(at the top of, necessarily)
            for (uint64_t i = 0; i < F.size(); ++i)
                intAtEnd[intAtTop[i]] = i;
            sdsl::serialize(intAtEnd, tempOutFile);
            Timer.stop(); //Computing and storing intAtEnd
        }
        /*
        {
            std::cout << "Checking if intAtTop is a permutation of [0,r-1]" << std::endl;
            std::vector<bool> a(F.size());
            for (uint64_t i = 0; i < a.size(); ++i) {
                if (intAtTop[i] >= F.size() || a[intAtTop[i]]) {
                    std::cerr << "ERROR: intAtTop is not a permutation of [0,r-1]" << std::endl;
                    exit(1);
                }
                a[intAtTop[i]] = true;
            }
        }

        {
            std::cout << "Psi SA order\ni\trunInd\tF\tintLen\tD_index\tD_offset\n";
            uint64_t k = 0;
            for (uint64_t i = 0; i < Psi.D_index.size(); ++i) {
                for (uint64_t j = 0; j < Flens[i]; ++j)
                    std::cout 
                        << k++ << '\t'
                        << i << '\t'
                        << F[i] << '\t'
                        << (*Psi.intLens)[i] << '\t'
                        << Psi.D_index[i] << '\t'
                        << Psi.D_offset[i] << '\n';
            }
        }
        {
            MoveStructure::IntervalPoint end{static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            curr = end;
            char con[] = "$ACGTN";
            uint64_t count = 0;
            do {
                curr = Psi.map(curr);
                std::cout << con[F[curr.interval]];
                ++count;
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }

        Timer.start("Verifying Psi");
        if (!Psi.permutationLengthN(totalLen)) {
            std::cerr << "ERROR: Psi is not a permutation of length n!" << std::endl;
            exit(1);
        }
        std::cout << "Psi is a permutation of length n\n";
        Timer.stop(); //Verifying Psi
        */

        uint64_t sampleInterval;
        sdsl::int_vector<> Psi_Index_Samples, Psi_Offset_Samples;
        {
            auto event = sdsl::memory_monitor::event("Construct Phi and Equidistant ISA Samples");
            ConstructPhiAndSamples(F, Psi, PhiIntLen, Psi.data.c, numTopRuns, seqLens, intAtTop, numSequences, maxPhiIntLen, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples);
        }

        /*
        {
            uint64_t count = 0, inInt = static_cast<uint64_t>(-1), outInt = static_cast<uint64_t>(-1);
            MoveStructure::IntervalPoint end = {static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            std::cout << "Phi text order\ni\tinInt\toutInt\tPsiInt\tPsiOff\tPhiInt\tPhiOff\tPhiLen\tPhiDind\tPhiDoff\n";
            std::cout << std::endl;
            curr = end;
            MoveStructure::IntervalPoint pPoint = {static_cast<uint64_t>(-1),0,0}, temp;
            do {
                inInt += (curr.offset == 0);
                outInt += (curr.offset == (*Psi.intLens)[curr.interval] - 1);
                curr = Psi.map(curr);
                temp = Phi.map(pPoint);
                std::cout << count++ << '\t'
                    << inInt << '\t'
                    << outInt << '\t'
                    << curr.interval << '\t'
                    << curr.offset << '\t'
                    << pPoint.interval << '\t'
                    << pPoint.offset << '\t'
                    << (*Phi.intLens)[pPoint.interval] << '\t' 
                    << Phi.D_index[pPoint.interval] << '\t'
                    << Phi.D_offset[pPoint.interval]+pPoint.offset << '\t'
                    << temp.interval << '\t'
                    << temp.offset << std::endl;
                ++pPoint.offset;
                if (pPoint.offset == (*Phi.intLens)[pPoint.interval]){
                    pPoint.offset = 0;
                    ++pPoint.interval;
                }
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }

        for (uint64_t i = 0; i < Phi.D_index.size(); ++i) {
            std::cout << Phi.D_index[i] << '\t'
                << Phi.D_offset[i] << '\t'
                << (*Phi.intLens)[i] << '\n';
        }

        {
            MoveStructure::IntervalPoint end = {static_cast<uint64_t>(-1), intAtTop[0], 0}, curr;
            curr = end;
            std::cout << "intAtTop\n";
            std::cout << curr.interval << ' ' << curr.offset << '\n';
            std::cout << "etc\n";
            do {
                curr = Phi.map(curr);
                std::cout << curr.interval << ' ' << curr.offset << '\n';
            } while (curr != end);
        }
        */

        /*
        Verifying Phi is VERY slow compared to verifying Psi ~300 seconds on mtb152 with 64 cores vs ~30 seconds for Psi on coombs c0-4. Why?
        Timer.start("Verifying Phi");
        if (!PhiNEWONEHE.permutationLengthN<EXPONENTIAL>(totalLen)) {
            std::cerr << "ERROR: Phi is not a permutation of length n!" << std::endl;
            exit(1);
        }
        std::cout << "Phi is a permutation of length n\n";
        Timer.stop(); //Verifying Phi
        */
        sdsl::int_vector<> intAtEnd;
        {
            auto event = sdsl::memory_monitor::event("Recover intAtEnd from disk");
            Timer.start("Recover intAtEnd from disk");
            intAtTop = sdsl::int_vector<>();
            tempOutFile.close();
            std::ifstream tempInFile(safeTempName);
            if (!tempInFile.is_open()) {
                std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' failed to open for reading!" << std::endl;
                exit(1);
            }
            sdsl::load(intAtEnd, tempInFile);
            tempInFile.close();
            Timer.stop(); //Recover intAtEnd from disk
        }


        {
            auto event = sdsl::memory_monitor::event("Compute PLCP Samples");
            ComputePLCPSamples(intAtEnd, numSequences, F, Psi, numTopRuns, seqLens, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples);
        }

        /*
        {
            Timer.start("Verifying Psi");
            {
                auto event = sdsl::memory_monitor::event("Verifying Psi");
                if (!Psi.permutationLengthN(totalLen)) {
                    std::cerr << "ERROR: Psi is not a permutation of length n!" << std::endl;
                    exit(1);
                }
                std::cout << "Psi is a permutation of length n\n";
            }
            Timer.stop(); //Verifying Psi
            //Verifying Phi is VERY slow compared to verifying Psi ~300 seconds on mtb152 with 64 cores vs ~30 seconds for Psi on coombs c0-4. Why?
            Timer.start("Verifying Phi");
            {
                auto event = sdsl::memory_monitor::event("Verifying Phi");
                if (!Phi.permutationLengthN(totalLen)) {
                    std::cerr << "ERROR: Phi is not a permutation of length n!" << std::endl;
                    exit(1);
                }
                std::cout << "Phi is a permutation of length n\n";
            }
            Timer.stop(); //Verifying Phi
        }

        //printPhiAndLCP(Phi, PLCPsamples);
        
        Timer.start("Computing minLCP per run");
        {
            auto event = sdsl::memory_monitor::event("Computing minLCP per run");
            //Psi = MoveStructure();
            Psi_Index_Samples = sdsl::int_vector<>();
            Psi_Offset_Samples = sdsl::int_vector<>();
            

            intAtTop = sdsl::int_vector<>(F.size(), 0, sdsl::bits::hi(F.size() - 1) + 1);
            //reinvert intAtEnd
            for (uint64_t i = 0; i < F.size(); ++i)
                intAtTop[intAtEnd[i]] = i;
            intAtEnd = sdsl::int_vector<>();
            
            ComputeMinLCPRun(intAtTop, F, Psi, lcpOut);
        }
        Timer.stop(); //Computing minLCP per run"
        */

        sdsl::memory_monitor::stop();
        std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;

        std::ofstream cstofs("construction.html");
        std::cout << "writing memory usage visualization to construction.html\n";
        sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
        cstofs.close();

        /*
        needs to be rewritten, affects bench
        std::ofstream outpsi("lcpcomputer.psi");
        //Flens.serialize(outpsi);
        Psi.serialize(outpsi);
        std::ofstream outphi("lcpcomputer.phi");
        (*Phi.intLens).serialize(outphi);
        Phi.serialize(outphi);
        */
    }

    void printRaw(const sdsl::int_vector<>& intAtTop) const {
        std::cout << "LCP\n";
        std::vector<uint64_t> lcp(totalLen);
        MoveStructureStartTable::IntervalPoint phiPoint{static_cast<uint64_t>(-1), intAtTop[0], 0};
        phiPoint.offset = PhiNEWONEHE.data.get<2>(phiPoint.interval) - 1;
        phiPoint = PhiNEWONEHE.map(phiPoint);
        for (uint64_t i = 0; i < totalLen; ++i) {
            lcp[totalLen - 1 - i] = PLCPsamples[phiPoint.interval] - phiPoint.offset;
            phiPoint = PhiNEWONEHE.map(phiPoint);
        }
        for (uint64_t i = 0; i < totalLen; ++i) {
            std::cout << lcp[i] << '\n';
        }

        std::cout << "\n\n";
        std::cout << "PLCP";
        phiPoint = {static_cast<uint64_t>(-1), 0, 0};
        for (uint64_t i = 0; i < totalLen; ++i) {
            std::cout << '\t' << PLCPsamples[phiPoint.interval] - phiPoint.offset;
            ++phiPoint.offset;
            phiPoint.offset %= PhiNEWONEHE.data.get<2>(phiPoint.interval);
            phiPoint.interval += (phiPoint.offset == 0);
        }
        std::cout << "\n";
    }

    void printPhiAndLCP(const sdsl::int_vector<>& PLCPsamples) const {
        std::cout << "runInd\tD_index\tD_offset\tintlen\tPLCPsamp\n";
        uint64_t numRuns = PLCPsamples.size();
        for (uint64_t i = 0; i < numRuns; ++i) {
            std::cout << i << '\t'
                << PhiNEWONEHE.data.get<0>(i) << '\t'
                << PhiNEWONEHE.data.get<1>(i) << '\t'
                << PhiNEWONEHE.data.get<2>(i) << '\t'
                << PLCPsamples[i] << '\n';
        }
    }

    static bool validateRB3(const rb3_fmi_t* rb3);

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;

        bytes += sdsl::serialize(totalLen, out, child, "totalLen");
        bytes += sdsl::serialize(F, out, child, "F");
        //bytes += sdsl::serialize(Flens, out, child, "Flens");
        bytes += sdsl::serialize(Psi, out, child, "Psi");
        bytes += sdsl::serialize(intAtTop, out, child, "intAtTop");
        //bytes += sdsl::serialize(PhiIntLen, out, child, "PhiIntLen");
        bytes += sdsl::serialize(PhiNEWONEHE, out, child, "Phi");
        bytes += sdsl::serialize(PLCPsamples, out, child, "PLCPsamples");

        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(totalLen, in);
        sdsl::load(F, in);
        //sdsl::load(Flens, in);
        sdsl::load(Psi, in);
        sdsl::load(intAtTop, in);
        //sdsl::load(PhiIntLen, in);
        sdsl::load(PhiNEWONEHE, in);
        sdsl::load(PLCPsamples, in);
        //Phi.intLens = &PhiIntLen;
        //Psi.intLens = &Flens;
    }
};

bool LCPComputer::validateRB3(const rb3_fmi_t* rb3){
    if (!rb3->e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return false;
    }
    return true;
}
