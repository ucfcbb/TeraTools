#include"optbwtrl/util.h"
//#include"optbwtrl/moveStructure.h"
#include<sdsl/sd_vector.hpp>
#include<sdsl/select_support.hpp>

struct packedTripleVector {
    typedef uint64_t size_type;
    private:
    const uint8_t a,b,c,width;
    public:
    sdsl::bit_vector data;
    packedTripleVector(const uint64_t a, const uint64_t b, const uint64_t c, const uint64_t elements): a(a), b(b), c(c), width(a+b+c), data(width*elements,0) {
        if (a > 64 || b > 64 || c > 64) {
            std::cerr << "ELEMENTS MUST BE <= 64 BITS!" << std::endl;
            exit(1);
        }
    }

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

enum FFMethod { LINEAR, EXPONENTIAL };

enum LINEARFFAccessMethod { RANDOM_ACCESS, ITERATOR };

enum storageMethod { PACKEDINT, SPARSEBV, SPARSECOMPBV };

enum initialInputIntervalRecoveryMethodSPARSECOMPBV { NA, RANK_vCOMPBV, RANK_v5COMPBV, RANK_scanCOMPBV, SELECT_mclCOMPBV, SELECT_scanCOMPBV };

struct MoveStructure {
    typedef uint64_t size_type;
    sdsl::int_vector<> intLens;
    sdsl::int_vector<> D_index;
    sdsl::int_vector<> D_offset;

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

    //NOTE: interval points returned by mapLF don't have valid position fields, they are set to -1
    //assumptions: inputs are valid and correspond to runs and runlens
    template<LINEARFFAccessMethod ac>
    IntervalPoint map(const IntervalPoint& intPoint) const;

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(intLens, out, child, "intLens");
        bytes += sdsl::serialize(D_index, out, child, "D_index");
        bytes += sdsl::serialize(D_offset, out, child, "D_offset");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    //WARNING!!!!!!!! THIS LOAD FUNCTION IS INCORRECT FOR THIS CLASS IMPLEMENTATION, IT IS WRITTEN FOR THE CURRENT IMPLEMENTED OUTPUT OF LCPCOMPUTER AND BUILDER MOVE STRUCTURE
    void load(std::istream& in) {
        sdsl::int_vector<> *dummy;
        sdsl::load(dummy, in);
        sdsl::load(D_index, in);
        sdsl::load(D_offset, in);
    }


    //assumptions:
    //no runs of length 0
    //the same runlens vector is passed for every call and unmodified
    //distance provided added to current position doesn't result in an out-of-bounds IntervalPoint 
    //  (except for position = n, the first position after the range
    void AdvanceIntervalPoint_unsafe(IntervalPoint& intPoint, uint64_t distance) const {
        uint64_t remaining = intPoint.offset + distance;
        //if intPoint.interval == intLens->size(), this will result in a BUG
        //TODO: FIX?
        while (remaining && remaining >= intLens[intPoint.interval])
            remaining -= intLens[intPoint.interval++];
        intPoint.offset = remaining;
        intPoint.position += distance;
    }

    void MakeHistogram(std::vector<uint64_t>& hist) const {
        //std::cout << "In MakeHistogram" << std::endl;
        hist = std::vector<uint64_t>();
        uint64_t intervals = intLens.size(), numTraversed;

        uint64_t interval, offset;
        for (uint64_t i = 0; i < intervals; ++i) {
            //std::cout << "HERE" << std::endl;
            //position of first position after interval, so not counted
            interval = D_index[i];
            offset = D_offset[i] + intLens[i];
            numTraversed = 0;

            if (numTraversed >= hist.size())
                hist.resize(numTraversed+1);
            hist[numTraversed] += std::min(offset, uint64_t(intLens[interval])) - D_offset[i];
            offset -= std::min(offset, uint64_t(intLens[interval]));

            while (offset) {
                ++interval;
                ++numTraversed;
                if (numTraversed >= hist.size())
                    hist.resize(numTraversed+1);
                hist[numTraversed] += std::min(offset, uint64_t(intLens[interval]));
                offset -= std::min(offset, uint64_t(intLens[interval]));
            }
        }
    }

    //returns whether the move structure is a permutation of length N
    //uses numIntervals log numIntervals auxiliary bits
    template<LINEARFFAccessMethod ac>
    bool permutationLengthN(size_type N) const {
        size_type runs = D_index.size();
        sdsl::int_vector<> nextInt(runs, runs, sdsl::bits::hi(runs) + 1);
        size_type totalOps = 0;

        #pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t i = 0; i < runs; ++i) {
            IntervalPoint curr{static_cast<uint64_t>(-1), i, 0};
            uint64_t ops = 0;
            do {
                curr = map<ac>(curr);
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

    //returns whether the move structure is a permutation of length N with exactly one cycle
    template<LINEARFFAccessMethod ac>
    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{static_cast<uint64_t>(-1), 0, 0};

        do {
            curr = map<ac>(curr);
            ++totalOps;
        } while ((curr.interval || curr.offset) && totalOps < N + 1);

        return totalOps == N;
    }
};

template<>
MoveStructure::IntervalPoint MoveStructure::map<RANDOM_ACCESS>(const IntervalPoint& intPoint) const {
    IntervalPoint res;
    res.position = static_cast<uint64_t>(-1);
    res.interval = D_index[intPoint.interval];
    res.offset = D_offset[intPoint.interval] + intPoint.offset;
    while (intLens[res.interval] <= res.offset)
        res.offset -= intLens[res.interval++];
    /*
       auto it = intLens->begin();
       it += static_cast<uint64_t>(res.interval);
       while (*it <= res.offset) {
       res.offset -= *it;
       ++res.interval;
       ++it;
       }
     */
    return res;
}

template<>
MoveStructure::IntervalPoint MoveStructure::map<ITERATOR>(const IntervalPoint& intPoint) const {
    IntervalPoint res;
    res.position = static_cast<uint64_t>(-1);
    res.interval = D_index[intPoint.interval];
    res.offset = D_offset[intPoint.interval] + intPoint.offset;
    auto it = intLens.begin() + static_cast<uint64_t>(res.interval);
    while (*it <= res.offset) {
        res.offset -= *it;
        ++res.interval;
        ++it;
    }
    return res;
}

struct MoveStructureTable {
    typedef uint64_t size_type;
    //D_index, D_offset, intlens
    packedTripleVector data;

    MoveStructureTable(const MoveStructure& mv): data(mv.D_index.width(), mv.D_offset.width(), mv.intLens.width(), mv.D_index.size()) {
        //std::cout << "calc: " << ((mv.D_index.width()+mv.D_offset.width() + mv.intLens->width())*mv.D_index.size()+7)/8 << std::endl;
        //std::cout << "Input index size: " << sdsl::size_in_bytes(mv) << std::endl;
        //std::cout << "This index size: " << sdsl::size_in_bytes(*this) << std::endl;
        //std::cout << "In constructor" << std::endl;
        //std::cout << data.data.size() << std::endl;
        for (uint64_t i = 0; i < mv.D_index.size(); ++i) {
            //std::cout << "i " << i << std::endl;
            data.set<0>(i, mv.D_index[i]);
            data.set<1>(i, mv.D_offset[i]);
            data.set<2>(i, mv.intLens[i]);
        }
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

    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{static_cast<uint64_t>(-1), 0, 0};

        do {
            curr = map(curr);
            ++totalOps;
        } while ((curr.interval || curr.offset) && totalOps < N + 1);

        return totalOps == N;
    }
};

template<storageMethod s, initialInputIntervalRecoveryMethodSPARSECOMPBV = NA>
struct MoveStructureStart;

template<>
struct MoveStructureStart<PACKEDINT> {
    typedef uint64_t size_type;
    //position i is the start position of the i-th input interval
    //start position has r+1 values, not r.
    sdsl::int_vector<> startPos;
    sdsl::int_vector<> D_index;
    sdsl::int_vector<> D_offset;

    MoveStructureStart(const MoveStructure& mv): D_index(mv.D_index), D_offset(mv.D_offset), startPos(mv.intLens) {
        uint64_t totalLen = 0;
        for (uint64_t i = 1; i < startPos.size(); ++i)
            totalLen += startPos[i];
        sdsl::util::expand_width(startPos, sdsl::bits::hi(totalLen) + 1);
        startPos.resize(startPos.size() + 1);

        uint64_t prevTemp = startPos[0], temp = 0;
        startPos[0] = 0;
        for (uint64_t i = 1; i < startPos.size(); ++i) {
            temp = startPos[i];
            startPos[i] = startPos[i-1] + prevTemp;
            prevTemp = temp;
        }
    }

    MoveStructureStart(MoveStructure&& mv) {
        D_index.swap(mv.D_index);
        D_offset.swap(mv.D_offset);
        startPos.swap(mv.intLens);
        uint64_t totalLen = 0;
        for (uint64_t i = 1; i < startPos.size(); ++i)
            totalLen += startPos[i];
        sdsl::util::expand_width(startPos, sdsl::bits::hi(totalLen) + 1);
        startPos.resize(startPos.size() + 1);

        uint64_t prevTemp = startPos[0], temp = 0;
        startPos[0] = 0;
        for (uint64_t i = 1; i < startPos.size(); ++i) {
            temp = startPos[i];
            startPos[i] = startPos[i-1] + prevTemp;
            prevTemp = temp;
        }
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
            return position != rhs.position;
        }

        bool operator<(const IntervalPoint& rhs) const {
            return position < rhs.position;
        }
    };

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(startPos, out, child, "startPos");
        bytes += sdsl::serialize(D_index, out, child, "D_index");
        bytes += sdsl::serialize(D_offset, out, child, "D_offset");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    template <FFMethod M>
    IntervalPoint map(const IntervalPoint& intPoint) const;

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

    void MakeHistogram(std::vector<uint64_t>& hist) const {
        hist = std::vector<uint64_t>();
        uint64_t intervals = D_index.size(), numTraversed;

        uint64_t interval, position;
        for (uint64_t i = 0; i < intervals; ++i) {
            interval = D_index[i];
            position = startPos[interval] + D_offset[i] + (startPos[i+1] - startPos[i]);
            numTraversed = 0;

            if (numTraversed >= hist.size())
                hist.resize(numTraversed+1);
            hist[numTraversed] += std::min(position - startPos[interval], startPos[interval + 1] - startPos[interval]) - D_offset[i];

            while (position > startPos[interval + 1]) {
                ++interval;
                ++numTraversed;
                if (numTraversed >= hist.size())
                    hist.resize(numTraversed+1);
                hist[numTraversed] += std::min(position - startPos[interval], startPos[interval + 1] - startPos[interval]);
            }
        }
    }

    template<FFMethod M>
    uint64_t sumOpsFastForward() const;
    uint64_t sumOpsFastForward() const {
        std::vector<uint64_t> hist;
        MakeHistogram(hist);

        uint64_t sum = 0;
        for (uint64_t i = 0; i < hist.size(); ++i)
            sum += i * hist[i];
        return sum;
    }
};

template<>
uint64_t MoveStructureStart<PACKEDINT>::sumOpsFastForward<LINEAR>() const {
    std::vector<uint64_t> hist;
    MakeHistogram(hist);

    uint64_t sum = 0;
    for (uint64_t i = 0; i < hist.size(); ++i)
        sum += i * hist[i];
    return sum;
}

template<>
uint64_t MoveStructureStart<PACKEDINT>::sumOpsFastForward<EXPONENTIAL>() const {
    std::vector<uint64_t> hist;
    MakeHistogram(hist);

    uint64_t sum = 0;
    uint64_t curLog = 0;
    uint64_t curMax = 1;
    for (uint64_t i = 1; i < hist.size(); ++i) {
        sum += (curLog + 1) * hist[i];
        if (i == curMax) {
            curMax *= 2;
            curLog += 1;
        }
    }
    return sum;
}


template<> 
MoveStructureStart<PACKEDINT>::IntervalPoint MoveStructureStart<PACKEDINT>::map<LINEAR>(const IntervalPoint& intPoint) const {
    //std::cout << "enter map" << std::endl;
    IntervalPoint res;
    res.interval = D_index[intPoint.interval];
    res.offset = D_offset[intPoint.interval] + intPoint.offset;
    res.position = startPos[res.interval] + res.offset;
    if (startPos[res.interval + 1] > res.position)
        return res;
    //std::cout << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
    //res.interval should never be > startPos.size()
    while (startPos[res.interval + 1] <= res.position)
        ++res.interval;
    res.offset = res.position - startPos[res.interval];
    //std::cout << "exit map" << std::endl;
    return res;
}

template<> 
MoveStructureStart<PACKEDINT>::IntervalPoint MoveStructureStart<PACKEDINT>::map<EXPONENTIAL>(const IntervalPoint& intPoint) const {
    //std::cout << "enter map" << std::endl;
    IntervalPoint res;
    res.interval = D_index[intPoint.interval];
    res.offset = D_offset[intPoint.interval] + intPoint.offset;
    res.position = startPos[res.interval] + res.offset;
    if (startPos[res.interval + 1] > res.position)
        return res;
    uint64_t dist = 2;
    //std::cout << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
    //res.interval should never be > startPos.size()
    while (dist < D_index.size() - res.interval && startPos[res.interval + dist] <= res.position)
        dist *= 2;
    res.interval += dist/2;
    
    dist = dist/2 - 1;
    using myint = uint64_t;

    myint lo = res.interval;
    myint hi = std::min(lo + dist, D_index.size() - 1);
    myint ans = -1;
    while(lo <= hi) {
        myint mid = (lo + hi) >> 1;
        myint here = startPos[mid];
        if (here <= res.position){
            lo = (ans = mid) + 1;
        } else {
            hi = mid - 1;
        }
    }
    
    res.interval = ans;
    res.offset = res.position - startPos[res.interval];
    //std::cout << "exit map" << std::endl;
    return res;
}

struct MoveStructureStartTable {
    typedef uint64_t size_type;
    //D_index, D_offset, startPos
    packedTripleVector data;
    
    MoveStructureStartTable(const MoveStructureStart<PACKEDINT>& mv): data(mv.D_index.width(), mv.D_offset.width(), mv.startPos.width(), mv.D_index.size() + 1) {
        //std::cout << "In constructor" << std::endl;
        //std::cout << data.data.size() << std::endl;
        for (uint64_t i = 0; i < mv.D_index.size(); ++i) {
            //std::cout << "i " << i << std::endl;
            data.set<0>(i, mv.D_index[i]);
            data.set<1>(i, mv.D_offset[i]);
            data.set<2>(i, mv.startPos[i]);
        }
        data.set<2>(mv.D_index.size(), mv.startPos[mv.D_index.size()]);
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
    
    template <FFMethod M>
    IntervalPoint map(const IntervalPoint& intPoint) const;

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(data, out, child, "data");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
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

template<>
struct MoveStructureStart<SPARSEBV> {
    typedef uint64_t size_type;
    //r+1 set bits
    sdsl::sd_vector<> inputStarts;
    //sdsl::sd_vector<>::select_1_type inpSLS;
    sdsl::select_support_sd<> inpSLS;
    sdsl::rank_support_sd<> inpRS;
    //r+1 set bits
    sdsl::sd_vector<> outputStarts;
    //sdsl::sd_vector<>::select_1_type outSLS;
    sdsl::select_support_sd<> outSLS;
    sdsl::int_vector<> outputRank;

    MoveStructureStart(const MoveStructure& mv) {
        //std::cout << "In constructor" << std::endl;
        std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> outputInputPairs(mv.D_index.size());
        std::vector<uint64_t> list, list2;
        outputRank = sdsl::int_vector<>(mv.D_index.size(), 0, mv.D_index.width());

        uint64_t curPos = 0;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            outputInputPairs[i].first = {mv.D_index[i], mv.D_offset[i]};
            outputInputPairs[i].second = i;
            list.push_back(curPos);
            curPos += mv.intLens[i];
        }
        list.push_back(curPos);
        //std::cout << "Done list" << std::endl;

        inputStarts = sdsl::sd_vector<>(list.begin(), list.end());

        std::sort(outputInputPairs.begin(), outputInputPairs.end());

        //std::cout << "Start list2" << std::endl;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            auto &a = outputInputPairs[i];
            list2.push_back(list[a.first.first] + a.first.second);
            outputRank[a.second] = i;
        }
        list2.push_back(curPos);
        //std::cout << "Done list2" << std::endl;
        outputStarts = sdsl::sd_vector<>(list2.begin(), list2.end());
        //std::cout << "done constructor" << std::endl;

        //inpSLS = sdsl::sd_vector<>::select_1_type(&inputStarts);
        //outSLS = sdsl::sd_vector<>::select_1_type(&outputStarts);
        inpSLS = sdsl::select_support_sd<>(&inputStarts);
        inpRS = sdsl::rank_support_sd<>(&inputStarts);
        outSLS = sdsl::select_support_sd<>(&outputStarts);
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
            return position != rhs.position;
        }

        bool operator<(const IntervalPoint& rhs) const {
            return position < rhs.position;
        }
    };

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(inputStarts, out, child, "inputStarts");
        bytes += sdsl::serialize(inpSLS, out, child, "inpSLS");
        bytes += sdsl::serialize(outputStarts, out, child, "outputStarts");
        bytes += sdsl::serialize(outSLS, out, child, "outSLS");
        bytes += sdsl::serialize(outputRank, out, child, "outputRank");

        return bytes;
    }

    IntervalPoint map(const IntervalPoint& intPoint) const {
        //std::cout << "map " << intPoint.position << ' ' << intPoint.interval << ' ' << intPoint.offset << std::endl;
        IntervalPoint res;
        uint64_t t = outSLS.select(outputRank[intPoint.interval] + 1) + intPoint.offset;
        //std::cout << "t " << t << std::endl;
        res.position = t;
        res.interval = inpRS(t + 1) - 1;
        //std::cout << "res.interval " << res.interval << std::endl;
        res.offset = t - inpSLS(res.interval + 1);
        //std::cout << "res " << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
        return res;
    }

    bool equalTo(const MoveStructureStart<PACKEDINT>& mv) {
        //std::cout << "In equal to" << std::endl;
        //std::cout << "mv.startPos.size() " << mv.startPos.size() << std::endl;
        //std::cout << "outputRank.size() " << outputRank.size() << std::endl;
        if (mv.D_index.size() != outputRank.size())
            return false;
        for (uint64_t i = 0; i < mv.D_index.size(); ++i) {
            //std::cout << "In loop" << std::endl;
            if (mv.startPos[i] != inpSLS.select(i + 1))
                return false;
            //std::cout << "In loop2" << std::endl;
            uint64_t pos = mv.startPos[mv.D_index[i]] + mv.D_offset[i];
            uint64_t pos2 = outSLS.select(outputRank[i] + 1);
             
            //std::cout << "In loop3" << std::endl;
            if (pos != pos2)
                return false;
        }
        //std::cout << "Done for" << std::endl;
        return inpSLS.select(outputRank.size() + 1) == outSLS.select(outputRank.size() + 1) &&
            inpSLS.select(outputRank.size() + 1) == mv.startPos[outputRank.size()];
    }

    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{0, 0, 0};

        do {
            curr = map(curr);
            ++totalOps;
        } while (curr.position && totalOps < N + 1);

        return totalOps == N;
    }
};

template<>
struct MoveStructureStart<SPARSECOMPBV, NA>;

//2 options: store outputRank in output sparse bv (then use select on comp bv + select on input sparse bv) or 
//                 outputRank in input+output bv (then use rank on comp bv + select on input sparse bv)
template<>
struct MoveStructureStart<SPARSECOMPBV, SELECT_mclCOMPBV> {
    typedef uint64_t size_type;
    //r+1 set bits
    sdsl::sd_vector<> inputStarts;
    sdsl::select_support_sd<> inpSLS;
    //r+1 set bits
    sdsl::sd_vector<> outputStarts;
    sdsl::select_support_sd<> outSLS;
    sdsl::int_vector<> outputRank;
    sdsl::bit_vector interleavedInputOutputStarts;
    sdsl::select_support_mcl<> SLSIO;
    /*
    if constexpr (inp == RANK_vCOMPBV)
        sdsl::rank_support_v<> r1;
    if constexpr (inp == RANK_V5COMPBV) {
        sdsl::rank_support_v5<> r2;
    }
    if constexpr (inp = RANK_scanCOMPBV) {
        sdsl::rank_support_scan<> r3;
    }
    if constexpr (inp == SELECT_mclCOMPBV) {
    }
    if constexpr (inp == SELECT_scanCOMPBV) {
        sdsl::select_support_scan<> s2;
    }
    */
    
    MoveStructureStart(const MoveStructure& mv) {
        //std::cout << "In constructor" << std::endl;
        std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> outputInputPairs(mv.D_index.size());
        std::vector<uint64_t> list, list2;
        outputRank = sdsl::int_vector<>(mv.D_index.size(), 0, mv.D_index.width());

        uint64_t curPos = 0;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            outputInputPairs[i].first = {mv.D_index[i], mv.D_offset[i]};
            outputInputPairs[i].second = i;
            list.push_back(curPos);
            curPos += mv.intLens[i];
        }
        list.push_back(curPos);
        //std::cout << "Done list" << std::endl;

        inputStarts = sdsl::sd_vector<>(list.begin(), list.end());

        std::sort(outputInputPairs.begin(), outputInputPairs.end());

        //std::cout << "Start list2" << std::endl;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            auto &a = outputInputPairs[i];
            list2.push_back(list[a.first.first] + a.first.second);
            outputRank[a.second] = i;
        }
        list2.push_back(curPos);
        //std::cout << "Done list2" << std::endl;
        outputStarts = sdsl::sd_vector<>(list2.begin(), list2.end());
        //std::cout << "done constructor" << std::endl;

        //inpSLS = sdsl::sd_vector<>::select_1_type(&inputStarts);
        //outSLS = sdsl::sd_vector<>::select_1_type(&outputStarts);
        inpSLS = sdsl::select_support_sd<>(&inputStarts);
        outSLS = sdsl::select_support_sd<>(&outputStarts);
        interleavedInputOutputStarts = sdsl::bit_vector(2 * outputInputPairs.size(), 0);
        uint64_t inIndex = 0, outIndex = 0;
        for (uint64_t i = 0; i < 2*outputInputPairs.size(); ++i) {
            if (inpSLS.select(inIndex+1) <= outSLS.select(outIndex+1)) {
                interleavedInputOutputStarts[i] = 0;
                ++inIndex;
            }
            else {
                interleavedInputOutputStarts[i] = 1;
                ++outIndex;
            }
        }
        SLSIO = sdsl::select_support_mcl<>(&interleavedInputOutputStarts);
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
            return position != rhs.position;
        }

        bool operator<(const IntervalPoint& rhs) const {
            return position < rhs.position;
        }
    };

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(inputStarts, out, child, "inputStarts");
        bytes += sdsl::serialize(inpSLS, out, child, "inpSLS");
        bytes += sdsl::serialize(outputStarts, out, child, "outputStarts");
        bytes += sdsl::serialize(outSLS, out, child, "outSLS");
        bytes += sdsl::serialize(outputRank, out, child, "outputRank");
        bytes += sdsl::serialize(interleavedInputOutputStarts, out, child, "interleavedInputOutputStarts");
        bytes += sdsl::serialize(SLSIO, out, child, "SLSIO");

        return bytes;
    }


    IntervalPoint map(const IntervalPoint& intPoint) const {
        //std::cout << "map " << intPoint.position << ' ' << intPoint.interval << ' ' << intPoint.offset << std::endl;
        IntervalPoint res;
        //std::cout << "outputRank[intPoint.interval] " << outputRank[intPoint.interval] << std::endl;
        uint64_t inRank = SLSIO.select(outputRank[intPoint.interval] + 1) - outputRank[intPoint.interval] - 1;
        //std::cout << "inRank " << inRank << std::endl;
        res.position = outSLS.select(outputRank[intPoint.interval] + 1) + intPoint.offset;
        //std::cout << "res.position " << res.position << std::endl;
        while (inRank < outputRank.size() - 1 && inpSLS.select(inRank + 2) <= res.position)
            ++inRank;
        //std::cout << "t " << t << std::endl;
        //std::cout << "res.interval " << res.interval << std::endl;
        res.interval = inRank;
        //std::cout << "res.interval " << res.interval << std::endl;
        res.offset = res.position - inpSLS(res.interval + 1);
        //std::cout << "res.offset " << res.offset << std::endl;
        //std::cout << "res " << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
        return res;
    }

    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{0, 0, 0};

        do {
            curr = map(curr);
            ++totalOps;
        } while (curr.position && totalOps < N + 1);

        return totalOps == N;
    }
};

template<>
struct MoveStructureStart<SPARSECOMPBV, SELECT_scanCOMPBV> {
    typedef uint64_t size_type;
    //r+1 set bits
    sdsl::sd_vector<> inputStarts;
    sdsl::select_support_sd<> inpSLS;
    //r+1 set bits
    sdsl::sd_vector<> outputStarts;
    sdsl::select_support_sd<> outSLS;
    sdsl::int_vector<> outputRank;
    sdsl::bit_vector interleavedInputOutputStarts;
    sdsl::select_support_mcl<> SLSIO;
    /*
    if constexpr (inp == RANK_vCOMPBV)
        sdsl::rank_support_v<> r1;
    if constexpr (inp == RANK_V5COMPBV) {
        sdsl::rank_support_v5<> r2;
    }
    if constexpr (inp = RANK_scanCOMPBV) {
        sdsl::rank_support_scan<> r3;
    }
    if constexpr (inp == SELECT_mclCOMPBV) {
    }
    if constexpr (inp == SELECT_scanCOMPBV) {
        sdsl::select_support_scan<> s2;
    }
    */
    
    MoveStructureStart(const MoveStructure& mv) {
        //std::cout << "In constructor" << std::endl;
        std::vector<std::pair<std::pair<uint64_t, uint64_t>, uint64_t>> outputInputPairs(mv.D_index.size());
        std::vector<uint64_t> list, list2;
        outputRank = sdsl::int_vector<>(mv.D_index.size(), 0, mv.D_index.width());

        uint64_t curPos = 0;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            outputInputPairs[i].first = {mv.D_index[i], mv.D_offset[i]};
            outputInputPairs[i].second = i;
            list.push_back(curPos);
            curPos += mv.intLens[i];
        }
        list.push_back(curPos);
        //std::cout << "Done list" << std::endl;

        inputStarts = sdsl::sd_vector<>(list.begin(), list.end());

        std::sort(outputInputPairs.begin(), outputInputPairs.end());

        //std::cout << "Start list2" << std::endl;
        for (uint64_t i = 0; i < outputInputPairs.size(); ++i) {
            auto &a = outputInputPairs[i];
            list2.push_back(list[a.first.first] + a.first.second);
            outputRank[a.second] = i;
        }
        list2.push_back(curPos);
        //std::cout << "Done list2" << std::endl;
        outputStarts = sdsl::sd_vector<>(list2.begin(), list2.end());
        //std::cout << "done constructor" << std::endl;

        //inpSLS = sdsl::sd_vector<>::select_1_type(&inputStarts);
        //outSLS = sdsl::sd_vector<>::select_1_type(&outputStarts);
        inpSLS = sdsl::select_support_sd<>(&inputStarts);
        outSLS = sdsl::select_support_sd<>(&outputStarts);
        interleavedInputOutputStarts = sdsl::bit_vector(2 * outputInputPairs.size(), 0);
        uint64_t inIndex = 0, outIndex = 0;
        for (uint64_t i = 0; i < 2*outputInputPairs.size(); ++i) {
            if (inpSLS.select(inIndex+1) <= outSLS.select(outIndex+1)) {
                interleavedInputOutputStarts[i] = 0;
                ++inIndex;
            }
            else {
                interleavedInputOutputStarts[i] = 1;
                ++outIndex;
            }
        }
        SLSIO = sdsl::select_support_mcl<>(&interleavedInputOutputStarts);
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
            return position != rhs.position;
        }

        bool operator<(const IntervalPoint& rhs) const {
            return position < rhs.position;
        }
    };

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(inputStarts, out, child, "inputStarts");
        bytes += sdsl::serialize(inpSLS, out, child, "inpSLS");
        bytes += sdsl::serialize(outputStarts, out, child, "outputStarts");
        bytes += sdsl::serialize(outSLS, out, child, "outSLS");
        bytes += sdsl::serialize(outputRank, out, child, "outputRank");
        bytes += sdsl::serialize(interleavedInputOutputStarts, out, child, "interleavedInputOutputStarts");
        bytes += sdsl::serialize(SLSIO, out, child, "SLSIO");

        return bytes;
    }


    IntervalPoint map(const IntervalPoint& intPoint) const {
        //std::cout << "map " << intPoint.position << ' ' << intPoint.interval << ' ' << intPoint.offset << std::endl;
        IntervalPoint res;
        //std::cout << "outputRank[intPoint.interval] " << outputRank[intPoint.interval] << std::endl;
        uint64_t inRank = SLSIO.select(outputRank[intPoint.interval] + 1) - outputRank[intPoint.interval] - 1;
        //std::cout << "inRank " << inRank << std::endl;
        res.position = outSLS.select(outputRank[intPoint.interval] + 1) + intPoint.offset;
        //std::cout << "res.position " << res.position << std::endl;
        while (inRank < outputRank.size() - 1 && inpSLS.select(inRank + 2) <= res.position)
            ++inRank;
        //std::cout << "t " << t << std::endl;
        //std::cout << "res.interval " << res.interval << std::endl;
        res.interval = inRank;
        //std::cout << "res.interval " << res.interval << std::endl;
        res.offset = res.position - inpSLS(res.interval + 1);
        //std::cout << "res.offset " << res.offset << std::endl;
        //std::cout << "res " << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
        return res;
    }

    bool permutationLengthNOneCycleSequential(size_type N) const {
        size_type totalOps = 0;

        IntervalPoint curr{0, 0, 0};

        do {
            curr = map(curr);
            ++totalOps;
        } while (curr.position && totalOps < N + 1);

        return totalOps == N;
    }
};

MoveStructure generateBalanced(const MoveStructure& mv, const uint64_t N, const uint64_t d) {
    MoveStructure res;

    res = mv;
    return res;
}

int main(int argc, char* argv[]) {
    Timer.start("bench");
    if (argc != 2 && argc != 3) {
        std::cerr << "1 or 2 parameter needs to be passed: the name of the file containing the move structure"
            << " to be benchmarked (and optionally the outputPref of the outputFiles)! " << argc - 1 << " passed!\n";
        exit(1);
    }
    std::string outputPref(argv[argc-1]);
    std::cout << "Timing move structures for move structure found at: " << argv[1] << std::endl;
    sdsl::memory_monitor::start(); 

    MoveStructure mvOG;
    sdsl::int_vector<> Lens;
    {
        Timer.start("Reading move structure");
        auto event = sdsl::memory_monitor::event("Reading move structure");
        std::ifstream mvIn(argv[1]);
        if (!mvIn.is_open()) {
            std::cerr << "Move structure file '" << argv[1] << "' failed to open!" << std::endl;
            exit(1);
        }
        Lens.load(mvIn);
        mvOG.load(mvIn);
        mvOG.intLens = Lens;
        Timer.stop(); //Reading move structure
    }


    //get maximum fast forwards
    uint64_t maxFF = 0;
    {
        std::vector<uint64_t> t;
        mvOG.MakeHistogram(t);
        maxFF = t.size() - 1;
    }
    uint64_t firstD = 2;
    while (2*firstD <= maxFF)
        firstD *= 2;
    std::cout << "Maximum fast forwards from passed move structure: " << maxFF << std::endl;
    std::cout << "Minimum balancing parameter d (that is a power of two) s.t. input move structure is d-balanced: " << firstD << std::endl; 

    std::vector<std::pair<uint64_t,std::vector<double>>> allResults;
    uint64_t currD = firstD;
    while (currD != 1) {
        Timer.start("Timing " + std::string((currD == firstD)? "unbalanced" : "balanced" ) + " move structures, d = " + std::to_string(std::min(currD, (maxFF+1)/2)));



        uint64_t N = 0;
        for (uint64_t i = 0; i < mvOG.D_index.size(); ++i)
            N += mvOG.intLens[i];

        MoveStructure mv;
        if (currD == firstD)
            mv = mvOG;
        else 
            mv = generateBalanced(mvOG, N, currD);

        //std::vector<uint64_t> t, v;
        //mvOG.MakeHistogram(t);
        MoveStructureStart<PACKEDINT> mvS(mvOG);
        //mvS.MakeHistogram(v);
        //if (t != v) {
        //std::cerr << "ERROR: Histograms not equal!" << std::endl;
        //std::cout << t.size() << ' ' << v.size() << std::endl;; 
        //exit(1);
        //}
        MoveStructureStart<SPARSEBV> mvSP(mvOG);
        if (!mvSP.equalTo(mvS)) {
            std::cerr << "MOVE STRUCTURE START SPARSEBV NOT EQUAL TO PACKEDINT!" << std::endl;
            exit(1);
        }

        uint64_t totLin, totExp;
        std::cout << "total operations for linear fast forwards: " << (totLin = mvS.sumOpsFastForward<LINEAR>()) << std::endl;
        std::cout << "total operations for exponential fast forwards: " << (totExp = mvS.sumOpsFastForward<EXPONENTIAL>()) << std::endl;

        std::cout << "n: " << N << " r: " << mvOG.D_index.size() 
            << " AvgFFLinear: " << static_cast<double>(totLin)/N 
            << " AvgFFExpon: " << static_cast<double>(totExp)/N << std::endl;

        std::vector<double> results;

        {
            Timer.start("run length, int, linear scan with random access");
            if (!mv.permutationLengthNOneCycleSequential<RANDOM_ACCESS>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //run length, int, linear scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mv);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            MoveStructureTable mvT(mv);
            Timer.start("run length, int, linear scan with random access, tablefied");
            if (!mvT.permutationLengthNOneCycleSequential(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //run length, int, linear scan, tablefied
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvT);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            Timer.start("run length, int, linear scan with iterator");
            if (!mv.permutationLengthNOneCycleSequential<ITERATOR>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //run length, int, linear scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mv);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            Timer.start("start pos, int, linear scan");
            if (!mvS.permutationLengthNOneCycleSequential<LINEAR>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, linear scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvS);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            MoveStructureStartTable mvST = MoveStructureStartTable(MoveStructureStart<PACKEDINT>(mv));
            Timer.start("start pos, int, linear scan, tablefied");
            if (!mvST.permutationLengthNOneCycleSequential<LINEAR>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, linear scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvST);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            Timer.start("start pos, int, exponential scan");
            if (!mvS.permutationLengthNOneCycleSequential<EXPONENTIAL>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, exponential scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvS);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            MoveStructureStartTable mvST = MoveStructureStartTable(MoveStructureStart<PACKEDINT>(mv));
            Timer.start("start pos, int, exponential scan, tablefied");
            if (!mvST.permutationLengthNOneCycleSequential<EXPONENTIAL>(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, exponential scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvST);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            Timer.start("start pos, sparse bv, rank/select on sparse bv");
            if (!mvSP.permutationLengthNOneCycleSequential(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, exponential scan
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvSP);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            MoveStructureStart<SPARSECOMPBV, SELECT_mclCOMPBV> mvSP(mv);
            Timer.start("start pos, sparse bv, LINEAR, select_mcl on comp bv");
            if (!mvSP.permutationLengthNOneCycleSequential(N)) {
                std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
                exit(1);
            }
            double time = Timer.stop(); //start pos, int, select_mcl on comp bv
            std::cout << "Per operation: " << time/N << " seconds." << std::endl;
            uint64_t bytes = sdsl::size_in_bytes(mvSP);
            std::cout << "Struct size: " << bytes << " bytes.\n";
            results.push_back(time/N*(1e9));
            std::cout << "ns\tMiB\n" << results.back(); 
            results.push_back(static_cast<double>(bytes)/1024/1024); 
            std::cout << '\t' << results.back()  << std::endl;
        }

        {
            //MoveStructureStart<SPARSECOMPBV, RANK_vCOMPBV> mvSP(mv);
            /*
               Timer.start("start pos, sparse bv, LINEAR, rank_v on comp bv");
               if (!mvSP.permutationLengthNOneCycleSequential(N)) {
               std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
               exit(1);
               }
               double time = Timer.stop(); //start pos, int, exponential scan
               std::cout << "Per operation: " << time/N << " seconds." << std::endl;
               uint64_t bytes = sdsl::size_in_bytes(mvSP);
               std::cout << "Struct size: " << bytes << " bytes.\n";
               std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
             */
        }

        {
            //MoveStructureStart<SPARSECOMPBV, RANK_v5COMPBV> mvSP(mv);
            /*
               Timer.start("start pos, sparse bv, LINEAR, rank_v5 on comp bv");
               if (!mvSP.permutationLengthNOneCycleSequential(N)) {
               std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
               exit(1);
               }
               double time = Timer.stop(); //start pos, int, exponential scan
               std::cout << "Per operation: " << time/N << " seconds." << std::endl;
               uint64_t bytes = sdsl::size_in_bytes(mvSP);
               std::cout << "Struct size: " << bytes << " bytes.\n";
               std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
             */
        }

        {
            //MoveStructureStart<SPARSECOMPBV, RANK_scanCOMPBV> mvSP(mv);
            /*
               Timer.start("start pos, sparse bv, LINEAR, rank_scan on comp bv");
               if (!mvSP.permutationLengthNOneCycleSequential(N)) {
               std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
               exit(1);
               }
               double time = Timer.stop(); //start pos, int, exponential scan
               std::cout << "Per operation: " << time/N << " seconds." << std::endl;
               uint64_t bytes = sdsl::size_in_bytes(mvSP);
               std::cout << "Struct size: " << bytes << " bytes.\n";
               std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
             */
        }

        {
            //MoveStructureStart<SPARSECOMPBV, SELECT_scanCOMPBV> mvSP(mv);
            /*
               Timer.start("start pos, sparse bv, LINEAR, select_scan on comp bv");
               if (!mvSP.permutationLengthNOneCycleSequential(N)) {
               std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
               exit(1);
               }
               double time = Timer.stop(); //start pos, int, exponential scan
               std::cout << "Per operation: " << time/N << " seconds." << std::endl;
               uint64_t bytes = sdsl::size_in_bytes(mvSP);
               std::cout << "Struct size: " << bytes << " bytes.\n";
               std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
             */
        }
        /*
           sdsl::memory_monitor::stop();
           std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;

           std::ofstream cstofs(outputPref + ".construction.html");
           std::cout << "writing memory usage visualization to construction.html\n";
           sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
           cstofs.close();
         */

        for (auto d: results) 
            std::cout << '\t' << d;
        std::cout << std::endl;

        allResults.emplace_back(std::min(currD, (maxFF+1)/2), results);

        currD /= 2;
        Timer.stop(); //this d
    }

    std::cout << maxFF << std::endl;
    for (auto [b, r] : allResults) {
        std::cout << b;
        for (auto d: r) 
            std::cout << '\t' << d;
        std::cout << std::endl;
    }

    Timer.stop(); //bench
}
