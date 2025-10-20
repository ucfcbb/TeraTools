#ifndef R_SA_LCP_MOVESTRUCTURE_H
#define R_SA_LCP_MOVESTRUCTURE_H
#include "util/packedTripleVector.h"

enum FFMethod { LINEAR, EXPONENTIAL };

enum PositionStorageMethod { ABSOLUTE, RELATIVE };

//copied from bench.cpp, clean up structure of code later
struct MoveStructureTable {
    typedef uint64_t size_type;
    //D_index, D_offset, intlens
    packedTripleVector data;

    MoveStructureTable() = default;

    /*
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
    */

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

    /*
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
    */

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

#endif //#ifndef R_SA_LCP_MOVESTRUCTURE_H
