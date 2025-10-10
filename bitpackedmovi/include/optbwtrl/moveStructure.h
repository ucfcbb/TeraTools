#include<sdsl/int_vector.hpp>
#include<cstdint>

struct MoveStructure {
    typedef uint64_t size_type;
    sdsl::int_vector<>* intLens;
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
    IntervalPoint map(const IntervalPoint& intPoint) const {
        IntervalPoint res;
        res.position = static_cast<uint64_t>(-1);
        res.interval = D_index[intPoint.interval];
        res.offset = D_offset[intPoint.interval] + intPoint.offset;
        while ((*intLens)[res.interval] <= res.offset)
            res.offset -= (*intLens)[res.interval++];
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

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(intLens, out, child, "intLens");
        bytes += sdsl::serialize(D_index, out, child, "D_index");
        bytes += sdsl::serialize(D_offset, out, child, "D_offset");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(intLens, in);
        intLens = nullptr;
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
        while (remaining && remaining >= (*intLens)[intPoint.interval])
            remaining -= (*intLens)[intPoint.interval++];
        intPoint.offset = remaining;
        intPoint.position += distance;
    }

    void MakeHistogram(std::vector<uint64_t>& hist) const {
        //std::cout << "In MakeHistogram" << std::endl;
        hist = std::vector<uint64_t>();
        uint64_t intervals = intLens->size(), numTraversed;

        uint64_t interval, offset;
        for (uint64_t i = 0; i < intervals; ++i) {
            //std::cout << "HERE" << std::endl;
            //position of first position after interval, so not counted
            interval = D_index[i];
            offset = D_offset[i] + (*intLens)[i];
            numTraversed = 0;

            if (numTraversed >= hist.size())
                hist.resize(numTraversed+1);
            hist[numTraversed] += std::min(offset, uint64_t((*intLens)[interval])) - D_offset[i];
            offset -= std::min(offset, uint64_t((*intLens)[interval]));

            while (offset) {
                ++interval;
                ++numTraversed;
                if (numTraversed >= hist.size())
                    hist.resize(numTraversed+1);
                hist[numTraversed] += std::min(offset, uint64_t((*intLens)[interval]));
                offset -= std::min(offset, uint64_t((*intLens)[interval]));
            }
        }
    }

    //returns whether the move structure is a permutation of length N
    //uses numIntervals log numIntervals auxiliary bits
    bool permutationLengthN(size_type N) const {
        size_type runs = D_index.size();
        sdsl::int_vector<> nextInt(runs, runs, sdsl::bits::hi(runs) + 1);
        size_type totalOps = 0;

        #pragma omp parallel for schedule(dynamic, 1)
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

    //returns whether the move structure is a permutation of length N with exactly one cycle
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
