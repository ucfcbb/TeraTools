#include"optbwtrl/util.h"
#include"optbwtrl/moveStructure.h"

enum FFMethod { LINEAR, EXPONENTIAL };

struct MoveStructureStart {
    typedef uint64_t size_type;
    //position i is the start position of the i-th input interval
    //start position has r+1 values, not r.
    sdsl::int_vector<> startPos;
    sdsl::int_vector<> D_index;
    sdsl::int_vector<> D_offset;

    MoveStructureStart(const MoveStructure& mv): D_index(mv.D_index), D_offset(mv.D_offset), startPos(*(mv.intLens)) {
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
        startPos.swap(*(mv.intLens));
        mv.intLens = nullptr;
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
};

template<> 
MoveStructureStart::IntervalPoint MoveStructureStart::map<LINEAR>(const IntervalPoint& intPoint) const {
    //std::cout << "enter map" << std::endl;
    IntervalPoint res;
    res.interval = D_index[intPoint.interval];
    res.offset = D_offset[intPoint.interval] + intPoint.offset;
    res.position = startPos[res.interval] + res.offset;
    //std::cout << res.position << ' ' << res.interval << ' ' << res.offset << std::endl;
    //res.interval should never be > startPos.size()
    while (startPos[res.interval + 1] <= res.position)
        ++res.interval;
    res.offset = res.position - startPos[res.interval];
    //std::cout << "exit map" << std::endl;
    return res;
}

template<> 
MoveStructureStart::IntervalPoint MoveStructureStart::map<EXPONENTIAL>(const IntervalPoint& intPoint) const {
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

    MoveStructure mv;
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
        mv.load(mvIn);
        mv.intLens = &Lens;
        Timer.stop(); //Reading move structure
    }

    uint64_t N = 0;
    for (uint64_t i = 0; i < mv.D_index.size(); ++i)
        N += (*mv.intLens)[i];

    std::cout << "n: " << N << " r: " << mv.D_index.size() << std::endl;

    {
        Timer.start("Timing move structure with run lengths");
        auto event = sdsl::memory_monitor::event("Run length packed integer move structure");
        if (!mv.permutationLengthNOneCycleSequential(N)) {
            std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
            exit(1);
        }
        double time = Timer.stop(); //Timing move structure with run lengths
        std::cout << "Per operation: " << time/N << " seconds." << std::endl;
    }

    MoveStructureStart mvS(std::move(mv));

    {
        Timer.start("Timing move structure with starting positions");
        auto event = sdsl::memory_monitor::event("Starting position packed integer move structure");
        if (!mvS.permutationLengthNOneCycleSequential<LINEAR>(N)) {
            std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
            exit(1);
        }
        double time = Timer.stop(); //Timing move structure with starting positions
        std::cout << "Per operations: " << time/N << " seconds." << std::endl;
    }

    {
        Timer.start("Timing move structure with starting positions with exponential ff");
        auto event = sdsl::memory_monitor::event("Starting position packed integer move structure exponential ff");
        if (!mvS.permutationLengthNOneCycleSequential<EXPONENTIAL>(N)) {
            std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
            exit(1);
        }
        double time = Timer.stop(); //Timing move structure with starting positions
        std::cout << "Per operations: " << time/N << " seconds." << std::endl;
    }
    sdsl::memory_monitor::stop();
    std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;

    std::ofstream cstofs(outputPref + ".construction.html");
    std::cout << "writing memory usage visualization to construction.html\n";
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
    cstofs.close();

    Timer.stop(); //bench
}
