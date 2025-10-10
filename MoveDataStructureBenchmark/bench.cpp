#include"optbwtrl/util.h"
#include"optbwtrl/moveStructure.h"
#include<sdsl/sd_vector.hpp>
#include<sdsl/select_support.hpp>

enum FFMethod { LINEAR, EXPONENTIAL };

enum storageMethod { PACKEDINT, SPARSEBV };

template<storageMethod s>
struct MoveStructureStart;

template<>
struct MoveStructureStart<PACKEDINT> {
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
            curPos += (*mv.intLens)[i];
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

    //std::vector<uint64_t> t, v;
    //mv.MakeHistogram(t);
    MoveStructureStart<PACKEDINT> mvS(mv);
    //mvS.MakeHistogram(v);
    //if (t != v) {
        //std::cerr << "ERROR: Histograms not equal!" << std::endl;
        //std::cout << t.size() << ' ' << v.size() << std::endl;; 
        //exit(1);
    //}
    MoveStructureStart<SPARSEBV> mvSP(mv);
    if (!mvSP.equalTo(mvS)) {
        std::cerr << "MOVE STRUCTURE START SPARSEBV NOT EQUAL TO PACKEDINT!" << std::endl;
        exit(1);
    }
    
    uint64_t totLin, totExp;
    std::cout << "total operations for linear fast forwards: " << (totLin = mvS.sumOpsFastForward<LINEAR>()) << std::endl;
    std::cout << "total operations for exponential fast forwards: " << (totExp = mvS.sumOpsFastForward<EXPONENTIAL>()) << std::endl;

    std::cout << "n: " << N << " r: " << mv.D_index.size() 
        << " AvgFFLinear: " << static_cast<double>(totLin)/N 
        << " AvgFFExpon: " << static_cast<double>(totExp)/N << std::endl;

    /*
    {
        Timer.start("run length, int, linear scan");
        if (!mv.permutationLengthNOneCycleSequential(N)) {
            std::cerr << "ERROR: move data structure is not a permutation of length N with one cycle!" << std::endl;
            exit(1);
        }
        double time = Timer.stop(); //run length, int, linear scan
        std::cout << "Per operation: " << time/N << " seconds." << std::endl;
        uint64_t bytes = sdsl::size_in_bytes(mv);
        std::cout << "Struct size: " << bytes << " bytes.\n";
        std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
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
        std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
    }
    */

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
        std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
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
        std::cout << "ns\tMiB\n" << time/N*(1e9) << '\t' << static_cast<double>(bytes)/1024/1024 << std::endl;
    }
    /*
    sdsl::memory_monitor::stop();
    std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;

    std::ofstream cstofs(outputPref + ".construction.html");
    std::cout << "writing memory usage visualization to construction.html\n";
    sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
    cstofs.close();
    */

    Timer.stop(); //bench
}
