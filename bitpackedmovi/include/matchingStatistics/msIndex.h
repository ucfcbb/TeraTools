#ifndef R_SA_LCP_MSINDEX_H
#define R_SA_LCP_MSINDEX_H
#include <optional>
#include<sdsl/int_vector.hpp>
#include"moveStructure/moveStructure.h"

class MSIndex {
    uint64_t totalLen;
    sdsl::int_vector<> F, rlbwt;

    MoveStructureTable Psi, LF;

    sdsl::int_vector<> intAtTop;

    MoveStructureStartTable Phi, InvPhi;

    sdsl::int_vector<> PLCPsamples;

    void generateRLBWTfromLFPsiandF() {
        const uint64_t numRuns = LF.data.size();
        auto pointToInt = [numRuns] (MoveStructureTable::IntervalPoint a) -> uint64_t {
            return a.interval + a.offset*numRuns;
        };

        uint64_t numSequences = 0;
        while (numSequences < numRuns && F[numSequences] == 0) 
            ++numSequences;
        assert(numSequences);
        assert(numSequences != numRuns);
        rlbwt = F;
        //map from LF position to character
        std::unordered_map<uint64_t, uint64_t> currentStarts;
        //initialize currentStarts
        MoveStructureTable::IntervalPoint p{static_cast<uint64_t>(-1), 0, numSequences}, start;
        //std::cout << "starting initialization loop" << std::endl;
        for (uint64_t i = numSequences; i < numRuns; ++i) {
            while (p.interval < numRuns && p.offset >= LF.data.get<2>(p.interval))
                p.offset -= LF.data.get<2>(p.interval++);
            if (i == numSequences) { start = p; }
            if (F[i] != F[i - 1])
                currentStarts[pointToInt(p)] = F[i];
            p.offset += Psi.data.get<2>(i);
        }

        //std::cout << "starting rlbwt generation loop" << std::endl;
        //generate rlbwt
        for (uint64_t i = 0; i < numRuns; ++i) {
            auto p = LF.map({static_cast<uint64_t>(-1), i, 0});
            if (p.interval < start.interval || (p.interval == start.interval && p.offset < start.offset))
                rlbwt[i] = 0;
            else {
                assert(currentStarts.count(pointToInt(p)));
                rlbwt[i] = currentStarts[pointToInt(p)];
                currentStarts.erase(pointToInt(p));
            }
            p.offset += LF.data.get<2>(i);
            while (p.interval < numRuns && p.offset >= LF.data.get<2>(p.interval))
                p.offset -= LF.data.get<2>(p.interval++);
            assert(p.interval != numRuns || p.offset == 0);
            currentStarts[pointToInt(p)] = rlbwt[i];
        }
    }

    public:
    void constructFromLCPIndexFileWriteAndClear(std::ifstream& lcpIn, std::ofstream& MSIout, verbosity v = TIME,
            bool vLF = false,
            bool vPsi = false,
            bool vText = false,
            bool vPhi = false,
            bool vInvPhi = false
            ) {
        auto test = [&] (bool test, const MoveStructureTable& a, std::string name) {
            if (!test) return;
            if (v >= TIME) { Timer.start("Verifying " + name + " is a permutation of length " + std::to_string(totalLen)); }
                //a.print();
            if (a.permutationLengthN(totalLen))
                std::cout << name << " is a permutation of length " << totalLen << std::endl;
            else {
                std::cerr << name << " is not a permutation of length " << totalLen << "!" << std::endl;
                exit(1);
            }
            if (v >= TIME) { Timer.stop(); }
        };
        auto test2 = [&] (bool test, const MoveStructureStartTable& a, std::string name) {
            if (!test) return;
            if (v >= TIME) { Timer.start("Verifying " + name + " is a permutation of length " + std::to_string(totalLen)); }
            if (a.permutationLengthN<EXPONENTIAL>(totalLen))
                std::cout << name << " is a permutation of length " << totalLen << std::endl;
            else {
                std::cerr << name << " is not a permutation of length " << totalLen << "!" << std::endl;
                exit(1);
            }
            if (v >= TIME) { Timer.stop(); }
        };
        if (v >= TIME) { Timer.start("constructFromLCPIndexFileWriteAndClear"); }
        //load
        sdsl::load(totalLen, lcpIn);
        sdsl::load(F, lcpIn);
        sdsl::load(Psi, lcpIn);
        test(vPsi, Psi, "Psi");
        LF = Psi.invert();
        test(vLF, LF, "LF");
        generateRLBWTfromLFPsiandF();
        if (vText) {
            if (v >= TIME) { Timer.start("Recovering texts from LF and Psi"); }
            if (recoverTextLF() == recoverTextPsi()){
                std::cout << "LF and Psi generate equivalent texts." << std::endl;
            }
            else {
                std::cerr << "LF and Psi generate unequal texts!" << std::endl;
                //std::cout << recoverTextLF() << std::endl;
                //std::cout << recoverTextPsi() << std::endl;
                exit(1);
            }
            if (v >= TIME) { Timer.stop(); }
        }
        //write
        sdsl::serialize(totalLen, MSIout);
        sdsl::serialize(F, MSIout);
        sdsl::serialize(rlbwt, MSIout);
        sdsl::serialize(Psi, MSIout);
        sdsl::serialize(LF, MSIout);
        //clear
        rlbwt = F = sdsl::int_vector<>();
        LF = Psi = MoveStructureTable();

        //intAtTop
        sdsl::load(intAtTop, lcpIn);
        sdsl::serialize(intAtTop, MSIout);
        intAtTop = sdsl::int_vector<>();
        
        //load
        sdsl::load(Phi, lcpIn);
        //Phi.print();
        test2(vPhi, Phi, "Phi");
        InvPhi = Phi.invert();
        //InvPhi.print();
        test2(vInvPhi, InvPhi, "InvPhi");
        //write
        sdsl::serialize(Phi, MSIout);
        sdsl::serialize(InvPhi, MSIout);
        //clear
        Phi = InvPhi = MoveStructureStartTable();

        //PLCPsamples
        sdsl::load(PLCPsamples,lcpIn);
        sdsl::serialize(PLCPsamples, MSIout);
        PLCPsamples = sdsl::int_vector<>();
        if (v >= TIME) { Timer.stop(); }
    }

    std::string recoverTextLF() const {
        uint64_t numSequences = 0;
        while (numSequences < LF.data.size() && F[numSequences] == 0) 
            ++numSequences;
        std::string text = "", converter="$ACGTN";
        MoveStructureTable::IntervalPoint p = {static_cast<uint64_t>(-1), 0, numSequences - 1}, start;
        while (LF.data.get<2>(p.interval) <= p.offset)
            p.offset -= LF.data.get<2>(p.interval++);
        start = p;
        uint64_t prevChar = 0;
        do {
            text.push_back(converter[prevChar]);
            prevChar = rlbwt[p.interval];
            p = LF.map(p);
        } while (p != start);
        std::reverse(text.begin(), text.end());
        return text;
    }

    std::string recoverTextPsi() const {
        uint64_t numSequences = 0;
        while (numSequences < Psi.data.size() && F[numSequences] == 0) 
            ++numSequences;
        std::string text = "", converter="$ACGTN";
        MoveStructureTable::IntervalPoint p = {static_cast<uint64_t>(-1), numSequences - 1, 0}, start;
        start = p;
        do {
            p = Psi.map(p);
            text.push_back(converter[F[p.interval]]);
        } while (p != start);
        return text;
    }

using LF_IntervalPoint = MoveStructureTable::IntervalPoint;
using Phi_IntervalPoint = MoveStructureStartTable::IntervalPoint;
    
    // Make static contexpr, but this is probably fine
    static uint8_t charToBits(const char c) {
        switch (c) {
            case '\0': return 0;
            case  'A': return 1;
            case  'C': return 2;
            case  'G': return 3;
            case  'T': return 4;
            default: throw std::invalid_argument("Invalid character: " + std::string(1, c));
        };
    }


    LF_IntervalPoint start_bwt_pos() const {
        return LF_IntervalPoint{static_cast<uint64_t>(-1), 0, 0};
    }

    LF_IntervalPoint end_bwt_pos() const {
        return LF_IntervalPoint{static_cast<uint64_t>(-1), LF.num_intervals() - 1, LF.get_length(LF.num_intervals() - 1) - 1};
    }
    
    std::optional<LF_IntervalPoint> pred_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        uint64_t interval = pos.interval;
        while (rlbwt[interval] != c) {
            if (interval == 0) return std::nullopt;
            --interval;
        }
        return LF_IntervalPoint{static_cast<uint64_t>(-1), interval, LF.get_length(interval) - 1};
    }

    LF_IntervalPoint pred_char_cyclic(const LF_IntervalPoint& pos, const uint8_t c) const {
        auto pred = pred_char(pos, c);
        if (pred.has_value()) {
            return pred.value();
        }
        else {
            pred = pred_char(end_bwt_pos(), c);
            if (pred.has_value()) {
                return pred.value();
            }
            else {
                throw std::runtime_error("Character not found in BWT: " + std::string(1, c));
            }
        }
    }
    
    std::optional<LF_IntervalPoint> succ_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        uint64_t interval = pos.interval;
        while (rlbwt[interval] != c) {
            if (interval == LF.data.size() - 1) return std::nullopt;
            ++interval;
        }
        return LF_IntervalPoint{static_cast<uint64_t>(-1), interval, 0};
    }

    LF_IntervalPoint succ_char_cyclic(const LF_IntervalPoint& pos, const uint8_t c) const {
        auto succ = succ_char(pos, c);
        if (succ.has_value()) {
            return succ.value();
        }
        else {
            succ = succ_char(start_bwt_pos(), c);
            if (succ.has_value()) {
                return succ.value();
            }
            else {
                throw std::runtime_error("Character not found in BWT: " + std::string(1, c));
            }
        }
    }

    // Returns vector (row[0..m-1]) where:
    // row[i] = IntervalPoint of suffix which matches lexicographically smallest suffix with maximum match to P[i..m-1] 
    std::vector<LF_IntervalPoint> pred_row(const char* pattern, const uint64_t m) const {
        std::vector<LF_IntervalPoint> pred_row(m);
        LF_IntervalPoint rlbwt_pos = end_bwt_pos();

        for (uint64_t i = 0; i < m; ++i) {
            uint8_t c = charToBits(pattern[m - i - 1]);
            if (c != rlbwt[rlbwt_pos.interval]) {
                rlbwt_pos = pred_char_cyclic(rlbwt_pos, c);
            }
            rlbwt_pos = LF.map(rlbwt_pos);
            pred_row[m - i - 1] = rlbwt_pos;
        }
        return pred_row;
    }

    void serialize() {

    }

    void load() {

    }

// std::pair<std::vector<uint64_t>, std::vector<uint64_t>> _query_ms(const char* query, const uint64_t m) {
//     // std::vector<uint64_t> ms_len;
//     // std::vector<uint64_t> ms_pos;

//     // position, interval, offset
//     MoveStructureTable::IntervalPoint rlbwt_pos = {static_cast<uint64_t>(-1), rlbwt.size() - 1, LF.data.get<2>(rlbwt.size() - 1) - 1};
//     MoveStructureStartTable::IntervalPoint phi_pos = {0, intAtTop[0], 0};
//     phi_pos = Phi.map(phi_pos);

//     uint64_t current_ms_len = 0;
//     for (uint64_t i = 0; i < m; ++i) {
//         const char c = charToBits[query[m - i - 1]];

//         if (c == rlbwt[rlbwt_pos.interval]) {
//             ms_len[m - i - 1] = current_ms_len + 1;
//             ms_pos[m - i - 1] = phi_pos.position;
//         } else {
//             // search for pred
//             // iterate the lcp
//             // forward extension (Psi/FL)

//             // search for succ
//             // iterate the lcp
//             // forward extension (Psi/FL)
//         }
//         // LF step
//         // Update Phi position
//     }
//     return std::make_pair(ms_len, ms_pos);
// }
};
#endif //#ifndef R_SA_LCP_MSINDEX_H
