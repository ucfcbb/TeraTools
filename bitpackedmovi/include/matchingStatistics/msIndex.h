#ifndef R_SA_LCP_MSINDEX_H
#define R_SA_LCP_MSINDEX_H
#include <optional>
#include<sdsl/int_vector.hpp>
#include"moveStructure/moveStructure.h"

#define STATS

class MSIndex {
    uint64_t totalLen;
    sdsl::int_vector<> F, rlbwt;

    sdsl::int_vector<> PsiIntAtTop, PsiOffAtTop;
    MoveStructureTable Psi, LF;

    sdsl::int_vector<> intAtTop, intAtBot;
    MoveStructureStartTable Phi, InvPhi;
    sdsl::int_vector<> PLCPsamples, PLCPBelowsamples;

    using LF_IntervalPoint = MoveStructureTable::IntervalPoint;
    using Psi_IntervalPoint = MoveStructureTable::IntervalPoint;
    using Phi_IntervalPoint = MoveStructureStartTable::IntervalPoint;

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
        while (p.interval < numRuns && p.offset >= LF.data.get<2>(p.interval))
            p.offset -= LF.data.get<2>(p.interval++);
        start = p;
        //std::cout << "starting initialization loop" << std::endl;
        for (uint64_t i = numSequences; i < numRuns; ++i) {
            while (p.interval < numRuns && p.offset >= LF.data.get<2>(p.interval))
                p.offset -= LF.data.get<2>(p.interval++);
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

    void generateInvPhiintAtBotPLCPBelow() {
        sdsl::int_vector<> pi, invPi; 
        std::tie(InvPhi, pi, invPi) = Phi.invertAndRetPiInvPi();
        PLCPBelowsamples = sdsl::int_vector<>(PLCPsamples.size(), 0, PLCPsamples.width());
        for (uint64_t i = 0; i < PLCPBelowsamples.size(); ++i)
            PLCPBelowsamples[i] = PLCPsamples[pi[i]];
        pi = sdsl::int_vector<>();
        intAtBot = sdsl::int_vector<>(intAtTop.size(), 0, intAtTop.width());
        for (uint64_t i = 0; i < intAtBot.size(); ++i)
            intAtBot[i] = invPi[intAtTop[(i+1)%intAtTop.size()]];
    }

    void generatePsiAtTop() {
        PsiIntAtTop = sdsl::int_vector<>(LF.num_intervals(), 0, LF.data.a);
        PsiOffAtTop = sdsl::int_vector<>(LF.num_intervals(), 0, LF.data.b);
        uint64_t L_pos = 0;
        uint64_t F_pos = 0;
        uint64_t F_int = 0;
        for (uint64_t L_int = 0; L_int < LF.num_intervals(); ++L_int) {
            PsiIntAtTop[L_int] = F_int;
            PsiOffAtTop[L_int] = L_pos - F_pos;
            
            L_pos += LF.get_length(L_int);
            while (F_int < Psi.num_intervals() && F_pos + Psi.get_length(F_int) <= L_pos) {
                F_pos += Psi.get_length(F_int);
                ++F_int;
            }
        }
    }

    public:
    typedef uint64_t size_type;

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
        if (v >= TIME) { Timer.start("Constructing LF, Psi, rlbwt, F, PsiIntAtTop, and PsiOffAtTop"); }
        //load
        sdsl::load(totalLen, lcpIn);
        sdsl::load(F, lcpIn);
        sdsl::load(Psi, lcpIn);
        test(vPsi, Psi, "Psi");
        sdsl::int_vector<> pi;
        std::tie(LF, pi, std::ignore) = Psi.invertAndRetPiInvPi();
        test(vLF, LF, "LF");
        generateRLBWTfromLFPsiandF();
        generatePsiAtTop();
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
        sdsl::serialize(PsiIntAtTop, MSIout);
        sdsl::serialize(PsiOffAtTop, MSIout);
        //clear
        rlbwt = F = sdsl::int_vector<>();
        LF = Psi = MoveStructureTable();
        PsiIntAtTop = PsiOffAtTop = sdsl::int_vector<>();
        if (v >= TIME) { Timer.stop(); } //Constructing LF, Psi, rlbwt, F, PsiIntAtTop, and PsiOffAtTop

        //intAtTop
        //PLCPsamples
        
        if (v >= TIME) { Timer.start("Constructing intAtTop, intAtBot, Phi, invPhi, PLCPsamples, and PLCPBelowsamples"); }
        //load
        sdsl::load(intAtTop, lcpIn);
        assert(pi.size() == intAtTop.size());
        for (uint64_t i = 0; i < pi.size(); ++i)
            pi[i] = intAtTop[pi[i]];
        intAtTop = pi;
        pi = sdsl::int_vector<>();
        sdsl::load(Phi, lcpIn);
        sdsl::load(PLCPsamples,lcpIn);
        generateInvPhiintAtBotPLCPBelow();
        test2(vPhi, Phi, "Phi");
        test2(vInvPhi, InvPhi, "InvPhi");
        //write
        sdsl::serialize(intAtTop, MSIout);
        sdsl::serialize(intAtBot, MSIout);
        sdsl::serialize(Phi, MSIout);
        sdsl::serialize(InvPhi, MSIout);
        sdsl::serialize(PLCPsamples, MSIout);
        sdsl::serialize(PLCPBelowsamples, MSIout);
        //clear
        intAtTop = sdsl::int_vector<>();
        intAtBot = sdsl::int_vector<>();
        Phi = InvPhi = MoveStructureStartTable();
        PLCPsamples = sdsl::int_vector<>();
        PLCPBelowsamples = sdsl::int_vector<>();
        if (v >= TIME) { Timer.stop(); } //Constructing intAtTop, intAtBot, Phi, invPhi, PLCPsamples, and PLCPBelowsamples
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

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;

        bytes += sdsl::serialize(totalLen, out, child, "totalLen");
        bytes += sdsl::serialize(F, out, child, "F");
        bytes += sdsl::serialize(rlbwt, out, child, "rlbwt");
        bytes += sdsl::serialize(Psi, out, child, "Psi");
        bytes += sdsl::serialize(LF, out, child, "LF");
        bytes += sdsl::serialize(PsiIntAtTop, out, child, "PsiIntAtTop");
        bytes += sdsl::serialize(PsiOffAtTop, out, child, "PsiOffAtTop");
        bytes += sdsl::serialize(intAtTop, out, child, "intAtTop");
        bytes += sdsl::serialize(intAtBot, out, child, "intAtBot");
        bytes += sdsl::serialize(Phi, out, child, "Phi");
        bytes += sdsl::serialize(InvPhi, out, child, "InvPhi");
        bytes += sdsl::serialize(PLCPsamples, out, child, "PLCPsamples");
        bytes += sdsl::serialize(PLCPBelowsamples, out, child, "PLCPBelowsamples");

        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(totalLen, in);
        sdsl::load(F, in);
        sdsl::load(rlbwt, in);
        sdsl::load(Psi, in);
        sdsl::load(LF, in);
        sdsl::load(PsiIntAtTop, in);
        sdsl::load(PsiOffAtTop, in);
        sdsl::load(intAtTop, in);
        sdsl::load(intAtBot, in);
        sdsl::load(Phi, in);
        sdsl::load(InvPhi, in);
        sdsl::load(PLCPsamples, in);
        sdsl::load(PLCPBelowsamples, in);
    }

    //matching algorithms-----------------------

    //WARNING, THIS IMPLEMENTATION ASSUMES NO RUN SPLITTING, MAY BE BUGGY IF RUN SPLITTING IS PERFORMED
    //NOTE: THIS ALGORITHM MIGHT BE FASTER IN PRACTICE IF WE COMPUTE IT IN THE TEXT ORDER,
    //  I.E. BY PLCP instead of BWT order (faster due to locality of reference)
    void superMaximalRepeats(std::ostream& out, const uint64_t lengthThreshold = 1) {
        //a supermaximal repeat is a substring of the text T[i,i+l) s.t.
        //  a. occ(T[i,i+l)) > 1
        //  b. occ(T[i-1,i+l)) = 1
        //  c. occ(T[i,i+l+1)) = 1
        //T[i,i+l) is a super maximal repeat iff
        //  a. i occurs in SA at the top or bottom of a run
        //  b. max(PLCP[i], PLCP[invphi[i]]) >= max(PLCP[i-1], PLCP[invphi[i-1]])
        //  c. l = max(PLCP[i], PLCP[invphi[i]])
        //
        //This function outputs all supermaximal repeats in the text and all of their occurrences in O(r + occ) time
        //It could also be easily modified to only output the supermaximal repeats 
        //in O(r + c) time where c is the number of repeats outputted
        //
        //If a lengthThreshold is provided, it only outputs supermaximal repeats of length at least the threshold (thresholds 0 and 1 have the same behavior)

        out << "seq\tpos\tlen\tocc\n";

        uint64_t runs = rlbwt.size();
        for (uint64_t run = 0; run < runs; ++run) {
            //check run boundary run (i.e. top of run [run] and bottom of run [run-1 mod runs]
            uint64_t runs = rlbwt.size();
            uint64_t matchingLengthBetweenRuns = PLCPsamples[intAtTop[run]];
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
                auto pred_pos_result = pred_char(rlbwt_pos, c);
                if (pred_pos_result.has_value()) {
                    rlbwt_pos = pred_pos_result.value().first;
                }
                else {
                    auto succ_pos_result = succ_char(rlbwt_pos, c);
                    if (succ_pos_result.has_value()) {
                        rlbwt_pos = succ_pos_result.value().first;
                    }
                    else {
                        throw std::runtime_error("Character not found in BWT: " + std::string(1, c));
                    }
                }
            }
            rlbwt_pos = LF.map(rlbwt_pos);
            pred_row[m - i - 1] = rlbwt_pos;
        }
        return pred_row;
    }

    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_phi(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length) {
            reposition_explicit(c, pos, phi_pos, length, [this](const Phi_IntervalPoint& phi_pos, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t dist, const uint64_t lower_lim, const uint64_t upper_lim) {
                return phi_lce(phi_pos, start, end, dist, lower_lim, upper_lim);
            });
        });
    }
    
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_psi(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length) {
            reposition_explicit(c, pos, phi_pos, length, [this](const Phi_IntervalPoint& phi_pos, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t dist, const uint64_t lower_lim, const uint64_t upper_lim) {
                return psi_lce(phi_pos, start, end, dist, lower_lim, upper_lim);
            });
        });
    }

    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_dual(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length) {
            reposition_explicit(c, pos, phi_pos, length, [this](const Phi_IntervalPoint& phi_pos, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t dist, const uint64_t lower_lim, const uint64_t upper_lim) {
                return dual_lce(phi_pos, start, end, dist, lower_lim, upper_lim);
            });
        });
    }

private:
    // ================================ MS Stats ================================
    #ifdef STATS
    size_t phi_steps;
    size_t psi_steps;
    size_t iteration;
    std::vector<std::pair<uint64_t, uint64_t>> bwt_row;
    std::vector<uint64_t> pred_lces;
    std::vector<uint64_t> succ_lces;    

    void init_ms_stats(const uint64_t m) {
        phi_steps = 0;
        psi_steps = 0;
        iteration = m - 1;
        bwt_row = std::vector<std::pair<uint64_t, uint64_t>>(m);
        pred_lces = std::vector<uint64_t>(m);
        succ_lces = std::vector<uint64_t>(m);
    }

    void print_ms_stats() const {
        std::cout << "  Phi Steps: " << phi_steps << std::endl;
        std::cout << "  Psi Steps: " << psi_steps << std::endl;
        std::cout << "Total Steps: " << phi_steps + psi_steps << std::endl;
        std::cout << "BWT Row: " << std::endl;
        for (auto& [interval, offset] : bwt_row) {
            std::cout << LF.get_start(interval) + offset << " ";
        }
        std::cout << std::endl;
        std::cout << "Pred.LCE: " << std::endl;
        for (auto lce : pred_lces) {
            std::cout << lce << " ";
        }
        std::cout << std::endl;
        std::cout << "Succ.LCE: " << std::endl;
        for (auto lce : succ_lces) {
            std::cout << lce << " ";
        }
        std::cout << std::endl;
    }
    #endif

    // ================================ General helper functions ================================
    // Make static contexpr, but this is probably fine
    static uint8_t charToBits(const char c) {
        switch (c) {
            case '\0': return 0;
            case  '$': return 0;
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
    
    Phi_IntervalPoint rlbwt_to_phi(const uint64_t interval) const {
        return Phi_IntervalPoint{Phi.get_start(intAtTop[interval]), intAtTop[interval], 0};
    }
    
    Psi_IntervalPoint rlbwt_to_psi(const uint64_t interval) const {
        return Psi_IntervalPoint{static_cast<uint64_t>(-1), PsiIntAtTop[interval], PsiOffAtTop[interval]};
    }

    // Might need to "fast forward" to the correct interval
    Psi_IntervalPoint rlbwt_to_psi(LF_IntervalPoint pos) const {
        Psi_IntervalPoint psi_pos = rlbwt_to_psi(pos.interval);
        psi_pos.offset += pos.offset;
        while (psi_pos.offset >= Psi.get_length(psi_pos.interval)) {
            psi_pos.offset -= Psi.get_length(psi_pos.interval);
            ++psi_pos.interval;
        }
        return psi_pos;
    }
    
    // Perform one LF step and update Phi accordingly
    void backward_step(LF_IntervalPoint& rlbwt_pos, Phi_IntervalPoint& phi_pos) {
        rlbwt_pos = LF.map(rlbwt_pos);
        if (phi_pos.offset == 0) {
            if (phi_pos.interval == 0) { phi_pos.interval = Phi.data.size() - 1; }
            else { phi_pos.interval--; }
            phi_pos.offset = Phi.get_length(phi_pos.interval) - 1;
            phi_pos.position = Phi.get_start(phi_pos.interval) + phi_pos.offset;
        }
        else {
            --phi_pos.position;
            --phi_pos.offset;
        }
    }
    
    // Returns LF position and distance to the previous occurrence of character c, if the predecessor exists. Returns passed position if already at the predecessor.
    std::optional<std::pair<LF_IntervalPoint, uint64_t>> pred_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        if (rlbwt[pos.interval] == c) {
            return std::make_pair(pos, 0);
        }
        uint64_t distance = pos.offset + 1;

        if (pos.interval == 0) return std::nullopt;
        uint64_t interval = pos.interval - 1;

        while (rlbwt[interval] != c) {
            distance += LF.get_length(interval);
            if (interval == 0) return std::nullopt;
            --interval;
        }
        return std::make_pair(LF_IntervalPoint{static_cast<uint64_t>(-1), interval, LF.get_length(interval) - 1}, distance);
    }
    
    // Returns LF position and distance to the next occurrence of character c, if the successor exists. Returns passed position if already at the successor.
    std::optional<std::pair<LF_IntervalPoint, uint64_t>> succ_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        if (rlbwt[pos.interval] == c) {
            return std::make_pair(pos, 0);
        }
        uint64_t distance = LF.get_length(pos.interval) - pos.offset;

        if (pos.interval == LF.data.size() - 1) return std::nullopt;
        uint64_t interval = pos.interval + 1;

        while (rlbwt[interval] != c) {
            distance += LF.get_length(interval);
            if (interval == LF.data.size() - 1) return std::nullopt;
            ++interval;
        }
        return std::make_pair(LF_IntervalPoint{static_cast<uint64_t>(-1), interval, 0}, distance);
    }

    // ================================ MS Loop (Main MS Function) ================================
    // Used to define different MS algorithms by passing a different reposition function
    using RepositionFunction = std::function<void(const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length)>;
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_loop(const char* pattern, const uint64_t m, RepositionFunction reposition) {
        #ifdef STATS
        init_ms_stats(m);
        #endif

        std::vector<uint64_t> ms_len(m);
        std::vector<uint64_t> ms_pos(m);

        // initialize to last position in BWT (for both rlbwt and corresponding Phi position)
        LF_IntervalPoint rlbwt_pos = end_bwt_pos();
        Phi_IntervalPoint phi_pos = rlbwt_to_phi(0); // Used to find MS positions
        phi_pos = Phi.map(phi_pos);

        uint64_t length = 0; // MS length
        for (uint64_t i = 0; i < m; ++i) {
            uint8_t c = charToBits(pattern[m - i - 1]);
            if (c != rlbwt[rlbwt_pos.interval]) {
                reposition(c, rlbwt_pos, phi_pos, length);
            }
            ++length;
            backward_step(rlbwt_pos, phi_pos);

            ms_len[m - i - 1] = length;
            ms_pos[m - i - 1] = phi_pos.position;
            #ifdef STATS
            --iteration;
            bwt_row[m - i - 1] = std::make_pair(rlbwt_pos.interval, rlbwt_pos.offset);
            #endif
        }
        
        #ifdef STATS
        print_ms_stats();
        #endif
        return std::make_pair(ms_len, ms_pos);
    }

    // ================================ Reposition Functions ================================
    // Doesn't use thresholds, so we look for both the predecessor and successor
    using LCEFunction = std::function<uint64_t(const Phi_IntervalPoint& phi_position, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t distance, const uint64_t lower_lim, const uint64_t upper_lim)>;
    void reposition_explicit(const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length, LCEFunction lce) {
        auto pred_pos_result = pred_char(pos, c);
        auto succ_pos_result = succ_char(pos, c);
        if (!pred_pos_result.has_value() && !succ_pos_result.has_value())
            throw std::runtime_error("No valid positions found for repositioning, character not found in BWT: " + std::string(1, c));

        LF_IntervalPoint pred_pos;
        Phi_IntervalPoint pred_phi_pos;
        uint64_t pred_lce = 0;
        if (pred_pos_result.has_value()) {
            pred_pos = pred_pos_result.value().first;
            uint64_t pred_distance = pred_pos_result.value().second;
            pred_lce = lce(phi_pos, pos, pred_pos, pred_distance, 0, length);
            // New phi position is a run tail, so map from the next head
            pred_phi_pos = rlbwt_to_phi(pred_pos.interval + 1);
            pred_phi_pos = Phi.map(pred_phi_pos);
            #ifdef STATS
            pred_lces[iteration] = pred_lce;
            #endif
        }

        LF_IntervalPoint succ_pos;
        Phi_IntervalPoint succ_phi_pos;
        uint64_t succ_lce = 0;
        if (succ_pos_result.has_value()) {
            succ_pos = succ_pos_result.value().first;
            uint64_t succ_distance = succ_pos_result.value().second;
            succ_lce = lce(phi_pos, pos, succ_pos, succ_distance, pred_lce, length);
            // New phi position is a run head
            succ_phi_pos = rlbwt_to_phi(succ_pos.interval);
            #ifdef STATS
            succ_lces[iteration] = succ_lce;
            #endif
        }

        // Need to check pred_pos_result in case both LCEs are 0 (could mean either no predecessor or that the LCE really was 0)
        if (pred_lce >= succ_lce && pred_pos_result.has_value()) {
            pos = pred_pos;
            phi_pos = pred_phi_pos;
            length = std::min(pred_lce, length);
        }
        else {
            pos = succ_pos;
            phi_pos = succ_phi_pos;
            length = std::min(succ_lce, length);
        }
    }

    // void reposition_threshold(const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length, LCEFunction lce) {
    // void reposition_pred()

    // ================================ LCE Functions ================================
    /**
    * @brief LCE between start and end position, where end is a run head or tail, by enumerating the LCP interval using Phi and PLCP samples
    * 
    * @param phi_position Current Phi position
    * @param start Start position (current MS position)
    * @param end End position (predecessor or successor position)
    * @param distance Distance to the end position (distance to the predecessor or successor)
    * @param lower_lim If the minimum LCE goes beneath this value, stop and return 0
    * @return uint64_t The LCE between the start and end position
    */
    uint64_t phi_lce(const Phi_IntervalPoint& phi_position, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t distance, const uint64_t lower_lim = 0) {
        if (start == end) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_down = start < end;

        uint64_t lce = std::numeric_limits<uint64_t>::max();
        Phi_IntervalPoint phi_extension_position = phi_position;
        // Start extending from the end position when extending down to run head
        if (extend_down) {
            phi_extension_position = rlbwt_to_phi(end.interval);
        }

        for (size_t i = 0; i < distance; ++i) {
            uint64_t lcp = PLCPsamples[phi_extension_position.interval] - phi_extension_position.offset;
            lce = std::min(lce, lcp);
            if (lce < lower_lim) { lce = 0; break; }
            if (i == distance - 1) { ++i; break; } // Increment i and break to avoid redundant mapping step during last iteration
            phi_extension_position = Phi.map(phi_extension_position);

            #ifdef STATS
            ++phi_steps;
            #endif
        }

        return lce;
    }
    // Wrapper for reposition methods
    uint64_t phi_lce(const Phi_IntervalPoint& phi_position, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t distance, const uint64_t lower_lim, const uint64_t /*upper_lim*/) {
        return phi_lce(phi_position, start, end, distance, lower_lim);
    }

    /**
    * @brief LCE between start and end position, where end is a run head or tail, by forward character comparisons using Psi
    * 
    * @param start Start position (current MS position)
    * @param end End position (predecessor or successor position)
    * @param upper_lim If the LCE reaches this value, stop and return the value
    * @return uint64_t The LCE between the start and end position
    */
    uint64_t psi_lce(const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t upper_lim = 0) {
        if (start == end) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_up = start > end;

        Psi_IntervalPoint start_psi_pos = rlbwt_to_psi(start);
        Psi_IntervalPoint end_psi_pos;
        // Extend to run tail
        if (extend_up) {
            end_psi_pos = rlbwt_to_psi(end.interval + 1);
            --end_psi_pos.interval;
            end_psi_pos.offset = Psi.get_length(end_psi_pos.interval) - 1;
        }
        // Extend to run head
        else {
            end_psi_pos = rlbwt_to_psi(end.interval);
        }

        // Move past BWT character (want to compare suffixes)
        start_psi_pos = Psi.map(start_psi_pos);
        end_psi_pos = Psi.map(end_psi_pos);

        uint64_t lce = 0;
        // Use true loops to avoid redundant mapping steps during last iteration
        while (true) {
            if (F[start_psi_pos.interval] != F[end_psi_pos.interval]) { break; }
            ++lce;
            if (lce >= upper_lim) { break; }
            start_psi_pos = Psi.map(start_psi_pos);
            end_psi_pos = Psi.map(end_psi_pos);

            #ifdef STATS
            psi_steps += 2;
            #endif
        }
        return lce;
    }
    uint64_t psi_lce(const Phi_IntervalPoint& /*_phi_position*/, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t /*_distance*/, const uint64_t /*lower_lim*/, const uint64_t upper_lim) {
        return psi_lce(start, end, upper_lim);
    }
    
    // Number of consecutive steps to take before switching to the other extension method, default is 10
    template<size_t consecutive_steps = 10>
    uint64_t dual_lce(const Phi_IntervalPoint& phi_position, const LF_IntervalPoint& start, const LF_IntervalPoint& end, const uint64_t distance, const uint64_t lower_lim = 0, const uint64_t upper_lim = std::numeric_limits<uint64_t>::max()) {
        if (start == end) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_up = start > end;
        
        // Phi LCE
        uint64_t phi_lce = std::numeric_limits<uint64_t>::max();
        size_t i = 0; // counter for the phi extension
        Phi_IntervalPoint phi_extension_position = phi_position;
        // Start extending from the end position when extending down
        if (!extend_up) {
            phi_extension_position = rlbwt_to_phi(end.interval);
        }
        
        // Psi LCE
        uint64_t psi_lce = 0;
        Psi_IntervalPoint start_psi_pos = rlbwt_to_psi(start);
        Psi_IntervalPoint end_psi_pos;
        // Extend to run tail
        if (extend_up) {
            // TODO no bounds check?
            end_psi_pos = rlbwt_to_psi(end.interval + 1);
            --end_psi_pos.interval;
            end_psi_pos.offset = Psi.get_length(end_psi_pos.interval) - 1;
        }
        // Extend to run head
        else {
            end_psi_pos = rlbwt_to_psi(end.interval);
        }
        // Move past BWT character (want to compare suffixes)
        start_psi_pos = Psi.map(start_psi_pos);
        end_psi_pos = Psi.map(end_psi_pos);

        size_t counter = 0; // Counter for mapping steps taken
        auto phi_turn = [&]() { return (counter < consecutive_steps); };
        auto phi_condition = [&]() { return (i < distance) && (phi_lce >= lower_lim); };
        auto psi_condition = [&]() { return (F[start_psi_pos.interval] == F[end_psi_pos.interval]) && (psi_lce < upper_lim); };
        // In this case, we know the psi will finish first, so we can skip the phi computation
        // Multiply by 2 because we take two mapping steps per iteration
        auto skip_phi = [&]() { return (2*upper_lim < distance); };
        while (true) {
            if (!skip_phi() && phi_turn()) {
                if (!phi_condition()) { break; } // If distance is reached, break and return the current LCE
                uint64_t lcp = PLCPsamples[phi_extension_position.interval] - phi_extension_position.offset;
                phi_lce = std::min(phi_lce, lcp);
                if (!phi_condition()) { phi_lce = 0; break; } // If the minimum LCE goes beneath the lower limit, stop and return 0
                if (i == distance - 1) { ++i; break; } // Increment i and break to avoid redundant mapping step during last iteration
                phi_extension_position = Phi.map(phi_extension_position);
                ++i;

                ++counter;
                #ifdef STATS
                ++phi_steps;
                #endif
            }
            else {
                if (!psi_condition()) { break; }
                ++psi_lce;
                if (!psi_condition()) { break; }
                start_psi_pos = Psi.map(start_psi_pos);
                end_psi_pos = Psi.map(end_psi_pos);

                counter += 2;
                #ifdef STATS
                psi_steps += 2;
                #endif
            }
            
            if (counter >= 2*consecutive_steps) {
                counter = 0;
            }
        }
        
        if (!phi_condition() && !psi_condition()) {
            return std::min(phi_lce, psi_lce);
        }
        else if (!phi_condition()) {
            return phi_lce;
        }
        else {
            return psi_lce;
        }
    }
};
#endif //#ifndef R_SA_LCP_MSINDEX_H
