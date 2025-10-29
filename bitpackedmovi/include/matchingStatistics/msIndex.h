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
            pi[i] = (intAtTop[pi[i]] + 1)%pi.size();
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

    #ifdef STATS
    void print_ms_stats() const {
        std::cout << "\t  Phi Steps: " << phi_steps << std::endl;
        std::cout << "\t  Psi Steps: " << psi_steps << std::endl;
        std::cout << "\tTotal Steps: " << phi_steps + psi_steps << std::endl;
    }
    #endif

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

    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_phi(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](MSState& state, const uint8_t c) {
            reposition_explicit(state, c, [this](const MSState& state, const PosDist& end, const uint64_t lower_lim) {
                return phi_lce(state, end, lower_lim);
            });
        });
    }
    
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_psi(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](MSState& state, const uint8_t c) {
            reposition_explicit(state, c, [this](const MSState& state, const PosDist& end, const uint64_t lower_lim) {
                return psi_lce(state, end, lower_lim);
            });
        });
    }

    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_dual(const char* pattern, const uint64_t m) {
        return ms_loop(pattern, m, [this](MSState& state, const uint8_t c) {
            reposition_explicit(state, c, [this](const MSState& state, const PosDist& end, const uint64_t lower_lim) {
                return dual_lce(state, end, lower_lim);
            });
        });
    }

private:
    // ================================ MS Stats ================================
    #ifdef STATS
    size_t phi_steps;
    size_t psi_steps;

    void init_ms_stats(const uint64_t m) {
        phi_steps = 0;
        psi_steps = 0;
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

    uint64_t get_PLCP(const Phi_IntervalPoint& phi_pos) const {
        return PLCPsamples[phi_pos.interval] - phi_pos.offset;
    }

    // ================================ BWT Position Conversion Functions ================================
    LF_IntervalPoint start_bwt_pos() const {
        return LF_IntervalPoint(0, 0);
    }

    LF_IntervalPoint end_bwt_pos() const {
        return LF_IntervalPoint(LF.num_intervals() - 1, LF.get_length(LF.num_intervals() - 1) - 1);
    }
    
    Phi_IntervalPoint rlbwt_head_to_phi(const uint64_t interval) const {
        return Phi.get_interval_point(intAtTop[interval]);
    }
    Phi_IntervalPoint rlbwt_head_to_phi(LF_IntervalPoint pos) const {
        assert(pos.offset == 0);
        return rlbwt_head_to_phi(pos.interval);
    }
    
    // Should only be called if pos is a run tail --> is quicker than using the general method above
    Phi_IntervalPoint rlbwt_tail_to_phi(LF_IntervalPoint pos) const {
        assert(pos.offset == LF.get_length(pos) - 1);
        uint64_t next_interval = pos.interval + 1;
        if (pos.interval == LF.num_intervals() - 1) { next_interval = 0; }
        Phi_IntervalPoint phi_pos = rlbwt_head_to_phi(next_interval);
        return Phi.map(phi_pos);
    }

    Psi_IntervalPoint rlbwt_head_to_psi(const uint64_t interval) const {
        return Psi_IntervalPoint(PsiIntAtTop[interval], PsiOffAtTop[interval]);
    }
    Psi_IntervalPoint rlbwt_head_to_psi(LF_IntervalPoint pos) const {
        assert(pos.offset == 0);
        return rlbwt_head_to_psi(pos.interval);
    }

    // Might need to "fast forward" to the correct interval
    Psi_IntervalPoint rlbwt_to_psi(LF_IntervalPoint pos) const {
        Psi_IntervalPoint psi_pos = rlbwt_head_to_psi(pos.interval);
        psi_pos.offset += pos.offset;
        while (psi_pos.offset >= Psi.get_length(psi_pos)) {
            psi_pos.offset -= Psi.get_length(psi_pos);
            ++psi_pos.interval;
        }
        return psi_pos;
    }

    // Should only be called if pos is a run tail --> is quicker than using the general method above
    // Doesn't work if called on the last interval
    Psi_IntervalPoint rlbwt_tail_to_psi(LF_IntervalPoint pos) const {
        assert(pos.offset == LF.get_length(pos) - 1);
        assert(pos.interval < LF.num_intervals() - 1);
        Psi_IntervalPoint psi_pos = rlbwt_head_to_psi(pos.interval + 1);
        if (psi_pos.offset == 0) { 
            --psi_pos.interval; 
            psi_pos.offset = Psi.get_length(psi_pos) - 1;
        }
        else {
            --psi_pos.offset;
        }
        return psi_pos;
    }
    
    // Perform one LF step and update Phi accordingly
    void backward_step(LF_IntervalPoint& rlbwt_pos, Phi_IntervalPoint& phi_pos) {
        rlbwt_pos = LF.map(rlbwt_pos);
        if (phi_pos.offset == 0) {
            if (phi_pos.interval == 0) { phi_pos.interval = Phi.data.size() - 1; }
            else { phi_pos.interval--; }
            phi_pos.offset = Phi.get_length(phi_pos) - 1;
            phi_pos.position = Phi.get_start(phi_pos) + phi_pos.offset;
        }
        else {
            --phi_pos.position;
            --phi_pos.offset;
        }
    }
    
    // ================================ Predecessor and Successor Functions ================================
    struct PosDist {
        LF_IntervalPoint pos;
        uint64_t dist;
    };

    // Returns LF position and distance to the previous occurrence of character c, if the predecessor exists. Returns passed position if already at the predecessor.
    std::optional<PosDist> pred_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        if (rlbwt[pos.interval] == c) {
            return PosDist{pos, 0};
        }
        uint64_t distance = pos.offset + 1;

        if (pos.interval == 0) return std::nullopt;
        uint64_t interval = pos.interval - 1;

        while (rlbwt[interval] != c) {
            distance += LF.get_length(interval);
            if (interval == 0) return std::nullopt;
            --interval;
        }
        return PosDist{LF_IntervalPoint(interval, LF.get_length(interval) - 1), distance};
    }
    
    // Returns LF position and distance to the next occurrence of character c, if the successor exists. Returns passed position if already at the successor.
    std::optional<PosDist> succ_char(const LF_IntervalPoint& pos, const uint8_t c) const {
        if (rlbwt[pos.interval] == c) {
            return PosDist{pos, 0};
        }
        uint64_t distance = LF.get_length(pos.interval) - pos.offset;

        if (pos.interval == LF.data.size() - 1) return std::nullopt;
        uint64_t interval = pos.interval + 1;

        while (rlbwt[interval] != c) {
            distance += LF.get_length(interval);
            if (interval == LF.data.size() - 1) return std::nullopt;
            ++interval;
        }
        return PosDist{LF_IntervalPoint(interval, 0), distance};
    }

    // ================================ MS Loop (Main MS Function) ================================
    struct MSState {
        const char* pattern;
        const uint64_t m;
        uint64_t i; // position on pattern

        LF_IntervalPoint rlbwt_pos;
        Phi_IntervalPoint phi_pos;
        uint64_t length;

        MSState(const char* pattern, const uint64_t m, const LF_IntervalPoint& rlbwt_pos, const Phi_IntervalPoint& phi_pos) 
        : pattern(pattern), m(m), i(0), rlbwt_pos(rlbwt_pos), phi_pos(phi_pos), length(0) {}
    };

    // Used to define different MS algorithms by passing a different reposition function
    using RepositionFunction = std::function<void(MSState& state, const uint8_t c)>;
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> ms_loop(const char* pattern, const uint64_t m, RepositionFunction reposition) {
        #ifdef STATS
        init_ms_stats(m);
        #endif

        // Initial state is the end of the BWT, end of pattern, length of 0
        MSState state(pattern, m, end_bwt_pos(), rlbwt_tail_to_phi(end_bwt_pos()));

        std::vector<uint64_t> ms_len(m);
        std::vector<uint64_t> ms_pos(m);

        for (state.i = 0; state.i < m; ++state.i) {
            uint8_t c = charToBits(state.pattern[m - state.i - 1]);
            if (c != rlbwt[state.rlbwt_pos.interval]) {
                reposition(state, c);
            }
            ++state.length;
            backward_step(state.rlbwt_pos, state.phi_pos);

            ms_len[state.m - state.i - 1] = state.length;
            ms_pos[state.m - state.i - 1] = state.phi_pos.position;
        }
        return std::make_pair(ms_len, ms_pos);
    }

    // ================================ Reposition Functions ================================
    // Doesn't use thresholds, so we look for both the predecessor and successor
    using LCEFunction = std::function<uint64_t(const MSState& state, const PosDist& end, const uint64_t lower_lim)>;
    void reposition_explicit(MSState& state, const uint8_t c, LCEFunction lce) {
        auto pred_pos_result = pred_char(state.rlbwt_pos, c);
        auto succ_pos_result = succ_char(state.rlbwt_pos, c);
        if (!pred_pos_result.has_value() && !succ_pos_result.has_value())
            throw std::runtime_error("No valid positions found for repositioning, character not found in BWT: " + std::string(1, c));

        LF_IntervalPoint pred_pos;
        Phi_IntervalPoint pred_phi_pos;
        uint64_t pred_lce = 0;
        if (pred_pos_result.has_value()) {
            pred_lce = lce(state, pred_pos_result.value(), 0);
            // New phi position is a run tail
            pred_pos = pred_pos_result.value().pos;
            pred_phi_pos = rlbwt_tail_to_phi(pred_pos);
        }

        LF_IntervalPoint succ_pos;
        Phi_IntervalPoint succ_phi_pos;
        uint64_t succ_lce = 0;
        if (succ_pos_result.has_value()) {
            succ_lce = lce(state, succ_pos_result.value(), pred_lce);
            // New phi position is a run head
            succ_pos = succ_pos_result.value().pos;
            succ_phi_pos = rlbwt_head_to_phi(succ_pos);
        }

        // Need to check pred_pos_result in case both LCEs are 0 (could mean either no predecessor or that the LCE really was 0)
        if (pred_lce >= succ_lce && pred_pos_result.has_value()) {
            state.rlbwt_pos = pred_pos;
            state.phi_pos = pred_phi_pos;
            state.length = std::min(pred_lce, state.length);
        }
        else {
            state.rlbwt_pos = succ_pos;
            state.phi_pos = succ_phi_pos;
            state.length = std::min(succ_lce, state.length);
        }
    }

    // void reposition_threshold(const uint8_t c, LF_IntervalPoint& pos, Phi_IntervalPoint& phi_pos, uint64_t& length, LCEFunction lce) {

    // ================================ LCE Functions ================================
    /**
    * @brief LCE between start and end position, where end is a run head or tail, by enumerating the LCP interval using Phi and PLCP samples
    * 
    * @param state Current MS state
    * @param end End position (predecessor or successor position)
    * @param lower_lim If the minimum LCE goes beneath this value, stop and return 0
    * @return uint64_t The LCE between the start and end position (or 0 if it goes beneath the lower limit)
    */
    uint64_t phi_lce(const MSState& state, const PosDist& end, const uint64_t lower_lim = 0) {
        if (state.rlbwt_pos == end.pos) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_down = state.rlbwt_pos < end.pos;

        uint64_t lce = std::numeric_limits<uint64_t>::max();
        Phi_IntervalPoint phi_extension_position = state.phi_pos;
        // Start extending from the end position when extending down to run head
        if (extend_down) {
            phi_extension_position = rlbwt_head_to_phi(end.pos);
        }

        for (size_t i = 0; i < end.dist; ++i) {
            uint64_t lcp = get_PLCP(phi_extension_position);
            lce = std::min(lce, lcp);
            if (lce < lower_lim) { lce = 0; break; }
            // Break to avoid redundant mapping step during last iteration
            if (i == end.dist - 1) { break; }
            phi_extension_position = Phi.map(phi_extension_position);

            #ifdef STATS
            ++phi_steps;
            #endif
        }

        return lce;
    }

    /**
    * @brief LCE between start and end position, where end is a run head or tail, by forward character comparisons using Psi
    * 
    * @param state Current MS state
    * @param end End position (predecessor or successor position)
    * @return uint64_t The LCE between the start and end position (or the current matched length if it is reached)
    */
    uint64_t psi_lce(const MSState& state, const PosDist& end) {
        if (state.rlbwt_pos == end.pos) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_up = state.rlbwt_pos > end.pos;

        size_t i = state.m - state.i; // position on pattern for psi extension
        Psi_IntervalPoint end_psi_pos;
        // Extend to run tail
        if (extend_up) {
            end_psi_pos = rlbwt_tail_to_psi(end.pos);
        }
        // Extend to run head
        else {
            end_psi_pos = rlbwt_head_to_psi(end.pos);
        }

        uint64_t lce = 0;
        // Use true loops to avoid redundant mapping steps during last iteration
        while (true) {
            if (i >= state.m 
                || charToBits(state.pattern[i]) != F[end_psi_pos.interval]) { break; }
            ++lce;
            if (lce >= state.length) { break; }
            end_psi_pos = Psi.map(end_psi_pos);
            
            ++i;
            #ifdef STATS
            ++psi_steps;
            #endif
        }
        return lce;
    }
    // Wrapper for psi_lce method
    uint64_t psi_lce(const MSState& state, const PosDist& end, const uint64_t /*lower_lim*/) {
        return psi_lce(state, end);
    }

    /**
    * @brief LCE between start and end position, where end is a run head or tail, by using both Phi and Psi, with early stopping conditions
    * 
    * @tparam consecutive_steps Number of consecutive steps to take before switching to the other extension method, best if multiple of 2 since Psi uses 2 mapping steps per iteration
    * @param phi_position Current Phi position
    * @param start Start position (current MS position)
    * @param end End position (predecessor or successor position)
    * @param distance Distance to the end position (distance to the predecessor or successor)
    * @param lower_lim If the minimum LCE goes beneath this value, stop and return 0
    * @param upper_lim If the LCE reaches this value, stop and return the value
    * @return uint64_t The LCE between the start and end position
    */
    template<size_t consecutive_steps = 10>
    uint64_t dual_lce(const MSState& state, const PosDist& end, const uint64_t lower_lim = 0) {
        if (state.rlbwt_pos == end.pos) { throw std::runtime_error("Calling LCE on the same position!"); }
        bool extend_up = state.rlbwt_pos > end.pos;
        
        // Phi LCE
        uint64_t phi_lce = std::numeric_limits<uint64_t>::max();
        size_t phi_i = 0; // counter for the phi extension
        Phi_IntervalPoint phi_extension_position = state.phi_pos;
        // Start extending from the end position when extending down to run head
        if (!extend_up) {
            phi_extension_position = rlbwt_head_to_phi(end.pos);
        }
        
        // Psi LCE
        uint64_t psi_lce = 0;
        size_t psi_i = state.m - state.i; // position on pattern for Psi extension
        Psi_IntervalPoint end_psi_pos;
        // Extend to run tail
        if (extend_up) {
            end_psi_pos = rlbwt_tail_to_psi(end.pos);
        }
        // Extend to run head
        else {
            end_psi_pos = rlbwt_head_to_psi(end.pos);
        }

        size_t counter = 0; // Counter for mapping steps taken
        auto phi_turn = [&]() { return (counter < consecutive_steps); };
        // Determines regular condition and whether early stopping condition not met for both phi and psi
        auto phi_condition = [&]() { return (phi_i < end.dist) && !(phi_lce < lower_lim); };
        auto psi_condition = [&]() { 
            return (psi_i < state.m) && (charToBits(state.pattern[psi_i]) == F[end_psi_pos.interval]) && !(psi_lce >= state.length); 
        };
        // In this case, we know the psi will finish first with its early stopping condition met, so we can skip the phi computation
        bool skip_phi = (state.length < end.dist);
        while (true) {
            if (!skip_phi && phi_turn()) {
                if (!phi_condition()) { break; } // If distance is reached, break and return the current LCE
                uint64_t lcp = get_PLCP(phi_extension_position);
                phi_lce = std::min(phi_lce, lcp);
                // Will only change from above if early stopping condition is met
                if (!phi_condition()) { phi_lce = 0; break; }
                // Increment i and break to avoid redundant mapping step during last iteration
                if (phi_i == end.dist - 1) { ++phi_i; break; }
                phi_extension_position = Phi.map(phi_extension_position);
                
                ++phi_i;
                ++counter;
                #ifdef STATS
                ++phi_steps;
                #endif
            }
            else {
                if (!psi_condition()) { break; }
                ++psi_lce;
                // Will only change from above if early stopping condition is met
                if (!psi_condition()) { break; }
                end_psi_pos = Psi.map(end_psi_pos);

                ++psi_i;
                ++counter;
                #ifdef STATS
                ++psi_steps;
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

static constexpr const char* ms_index_extension = ".ms_index";

#endif //#ifndef R_SA_LCP_MSINDEX_H
