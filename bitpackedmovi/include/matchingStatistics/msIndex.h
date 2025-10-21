#ifndef R_SA_LCP_MSINDEX_H
#define R_SA_LCP_MSINDEX_H
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

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;

        bytes += sdsl::serialize(F, out, child, "F");
        bytes += sdsl::serialize(rlbwt, out, child, "rlbwt");
        bytes += sdsl::serialize(Psi, out, child, "Psi");
        bytes += sdsl::serialize(LF, out, child, "LF");
        bytes += sdsl::serialize(intAtTop, out, child, "intAtTop");
        bytes += sdsl::serialize(Phi, out, child, "Phi");
        bytes += sdsl::serialize(InvPhi, out, child, "InvPhi");
        bytes += sdsl::serialize(PLCPsamples, out, child, "PLCPsamples");

        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(F, in);
        sdsl::load(rlbwt, in);
        sdsl::load(Psi, in);
        sdsl::load(LF, in);
        sdsl::load(intAtTop, in);
        sdsl::load(Phi, in);
        sdsl::load(InvPhi, in);
        sdsl::load(PLCPsamples, in);
    }
};
#endif //#ifndef R_SA_LCP_MSINDEX_H
