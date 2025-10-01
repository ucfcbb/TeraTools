#include"util.h"
#include"moveStructure.h"
#include"fm-index.h"

class LCPComputer {
    uint64_t totalLen;

    void ConstructPsi(const rb3_fmi_t* rb3, sdsl::int_vector<> & F, MoveStructure &Psi, uint64_t & numSequences) {
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

        Timer.start("Reading fmd to create Flens");
        {
            *Psi.intLens = sdsl::int_vector<>(runs, 0, lenbits);

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
                        (*Psi.intLens)[alphFRunStarts[c]++] = 1;
                }
                else {
                    (*Psi.intLens)[alphFRunStarts[c]++] = l;
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
            Psi.D_index = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(runs - 1) + 1);
            Psi.D_offset = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(lenRange.max - 1) + 1);

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
                        Psi.D_index[alphFRunStarts[c]] = currentRun;
                        Psi.D_offset[alphFRunStarts[c]] = currentOffset;

                        ++alphFRunStarts[c];

                        ++currentOffset;
                        currentOffset %= (*Psi.intLens)[currentRun];
                        currentRun += (currentOffset == 0);
                    }
                }
                else {
                    Psi.D_index[alphFRunStarts[c]] = currentRun;
                    Psi.D_offset[alphFRunStarts[c]] = currentOffset;

                    ++alphFRunStarts[c];

                    currentOffset += l;
                    while (currentOffset && currentOffset >= (*Psi.intLens)[currentRun])
                        currentOffset -= (*Psi.intLens)[currentRun++];
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
    }

    //O(n) time, not needed for LCP computation
    //LCP computation is O(n) anyways and we can repair while we do it
    //although we don't use psi after so may as well not repair it
    //this function fixes the psi mappings for the endmarker runs (of F)
    void RepairPsi(const sdsl::int_vector<>& F, MoveStructure& Psi) {
        Timer.start("Repairing Psi of endmarkers in F");
        uint64_t numSequences = 0;

        while (numSequences < F.size() && F[numSequences] == 0)
            ++numSequences;

        //#pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t seq = 0; seq < numSequences; ++seq) {
            MoveStructure::IntervalPoint start = {static_cast<uint64_t>(-1), seq, 0}, curr;
            start = Psi.map(start);
            curr = start;

            while (curr.interval >= numSequences)
                curr = Psi.map(curr);

            if (curr.offset) {
                std::cerr << "ERROR: Run of endmarkers in F of length more than 1!" << std::endl;
                exit(1);
            }

            uint64_t seqStartingAtStart = (curr.interval)? curr.interval - 1 : numSequences - 1;
            //#pragma omp critical
            {
                Psi.D_index[seqStartingAtStart] = start.interval;
                Psi.D_offset[seqStartingAtStart] = start.offset;
            }
        }
        Timer.stop(); //Repairing Psi of endmarkers in F
    }

    public:
    typedef uint64_t size_type;

    //input: a run length encoding of a multidollar BWT where all dollars are represented by 0
    //All characters between (and including) 0 and max_char are assumed to have more than 0 occurrences
    //in the text. max_char is the maximum character in the text
    LCPComputer(rb3_fmi_t* rb3) {

        MoveStructure Psi;
        std::vector<uint64_t> alphStarts;

        sdsl::int_vector<> F;
        sdsl::int_vector<> Flens;
        Psi.intLens = &Flens;
        uint64_t numSequences;
        ConstructPsi(rb3, F, Psi, numSequences);

        /*
        for (uint64_t i = 0; i < Psi.D_index.size(); ++i) {
            std::cout 
                << F[i] << '\t'
                << (*Psi.intLens)[i] << '\t'
                << Psi.D_index[i] << '\t'
                << Psi.D_offset[i] << '\n';
        }
        */

        RepairPsi(F, Psi);

        MoveStructure::IntervalPoint end{static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
        curr = end;
        //char con[] = "$ACGTN";
        //do {
            //curr = Psi.map(curr);
            //std::cout << con[F[curr.interval]];
        //} while (curr != end);
        //std::cout << std::endl;


        Timer.start("Verifying Psi");
        if (!Psi.permutationLengthN(totalLen)) {
            std::cerr << "ERROR: Psi is not a permutation of length n!" << std::endl;
            exit(1);
        }
        std::cout << "Psi is a permutation of length n\n";
        Timer.stop(); //Verifying Psi

        rb3_fmi_free(rb3);
    }

    static bool validateRB3(const rb3_fmi_t* rb3);
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;

        bytes += sdsl::serialize(totalLen, out, child, "totalLen");

        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(totalLen, in);
    }
};

bool LCPComputer::validateRB3(const rb3_fmi_t* rb3){
    if (!rb3->e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return false;
    }
    return true;
}
