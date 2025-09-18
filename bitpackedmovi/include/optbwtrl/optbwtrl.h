#include"util.h"
#include<sdsl/int_vector.hpp>
#include<vector>
#include"fm-index.h"

class OptBWTRL {
    public:
        typedef uint64_t size_type;
    private:
    uint64_t totalLen = 0;
    sdsl::int_vector<> rlbwt, runlens;
    sdsl::int_vector<> SATopRunInt;

    struct MoveStructure {
        typedef OptBWTRL::size_type size_type;
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
        IntervalPoint map(const IntervalPoint& intPoint) {
            IntervalPoint res;
            res.position = static_cast<uint64_t>(-1);
            res.interval = D_index[intPoint.interval];
            res.offset = D_offset[intPoint.interval] + intPoint.offset;
            while ((*intLens)[res.interval] <= res.offset)
                res.offset -= (*intLens)[res.interval++];
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
        void AdvanceIntervalPoint_unsafe(IntervalPoint& intPoint, uint64_t distance) {
            uint64_t remaining = intPoint.offset + distance;
            //if intPoint.interval == intLens->size(), this will result in a BUG
            //TODO: FIX?
            while (remaining && remaining >= (*intLens)[intPoint.interval])
                remaining -= (*intLens)[intPoint.interval++];
            intPoint.offset = remaining;
            intPoint.position += distance;
        }

        void MakeHistogram(std::vector<uint64_t>& hist) {
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
    } LF;

    struct InvertibleMoveStructure {
        typedef OptBWTRL::size_type size_type;
        //suffix array samples at
        sdsl::int_vector<> SeqAt;
        sdsl::int_vector<> PosAt;

        //characters between this sample and the next (including this sample)
        //unnecessary, but makes coding cleaner
        sdsl::int_vector<> IntLen;

        MoveStructure phi;
        MoveStructure invPhi;
        //sdsl::int_vector<> AboveToInterval;
        //sdsl::int_vector<> AboveToOffset;
        //sdsl::int_vector<> BelowToInterval;
        //sdsl::int_vector<> BelowToOffset;

        sdsl::int_vector<> AboveLCP;

        /*
        IntervalPoint mapPhi(const IntervalPoint& intPoint) {
            IntervalPoint res;
            res.position = (uint64_t) -1;
            res.interval = AboveToInterval[intPoint.interval];
            res.offset = AboveToOffset[intPoint.interval] + intPoint.offset;
            while (IntLen[res.interval] <= res.offset)
                res.offset -= IntLen[res.interval++];
            return res;
        }

        IntervalPoint mapInvPhi(const IntervalPoint& intPoint) {
            IntervalPoint res;
            res.position = (uint64_t) -1;
            res.interval = BelowToInterval[intPoint.interval];
            res.offset = BelowToOffset[intPoint.interval] + intPoint.offset;
            while(IntLen[res.interval] <= res.offset)
                res.offset -= IntLen[res.interval++];
            return res;
        }
        */
        size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

            size_type bytes = 0;

            bytes += sdsl::serialize(SeqAt, out, child, "SeqAt");
            bytes += sdsl::serialize(PosAt, out, child, "PosAt");
            bytes += sdsl::serialize(IntLen, out, child, "IntLen");
            bytes += sdsl::serialize(phi, out, child, "phi");
            bytes += sdsl::serialize(invPhi, out, child, "invPhi");
            bytes += sdsl::serialize(AboveLCP, out, child, "AboveLCP");

            sdsl::structure_tree::add_size(child, bytes);
            
            return bytes;
        }

        void load(std::istream& in) {
            sdsl::load(SeqAt, in);
            sdsl::load(PosAt, in);
            sdsl::load(IntLen, in);
            sdsl::load(phi, in);
            sdsl::load(invPhi, in);
            sdsl::load(AboveLCP, in);
        }

        uint64_t LCP(const MoveStructure::IntervalPoint plPoint) {
//            std::cout << "InvMovStr::LCP( " << plPoint.interval 
//                << ", " << plPoint.offset << ") called."
//                << " AboveLCP[] = " << AboveLCP[plPoint.interval] << std::endl;
            return AboveLCP[plPoint.interval] - plPoint.offset;
        }
    } PL;

    void RLBWTconstruction(const rb3_fmi_t* rb3, uRange & alphRange, std::vector<uint64_t>& alphCounts) {
        Timer.start("Constructing RLBWT from FMD");
        //'repaired' values. These are the values of our constructed index.
        //They may differ from ropebwt3's values because we split endmarkers into separate runs
        uint64_t runs = 0, alphbits, lenbits;
        totalLen = 0;
        uRange lenRange;
        std::vector<uint64_t> alphRuns;

        uRange RB3_lenRange;
        Timer.start("Reading fmd for parameters");
        {
            //original ropebwt3 values
            uint64_t RB3_runs = 0, RB3_lenbits;



            rlditr_t itr1;
            rld_itr_init(rb3->e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;


            if ((l = rld_dec(rb3->e, &itr1, &c, 0)) > 0) {
                alphRange = {static_cast<uint64_t>(c),static_cast<uint64_t>(c)};

                RB3_lenRange = {static_cast<uint64_t>(l),static_cast<uint64_t>(l)};
                lenRange = (c == 0)? uRange{static_cast<uint64_t>(1), static_cast<uint64_t>(1)} : RB3_lenRange;

                runs += (c == 0)? static_cast<uint64_t>(l) : 1;
                ++RB3_runs;

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

                totalLen += static_cast<uint64_t>(l);
            }

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
        Timer.start("Reading fmd into sdsl");
        {
            alphCounts.resize(alphRange.max+1, 0);
            alphRuns.resize(alphRange.max+1, 0);
            std::vector<uint64_t> RB3_alphRuns(alphRange.max+1, 0);

            uint64_t alph = 0, len = 0, newTotalLen = 0, newRuns = 0;


            rlditr_t itr;
            rld_itr_init(rb3->e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?

            rlbwt = sdsl::int_vector<>(runs, 0, alphbits);
            runlens = sdsl::int_vector<>(runs, 0, lenbits);

            int64_t l;
            int c = 0;
            while ((l = rld_dec(rb3->e, &itr, &c, 0)) > 0) {
                alph = static_cast<uint64_t>(c);
                len = static_cast<uint64_t>(l);
                if (alph < alphRange.min || alph > alphRange.max) {
                    std::cerr << "ERROR: Run symbol outside of previously found range, symbol: " << alph << ", Range: " << alphRange << std::endl;
                    exit(1);
                }
                if (len < RB3_lenRange.min || len > RB3_lenRange.max) {
                    std::cerr << "ERROR: Run length outside of previously found range, length: " << len << ", Range: " << lenRange << std::endl;
                    exit(1);
                }
                newTotalLen += len;
                alphCounts[alph] += len;

                ++RB3_alphRuns[alph];
                alphRuns[alph] += (alph == 0)? len : 1;

                if (alph == 0) {
                    while (len) {
                        rlbwt[newRuns] = alph;
                        runlens[newRuns] = 1;
                        ++newRuns;
                        --len;
                    }
                }
                else {
                    rlbwt[newRuns] = alph;
                    runlens[newRuns] = len;
                    ++newRuns;
                }
            }

            if (newRuns != runs) {
                std::cerr << "ERROR: Number of runs found is different in first and second reads. First: " << runs << ", second: " << newRuns << std::endl;
                exit(1);
            }
            if (newTotalLen != totalLen) {
                std::cerr << "ERROR: Total bwt length found is different in first and second reads. First: " << totalLen << ", second: " << newTotalLen << std::endl;
                exit(1);
            }

            std::cout << "Number of runs and number of total occurrences in text for each character is printed below. Format:\n"
                << "\tsymbol\truns in this index\truns in input\toccurrences\n";
            for (uint64_t i = 0; i <= alphRange.max; ++i) 
                std::cout << '\t' << i << '\t' << alphRuns[i] << '\t' << RB3_alphRuns[i] << '\t' << alphCounts[i] << '\n';
        }
        Timer.stop(); //Reading fmd into sdsl

        /*
           Timer.start("Shrinking sdsl");
           std::cout << "Shrinking rlbwt to size" << std::endl;
           std::cout << "Previous element width in bits: " << (int)rlbwt.width() << std::endl;
           sdsl::util::bit_compress(rlbwt);
           std::cout << "New element width in bits: " << (int)rlbwt.width() << std::endl;
           std::cout << "Previous length in number of elements: " << rlbwt.size() << std::endl;
           std::cout << "Previous capacity by number of bits: " << rlbwt.capacity() << std::endl;
           rlbwt.resize(runs);
           std::cout << "New length in number of elements: " << rlbwt.size() << std::endl;
           std::cout << "New capacity by number of bits: " << rlbwt.capacity() << std::endl;
           std::cout << std::endl;
           std::cout << "Shrinking runlens to size" << std::endl;
           std::cout << "Previous element width in bits: " << (int)runlens.width() << std::endl;
           sdsl::util::bit_compress(runlens);
           std::cout << "New element width in bits: " << (int)runlens.width() << std::endl;
           std::cout << "Previous length in number of elements: " << runlens.size() << std::endl;
           std::cout << "Previous capacity by number of bits: " << runlens.capacity() << std::endl;
           runlens.resize(runs);
           std::cout << "New length in number of elements: " << runlens.size() << std::endl;
           std::cout << "New capacity by number of bits: " << runlens.capacity() << std::endl;
           Timer.stop(); //Shrinking sdsl
           */

        Timer.start("Computing RLBWT size");
        std::cout << "Final rlbwt size in bytes: " << sdsl::size_in_bytes(rlbwt) << std::endl;
        std::cout << "Final runlens size in bytes: " << sdsl::size_in_bytes(runlens) << std::endl;
        Timer.stop(); //Computing sdsl size

        /*
        Timer.start("Verifying correctness of constructed RLBWT");

        std::cout << "Built RLBWT and input FMD are " 
            << ((equalToFmi(rlbwt, runlens, *rb3))? "equal" : "not equal!") << std::endl;
        Timer.stop(); //Verifying correctness of constructed RLBWT
        */


        Timer.stop(); //Constructing RLBWT from FMD
    }

    void LFconstruction(uRange alphRange, const std::vector<uint64_t> & alphCounts, std::vector<MoveStructure::IntervalPoint> & alphStarts) {
        Timer.start("Constructing LF from RLBWT");
        uint64_t runs = rlbwt.size();
        LF.intLens = &runlens;
        LF.D_index = sdsl::int_vector<>(rlbwt.size(), 0, sdsl::bits::hi(runs) + 1);
        LF.D_offset= sdsl::int_vector<>(rlbwt.size(), 0, runlens.width());
        alphStarts.resize(alphRange.max+2);
        {
            Timer.start("Computing alphStarts");
            {
                alphStarts[0] = {0, 0, 0};
                uint64_t thisAlph = 0, thisAlphLeft = alphCounts[0];
                uint64_t thisRunStart = 0;

                //at beginning of one iteration of this loop,
                //we are at the beginning of a run in the bwt with thisAlphLeft
                //characters of thisAlph to traverse
                for (uint64_t i = 0; i < runs; ++i) {
                    uint64_t l = runlens[i];
                    //at the beginning of one iteration of this loop, we are at a position
                    //within a run (possibly the first) such that the last position in the 
                    //F array of thisAlph occurs in this run (run i).
                    //!!!!!!!!!!!!!
                    //thisAlphLeft is the number of thisAlph characters left to traverse + the
                    //number of non thisAlph characters before the thisAlph characters in the 
                    //current run
                    //!!!!!!!!!!!!!
                    while (thisAlphLeft < l) {
                        ++thisAlph;
                        alphStarts[thisAlph] = { thisRunStart + thisAlphLeft, i, thisAlphLeft };
                        thisAlphLeft += alphCounts[thisAlph];
                    }

                    thisAlphLeft -= l;
                    thisRunStart += l;
                }

                if (thisAlph != alphRange.max) {
                    std::cerr << "ERROR: When computing alphStarts, didn't end up on largest alphabet value! Largest value: "
                        << alphRange.max << ". Ended up on: " << thisAlph << std::endl;
                    exit(1);
                }
                if (thisAlphLeft) {
                    std::cerr << "ERROR: When computing alphStarts, there are a nonzero amount of "
                        << alphRange.max << " characters left to look for, but the end of the rlbwt has been reached." << std::endl;
                    exit(1);
                }
                if (thisRunStart != totalLen) {
                    std::cerr << "Final position is not total length of bwt in alphStart computation! Final position: " 
                        << thisRunStart << ". Total length: " << totalLen;
                    exit(1);
                }

                //equivalent to alphStarts[thisAlph] = { thisRunStart + thisAlphLeft, i, thisAlphLeft };
                //since thisAlphLeft = 0
                alphStarts[++thisAlph] = { thisRunStart, runs, 0 };
            }
            Timer.stop(); //Computing alphStarts

            std::cout << "Computed alphStarts printed below."
                << " Note: ENDBWT is not a symbol, but is stored as a sentinel in the alphStarts array at the endFormat:\n"
                << "\tsymbol\tposition\trun\toffset\n";
            for (uint64_t i = 0; i <= alphRange.max; ++i)
                std::cout << '\t' << i << '\t' << alphStarts[i].position
                    << '\t' << alphStarts[i].interval
                    << '\t' << alphStarts[i].offset << '\n';
            std::cout << '\t' << "ENDBWT" << '\t' << alphStarts[alphRange.max+1].position
                << '\t' << alphStarts[alphRange.max+1].interval
                << '\t' << alphStarts[alphRange.max+1].offset << '\n';

            Timer.start("Computing run and offsets for LF from run starts");
            {
                std::vector<MoveStructure::IntervalPoint> currentAlphLFs(alphRange.max+1);
                for (uint64_t i = 0; i <= alphRange.max; ++i) 
                    currentAlphLFs[i] = alphStarts[i];

                for (uint64_t i = 0; i < runs; ++i) {
                    uint64_t l = runlens[i], c = rlbwt[i];

                    LF.D_index[i] = currentAlphLFs[c].interval;
                    LF.D_offset[i] = currentAlphLFs[c].offset;

                    LF.AdvanceIntervalPoint_unsafe(currentAlphLFs[c], l);
                }

                //verify all currentAlphLFs ended up at the start of the next alph
                bool allGood = true;
                for (uint64_t i = 0; i <= alphRange.max; ++i) {
                    if (currentAlphLFs[i] != alphStarts[i+1]) {
                        allGood = false;
                        std::cerr << "Symbol " << i << " LF didn't end up at the start of "
                            << " symbol " << i + 1 << " in F array.\nIt ended up at: "
                            << "{ position: " << currentAlphLFs[i].position
                            << ", interval: " << currentAlphLFs[i].interval
                            << ", offset: " << currentAlphLFs[i].offset << " }.\n"
                            << "The next symbol starts at " 
                            << "{ position: " << alphStarts[i+1].position
                            << ", interval: " << alphStarts[i+1].interval
                            << ", offset: " << alphStarts[i+1].offset << " }."
                            << std::endl;
                    }
                }

                if (allGood)
                    std::cout << "The LFs of all symbols ended up at the start of the next symbol.\n";
                else {
                    std::cerr << "ERROR: The LFs of some symbols didn't end up at the start of the next symbol. (See above)." << std::endl;
                    exit(1);
                }
            }
            Timer.stop(); //Computing run and offsets for LF from run starts

            Timer.start("Computing LF size");
            std::cout << "Final LF.D_index size in bytes: " << sdsl::size_in_bytes(LF.D_index) << std::endl;
            std::cout << "Final LF.D_offset size in bytes: " << sdsl::size_in_bytes(LF.D_offset) << std::endl;
            Timer.stop(); //Computing LF size

//            /*
//            Timer.start("Verifying computed LF by sum of sequence lengths");
//            {
//                uint64_t sumSeqLengths = 0;
//                //This below comment is outdated. It is kept for informational purposes and posterity.
//                //We now keep each dollar in a separate run. We still verify by sequence
//                //because it is easily parallelizable and complete LF traversal is by far the slowest part of construction.
//                //(Naturally, since it is the only O(n) instead of O(r) part. It's not even O(n) because I haven't
//                //split the intervals of the move data structure yet.)
//                //-----------------------------------------------------------------------------------------------------------
//                //In theory, this is a multi-dollar BWT. I.E. for strings, s_0, s_1, s_2,..., s_k,
//                //the text is s_0,$_0,s_1,$_1,s_2,$_2,...,s_k,$_k concatenated in that order,
//                //where $_i < $_j for i < j and $_i < a for all non $_* characters in the text
//                //HOWEVER: In practice, the multi-dollars incur a big cost by increasing the alphabet size
//                //therefore they are usually all treated as $. This costs us something: we are no longer able
//                //to perform LF on the BWT locations with $. If desired, the value i in BWT[x] = $_i must be recovered.
//                //then LF[x] = i. We forgo this ability in implementation because we don't care about it and costs index size
//
//                //However this does make verifying the LF function more difficult. (Typically, we could just LF n times
//                //and verify that we return to the origin, then LF is a permutation with one cycle and likely correct.)
//                //Here instead, we indpentently LF through each string contained in the text, sum up their lengths, 
//                //and claim that if the sum of the lengths is accurate then the LF is likely correct. We could
//                //also keep a bool array and make sure every position of the BWT is traversed exactly once,
//                //but this would be costly in terms of memory for large compressed BWTs (150 GB for human472)
//                //compute starts
//                //-----------------------------------------------------------------------------------------------------------
//                std::vector<IntervalPoint> starts(alphCounts[0]);
//                IntervalPoint start{ (uint64_t)-1, 0, 0};
//                starts[0] = start;
//                for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
//                    ++start.offset;
//                    if (start.offset == runlens[start.interval]) {
//                        start.offset = 0;
//                        ++start.interval;
//                    }
//                    starts[seq] = start;
//                }
//
//                //count the length of sequence seq
//                #pragma omp parallel for schedule(guided)
//                for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
//                    IntervalPoint current{starts[seq]};
//                    uint64_t seqLen = 1; //counting LF to endmarker, which isn't performed in actuality (simulated by incrementing start.offset...)
//                    while (rlbwt[current.interval] != 0) {
//                        ++seqLen;
//                        current = );
//                    }
//                    #pragma omp critical 
//                    {
//                        sumSeqLengths += seqLen;
//                        stringStarts[seq] = current;
//                    }
//                }
//
//                if (sumSeqLengths != totalLen) {
//                    std::cerr << "ERROR: Sum of sequence length not equal to the length of the BWT! Sequence length is computed by LF.\n"
//                        << "ERROR: A total of " << sumSeqLengths << " LFs were computed. This is not equal to the length of the BWT, which is "
//                        << totalLen << "." << std::endl;
//                    return 1;
//                }
//                std::cout << "LF is likely correct, the sum of sequence lengths computed by it is " << sumSeqLengths << ". " 
//                    << sumSeqLengths - alphCounts[0] << " LFs were computed in order to verify this.\n";
//            }
//            Timer.stop(); //Verifying computed LF by sum of sequence lengths
//            */
        }
        Timer.stop(); //Constructing LF from RLBWT
    }

    void SAconstruction(const std::vector<uint64_t> & alphCounts, std::vector<uint64_t> & seqNumsTopOrBotRun, std::vector<uint64_t> & seqLens, std::vector<MoveStructure::IntervalPoint>& stringStarts) {
        Timer.start("SA sampling");
        stringStarts.resize(alphCounts[0]);

        //suffixes are 0-indexed
        //sequences are 0-indexed

        //suffix array samples at the top of runs
        //the suffix at the top of run i is the SATopRunInt[i] interval in the text order
        //sdsl::int_vector<> SATopRunInt;
        //SABotRunInt is not needed after the data structure is built but is needed to build the data structre (for sampling)
        sdsl::int_vector<> SABotRunInt;

        sdsl::int_vector<> runSampledAt;

        //Samples in Text order
        //InvertibleMoveStructure PhiInvPhi;

        uint64_t maxSeqLen = 0;
        seqLens.resize(alphCounts[0]); //seq lengths, counts endmarker so the empty string has length 1
                                       //number of times each string is at the top (and bottom respectively) of a run, 
                                       //counts endmarker, so value for the empty string would be 1
        std::vector<uint64_t> seqNumsTopRun(alphCounts[0]), seqNumsBotRun(alphCounts[0]); 
        seqNumsTopOrBotRun.resize(alphCounts[0]); 
        {
            uint64_t sumSeqLengths = 0;
            uint64_t maxIntLen = 0;
            Timer.start("Auxiliary info computation (seqLens, seqNumsTopRun, seqNumsBotRun, seqNumsTopOrBotRun)");
            {
                std::vector<MoveStructure::IntervalPoint> starts(alphCounts[0]);
                MoveStructure::IntervalPoint start{ static_cast<uint64_t>(-1), 0, 0};
                starts[0] = start;
                for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                    ++start.offset;
                    if (start.offset == runlens[start.interval]) {
                        start.offset = 0;
                        ++start.interval;
                    }
                    starts[seq] = start;
                }

                #pragma omp parallel for schedule(guided)
                for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                    MoveStructure::IntervalPoint current{starts[seq]};
                    uint64_t seqLen = 0, seqNumTopRun = 0, seqNumBotRun = 0, seqNumTopOrBotRun = 0, seqNumOneRun = 0;
                    uint64_t currentTopOrBotIntervalLen = 0; //number of characters in the current interval on the sequence since the last time the sequence was at the top or bottom of a run
                    uint64_t maxTopOrBotIntervalLen = 0;
                    while (rlbwt[current.interval] != 0) {
                        //add suffix at current in SA
                        ++seqLen;
                        ++currentTopOrBotIntervalLen;
                        seqNumOneRun += runlens[current.interval] == 1;
                        seqNumTopRun += (current.offset == 0);
                        seqNumBotRun += (current.offset == runlens[current.interval] - 1);
                        if ((current.offset == 0) || (current.offset == runlens[current.interval] - 1)) {
                            seqNumTopOrBotRun++;
                            maxTopOrBotIntervalLen = std::max(maxTopOrBotIntervalLen, currentTopOrBotIntervalLen);
                            currentTopOrBotIntervalLen = 0;
                        }

                        current = LF.map(current);
                    } 
                    //add suffix 0 of seq
                    ++seqLen;
                    ++currentTopOrBotIntervalLen;
                    seqNumOneRun += runlens[current.interval] == 1;
                    seqNumTopRun += (current.offset == 0);
                    seqNumBotRun += (current.offset == runlens[current.interval] - 1);
                    if ((current.offset == 0) || (current.offset == runlens[current.interval] - 1)) {
                        seqNumTopOrBotRun++;
                        maxTopOrBotIntervalLen = std::max(maxTopOrBotIntervalLen, currentTopOrBotIntervalLen);
                        currentTopOrBotIntervalLen = 0;
                    }
                    else {
                        std::cerr << "ERROR: suffix 0 of seq " << seq << " not at top or bottom of run (should be both, since in a run of length 1!\n";
                        exit(1);
                    }

                    //save results
                    #pragma omp critical
                    {
                        sumSeqLengths += seqLen;
                        stringStarts[seq] = current;

                        seqLens[seq] = seqLen;
                        seqNumsTopRun[seq] = seqNumTopRun;
                        seqNumsBotRun[seq] = seqNumBotRun;
                        seqNumsTopOrBotRun[seq] = seqNumTopOrBotRun;

                        maxIntLen = std::max(maxIntLen, maxTopOrBotIntervalLen);

                        if (seqNumTopOrBotRun != seqNumTopRun + seqNumBotRun - seqNumOneRun) {
                            std::cerr << "Number of times at bottom, top, boundary, and number runs of length one not consistent for seq " 
                                << seq << "!\nTop: " << seqNumTopRun
                                << "\nBot: " << seqNumBotRun 
                                << "\nTop or Bot (Boundary): " << seqNumTopOrBotRun
                                << "\nRuns of length one: " << seqNumOneRun << "\n";
                            exit(1);
                        }
                    }
                }
            }
            Timer.stop(); //Auxiliary info computation (seqLens, seqNumsTopRun, seqNumsBotRun)

            if (sumSeqLengths != totalLen) {
                std::cerr << "ERROR: Sum of sequence length not equal to the length of the BWT! Sequence length is computed by LF.\n"
                    << "ERROR: A total of " << sumSeqLengths << " LFs were computed. This is not equal to the length of the BWT, which is "
                    << totalLen << "." << std::endl;
                exit(1);
            }
            std::cout << "LF is likely correct, the sum of sequence lengths computed by it is " << sumSeqLengths << ". " 
                << sumSeqLengths - alphCounts[0] << " LFs were computed in order to verify this.\n";

            Timer.start("Prefix summing auxiliary info");
            maxSeqLen = seqLens[0];
            for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                maxSeqLen = std::max(maxSeqLen, seqLens[seq]);

                seqNumsTopRun[seq] += seqNumsTopRun[seq-1];
                seqNumsBotRun[seq] += seqNumsBotRun[seq-1];
                seqNumsTopOrBotRun[seq] += seqNumsTopOrBotRun[seq-1];
            }
            Timer.stop(); //Prefix summing auxiliary info

            std::cout << "Maximum sequence length: " << maxSeqLen << '\n';
            std::cout << "Maximum interval length: " << maxIntLen << '\n';

            uint64_t runs = rlbwt.size();

            if (seqNumsTopRun.back() != runs) {
                std::cerr << "ERROR:seqNumsTopRun.back() != runs, " << seqNumsTopRun.back() << " != " << runs << '\n';
                exit(1);
            }
            if (seqNumsBotRun.back() != runs) {
                std::cerr << "ERROR: seqNumsBotRun.back() != runs, " << seqNumsBotRun.back() << " != " << runs << '\n';
                exit(1);
            }
            std::cout << "Total run tops: " << seqNumsTopRun.back() << '\n'
                << "Total run bots: " << seqNumsBotRun.back() << '\n'
                << "Total run top and bots: " << seqNumsTopOrBotRun.back() << '\n'
                << "Therefore, total runs of length 1 [(tops + bots) - (tops and bots)]: " << seqNumsTopRun.back() + seqNumsBotRun.back() - seqNumsTopOrBotRun.back() << '\n';

            SATopRunInt = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
            SABotRunInt = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);

            PL.SeqAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(alphCounts[0] - 1) + 1);
            PL.PosAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxSeqLen - 1) + 1);

            PL.IntLen = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen) + 1);

            PL.phi.D_index = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
            PL.phi.D_offset = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen - 1) + 1);
            PL.invPhi.D_index = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(seqNumsTopOrBotRun.back() - 1) + 1);
            PL.invPhi.D_offset = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(maxIntLen - 1) + 1);


            runSampledAt = sdsl::int_vector<>(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(runs - 1) + 1);

            Timer.start("Sampling");
            {
                std::vector<MoveStructure::IntervalPoint> starts(alphCounts[0]);
                MoveStructure::IntervalPoint start{ static_cast<uint64_t>(-1), 0, 0};
                starts[0] = start;
                for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
                    ++start.offset;
                    if (start.offset == runlens[start.interval]) {
                        start.offset = 0;
                        ++start.interval;
                    }
                    starts[seq] = start;
                }

                Timer.start("Sampling the SA order");
                #pragma omp parallel for schedule(guided)
                for(uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                    MoveStructure::IntervalPoint current{starts[seq]};
                    uint64_t pos = seqLens[seq] - 1;
                    uint64_t prevPos = seqLens[seq];
                    uint64_t posTopRun = seqNumsTopRun[seq] - 1;
                    uint64_t posBotRun = seqNumsBotRun[seq] - 1;
                    uint64_t posRun = seqNumsTopOrBotRun[seq] - 1;
                    uint64_t finalTopRun = (seq == 0)? 0 : seqNumsTopRun[seq-1];
                    uint64_t finalBotRun = (seq == 0)? 0 : seqNumsBotRun[seq-1];
                    uint64_t finalRun = (seq == 0)? 0 : seqNumsTopOrBotRun[seq-1];
                    while (rlbwt[current.interval] != 0) {
                        if (current.offset == 0 || current.offset == runlens[current.interval] - 1) {
                            if (current.offset == 0) {
                                #pragma omp critical
                                {
                                    if (SATopRunInt[current.interval] != 0) {
                                        std::cerr << "ERROR: this run's top sample has already been set!\n";
                                    }
                                    SATopRunInt[current.interval] = posRun;
                                    if (posTopRun == finalTopRun) {
                                        std::cerr << "ERROR: posTopRun to reach last position before endmarker found!\n";
                                    }
                                }
                                --posTopRun;
                            }
                            if (current.offset == runlens[current.interval] - 1) {
                                #pragma omp critical
                                {
                                    if (SABotRunInt[current.interval] != 0) {
                                        std::cerr << "ERROR: this run's bot sample has already been set!\n";
                                    }
                                    SABotRunInt[current.interval] = posRun;
                                    if (posBotRun == finalBotRun) {
                                        std::cerr << "ERROR: posBotRun to reach last position before endmarker found!\n";
                                    }
                                }
                                --posBotRun;
                            }

                            #pragma omp critical
                            {
                                if (PL.SeqAt[posRun] != 0 || PL.PosAt[posRun] != 0) {
                                    std::cerr << "ERROR: this interval's phi sample has already been set!\n";
                                }
                                PL.SeqAt[posRun] = seq;
                                PL.PosAt[posRun] = pos;
                                PL.IntLen[posRun] = prevPos - pos;
                                runSampledAt[posRun] = current.interval;
                                if (posRun == finalRun) {
                                    std::cerr << "ERROR: posRun to reach last position before endmarker found!\n";
                                }
                            }
                            prevPos = pos;
                            --posRun;

                        }
                        if (pos == 0) {
                            std::cerr << "ERROR: pos in sequence to reach -1 before endmarker found!\n";
                        }
                        --pos;
                        current = LF.map(current);
                    }
                    if (pos != 0) {
                        std::cerr << "ERROR: full sequence not traversed!\n";
                    }
                    if (posTopRun != finalTopRun) {
                        std::cerr << "ERROR: didn't traverse all top samples of this sequence!\n";
                    }
                    if (posBotRun != finalBotRun) {
                        std::cerr << "ERROR: didn't traverse all bottom samples of this sequence!\n";
                    }
                    if (posRun != finalRun) {
                        std::cerr << "ERROR: didn't traverse all interval samples of this sequence!\n";
                    }
                    if (current.offset == 0) {
                        #pragma omp critical
                        {
                            SATopRunInt[current.interval] = posRun;
                        }
                        --posTopRun;

                        PL.SeqAt[posRun] = seq;
                        PL.PosAt[posRun] = pos;
                        PL.IntLen[posRun] = prevPos - pos;
                        runSampledAt[posRun] = current.interval;
                    }
                    else {
                        std::cerr << "ERROR: offset != 0 in last run, endmarker run!\n";
                    }
                    if (current.offset == runlens[current.interval] - 1) {
                        #pragma omp critical
                        {
                            SABotRunInt[current.interval] = posRun;
                        }
                        --posBotRun;
                    }
                    else {
                        std::cerr << "ERROR: offset not at end of run in last run, endmarker run!\n";
                    }
                }
                Timer.stop(); //Sampling in SA order

                Timer.start("Sampling the Text order");
                //#pragma omp parallel for schedule(guided)
                for (uint64_t seq = 0; seq < alphCounts[0]; ++seq) {
                    //go from beginning of sequence to end, by interval in phi
                    uint64_t currentIntervalIndex = (seq == 0)? 0 :  seqNumsTopOrBotRun[seq-1];
                    uint64_t seqTraversed = 0;

                    //do first interval
                    if (runlens[runSampledAt[currentIntervalIndex]] != 1) {
                        std::cerr << "ERROR: first interval of sequences starts at a run of length not equal to 1"
                            << " (should be 1 since bwt value should be endmarker of previous sequence)!\n";
                    }
                    if (PL.SeqAt[currentIntervalIndex] != seq)
                        std::cerr << "ERROR: seq at interval doesn't match seq in Phi!\n";
                    if (SATopRunInt[runSampledAt[currentIntervalIndex]] != SABotRunInt[runSampledAt[currentIntervalIndex]])
                        std::cerr << "ERROR: top and bottom run sequence samples don't match in endmarker run (should be length 1)!\n";

                    PL.phi.D_index[currentIntervalIndex] = SABotRunInt[(runSampledAt[currentIntervalIndex] == 0)? SABotRunInt.size()-1 : runSampledAt[currentIntervalIndex]-1];
                    PL.phi.D_offset[currentIntervalIndex] = 0;
                    PL.phi.intLens = &PL.IntLen;
                    PL.invPhi.D_index[currentIntervalIndex] = SATopRunInt[(runSampledAt[currentIntervalIndex] == SATopRunInt.size() - 1)? 0 : runSampledAt[currentIntervalIndex]+1];
                    PL.invPhi.D_offset[currentIntervalIndex] = 0;
                    PL.invPhi.intLens = &PL.IntLen;

                    //std::cout << "PhiInvPhi.phi.D_index[currentIntervalIndex]: " << PhiInvPhi.phi.D_index[currentIntervalIndex] << std::endl;
                    //std::cout << "PhiInvPhi.invPhi.D_index[currentIntervalIndex]: " << PhiInvPhi.invPhi.D_index[currentIntervalIndex] << std::endl;

                    seqTraversed += PL.IntLen[currentIntervalIndex];
                    ++currentIntervalIndex;

                    //std::cout << "Here" << std::endl;
                    //std::cout << "seqTraversed: " << seqTraversed << ". currentIntervalIndex: " << currentIntervalIndex << std::endl;

                    while (seqTraversed < seqLens[seq]) {
                        //std::cout << "In here" << std::endl;
                        //add next interval
                        uint64_t runIndex = runSampledAt[currentIntervalIndex];
                        uint64_t topRunInt = SATopRunInt[runIndex];
                        uint64_t botRunInt = SABotRunInt[runIndex];
                        if (!((seq == PL.SeqAt[topRunInt] && seqTraversed == PL.PosAt[topRunInt]) ||
                                    (seq == PL.SeqAt[botRunInt] && seqTraversed == PL.PosAt[botRunInt])))
                            std::cerr << "ERROR: Beginning of run interval in sequence is not equal to the sample at"
                                << " the beginning or the end of the corresponding interval!\n";
                        //std::cout << "AFter this?" << std::endl;
                        //computing above sample
                        if (seq == PL.SeqAt[topRunInt] && seqTraversed == PL.PosAt[topRunInt]) {
                            uint64_t runAboveIndex = (runIndex == 0)? runs - 1 : runIndex - 1;
                            #pragma omp critical
                            {
                                PL.phi.D_index[currentIntervalIndex] = SABotRunInt[runAboveIndex];
                                PL.phi.D_offset[currentIntervalIndex] = 0;
                            }
                        }
                        else {
                            //use previous above sample
                            uint64_t prevInt = PL.phi.D_index[currentIntervalIndex - 1];
                            uint64_t prevOffset = PL.phi.D_offset[currentIntervalIndex - 1];
                            prevOffset += PL.IntLen[currentIntervalIndex - 1];
                            while (prevOffset >= PL.IntLen[prevInt]) {
                                prevOffset -= PL.IntLen[prevInt];
                                ++prevInt;
                            }
                            #pragma omp critical
                            {
                                PL.phi.D_index[currentIntervalIndex] = prevInt;
                                PL.phi.D_offset[currentIntervalIndex] = prevOffset;
                            }
                        }

                        //std::cout << "How about this?" << std::endl;

                        //computing below sample
                        if (seq == PL.SeqAt[botRunInt] && seqTraversed == PL.PosAt[botRunInt]) {
                            //std::cout << "In if" << std::endl;
                            uint64_t runBelowIndex = (runIndex == runs - 1)? 0 : runIndex + 1;
                            #pragma omp critical
                            {
                                PL.invPhi.D_index[currentIntervalIndex] = SATopRunInt[runBelowIndex];
                                PL.invPhi.D_offset[currentIntervalIndex] = 0;
                            }
                        }
                        else {
                            //std::cout << "In else" << std::endl;
                            //use previous below sample
                            uint64_t prevInt = PL.invPhi.D_index[currentIntervalIndex - 1];
                            uint64_t prevOffset = PL.invPhi.D_offset[currentIntervalIndex - 1];
                            prevOffset += PL.IntLen[currentIntervalIndex - 1];
                            //std::cout << "Got here" << std::endl;
                            while (prevOffset >= PL.IntLen[prevInt]) {
                                //std::cout << "prevInt: " << prevInt << ". prevOffset: " << prevOffset << std::endl;
                                prevOffset -= PL.IntLen[prevInt];
                                ++prevInt;
                            }
                            //std::cout << "prevInt: " << prevInt << ". prevOffset: " << prevOffset << std::endl;
                            //std::cout << "passed while loop" << std::endl;
                            #pragma omp critical
                            {
                                PL.invPhi.D_index[currentIntervalIndex] = prevInt;
                                PL.invPhi.D_offset[currentIntervalIndex] = prevOffset;
                            }
                        }

                        //std::cout << "How aboutttt this?" << std::endl;

                        seqTraversed += PL.IntLen[currentIntervalIndex];
                        ++currentIntervalIndex;
                    }

                    //std::cout << "Here after" << std::endl;

                    #pragma omp critical
                    if (seqTraversed != seqLens[seq]) {
                        std::cerr << "ERROR: Traversed a sequence some length not equal to it's actual length!\n";
                        std::cerr << "ERROR: Traversed seq " << seq << " " << seqTraversed << " characters. Actual length: " << seqLens[seq] << '\n';
                    }
                    #pragma omp critical
                    if (currentIntervalIndex != seqNumsTopOrBotRun[seq]) {
                        std::cerr << "ERROR: Traversed sequence but ended up on some interval in PL other than the first interval of the next sequence!\n";
                        std::cerr << "ERROR: Ended up on interval " << currentIntervalIndex << ". Should have ended up on " << seqNumsTopOrBotRun[seq] << '\n';
                    }
                }
                Timer.stop(); //Sampling in Textorder
            }
            Timer.stop(); //Sampling
        }
        Timer.stop(); //SA sampling
    }

    void LCPconstruction(const std::vector<MoveStructure::IntervalPoint> & alphStarts, const std::vector<uint64_t> & seqNumsTopOrBotRun, const std::vector<uint64_t> & seqLens) {
        Timer.start("LCP computation");
        uint64_t numStrings = alphStarts[1].position;

        std::vector<MoveStructure::IntervalPoint> starts(numStrings);
        //std::vector<MoveStructure::IntervalPoint> revEquivLF(seqNumsTopOrBotRun.back());
        sdsl::int_vector<> revEquivLF_index(seqNumsTopOrBotRun.back(), 0, sdsl::bits::hi(seqNumsTopOrBotRun.back()-1) + 1);
        sdsl::int_vector<> revEquivLF_offset(seqNumsTopOrBotRun.back(), 0, runlens.width());
        MoveStructure::IntervalPoint start{ static_cast<uint64_t>(-1), 0, 0};
        starts[0] = start;
        for (uint64_t seq = 1; seq < numStrings; ++seq) {
            ++start.offset;
            if (start.offset == runlens[start.interval]) {
                start.offset = 0;
                ++start.interval;
            }
            starts[seq] = start;
        }

        Timer.start("Reverse sampling");
        #pragma omp parallel for schedule(guided)
        for (uint64_t seq = 0; seq < numStrings; ++seq) {
            //traverse seq from end to beginning, storing samples of LF position of seq in the intervals of rev(seq) in the PhiInvPHi data structure
            uint64_t revSeq = (seq%2 == 0)? seq + 1 : seq - 1;
            MoveStructure::IntervalPoint current{starts[seq]};
            //uint64_t posSeq = seqLens[seq] - 1;
            uint64_t revSeqIntervalIndex = (revSeq == 0)? 0 : seqNumsTopOrBotRun[revSeq - 1];
            uint64_t revSeqIntervalOffset = 0;
            uint64_t finalRevSeqIntervalIndex = seqNumsTopOrBotRun[revSeq];

            //check if seq and revSeq lengths are equal
            #pragma omp critical
            if (seqLens[seq] != seqLens[revSeq]) {
                std::cerr << "ERROR: Length of sequence (" << seq << ") and its reverse (" << revSeq << ") are not equal! Respectively: " 
                    << seqLens[seq] << " and " << seqLens[revSeq] << "!\n";
            }

            while (revSeqIntervalIndex != finalRevSeqIntervalIndex) {
                if (revSeqIntervalOffset == 0) {
                    #pragma omp critical
                    {
                        if (revEquivLF_index[revSeqIntervalIndex] != 0 || revEquivLF_offset[revSeqIntervalIndex] != 0)
                            std::cerr << "ERROR: setting already set revEquivLF in reverse sampling!\n";
                        //revEquivLF[revSeqIntervalIndex] = current;
                        revEquivLF_index[revSeqIntervalIndex] = current.interval;
                        revEquivLF_offset[revSeqIntervalIndex] = current.offset;
                    }
                }

                current = LF.map(current);
                ++revSeqIntervalOffset;
                if (revSeqIntervalOffset == PL.IntLen[revSeqIntervalIndex]) {
                    revSeqIntervalOffset = 0;
                    ++revSeqIntervalIndex;
                }
            }
        }
        Timer.stop(); //Reverse sampling

        Timer.start("LCP Sampling");
        {
            PL.AboveLCP.resize(seqNumsTopOrBotRun.back());

//            Timer.start("locating alphStarts in PL");
//            uint64_t sampledAlphStarts = 0;
//            std::vector<MoveStructure::IntervalPoint> phiAlphStarts(alphStarts.size()-1);
//            for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
//                //std::cout << "Hre" << std::endl;
//                phiAlphStarts[al] = {uint64_t(-1), SATopRunInt[alphStarts[al].interval], 0};
//                for (uint64_t i = 0; i < alphStarts[al].offset; ++i)
//                    phiAlphStarts[al] = PL.invPhi.map(phiAlphStarts[al]);
//                sampledAlphStarts += phiAlphStarts[al].offset == 0;
//            }
//            Timer.stop(); //locating alphStarts in PL
            //std::cout << "Hredone" << std::endl;

//            std::cout << "Of " << phiAlphStarts.size() << ", " << sampledAlphStarts << " alphStarts are at the beginning or end of a run in the BWT\n";
//            std::cout << "Locations of alphStarts in PL intervals:\n\tsymbol\tinterval\toffset\n";
//            for (uint64_t al = 0; al < phiAlphStarts.size(); ++al)
//                std::cout << '\t' << al
//                    << '\t' << phiAlphStarts[al].interval
//                    << '\t' << phiAlphStarts[al].offset << '\n';
//
//            
//            std::vector<bool> alphStartFound(phiAlphStarts.size(), false);
            //compute only for those with phi.D_offset = 0
            #pragma omp parallel for schedule(guided)
            for (uint64_t seq = 0; seq < numStrings; ++seq) {
                uint64_t startInterval = (seq == 0)? 0 : seqNumsTopOrBotRun[seq - 1];
                for (uint64_t currentInterval = startInterval; currentInterval < seqNumsTopOrBotRun[seq]; ++currentInterval) {
                    if (PL.phi.D_offset[currentInterval] != 0) {
                        //#pragma omp critical
                        //{
                            PL.AboveLCP[currentInterval] = PL.AboveLCP[currentInterval -1] - PL.IntLen[currentInterval - 1];
                        //}
//                        //this whole if statement is just error checking that could (hopefully) be removed
//                        if (PL.AboveLCP[currentInterval] < PL.IntLen[currentInterval]){
//                            bool found = false;
//                            if (PL.AboveLCP[currentInterval] == PL.IntLen[currentInterval] - 1) {
//                                //possibly at the beginning of alphabet in F column, LCP of 0 is correct
//                                for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
//                                    if (phiAlphStarts[al].interval != currentInterval || phiAlphStarts[al].offset != PL.AboveLCP[currentInterval])
//                                        continue;
//                                    if (alphStartFound[al]) {
//                                        #pragma omp critical 
//                                        {
//                                            std::cerr << "ERROR: already found this alplhStart in phiinv data structure!\n";
//                                        }
//
//                                    }
//                                    #pragma omp critical 
//                                    {
//                                        alphStartFound[al] = true;
//                                    }
//                                    found = true;
//                                }
//                            }
//                            if (currentInterval == seqNumsTopOrBotRun[seq] - 1) {
//                                //LCP should be interval length - 1
//                                if (PL.AboveLCP[currentInterval] != PL.IntLen[currentInterval] - 1) {
//                                    #pragma omp critical 
//                                    {
//                                        std::cerr << "ERROR: Last interval of sequence has weird matching length, maybe including endmarker!\n";
//                                    }
//                                }
//                                else 
//                                    found = true;
//                            }
//                            if (!found) {
//                                #pragma omp critical 
//                                {
//                                    std::cerr << "ERROR: forwarded LCP computed smaller than interval length and is not one of alphStarts!\n"
//                                        << "ERROR: forwarded: " << PL.AboveLCP[currentInterval] << ". Interval length: " << PL.IntLen[currentInterval] << '\n';
//                                }
//                            }
//                        }
                        continue;
                    }

                    uint64_t matchingLength = 0;
                    MoveStructure::IntervalPoint revSeq{static_cast<uint64_t>(-1), revEquivLF_index[currentInterval], revEquivLF_offset[currentInterval]};
                    MoveStructure::IntervalPoint revSeqAbove{static_cast<uint64_t>(-1), revEquivLF_index[PL.phi.D_index[currentInterval]], revEquivLF_offset[PL.phi.D_index[currentInterval]]};

                    while (rlbwt[revSeq.interval] != 0 && rlbwt[revSeqAbove.interval] != 0 && rlbwt[revSeq.interval] == rlbwt[revSeqAbove.interval]) {
                        revSeq = LF.map(revSeq);
                        revSeqAbove = LF.map(revSeqAbove);
                        ++matchingLength;
                    }

                    //#pragma omp critical
                    //{
                        PL.AboveLCP[currentInterval] = matchingLength;
                    //}
                    //this whole if statement is just error checking that could (hopefully) be removed
//                    if (matchingLength < PL.IntLen[currentInterval]) {
//                        bool found = false;
//                        if (matchingLength == PL.IntLen[currentInterval] - 1) {
//                            //possibly at the beginning of alphabet in F column, LCP of 0 is correct
//                            for (uint64_t al = 0; al < phiAlphStarts.size(); ++al) {
//                                if (phiAlphStarts[al].interval != currentInterval || phiAlphStarts[al].offset != matchingLength)
//                                    continue;
//                                if (alphStartFound[al]) {
//                                    #pragma omp critical
//                                    {
//                                        std::cerr << "ERROR: already found this alphStart in the PL data structure!\n";
//                                    }
//                                }
//                                #pragma omp critical 
//                                {
//                                    alphStartFound[al] = true;
//                                }
//                                found = true;
//                            }
//                        }
//                        if (currentInterval == seqNumsTopOrBotRun[seq] - 1) {
//                            //LCP should be interval length - 1
//                            if (matchingLength !=  PL.IntLen[currentInterval] - 1) {
//                                #pragma omp critical 
//                                {
//                                    std::cerr << "ERROR: Last interval of sequence has weird matching that might include the endmarker!\n";
//                                }
//                            }
//                            else 
//                                found = true;
//                        }
//                        if (!found) {
//                            #pragma omp critical
//                            {
//                                std::cerr << "ERROR: Computed LCP smaller than interval length and is not one of alphStarts!\n"
//                                    << "ERROR: Computed: " << matchingLength <<". Interval length: " << PL.IntLen[currentInterval] << " for interval: " << currentInterval << '\n';
//                            }
//                        }
//                    }
                }
            }

//            bool re = true;
//            for (auto t : alphStartFound)
//                re = re && t;
//            if (re)
//                std::cout << "Found all " << sampledAlphStarts << " sampled alplh starts.\n";
//            else 
//                std::cerr << "ERROR: Didn't find all alph starts!\n";

            Timer.stop(); //LCP Sampling
            /*
            Timer.start("Printing Structures");
            printStructures(10, rlbwt, runlens, toRun, toOffset, SATopRunInt, SABotRunInt, PhiInvPhi, runSampledAt, revEquivLF);
            Timer.stop(); //Printing Structures
            */
        }
        Timer.stop(); //LCP computation

        Timer.start("Shrinking LCP");
        std::cout << "Shrinking AboveLCP to size\n";
        std::cout << "Previous element width in bits: " << static_cast<int>(PL.AboveLCP.width()) << '\n';
        sdsl::util::bit_compress(PL.AboveLCP);
        std::cout << "New element width in bits: " << static_cast<int>(PL.AboveLCP.width()) << '\n';
        Timer.stop(); //Shrinking LCP
    }

    void endmarkerRepair(const std::vector<uint64_t> & alphCounts, const std::vector<MoveStructure::IntervalPoint>& stringStarts) {
        Timer.start("RLBWT Repair");
        Timer.start("Detecting endmarkers in runs in RLBWT");
        //set of runIDs that are runs of multiple endmarkers in the RLBWT
        for (const auto& intPoint : stringStarts) {
            if (intPoint.offset) {
                std::cerr << "ERROR: There is a run of endmarkers in the RLBWT with length > 1!\n";
                exit(1);
            }
        }
        std::cout << "Every endmarker in the RLBWT is in a run of length 1.\n";
        Timer.stop(); //Detecting endmarkers in runs in RLBWT

        Timer.start("Correcting LFs of endmarkers");
        MoveStructure::IntervalPoint dollarSignF{ static_cast<uint64_t>(-1), 0, 0};
        for (uint64_t seq = 1; seq < alphCounts[0]; ++seq) {
            LF.D_index[stringStarts[seq].interval] = dollarSignF.interval;
            LF.D_offset[stringStarts[seq].interval] = dollarSignF.offset;
            LF.AdvanceIntervalPoint_unsafe(dollarSignF, 1);
        }
        LF.D_index[stringStarts[0].interval] = dollarSignF.interval;
        LF.D_offset[stringStarts[0].interval] = dollarSignF.offset;
        Timer.stop(); //Correcting LFs of endmarkers

        /*
        //This code is correct, but is commented out because it takes very long for BWTs with large n (O(n) non-constant LF calls)
        Timer.start("Testing all LFs by permutation with one cycle");
        IntervalPoint start{ (uint64_t)-1, 0, 0}, current{ (uint64_t)-1, 0, 0};
        uint64_t lfsDone = 0;
        uint64_t outputInterval = 1e8;

        Timer.start("Timing LFs 0 to " + std::to_string(std::min(totalLen, outputInterval - 1)));
        do {
            current = mapLF(current, runlens, toRun, toOffset);
            ++lfsDone;
            if (lfsDone % outputInterval == outputInterval - 1) {
                Timer.stop();
                Timer.start("Timing LFs " + std::to_string(lfsDone + 1) + " to " + std::to_string(std::min(totalLen, lfsDone + outputInterval)));
            }
        } while (lfsDone <= totalLen && current != start);
        Timer.stop();

        if (lfsDone != totalLen) {
            std::cerr << "ERROR: LF Invalid after endmarker repair. " << lfsDone 
                << " LFs done in an attempt to traverse one cycle of the BWT. Length is actually " << totalLen
                << ". NOTE: if sumSeqLengths = length + 1, LF (very likely) was terminated early and didn't compute "
                << "the full sequence lengths. This is because verification is automatically terminated after length + 1"
                << " LFs because the LF function is then known to be incorrect." << std::endl;
            return 1;
        }

        std::cout << "LF is a permutation with one cycle, therefore it is very likely correct.\n";
        Timer.stop(); //Testing all LFs by permutation with one cycle
        */
        Timer.stop(); //RLBWT Repair
    }

    public:
    struct LFPhiCoordinate {
        MoveStructure* LF;
        InvertibleMoveStructure* PL;
        sdsl::int_vector<>* SATopRunInt;

        MoveStructure::IntervalPoint LFpoint;
        MoveStructure::IntervalPoint phiPoint;

        void doPhi() {
            if (LFpoint.offset)
                --LFpoint.offset;
            else {
                LFpoint.interval = (LFpoint.interval == 0)? LF->intLens->size() - 1 : LFpoint.interval - 1;
                LFpoint.offset = (*LF->intLens)[LFpoint.interval] - 1;
            }

            phiPoint = PL->phi.map(phiPoint);
        }

        uint64_t LCP() {
            return PL->LCP(phiPoint);
        }

        void doInvPhi() {
            if (++LFpoint.offset == (*LF->intLens)[LFpoint.interval]) {
                ++LFpoint.interval;
                LFpoint.interval = LFpoint.interval % LF->intLens->size();
                LFpoint.offset = 0;
            }

            phiPoint = PL->invPhi.map(phiPoint);
        }

        void doLF() {
            if (phiPoint.offset)
                --phiPoint.offset;
            else {
                phiPoint.interval = (phiPoint.interval == 0)? PL->IntLen.size() - 1 : phiPoint.interval - 1;
                phiPoint.offset = PL->IntLen[phiPoint.interval] - 1;
            }

            LFpoint = LF->map(LFpoint);
        }

        //takes as input top of run to initialize to, defaults to top of BWT
        LFPhiCoordinate(MoveStructure* lf, InvertibleMoveStructure* pl, sdsl::int_vector<>* satoprunint, uint64_t run = 0): LF(lf), PL(pl), SATopRunInt(satoprunint) {
            setToTop(run);
        }

        void setToTop(uint64_t run = 0) {
            LFpoint.position = phiPoint.position = -1;
            LFpoint.interval = run;
            phiPoint.interval = (*SATopRunInt)[run];
            LFpoint.offset = phiPoint.offset = 0;
        }
    };

    OptBWTRL(rb3_fmi_t* rb3){
        validateRB3(rb3);

        uRange alphRange; 
        std::vector<uint64_t> alphCounts;
        RLBWTconstruction(rb3, alphRange, alphCounts);

        rb3_fmi_free(rb3);

        std::vector<MoveStructure::IntervalPoint> alphStarts;
        LFconstruction(alphRange, alphCounts, alphStarts);

        std::vector<uint64_t> seqNumsTopOrBotRun, seqLens;
        std::vector<MoveStructure::IntervalPoint> stringStarts;
        SAconstruction(alphCounts, seqNumsTopOrBotRun, seqLens, stringStarts);

        
        //if (!verifyPhi() || !verifyInvPhi()) { exit(1); } 

        LCPconstruction(alphStarts, seqNumsTopOrBotRun, seqLens);

        endmarkerRepair(alphCounts, stringStarts);
    }

    //read OptBWTRL from file
    OptBWTRL(const char* fileName) {
        std::ifstream in(fileName);
        if (!in.is_open()) {
            std::cerr << "ERROR: File provided as input for OptBWTRL loading, '" 
                << fileName << "' failed to open!\n";
            exit(1);
        }
        
        load(in);

        in.close();
    }

    void printRaw() { 
        std::cout << "i\tSA_S\tSA_O\tLCP\tLF\tBWT\n";
        std::vector<uint64_t> runlenPrefSum(runlens.size());
        for (uint64_t i = 1; i < runlens.size(); ++i)
            runlenPrefSum[i] = runlenPrefSum[i-1] + runlens[i-1];
        
        //saOrder traversal
        LFPhiCoordinate saOrder(&LF, &PL, &SATopRunInt);
        uint64_t ind = 0;
        do {
            MoveStructure::IntervalPoint LFto = LF.map(saOrder.LFpoint);
            std::cout << ind++ << '\t' 
                << PL.SeqAt[saOrder.phiPoint.interval] << '\t'
                << PL.PosAt[saOrder.phiPoint.interval] + saOrder.phiPoint.offset << '\t'
                << PL.LCP(saOrder.phiPoint) << '\t'
                << runlenPrefSum[LFto.interval] + LFto.offset << '\t'
                << rlbwt[saOrder.LFpoint.interval] << '\n';

            saOrder.doInvPhi();
        } while (saOrder.LFpoint.interval != 0 || saOrder.LFpoint.offset != 0);

        std::cout << "\n\n";


        //text order traversal
        LFPhiCoordinate tOrder(&LF, &PL, &SATopRunInt);
        while (rlbwt[tOrder.LFpoint.interval] != 0)
            tOrder.doLF();
        char FofcurrSuff = 0;
        tOrder.doLF();
        LFPhiCoordinate end = tOrder;
        std::vector<uint64_t> text, isa, plcp, ph_s, ph_o, iph_s, iph_o;
        do {
            MoveStructure::IntervalPoint phto = PL.phi.map(tOrder.phiPoint);
            MoveStructure::IntervalPoint iphto = PL.invPhi.map(tOrder.phiPoint);

            text.push_back(FofcurrSuff);
            isa.push_back(runlenPrefSum[tOrder.LFpoint.interval] + tOrder.LFpoint.offset);
            plcp.push_back(PL.LCP(tOrder.phiPoint));
            ph_s.push_back(PL.SeqAt[phto.interval]);
            ph_o.push_back(PL.PosAt[phto.interval] + phto.offset);
            iph_s.push_back(PL.SeqAt[iphto.interval]);
            iph_o.push_back(PL.PosAt[iphto.interval] + iphto.offset);

            FofcurrSuff = rlbwt[tOrder.LFpoint.interval];
            tOrder.doLF();
        } while (tOrder.LFpoint.interval != end.LFpoint.interval || tOrder.LFpoint.offset != end.LFpoint.offset);

        std::vector<std::pair<const char*,std::vector<uint64_t>*>> outs = { 
            std::pair<const char*,std::vector<uint64_t>*>("TEXT"  , &text), 
            std::pair<const char*,std::vector<uint64_t>*>("ISA"   , &isa),
            std::pair<const char*,std::vector<uint64_t>*>("PLCP"  , &plcp),
            std::pair<const char*,std::vector<uint64_t>*>("PHI_S" , &ph_s),
            std::pair<const char*,std::vector<uint64_t>*>("PHI_O" , &ph_o),
            std::pair<const char*,std::vector<uint64_t>*>("IPHI_S", &iph_s),
            std::pair<const char*,std::vector<uint64_t>*>("IPHI_O", &iph_o)
        };
        std::cout << "i";
        for (uint64_t i = 0; i < text.size(); ++i) 
            std::cout << '\t' << i;
        std::cout << '\n';
        for (auto a : outs) {
            std::cout << a.first;
            for (auto rit = a.second->rbegin(); rit != a.second->rend(); ++rit) 
                std::cout << '\t' << *rit;
            std::cout << '\n';
        }
    }

    static bool validateRB3(const rb3_fmi_t* rb3);

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));

        size_type bytes = 0;

        bytes += sdsl::serialize(totalLen, out, child, "totalLen");
        bytes += sdsl::serialize(rlbwt, out, child, "rlbwt");
        bytes += sdsl::serialize(runlens, out, child, "runlens");
        bytes += sdsl::serialize(SATopRunInt, out, child, "SATopRunInt");
        bytes += sdsl::serialize(LF, out, child, "LF");
        bytes += sdsl::serialize(PL, out, child, "PL");

        sdsl::structure_tree::add_size(child, bytes);

        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(totalLen, in);
        sdsl::load(rlbwt, in);
        sdsl::load(runlens, in);
        sdsl::load(SATopRunInt, in);
        sdsl::load(LF, in);
        LF.intLens = &runlens;
        sdsl::load(PL, in);
        PL.phi.intLens = PL.invPhi.intLens = &PL.IntLen;
    }

    //returns whether the RLBWT is equivalent to an fmi
    bool equalToFmi(const rb3_fmi_t& fmi) {
        if (rlbwt.size() != runlens.size()) {
            std::cerr << "areEqual called for a symbol array and length array of different sizes" << std::endl;
            exit(1);
        }

        validateRB3(&fmi);

        rlditr_t itr;
        rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
        uint64_t currentRun = 0;
        int64_t l = 0;
        int c = 0;
        while (currentRun < rlbwt.size() && (l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
            uint64_t thisRunLength = runlens[currentRun], thisRunChar = rlbwt[currentRun];
            //reconcatenating endmarker runs for comparison
            while (thisRunChar == 0 && currentRun + 1 < rlbwt.size() && rlbwt[currentRun + 1] == 0) {
                ++currentRun;
                thisRunLength += runlens[currentRun];
            }

            if (static_cast<uint64_t>(l) != thisRunLength || static_cast<uint64_t>(c) != thisRunChar)
                return false;
            ++currentRun;
        }
        return currentRun == rlbwt.size();
    }

    bool verifyPhi() {
        Timer.start("Verifying Phi");

        uint64_t runs = rlbwt.size();
        bool pass = true;
        #pragma omp parallel for schedule(guided)
        for (uint64_t run = 0; run < runs; ++run) {
            //check if top of run phis to top of previous run
            MoveStructure::IntervalPoint start{uint64_t(-1), SATopRunInt[(run == SATopRunInt.size()-1)? 0 : (run + 1)], 0},
                          end{uint64_t(-1), SATopRunInt[run], 0};
            for(uint64_t ops = 0; ops < runlens[run]; ++ops)
                start = PL.phi.map(start);
            if (start != end) {
                #pragma omp critical
                {
                    pass = false;
                    std::cerr << "ERROR: suffix did not end up at the beginning of run in verification of phi when phiing runlen times from the top of next run!\n";
                }
            }
        }
        if (pass)
            std::cout << "Phi move data structure is likely correct, it is a permutation with one cycle\n";

        Timer.stop(); //Verifying Phi
        return pass;
    }

    bool verifyInvPhi() {
        Timer.start("Verifying InvPhi");

        uint64_t runs = rlbwt.size();
        bool pass = true;
        #pragma omp parallel for schedule(guided)
        for (uint64_t run = 0; run < runs; ++run) {
            //check if top of run phis to top of next run
            MoveStructure::IntervalPoint start{uint64_t(-1), SATopRunInt[run], 0},
                          end{uint64_t(-1), SATopRunInt[(run == SATopRunInt.size()-1)? 0 : (run + 1)], 0};
            for(uint64_t ops = 0; ops < runlens[run]; ++ops)
                start = PL.invPhi.map(start);
            if (start != end) {
                #pragma omp critical
                {
                    pass = false;
                    std::cerr << "ERROR: suffix did not end up at the beginning of next run in verification of inverse phi when invPhiing runlen times from the top of a run!\n";
                }
            }
        }
        if (pass)
            std::cout << "InvPhi move data structure is likely correct, it is a permutation with one cycle\n";

        Timer.stop(); //Verifying InvPhi
        return pass;
    }

    bool validateAllExceptRLBWT() {
        return verifyPhi() && verifyInvPhi();
    }

    bool validateAll(const rb3_fmi_t& fmi) {
        return equalToFmi(fmi) && validateAllExceptRLBWT();
    }

    //stats-------------------------------------
    uint64_t sumLCPTopRun() {
        uint64_t numIntervals = SATopRunInt.size();
        uint64_t sum = 0;
        for (uint64_t i = 0; i < numIntervals; ++i) {
            sum += PL.AboveLCP[SATopRunInt[i]];
        }
        return sum;
    }

    struct Histograms {
        std::vector<uint64_t> LF;
        std::vector<uint64_t> Phi;
        std::vector<uint64_t> InvPhi;
        struct InvertibleMoveStructure {
            std::vector<uint64_t> Phi;
            std::vector<uint64_t> InvPhi;
        } INVMOVE;
    };

    Histograms IntervalTraversals() {
        Histograms a;
        
        LF.MakeHistogram(a.LF);
        PL.phi.MakeHistogram(a.INVMOVE.Phi);
        PL.invPhi.MakeHistogram(a.INVMOVE.InvPhi);

        /*
        uint64_t runs = rlbwt.size();
        std::cout << "finished" << std::endl;

        for (uint64_t i = 0; i < runs; ++i) {
            LFPhiCoordinate b(&LF, &PL, &SATopRunInt, i);


            


        }
        */
        return a;
    }

    //extract-----------------------------------

    //extracts up to len characters of sequence seq, starting at position pos
    //i.e. returns S_seq[pos, min(pos+len, |S_seq|) - 1]
    //where S_seq is the seq-th string in the text (0-indexed)
    std::string extract(uint64_t seq, uint64_t pos = 0, uint64_t len = static_cast<uint64_t>(-1)) {
        MoveStructure::IntervalPoint lfPoint{0, 0, 0};
        uint64_t compSeq = seq ^ 0x1;
        LF.AdvanceIntervalPoint_unsafe(lfPoint, compSeq);

        for(uint64_t i = 0; rlbwt[lfPoint.interval] && i < pos; ++i)
            lfPoint = LF.map(lfPoint);

        std::string res, convertComp = "$TGCAN";
        while (res.size() < len && rlbwt[lfPoint.interval] != 0) {
            res.push_back(convertComp[rlbwt[lfPoint.interval]]);
            lfPoint = LF.map(lfPoint);
        }
        return res;
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
            
            auto check = [&out,this,lengthThreshold] (MoveStructure::IntervalPoint plPoint) {
                MoveStructure::IntervalPoint plPointPrev{ static_cast<uint64_t>(-1), ((plPoint.interval)? plPoint.interval-1: PL.IntLen.size() - 1), 0 };
                plPointPrev.offset = PL.IntLen[plPointPrev.interval] - 1;

                uint64_t len, lenLF;

                len = std::max(PL.LCP(plPoint), PL.LCP(PL.invPhi.map(plPoint)));
                lenLF = std::max(PL.LCP(plPointPrev), PL.LCP(PL.invPhi.map(plPointPrev)));
                if (len >= lenLF && len >= lengthThreshold) {
                    //out << "LCP " << PL.LCP(plPoint) << " LCPBelow " << PL.LCP(PL.invPhi.map(plPoint)) << '\n';
                    //out << "LCPPrev " << PL.LCP(plPointPrev) << " LCPBelowPrev " << PL.LCP(PL.invPhi.map(plPointPrev)) << '\n';
                    //out << len << '\t' << lenLF << '\n';
                    //supermaximal repeat
                    //else, len = lenLF - 1
                    out << PL.SeqAt[plPoint.interval] << '\t'
                        << PL.PosAt[plPoint.interval] + plPoint.offset << '\t'
                        << len << '\t';

                    std::stack<std::pair<uint64_t,uint64_t>> occ;

                    //traverse up
                    MoveStructure::IntervalPoint plp = plPoint;
                    while (PL.LCP(plp) >= len) {
                        plp = PL.phi.map(plp);
                        occ.emplace(PL.SeqAt[plp.interval], PL.PosAt[plp.interval] + plp.offset);
                    }

                    //traverse down
                    plp = PL.invPhi.map(plPoint);
                    while(PL.LCP(plp) >= len) {
                        occ.emplace(PL.SeqAt[plp.interval], PL.PosAt[plp.interval] + plp.offset);
                        plp = PL.invPhi.map(plp);
                    }

                    //assert(occ.size());
                    out << 1+occ.size();
                    while (occ.size()) {
                        out << '\t' << occ.top().first
                            << '\t' << occ.top().second;
                        occ.pop();
                    }
                    out << '\n';
                    //out.flush();
                }
            };

            MoveStructure::IntervalPoint topPlPoint{ static_cast<uint64_t>(-1), SATopRunInt[run], 0};
            //std::cout << "run " << run
                //<< " topPlPoint.interval " << topPlPoint.interval
                //<< " topPlPoint.offset " << topPlPoint.offset 
                //<< std::endl;
            //if the run has length 1, its top is also its bottom so the suffix at the top will be 
            //evaluated in the next iteration of the loop
            if (runlens[run] != 1)  
                check(topPlPoint);

            MoveStructure::IntervalPoint botPlPoint = PL.phi.map(topPlPoint);
            check(botPlPoint);
        }
    }

    //WARNING, THIS IMPLEMENTATION ASSUMES NO RUN SPLITTING, MAY BE BUGGY IF RUN SPLITTING IS PERFORMED
    //no default value for lengthThreshold because there will be many of length >= 1
    //occ[c]^2 per character c?
    void repeats(std::ostream& out, const uint64_t lengthThreshold) {
        //a repeat is a match T[i,i+l) = T[j, j+l) s.t.
        //  a. T[i-1] != T[j-1]
        //  b. T[i+l] != T[j+l]
        //T[i,i+l) = T[j,j+l) is a repeat iff
        //  a. There is a run boundary between ISA[i] and ISA[j]
        //  b. BWT[ISA[i]] != BWT[ISA[j]]
        //  c. l = min(LCP[k]) for k = min(ISA[i],ISA[j])+1 to max(ISA[i],ISA[j])
        //
        //This function outputs all repeats in the text of length at least lengthThreshold
        //It runs in O(r + occ) time
        
        out << "seq1\tpos1\tseq2\tpos2\tlen\n";

        uint64_t runs = rlbwt.size();
        //minLCPRun is a stack of runs above run where the min LCP in each run
        //is at least lengthThreshold. minLCPRun.back() stores the min LCP value
        //in run run-1, minLCPRun[minLCPRun.size()-2] stores the min LCP value
        //in run run-2, and so on. minLCP[i] >= lengthThreshold for all i
        //run x is only in minLCPRun if run x+1 is also in it or x=run-1 and
        //the LCP value at the top of run run is at least length threshold
        std::vector<uint64_t> minLCPRun;
        for (uint64_t run = 0; run < runs; ++run) {
            LFPhiCoordinate topRun(&LF, &PL, &SATopRunInt, run);
            LFPhiCoordinate coord = topRun;
            
            uint64_t l, minLCP = static_cast<uint64_t>(-1), ch = rlbwt[run];
            
            //returns lcp of SA[coord] and first position above it
            //where BWT != ch 
            //assumes minLCPRun valid
            //if lcp < length threshold, may return any value < length threshold (not necessarily the correct one)
            auto nextDiffLCP = [lengthThreshold, run, minLCPRun, this](LFPhiCoordinate& coord, const uint64_t ch) -> uint64_t {
//                std::cout << "coord.lfpoint " 
//                    << coord.LFpoint.interval << ',' << coord.LFpoint.offset
//                    << " coord.phiPoint " 
//                    << coord.phiPoint.interval << ',' << coord.phiPoint.offset 
//                    << std::endl;
//                std::cout << "coord.Seq " << PL.SeqAt[coord.phiPoint.interval]
//                    << " coord.Pos " << PL.PosAt[coord.phiPoint.interval] + coord.phiPoint.offset
//                    << std::endl;
                uint64_t l = coord.LCP();
                // std::cerr << "l " << l << std::endl;
                if (coord.LFpoint.offset || l < lengthThreshold || ch == 0 || rlbwt[coord.LFpoint.interval - 1] != ch)
                    return l;
//                std::cout << "skipping run of equal to ch " << ch << std::endl;
                //else skipping run
                coord.setToTop(coord.LFpoint.interval - 1);
                //if minLCP of run < length threshold, return 0
                if (run - coord.LFpoint.interval > minLCPRun.size())
                    return 0;
//                std::cout << "run has large enough minLCP: " << minLCPRun[minLCPRun.size() - (run - coord.LFpoint.interval)] << std::endl;
                //else
		l = std::min(l, PL.LCP(coord.phiPoint));
                return std::min(l, minLCPRun[minLCPRun.size() - (run - coord.LFpoint.interval)]);
            };

            //sequence, position, minLCP to top of run
            std::vector<std::tuple<uint64_t,uint64_t,uint64_t>> minLCPSuff;
//            std::cout << "Going up run " << run << std::endl;
            //go up
            while ((l = nextDiffLCP(coord,ch)) >= lengthThreshold) {
                minLCP = std::min(l, minLCP);
//                std::cout << "computed nextDiffLCP: " << l << std::endl;
//                std::cout << "new minLCP: " << minLCP << std::endl;
                coord.doPhi();
                // TODO: THINK ABOUT SPACE COMPLEXITY OF THIS STACK
                minLCPSuff.emplace_back(
                        PL.SeqAt[coord.phiPoint.interval],
                        PL.PosAt[coord.phiPoint.interval]+coord.phiPoint.offset,
                        minLCP);
            }
//            std::cout << "finished going up" << std::endl;
            
            //go down current run
            uint64_t rlen = runlens[run], currLCP = static_cast<uint64_t>(-1);
            coord = topRun;
            for (uint64_t i = 0; currLCP >= lengthThreshold && i < rlen; ++i) {
                uint64_t seq = PL.SeqAt[coord.phiPoint.interval];
                uint64_t pos = PL.PosAt[coord.phiPoint.interval] + coord.phiPoint.offset;
                for (const auto& a : minLCPSuff) {
                    if (seq < std::get<0>(a) || (seq == std::get<0>(a) && pos < std::get<1>(a)))
                        out << seq << '\t'
                            << pos << '\t'
                            << std::get<0>(a) << '\t'
                            << std::get<1>(a) << '\t'
                            << std::min(std::get<2>(a), currLCP) << '\n';
                    else 
                        out << std::get<0>(a) << '\t'
                            << std::get<1>(a) << '\t'
                            << seq << '\t'
                            << pos << '\t'
                            << std::min(std::get<2>(a), currLCP) << '\n';
                }
                // if (i != rlen - 1) {
                    coord.doInvPhi();
                    currLCP = std::min(currLCP, coord.LCP());
                //}
            }

            if (currLCP < lengthThreshold)
                minLCPRun.clear();
            else 
                minLCPRun.push_back(currLCP);
        }
    }
};

bool OptBWTRL::validateRB3(const rb3_fmi_t* rb3){
    if (!rb3->e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return false;
    }

    if (!rb3_fmi_is_symmetric(rb3)) {
        std::cerr << "ERROR: fmd is not 'symmetric' (does this mean bidirectional?) fmd must be bidirectional" << std::endl;
        return false;
    }
    return true;
}

