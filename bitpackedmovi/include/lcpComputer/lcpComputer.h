#ifndef R_SA_LCP_LCPCOMPUTER_H
#define R_SA_LCP_LCPCOMPUTER_H

#include"util/util.h"
#include"moveStructure/moveStructure.h"
#include"fm-index.h"
#include<queue>
#include<omp.h>
#include<atomic>
#include<mutex>

static constexpr const char* lcp_index_extension = ".lcp_index";

class LCPComputer {
    uint64_t totalLen;

    sdsl::int_vector<> F;
    
    //Flens is now part of Psi
    //sdsl::int_vector<> Flens;
    MoveStructureTable Psi;

    sdsl::int_vector<> intAtTop;

    //sdsl::int_vector<> PhiIntLen;
    MoveStructureStartTable Phi;

    sdsl::int_vector<> PLCPsamples;

    void ConstructPsi(const rb3_fmi_t* rb3, sdsl::int_vector<> & F, MoveStructureTable &Psi, uint64_t & numSequences
            #ifndef BENCHFASTONLY
            , const verbosity v
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Constructing Psi from FMD"); }
        #endif

        uint64_t runs = 0, alphbits, lenbits;
        totalLen = 0;
        uRange lenRange;
        std::vector<uint64_t> alphRuns;

        uRange RB3_lenRange;
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Reading fmd for parameters"); }
        #endif
        {
            //original ropebwt3 values
            uint64_t RB3_runs = 0; 
            #ifndef BENCHFASTONLY
            uint64_t RB3_lenbits;
            #endif

            uRange alphRange;



            rlditr_t itr1;
            rld_itr_init(rb3->e, &itr1, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
            int64_t l;
            int c = 0;
            //int prevc;


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
            //prevc = c;

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

                //if (c == prevc) {
                //    std::cerr << "ERROR: Two runs of the same character follow each other from ropebwt3!" << std::endl;
                //    exit(1);
                //}
                //prevc = c;
            }

            numSequences = alphRuns[0];

            if (alphRange.max == static_cast<uint64_t>(-1)) {
                std::cerr << "Maximum alphabet symbol is 2^64 - 1. "
                    << "This program assumes this is not the case (it can only handle alphabet <= (2^64) - 2." << std::endl;
                exit(1);
            }

            #ifndef BENCHFASTONLY
            if (v >= VERB)
                std::cout << "INFO: The parameters for our constructed BWT (i.e. #runs, max length, etc.) may be "
                    << "different from those of the input (ropebwt3).\nINFO: This is because in our constructed BWT, "
                    << "each endmarker is contained in its own run.\n";
            #endif

            alphbits = sdsl::bits::hi(alphRange.max) + 1;
            lenbits = sdsl::bits::hi(lenRange.max) + 1;
            if (alphbits != static_cast<uint64_t>(rb3->e->abits)) 
                std::cout << "WARNING: computed bits per symbol not equal to bits used in fmd. Computed: " 
                    << alphbits << ", ropebwt3: " << static_cast<uint64_t>(rb3->e->abits) << std::endl;

            #ifndef BENCHFASTONLY
            RB3_lenbits = sdsl::bits::hi(RB3_lenRange.max) + 1;
            if (v >= VERB) {
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
            #endif
        }
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Reading fmd for parameters
        #endif

        Psi.data = packedTripleVector(sdsl::bits::hi(runs - 1) + 1, sdsl::bits::hi(lenRange.max - 1) + 1, lenbits, runs);
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Reading fmd to create Flens"); }
        #endif
        {

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
                    for (uint64_t i = 0; i < static_cast<uint64_t>(l); ++i)
                        Psi.data.set<2>(alphFRunStarts[c]++, 1);
                }
                else {
                    Psi.data.set<2>(alphFRunStarts[c]++, l);
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
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Reading fmd to create Flens

        //NOTE: Psi of endmarker runs will be incorrect, will fix in a later step
        if (v >= TIME) { Timer.start("Reading fmd to create D_index and D_offset"); }
        #endif
        {

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
                        Psi.data.set<0>(alphFRunStarts[c], currentRun);
                        Psi.data.set<1>(alphFRunStarts[c], currentOffset);

                        ++alphFRunStarts[c];

                        ++currentOffset;
                        currentOffset %= Psi.data.get<2>(currentRun);
                        currentRun += (currentOffset == 0);
                    }
                }
                else {
                    Psi.data.set<0>(alphFRunStarts[c], currentRun);
                    Psi.data.set<1>(alphFRunStarts[c], currentOffset);

                    ++alphFRunStarts[c];

                    currentOffset += l;
                    while (currentOffset && currentOffset >= Psi.data.get<2>(currentRun))
                        currentOffset -= Psi.data.get<2>(currentRun++);
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
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Reading fmd to create D_index and D_offset

        //obviously, we could just replace F with a bit vector of length r with sigma set bits
        //and use rank on the bitvector
        //for the same time complexity but r bits instead of r log sigma
        //I'm not sure how much slower that is. I haven't tried it yet.
        if (v >= TIME) { Timer.start("Constructing F"); }
        #endif
        {
            F = sdsl::int_vector<>(runs, 7, alphbits);
            uint64_t curr = 0;
            for (uint64_t alph = 0; alph < alphRuns.size(); ++alph)
                for (uint64_t i = 0; i < alphRuns[alph]; ++i)
                    F[curr++] = alph;
        }
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Constructing F

        if (v >= TIME) { Timer.stop(); } //Constructing Psi from FMD }
        #endif
    }

    //O(n) time, not needed for LCP computation
    //LCP computation is O(n) anyways and we can repair while we do it
    //although we don't use psi after so may as well not repair it
    //this function fixes the psi mappings for the endmarker runs (of F)
    //It also computes numTopRuns and seqLens
    //numTopRuns and seqLens are vectors of length equal to the number of sequences
    //the i+1-th value of numTopRuns is the number of times suffixes of sequence i
    //are at the (SA position corresponding to the) top of a run in the BWT. 
    //(includes the termination symbol at the end of sequence i, $_i)
    //the i+1-th value of seqLens is the length of seq i,
    //(includes the termination symbol) 
    //Actually, numTopRuns and seqLens are the prefix sums of the above definitions
    //
    //
    //for every suffix x at the top of a run in the BWT, 
    //there is an input interval of psi, j, where 
    //suffix x-1 is at the top of the input interval
    //intAtTop stores, for every input interval of psi
    //with suffix x-1 at the top of the input interval,
    //how many suffixes < x are at the top of a run in the BWT
    void ComputeAuxAndRepairPsi(uint64_t& maxPhiIntLen, std::vector<uint64_t> & numTopRuns, std::vector<uint64_t> & seqLens, sdsl::int_vector<> & intAtTop, 
            const sdsl::int_vector<>& F, MoveStructureTable& Psi, sdsl::int_vector<>& PhiIntLen, const uint64_t numSequences
            #ifndef BENCHFASTONLY
            , const verbosity v
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Computing numTopRuns, seqLens, and repairing Psi of endmarkers in F"); } 
        #endif

        numTopRuns.resize(numSequences + 1);
        seqLens = std::vector<uint64_t>(numSequences + 1);

        intAtTop = sdsl::int_vector<>(F.size(), -1, sdsl::bits::hi(F.size() - 1) + 1);

        //computing seqLens, numTopRuns, correctSeqPsis, and maxPhiIntLen
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Parallel seq traversal"); }
        #endif
        {
            std::vector<uint64_t> maxPhiIntLenPerSeq(numSequences);
            std::vector<MoveStructureTable::IntervalPoint> correctSeqPsis(numSequences);
            #pragma omp parallel for schedule(dynamic, 1)
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                MoveStructureTable::IntervalPoint start = {static_cast<uint64_t>(-1), seq, 0}, curr;
                start = Psi.map(start);
                curr = start;

                //suffix 0 of seq seq is at the top of a run in the bwt since its bwt value is a termination string
                uint64_t numTop = 1;
                uint64_t len = 1;
                uint64_t maxPhiIntLenThisSeq = 0;
                uint64_t currPhiIntLen = 1;

                while (curr.interval >= numSequences) {
                    //suffix x is at the top of a run in the bwt iff suffix x - 1 is at the start
                    //of an input interval of psi
                    numTop += (curr.offset == 0);
                    if (curr.offset == 0) {
                        maxPhiIntLenThisSeq = std::max(maxPhiIntLenThisSeq, currPhiIntLen);
                        currPhiIntLen = 0;
                    }
                    ++len;
                    ++currPhiIntLen;
                    curr = Psi.map(curr);
                }

                //for last interval in sequence, if condition is unnecessary, guaranteed to be true.
                //if (curr.offset == 0) {
                    maxPhiIntLenThisSeq = std::max(maxPhiIntLenThisSeq, currPhiIntLen);
                    currPhiIntLen = 0;
                //}

                //if (curr.offset) {
                    //std::cerr << "ERROR: Run of endmarkers in F of length more than 1!" << std::endl;
                    //exit(1);
                //}

                uint64_t seqStartingAtStart = curr.interval;

                //should not need omp critical
                correctSeqPsis[(seqStartingAtStart)? seqStartingAtStart - 1 : numSequences - 1] = start;
                seqLens[seqStartingAtStart + 1] = len;
                numTopRuns[seqStartingAtStart + 1] = numTop;
                maxPhiIntLenPerSeq[seqStartingAtStart] = maxPhiIntLenThisSeq;
            }

            maxPhiIntLen = 0;
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                Psi.data.set<0>(seq, correctSeqPsis[seq].interval);
                Psi.data.set<1>(seq, correctSeqPsis[seq].offset);
                maxPhiIntLen = std::max(maxPhiIntLen, maxPhiIntLenPerSeq[seq]);
            }
        }
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Parallel seq traversal 
        #endif

        //from here on, numTopRuns and seqLens are the exclusive prefix sums of their previous definition
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Prefix summing auxiliary data"); }
        #endif
        for (uint64_t i = 1; i < seqLens.size(); ++i) {
            seqLens[i] += seqLens[i-1];
            numTopRuns[i] += numTopRuns[i-1];
        }

        if (numTopRuns.back() != F.size()) {
            std::cerr << "ERROR: Number of runs in numTopRuns doesn't sum to total number of runs!" << std::endl;
            exit(1);
        }
        if (seqLens.back() != totalLen) {
            std::cerr << "ERROR: Lengths in seqLens doesn't sum to total length!" << std::endl;
            exit(1);
        }
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Prefix summing auxiliary data
        #endif

        #ifndef BENCHFASTONLY
        if (v >= VERB) { std::cout << "INFO: maximum interval length for phi data structure: " << maxPhiIntLen << std::endl; } 
        #endif

        PhiIntLen = sdsl::int_vector<>(F.size(), 0, sdsl::bits::hi(maxPhiIntLen) + 1);
        //computing intAtTop and PhiIntLen
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Second Parallel seq traversal"); }
        #endif
        {
            uint64_t dangerousInts = 64/std::min(PhiIntLen.width(), intAtTop.width()) + (64 % std::min(PhiIntLen.width(), intAtTop.width()) != 0);
            #pragma omp parallel for schedule(dynamic, 1)
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                uint64_t prevSeq = (seq)? seq - 1 : numSequences - 1;
                MoveStructureTable::IntervalPoint curr = {static_cast<uint64_t>(-1), prevSeq, 0};
                curr = Psi.map(curr);

                uint64_t currentInt = numTopRuns[seq];
                const uint64_t start = numTopRuns[seq];
                const uint64_t end = numTopRuns[seq + 1];
                const uint64_t safeStart = start + dangerousInts,
                         safeEnd = numTopRuns[seq + 1] - dangerousInts;
                uint64_t currIntLen = 1;
                sdsl::int_vector<> intAtTopIndex(end - start, 0, Psi.data.a);
                sdsl::int_vector<> intAtTopValue(end - start, 0, Psi.data.a);
                do {
                    if (curr.offset == 0) {
                        /*
                        #pragma omp critical
                        {
                            PhiIntLen[currentInt] = currIntLen;
                            currIntLen = 0;
                            intAtTop[curr.interval] = currentInt++;
                        }
                        */
                        if (currentInt >= safeStart && currentInt < safeEnd)
                            PhiIntLen[currentInt] = currIntLen;
                        else {
                            #pragma omp critical(philenwriting)
                            {
                                PhiIntLen[currentInt] = currIntLen;
                            }
                        }
                        intAtTopIndex[currentInt - start] = curr.interval;
                        intAtTopValue[currentInt - start] = currentInt;
                        currIntLen = 0;
                        currentInt++;
                    }
                    curr = Psi.map(curr);
                    ++currIntLen;
                    //std::cout << "curr.interval " << curr.interval << " curr.offset " << curr.offset << std::endl;
                } while (curr.interval >= numSequences);

                assert(curr.offset == 0);
                //if (curr.offset == 0) {
                #pragma omp critical(philenwriting)
                {
                    PhiIntLen[currentInt] = currIntLen;
                }
                intAtTopIndex[currentInt - start] = curr.interval;
                intAtTopValue[currentInt - start] = currentInt;
                currIntLen = 0;
                currentInt++;
                /*
                #pragma omp critical
                {
                    PhiIntLen[currentInt] = currIntLen;
                    currIntLen = 0;
                    intAtTop[curr.interval] = currentInt++;
                }
                */
                //}
                #pragma omp critical(intattopwriting)
                {
                    for (uint64_t i = 0; i < intAtTopIndex.size(); ++i) {
                        intAtTop[intAtTopIndex[i]] = intAtTopValue[i];
                    }
                }


                assert(currentInt == numTopRuns[seq + 1]);
                //if (currentInt != numTopRuns[seq + 1]) {
                //std::cerr << "ERROR: Didn't reach beginning of next sequence in numTopRuns!" << std::endl;
                //std::cout << currentInt << " " << numTopRuns[seq + 1] << std::endl;
                //for (uint64_t i = 0; i < F.size(); ++i) {
                //std::cout << intAtTop[i] << std::endl;
                //}
                //exit(1);
                //}
            }
        }
        /*
        uint64_t sum = 0;
        for (uint64_t i = 0; i < PhiIntLen.size(); ++i) {
            std::cout << "PhiIntLen[" << i << "]: " << PhiIntLen[i] << ", sum: " << sum << std::endl;
            sum += PhiIntLen[i];
        }
        std::cout << "total sum: " << sum << std::endl;
        */
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Second Parallel seq traversal }
        #endif
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Computing numTopRuns, seqLens, and repairing Psi of endmarkers in F
        #endif
    }

    void ConstructPhiAndSamples(const MoveStructureTable& Psi, sdsl::int_vector<>& PhiIntLen, const uint64_t FlensBits,
            const std::vector<uint64_t>& numTopRuns, const std::vector<uint64_t>& seqLens, const sdsl::int_vector<>& intAtTop, const uint64_t numSequences, const uint64_t maxPhiIntLen,
            uint64_t & sampleInterval, sdsl::int_vector<> &Psi_Index_Samples, sdsl::int_vector<> &Psi_Offset_Samples
            #ifndef BENCHFASTONLY
            , const verbosity v
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("Construct Phi and Samples"); }
        #endif
        //MoveStructure tempPhi;
        //tempPhi.intLens = &PhiIntLen;
        //computing, for each seq i, 
        // 1. whenever suffix j of i is at the bottom of a run, the input interval and offset of suffix j of i in Phi (D_index and D_offset of the top of the run below it)
        // 2. Psi_Index_Samples and Psi_Offset_Samples: ISA samples (in psi move data structure), for every suffix that is a multiple of n/r (of the original text)
        uint64_t numRuns = Psi.data.size();
        sampleInterval = totalLen/numRuns;
        uint64_t numSamples = (totalLen%sampleInterval != 0) + (totalLen/sampleInterval);
        //std::cout << "sampleInterval " << sampleInterval << std::endl;
        //std::cout << "numSamples " << numSamples << std::endl;
        //tempPhi.D_index = sdsl::int_vector<>(numRuns, 0, sdsl::bits::hi(numRuns - 1) + 1);
        //tempPhi.D_offset = sdsl::int_vector<>(numRuns, 0, sdsl::bits::hi(maxPhiIntLen - 1) + 1);
        {
            uint64_t sum = 0;
            for (uint64_t i = 0; i < PhiIntLen.size(); ++i)
                sum += PhiIntLen[i];
            Phi.data = packedTripleVector(sdsl::bits::hi(numRuns - 1) + 1, sdsl::bits::hi(maxPhiIntLen - 1) + 1, sdsl::bits::hi(sum) + 1, PhiIntLen.size() + 1);
            //std::cout << "Total sum: "  << sum << std::endl;
            sum = 0;
            for (uint64_t i = 0; i < PhiIntLen.size(); ++i){
                Phi.data.set<2>(i, sum);
                //std::cout << "current sum: " << sum << std::endl;
                sum += PhiIntLen[i];
            }
            Phi.data.set<2>(PhiIntLen.size(), sum);
            PhiIntLen = sdsl::int_vector<>();
            //std::cout << "Total sum: "  << sum << std::endl;
        }
        Psi_Index_Samples = sdsl::int_vector<>(numSamples, 0, sdsl::bits::hi(numRuns - 1) + 1);
        Psi_Offset_Samples = sdsl::int_vector<>(numSamples, 0, FlensBits);
        #ifndef BENCHFASTONLY
        if (v >= VERB) {
            std::cout << "FlensBits: " << FlensBits << std::endl;
            std::cout << "Psi_Offset_Samples.width(): " << static_cast<uint64_t>(Psi_Offset_Samples.width()) << std::endl;
        }
        #endif
        {
            /*
               Fixing the thrashing of the Phi computation method without 
               increasing the memory usage is not so simple. We have a few 
               options:

                    1. Separate Psi_Index_Samples and Psi_Index_Offsets 
                    computation into another O(n) work traversal after 
                    the current one. Use them to store PhiInfo temporarily,
                    then remap to Phi Order in O(r) time afterwards. This is 
                    worth it when thrashing would more than double the runtime 
                    of the following loop. I.E. when #threads is high. How to
                    tell when this is the case? Too complicated.

                    2. Store the D_index and D_offset values in Phi temporarily
                    here in output order. Later reorder them in O(r) time to 
                    input order using intAtEnd and intAtTop. This requires 
                    temporarily writing r log r bits to a file and reading
                    them later, since we need an extra r log r bits for 
                    intAtEnd. This can come from Psi_Index_Samples since it's not needed
                    for the reordering step.

                    3. Keep extra arrays for every thread, storing the temporary
                    results per sequence, writing in one step at the end of the
                    loop.
                    
                2. Is better than 1. when reading/writing + O(r) reordering
                is faster than a parallel O(n) traversal. We assume this to
                be the case. 
                
                We pick 3. it may not reduce the thrashing to 0 but it 
                should do enough. 1. Adds an O(n) traversal. This adds
                20-25% run time when there is no thrashing (since we do
                4 O(n) traversals for construction + 1 O(n) for minLCP).
                2. is too complicated to implement.
             */
            //the number of bits the buffer for 3. can use is at most r log max lcp len
            //max lcp len is not known ahead of time, but it is at least maxPhiIntLen - 1
            //if Ns match. In our case, they don't so we use at most r bits for the buffer.
            uint64_t threads = omp_get_max_threads();
            packedTripleVector buffer;
            uint64_t bufferElementsPerThread;
            {
                //uses 64 bytes at least
                const uint64_t bufferMaxBits = std::max(numRuns, static_cast<uint64_t>(512));
                const uint64_t bitsPerElement = (2*Phi.data.a + Phi.data.b);
                uint64_t bufferElements = (bufferMaxBits + bitsPerElement - 1)/bitsPerElement;
                //each thread should use at least 64 bytes:
                if (bufferElements*bitsPerElement < threads*512) {
                    //!!!!!!!!!!!!!!!!!!! should be rare, #threads is huge or r is very small
                    //8 is an arbitrary constant I chose so that the buffer is not useless
                    if (v >= TIME) {
                        std::cout << "WARNING: limiting number of threads in Phi computation to " << bufferElements/8
                            << " to save memory. Previously, would have been " << threads << std::endl;
                    }
                    threads = (bufferElements*bitsPerElement)/512;
                }
                bufferElementsPerThread = bufferElements/threads;
                bufferElements = bufferElementsPerThread*threads;
                buffer = packedTripleVector(Phi.data.a, Phi.data.a, Phi.data.b, bufferElements);
            }
            uint64_t dangerousInts = 64/std::min(Psi_Index_Samples.width(), Psi_Offset_Samples.width()) + (64 % std::min(Psi_Index_Samples.width(), Psi_Offset_Samples.width()) != 0);
            uint64_t dangerousBufferInts;
            if (buffer.width <= 64)
                dangerousBufferInts = 64/buffer.width + (64 % buffer.width != 0);
            else
                dangerousBufferInts = 1;
            #pragma omp parallel for num_threads(threads) schedule (dynamic, 1)
            for (uint64_t seq = 0; seq < numSequences; ++seq) {
                //curr is the interval point in the psi move data structure of suffix suff
                //if curr is at the top of a psi interval, then suff+1 is at the top of an rlbwt interval
                //if curr is at the bottom of a psi interval, then suff+1 is at the bottom of an rlbwt interval
                MoveStructureTable::IntervalPoint curr = {static_cast<uint64_t>(-1), (seq)? seq - 1 : numSequences - 1, 0};
                curr = Psi.map(curr);
                uint64_t suff = seqLens[seq];
                //phiPoint is the interval point in the move structure of suff
                MoveStructureStartTable::IntervalPoint phiPoint = {Phi.data.get<2>(numTopRuns[seq]), numTopRuns[seq], 0};
                MoveStructureStartTable::IntervalPoint phiPointAtPhiOutputIntervalStart = phiPoint;

                const uint64_t suffSampleSafeStart = (seqLens[seq]/sampleInterval) + (seqLens[seq] % sampleInterval != 0) + dangerousInts;
                const uint64_t suffSampleSafeEnd = (seqLens[seq + 1]/sampleInterval) + (seqLens[seq + 1] % sampleInterval != 0) - dangerousInts;

                const uint64_t bufferStart = omp_get_thread_num() * bufferElementsPerThread;
                const uint64_t bufferEnd = bufferStart + bufferElementsPerThread;
                const uint64_t bufferSafeStart = bufferStart + dangerousBufferInts;
                const uint64_t bufferSafeEnd = bufferEnd - dangerousBufferInts;
                uint64_t bufferInd = bufferStart;


                while (suff < seqLens[seq+1]) {
                    //2.
                    //store psi samples if needed
                    if (suff % sampleInterval == 0) {
                        if (suff / sampleInterval >= suffSampleSafeStart && suff / sampleInterval < suffSampleSafeEnd) {
                            Psi_Index_Samples[suff/sampleInterval] = curr.interval;
                            Psi_Offset_Samples[suff/sampleInterval] = curr.offset;
                        }
                        else {
                            #pragma omp critical
                            {
                                Psi_Index_Samples[suff/sampleInterval] = curr.interval;
                                Psi_Offset_Samples[suff/sampleInterval] = curr.offset;
                            }
                        }
                    }

                    //update phiPoint to suff + 1
                    ++phiPoint.offset;
                    ++phiPoint.position;
                    //phiPoint.offset == (*Phi.intLens)[phiPoint.interval] iff curr.offset == 0)
                    assert((phiPoint.position == Phi.data.get<2>(phiPoint.interval + 1)) ==
                            (curr.offset == 0));
                    //suff + 1 starts an input interval iff suff ends an input interval
                    if (curr.offset == 0) {
                        if (phiPoint.position != Phi.data.get<2>(phiPoint.interval + 1)) {
                            std::cerr << "Offset != length at the end of phi interval!" << std::endl;
                            exit(1);
                        }
                        phiPoint.offset = 0;
                        ++phiPoint.interval;
                    }

                    //if suff is the end of an output interval
                    if (curr.offset == Psi.data.get<2>(curr.interval) - 1) {
                        uint64_t runBelow = (curr.interval+1)%numRuns;
                        uint64_t phiInterval = intAtTop[runBelow];
                        if (bufferInd >= bufferSafeStart && bufferInd < bufferSafeEnd) {
                            buffer.set<0>(bufferInd, phiInterval);
                            buffer.set<1>(bufferInd, phiPointAtPhiOutputIntervalStart.interval);
                            buffer.set<2>(bufferInd, phiPointAtPhiOutputIntervalStart.offset);
                        }
                        else {
                            #pragma omp critical 
                            {
                                buffer.set<0>(bufferInd, phiInterval);
                                buffer.set<1>(bufferInd, phiPointAtPhiOutputIntervalStart.interval);
                                buffer.set<2>(bufferInd, phiPointAtPhiOutputIntervalStart.offset);
                            }
                        }
                        if (++bufferInd == bufferEnd) {
                            #pragma omp critical
                            {
                                for (uint64_t i = bufferStart; i < bufferEnd; ++i) {
                                    Phi.data.set<0>(buffer.get<0>(i), buffer.get<1>(i));
                                    Phi.data.set<1>(buffer.get<0>(i), buffer.get<2>(i));
                                }
                            }
                            bufferInd = bufferStart;
                        }
                        phiPointAtPhiOutputIntervalStart = phiPoint;
                    }
                    
                    //now, curr is the position of suff + 1
                    curr = Psi.map(curr);
                    //update suff
                    ++suff;
                    /*
                    std::cout << "suff " << suff << std::endl;
                    std::cout << "curr: " << curr.interval << ' ' << curr.offset << std::endl;
                    std::cout << "intLen: " << (*Psi.intLens)[curr.interval] << std::endl;
                    //1. 
                    //when curr at the bottom of a psi input interval, suff is the start of a phi output interval
                    //then suff-1 is the end of a phi output interval, then invPhi(suff-1) is the end of a phi input interval
                    //so, update the start of the previous input interval to phiPointAtPhiOutputIntervalStart
                    if (curr.offset == (*Psi.intLens)[curr.interval] - 1)
                        std::cout << "starts output interval" << std::endl;
                    if (curr.offset == (*Psi.intLens)[curr.interval] - 1) {
                        std::cout << "in if" << std::endl;
                        uint64_t runBelow = (curr.interval+1)%numRuns;
                        uint64_t phiInterval = intAtTop[runBelow];
                        if (phiInterval == 25) {
                            std::cout << "25: " << seq << " " << suff << std::endl;
                        }
                        #pragma omp critical 
                        {
                            Phi.D_index[phiInterval] = phiPointAtPhiOutputIntervalStart.interval;
                            Phi.D_offset[phiInterval] = phiPointAtPhiOutputIntervalStart.offset;
                        }
                        phiPointAtPhiOutputIntervalStart = phiPoint;
                    }

                    bool suffAtTop = (curr.offset == 0);
                    ++suff;
                    ++phiPoint.offset;

                    if (suff == 30) {
                        std::cout << suffAtTop << std::endl;
                    }
                    //2.
                    //start new interval in phi if suff is at the top of a run
                    if (suffAtTop) {
                        #pragma omp critical
                        {
                            ++phiPoint.interval;
                            phiPoint.offset = 0;
                        }
                    }
                    */
                }
                #pragma omp critical
                {
                    for (uint64_t i = bufferStart; i < bufferInd; ++i) {
                        Phi.data.set<0>(buffer.get<0>(i), buffer.get<1>(i));
                        Phi.data.set<1>(buffer.get<0>(i), buffer.get<2>(i));
                    }
                }
                /*
                std::cout << "done while" << std::endl;
                //condition in the following if statement is redundant, should never be false.
                //1. 
                //when at the bottom of a psi input interval, suff is the start of a phi output interval
                //so, update the start of the previous input interval to phiPointAtPhiOutputIntervalStart
                if (curr.offset == (*Psi.intLens)[curr.interval] - 1) {
                    uint64_t runBelow = (curr.interval+1)%numRuns;
                    uint64_t phiInterval = intAtTop[runBelow];
                    #pragma omp critical 
                    {
                        Phi.D_index[phiInterval] = phiPointAtPhiOutputIntervalStart.interval;
                        Phi.D_offset[phiInterval] = phiPointAtPhiOutputIntervalStart.offset;
                    }
                    //phiPointAtPhiOutputIntervalStart = phiPoint;
                }
                else {
                    std::cerr << "ERROR: didn't set phi values of last interval of seq!" << std::endl;
                    exit(1);
                }
                */

                if (suff != seqLens[seq+1]) {
                    std::cerr << "ERROR: suff didn't end up at the beginning of the next sequence!" << std::endl;
                    exit(1);
                }
                if (phiPoint.interval != numTopRuns[seq+1]) {
                    std::cerr << "ERROR: phiPoint didn't end up at the beginning of the next sequence!" << std::endl;
                    exit(1);
                }
            }
        }
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //Construct Phi and Samples
        #endif
    }

    void ComputePLCPSamples(sdsl::int_vector<>& intAtEnd, const uint64_t numSequences, const sdsl::int_vector<>& F, 
            const MoveStructureTable& Psi, const std::vector<uint64_t>& numTopRuns, const std::vector<uint64_t>& seqLens,
            const uint64_t sampleInterval, const sdsl::int_vector<>& Psi_Index_Samples, const sdsl::int_vector<>& Psi_Offset_Samples
            #ifndef BENCHFASTONLY
            , const verbosity v
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("LCP Computation"); }
        #endif
        //PLCPsamples = sdsl::int_vector<>(F.size(), 0, 1);

        /*
        sdsl::int_vector<> suffStartingPhiInt(F.size(), 0, sdsl::bits::hi(seqLens.back() - 1) + 1);
        for (uint64_t i = 1; i < F.size(); ++i)
            suffStartingPhiInt[i] = suffStartingPhiInt[i-1] + Phi.data.get<2>(i-1);
         */

        std::atomic<uint64_t> updateWidthsWaiting(0), plcpWritesOccurring(0);
        std::mutex updateWidthMutex;

        uint64_t dangerousInts = 64/intAtEnd.width() + (64 % intAtEnd.width() != 0);

        bool needCompress = true;
        std::vector<uint64_t> prevPsiIntSeqStart(numSequences);
        prevPsiIntSeqStart[0] = intAtEnd[F.size() - 1];
        for (uint64_t i = 1; i < numSequences; ++i)
            prevPsiIntSeqStart[i] = intAtEnd[numTopRuns[i] - 1];
        #pragma omp parallel for schedule (dynamic, 1)
        for (uint64_t seq = 0; seq < numSequences; ++seq) {
            uint64_t suffMatchEnd = seqLens[seq], currIntStart = seqLens[seq];
            MoveStructureTable::IntervalPoint suffMatchEndIntPoint = Psi.map({static_cast<uint64_t>(-1), ((seq)? seq - 1 : numSequences - 1), 0});
            const uint64_t start = numTopRuns[seq];
            const uint64_t end = numTopRuns[seq+1];
            uint64_t prevIntAtEnd = prevPsiIntSeqStart[seq];
            for (uint64_t phiInt = start; phiInt < end; ++phiInt) { 
                //get suffix above phiInt
                MoveStructureStartTable::IntervalPoint suffMatchingToPhiIntPoint = Phi.map({static_cast<uint64_t>(-1), phiInt, 0});
                uint64_t suffMatchingTo = Phi.data.get<2>(suffMatchingToPhiIntPoint.interval) + suffMatchingToPhiIntPoint.offset;
                uint64_t psiIntAbove = prevIntAtEnd;
                MoveStructureTable::IntervalPoint coordAbove = Psi.map({static_cast<uint64_t>(-1), psiIntAbove, 0});
                if (coordAbove.offset == 0) {
                    coordAbove.interval = (coordAbove.interval)? coordAbove.interval - 1 : F.size() - 1;
                    coordAbove.offset = Psi.data.get<2>(coordAbove.interval);
                }
                --coordAbove.offset;
                if (suffMatchEnd - currIntStart) {
                    uint64_t destinationSuff = suffMatchingTo + (suffMatchEnd - currIntStart);
                    uint64_t sampleNum = destinationSuff/sampleInterval;
                    uint64_t closestSampledSuff = sampleNum * sampleInterval;
                    //by default, use above run
                    uint64_t dist = suffMatchEnd - currIntStart;
                    //otherwise, if sample is closer
                    if (dist > destinationSuff - closestSampledSuff) {
                        coordAbove = {static_cast<uint64_t>(-1), Psi_Index_Samples[sampleNum], Psi_Offset_Samples[sampleNum]};
                        dist = destinationSuff - closestSampledSuff;
                    }

                    for (uint64_t i = 0; i < dist; ++i) 
                        coordAbove = Psi.map(coordAbove);
                }
                //compute lcp of interval phiInt
                //the LCP value at the start of this interval is at least the interval length - 1
                //(it is only the interval length - 1 when the suffix at the start of the interval starts with c
                //and is the lexicographically smallest suffix that starts with c.
                //happens by default for all endmarkers since multidollar bwt

                //In our application, LCP can be less than the interval length when Ns are present in the interval
                //since we don't allow matches to include Ns. If this wasn't the case, the currLCP could be initialized to at least (*Phi.intLens)[phiInt] - 1

                while (F[suffMatchEndIntPoint.interval] != 0 && F[suffMatchEndIntPoint.interval] != 5 &&
                        F[suffMatchEndIntPoint.interval] == F[coordAbove.interval]) {
                    suffMatchEndIntPoint = Psi.map(suffMatchEndIntPoint);
                    coordAbove = Psi.map(coordAbove);
                    ++suffMatchEnd;
                }

                //we could also handle the parallelization very easily (and probably faster) with an r bit 
                //buffer. But that increases peak memory by r bits (very minor of course, but an increase nonetheless
                //Maybe implement it as a separate function?

                //phiInt is the index in PLCPsamples to write lcpVal to
                uint64_t lcpVal = suffMatchEnd - currIntStart;
                uint64_t w = sdsl::bits::hi(lcpVal) + 1;

                //use intAtEnd before it's overwritten
                prevIntAtEnd = intAtEnd[phiInt];
                currIntStart = Phi.data.get<2>(phiInt+1);
                if (suffMatchEnd < currIntStart) {
                    suffMatchEnd = currIntStart;
                    suffMatchEndIntPoint = Psi.map({static_cast<uint64_t>(-1), intAtEnd[phiInt], 0});
                }

                //write lcpVal to intAtEnd[phiInt
                if (w > intAtEnd.width()) {
                    ++updateWidthsWaiting;
                    {
                        //(spinlock) busy wait for current writes to finish
                        //this should be fast enough, since writes should be fast
                        //and new ones don't start while updates are pending
                        //I would use shared_lock for this but it's not clear whether
                        //shared_lock or unique_lock gets preference. Is it fair? I suspect
                        //shared_lock can starve unique_lock, but we want unique_lock to take 
                        //precedence
                        //
                        //I can't think of a better way to do this other than spinlock right now.
                        while (plcpWritesOccurring);

                        std::lock_guard<std::mutex> lock(updateWidthMutex);
                        //this function checks if w <= width and if so exits early, so the non atomic check in the above if statement is fine.
                        sdsl::util::expand_width(intAtEnd, w);
                        needCompress = false;
                        dangerousInts = 64/intAtEnd.width() + (64 % intAtEnd.width() != 0);
                    }
                    --updateWidthsWaiting;
                }
                //write value
                bool written = false;
                while (!written) {
                    //busy wait
                    while (updateWidthsWaiting);
                    ++plcpWritesOccurring;
                    if (updateWidthsWaiting) {
                        --plcpWritesOccurring;
                        continue;
                    }
                    written = true;
                    if (phiInt >= start + dangerousInts && phiInt < end - dangerousInts)
                        intAtEnd[phiInt] = lcpVal;
                    else {
                        //mixing omp and std concurrency handling should be fine
                        //if my logic is right ... I think?
                        #pragma omp critical
                        {
                            intAtEnd[phiInt] = lcpVal;
                        }
                    }
                    --plcpWritesOccurring;
                }
            }
        }

        PLCPsamples = std::move(intAtEnd);

        if (needCompress)
            sdsl::util::bit_compress(PLCPsamples);
        #ifndef BENCHFASTONLY
        if (v >= VERB)
            std::cout << "PLCP width: " << static_cast<uint64_t>(PLCPsamples.width()) << std::endl;
        if (v >= TIME) { Timer.stop(); } //LCP Computation
        #endif
    }

    public:
    typedef uint64_t size_type;

	// Constructor to assist with matching statistic computation
	LCPComputer(const std::string& inFile) {
        Timer.start("Loading input lcp_index");
		std::ifstream in = safeOpenFile<std::ifstream>(inFile);
        load(in);
        in.close();
        Timer.stop(); //Loading input lcp_index
	}

    //input: a run length encoding of a multidollar BWT where all dollars are represented by 0
    //All characters between (and including) 0 and max_char are assumed to have more than 0 occurrences
    //in the text. max_char is ththe text
    LCPComputer(rb3_fmi_t * rb3, const std::string& safeTempName
            #ifndef BENCHFASTONLY
            , verbosity v = QUIET
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= VERB) { std::cout << "Number of threads: " << omp_get_max_threads() << "\n"; }
        #endif
        std::ofstream tempOutFile(safeTempName);
        if (!tempOutFile.is_open()) {
            std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' failed to open for writing!" << std::endl;
            exit(1);
        }

        #ifndef BENCHFASTONLY
        sdsl::memory_monitor::granularity(std::chrono::milliseconds(1));
        sdsl::memory_monitor::start();
        #endif

        //Psi.intLens = &Flens;
        uint64_t numSequences;
        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Construct Psi");
            #endif
            ConstructPsi(rb3, F, Psi, numSequences
                    #ifndef BENCHFASTONLY
                    , v
                    #endif
                    );
            rb3_fmi_free(rb3);
        }

        /*
        for (uint64_t i = 0; i < Psi.D_index.size(); ++i) {
            std::cout 
                << F[i] << '\t'
                << (*Psi.intLens)[i] << '\t'
                << Psi.D_index[i] << '\t'
                << Psi.D_offset[i] << '\n';
        }
        {
            MoveStructure::IntervalPoint end{static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            curr = end;
            char con[] = "$ACGTN";
            uint64_t count = 0;
            do {
                curr = Psi.map(curr);
                std::cout << con[F[curr.interval]];
                ++count;
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }
        */

        //numTopRuns and seqLens are vectors of length equal to the number of sequences+1
        //the i+1-th value of numTopRuns is the number of times suffixes of sequence i
        //are at the (SA position corresponding to the) top of a run in the BWT. 
        //(includes the termination symbol at the end of sequence i, $_i)
        //the i+1-th value of seqLens is the length of seq i,
        //(includes the termination symbol) 
        std::vector<uint64_t> numTopRuns, seqLens;
        //for every suffix x at the top of a run in the BWT, 
        //there is an input interval of psi, j, where 
        //suffix x-1 is at the top of the input interval
        //intAtTop stores, for every input interval of psi
        //with suffix x-1 at the top of the input interval,
        //how many suffixes < x - 1 are at the top of a run in the BWT
        uint64_t maxPhiIntLen;
        sdsl::int_vector<> PhiIntLen;
        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Compute Auxiliary Data and Repair Endmarker Psis");
            #endif
            ComputeAuxAndRepairPsi(maxPhiIntLen, numTopRuns, seqLens, intAtTop, F, Psi, PhiIntLen, numSequences
                    #ifndef BENCHFASTONLY
                    , v
                    #endif
                    );
        }

        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Computing and storing intAtEnd and F");
            if (v >= TIME) { Timer.start("Computing and storing intAtEnd"); }
            #endif
            sdsl::int_vector<> intAtEnd(F.size(), 0, sdsl::bits::hi(F.size() - 1) + 1);
            //intAtEnd[j] is the input interval of Psi where the suffix at the end of input interval j of Phi occurs
            //(at the top of, necessarily)
            for (uint64_t i = 0; i < F.size(); ++i)
                intAtEnd[intAtTop[i]] = i;
            sdsl::serialize(intAtEnd, tempOutFile);
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Computing and storing intAtEnd
            if (v >= TIME) { Timer.start("Storing F"); }
            #endif
            sdsl::serialize(F, tempOutFile);
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); }
            #endif
            F = sdsl::int_vector<>();
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.start("Storing intAtTop"); }
            #endif
            sdsl::serialize(intAtTop, tempOutFile);
            tempOutFile.close();
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); }
            #endif
        }
        /*
        {
            std::cout << "Checking if intAtTop is a permutation of [0,r-1]" << std::endl;
            std::vector<bool> a(F.size());
            for (uint64_t i = 0; i < a.size(); ++i) {
                if (intAtTop[i] >= F.size() || a[intAtTop[i]]) {
                    std::cerr << "ERROR: intAtTop is not a permutation of [0,r-1]" << std::endl;
                    exit(1);
                }
                a[intAtTop[i]] = true;
            }
        }

        {
            std::cout << "Psi SA order\ni\trunInd\tF\tintLen\tD_index\tD_offset\n";
            uint64_t k = 0;
            for (uint64_t i = 0; i < Psi.D_index.size(); ++i) {
                for (uint64_t j = 0; j < Flens[i]; ++j)
                    std::cout 
                        << k++ << '\t'
                        << i << '\t'
                        << F[i] << '\t'
                        << (*Psi.intLens)[i] << '\t'
                        << Psi.D_index[i] << '\t'
                        << Psi.D_offset[i] << '\n';
            }
        }
        {
            MoveStructure::IntervalPoint end{static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            curr = end;
            char con[] = "$ACGTN";
            uint64_t count = 0;
            do {
                curr = Psi.map(curr);
                std::cout << con[F[curr.interval]];
                ++count;
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }

        Timer.start("Verifying Psi");
        if (!Psi.permutationLengthN(totalLen)) {
            std::cerr << "ERROR: Psi is not a permutation of length n!" << std::endl;
            exit(1);
        }
        std::cout << "Psi is a permutation of length n\n";
        Timer.stop(); //Verifying Psi
        */

        uint64_t sampleInterval;
        sdsl::int_vector<> Psi_Index_Samples, Psi_Offset_Samples;
        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Construct Phi and Equidistant ISA Samples");
            #endif
            ConstructPhiAndSamples(Psi, PhiIntLen, Psi.data.c, numTopRuns, seqLens, intAtTop, numSequences, maxPhiIntLen, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples
                    #ifndef BENCHFASTONLY
                    , v
                    #endif
                    );
        }

        /*
        {
            uint64_t count = 0, inInt = static_cast<uint64_t>(-1), outInt = static_cast<uint64_t>(-1);
            MoveStructure::IntervalPoint end = {static_cast<uint64_t>(-1), numSequences-1, 0}, curr;
            std::cout << "Phi text order\ni\tinInt\toutInt\tPsiInt\tPsiOff\tPhiInt\tPhiOff\tPhiLen\tPhiDind\tPhiDoff\n";
            std::cout << std::endl;
            curr = end;
            MoveStructure::IntervalPoint pPoint = {static_cast<uint64_t>(-1),0,0}, temp;
            do {
                inInt += (curr.offset == 0);
                outInt += (curr.offset == (*Psi.intLens)[curr.interval] - 1);
                curr = Psi.map(curr);
                temp = Phi.map(pPoint);
                std::cout << count++ << '\t'
                    << inInt << '\t'
                    << outInt << '\t'
                    << curr.interval << '\t'
                    << curr.offset << '\t'
                    << pPoint.interval << '\t'
                    << pPoint.offset << '\t'
                    << (*Phi.intLens)[pPoint.interval] << '\t' 
                    << Phi.D_index[pPoint.interval] << '\t'
                    << Phi.D_offset[pPoint.interval]+pPoint.offset << '\t'
                    << temp.interval << '\t'
                    << temp.offset << std::endl;
                ++pPoint.offset;
                if (pPoint.offset == (*Phi.intLens)[pPoint.interval]){
                    pPoint.offset = 0;
                    ++pPoint.interval;
                }
            } while (curr != end);
            std::cout << std::endl;
            std::cout << count << " characters" << std::endl;
        }

        for (uint64_t i = 0; i < Phi.D_index.size(); ++i) {
            std::cout << Phi.D_index[i] << '\t'
                << Phi.D_offset[i] << '\t'
                << (*Phi.intLens)[i] << '\n';
        }

        {
            MoveStructure::IntervalPoint end = {static_cast<uint64_t>(-1), intAtTop[0], 0}, curr;
            curr = end;
            std::cout << "intAtTop\n";
            std::cout << curr.interval << ' ' << curr.offset << '\n';
            std::cout << "etc\n";
            do {
                curr = Phi.map(curr);
                std::cout << curr.interval << ' ' << curr.offset << '\n';
            } while (curr != end);
        }
        */

        /*
        Verifying Phi is VERY slow compared to verifying Psi ~300 seconds on mtb152 with 64 cores vs ~30 seconds for Psi on coombs c0-4. Why?
        Timer.start("Verifying Phi");
        if (!Phi.permutationLengthN<EXPONENTIAL>(totalLen)) {
            std::cerr << "ERROR: Phi is not a permutation of length n!" << std::endl;
            exit(1);
        }
        std::cout << "Phi is a permutation of length n\n";
        Timer.stop(); //Verifying Phi
        */
        sdsl::int_vector<> intAtEnd;
        std::ifstream tempInFile;
        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Recover intAtEnd and F from disk");
            if (v >= TIME) { Timer.start("Recover intAtEnd from disk"); }
            #endif
            intAtTop = sdsl::int_vector<>();
            tempInFile.open(safeTempName);
            if (!tempInFile.is_open()) {
                std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' failed to open for reading!" << std::endl;
                exit(1);
            }
            sdsl::load(intAtEnd, tempInFile);
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Recover intAtEnd from disk
            if (v >= TIME) { Timer.start("Recover F from disk"); } 
            #endif
            sdsl::load(F, tempInFile);
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Recover F from disk
            #endif
        }


        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Compute PLCP Samples");
            #endif
            ComputePLCPSamples(intAtEnd, numSequences, F, Psi, numTopRuns, seqLens, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples
                    #ifndef BENCHFASTONLY
                    , v
                    #endif
                    );
            Psi_Index_Samples = sdsl::int_vector<>();
            Psi_Offset_Samples = sdsl::int_vector<>();
        }

        /*
        {
            Timer.start("Verifying Psi");
            {
                auto event = sdsl::memory_monitor::event("Verifying Psi");
                if (!Psi.permutationLengthN(totalLen)) {
                    std::cerr << "ERROR: Psi is not a permutation of length n!" << std::endl;
                    exit(1);
                }
                std::cout << "Psi is a permutation of length n\n";
            }
            Timer.stop(); //Verifying Psi
            //Verifying Phi is VERY slow compared to verifying Psi ~300 seconds on mtb152 with 64 cores vs ~30 seconds for Psi on coombs c0-4. Why?
            Timer.start("Verifying Phi");
            {
                auto event = sdsl::memory_monitor::event("Verifying Phi");
                if (!Phi.permutationLengthN<EXPONENTIAL>(totalLen)) {
                    std::cerr << "ERROR: Phi is not a permutation of length n!" << std::endl;
                    exit(1);
                }
                std::cout << "Phi is a permutation of length n\n";
            }
            Timer.stop(); //Verifying Phi
        }

        */
        //printPhiAndLCP(PLCPsamples);

        //reload intAtEnd
        {
            #ifndef BENCHFASTONLY
            auto event = sdsl::memory_monitor::event("Recover intAtTop from disk");
            if (v >= TIME) { Timer.start("Recover intAtTop from disk"); } 
            #endif
            if (!tempInFile.is_open()) {
                std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' is no longer open for reading!" << std::endl;
                exit(1);
            }
            sdsl::load(intAtTop, tempInFile);
            tempInFile.close();
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Recover intAtTop from disk 
            #endif
        }

        /*
        Timer.start("Computing minLCP per run");
        {
            auto event = sdsl::memory_monitor::event("Computing minLCP per run");
            ComputeMinLCPRun(intAtTop, F, Psi, lcpOut);
        }
        Timer.stop(); //Computing minLCP per run"
        */

        #ifndef BENCHFASTONLY
        sdsl::memory_monitor::stop();
        if (v >= TIME) {
            std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;
            std::cout << "peak usage = " << static_cast<double>(sdsl::memory_monitor::peak())/1024 << " kibibytes" << std::endl;
        }

        std::ofstream cstofs("construction.html");
        if (v >= VERB) { std::cout << "writing memory usage visualization to construction.html\n"; }
        sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
        cstofs.close();
        #endif

        /*
        needs to be rewritten, affects bench
        std::ofstream outpsi("lcpcomputer.psi");
        //Flens.serialize(outpsi);
        Psi.serialize(outpsi);
        std::ofstream outphi("lcpcomputer.phi");
        (*Phi.intLens).serialize(outphi);
        Phi.serialize(outphi);
        */
    }

    void ComputeMinLCPRun(std::ofstream& out
            #ifndef BENCHFASTONLY
            , const verbosity v
            #endif
            ) const {
        //O(r log sigma) time, can be skipped if RLBWT is maintained or re-read
        if (v >= TIME) { Timer.start("Computing Min LCP per Run"); }
        struct MappedPositionRunPair {
            uint64_t psiInputInt;
            MoveStructureTable::IntervalPoint runStart;

            MappedPositionRunPair(uint64_t interval, const MoveStructureTable& psi): psiInputInt(interval) {
                runStart = psi.map({static_cast<uint64_t>(-1), psiInputInt, 0});
            }
        };

        struct ComparePositionPair {
            bool operator()(const MappedPositionRunPair& lhs, const MappedPositionRunPair& rhs) const {
                return lhs.runStart.interval > rhs.runStart.interval || (lhs.runStart.interval == rhs.runStart.interval && lhs.runStart.offset > rhs.runStart.offset);
            }
        };


        std::priority_queue<MappedPositionRunPair,std::vector<MappedPositionRunPair>,ComparePositionPair> nextAlph;
        std::vector<uint64_t> endmarkerOrder;
        {
            uint64_t seq = 0;
            while (F[seq] == 0)
                nextAlph.emplace(seq++,Psi);
            endmarkerOrder.reserve(nextAlph.size());
            while(nextAlph.size()) {
                endmarkerOrder.push_back(nextAlph.top().psiInputInt);
                nextAlph.pop();
            }
        }

        //get first of each character
        nextAlph.emplace(endmarkerOrder[0], Psi);
        for (uint64_t i = endmarkerOrder.size(); i < F.size(); ++i)
            if (F[i] != F[i-1])
                nextAlph.emplace(i, Psi);

        MappedPositionRunPair firstRun = nextAlph.top();

        uint64_t runs = 0, endmarkerPosition = 0;
        while (nextAlph.size()) {
            MappedPositionRunPair t = nextAlph.top();
            nextAlph.pop();
            if (F[t.psiInputInt] != 0 && t.psiInputInt != F.size() - 1 && F[t.psiInputInt] == F[t.psiInputInt+1])
                nextAlph.emplace(t.psiInputInt+1, Psi);
            if (F[t.psiInputInt] == 0 && endmarkerPosition != endmarkerOrder.size() - 1)
                nextAlph.emplace(endmarkerOrder[++endmarkerPosition], Psi);
            ++runs;

            uint64_t runLen = Psi.data.get<2>(t.psiInputInt);
            //out << F[t.psiInputInt] << ' ' << runLen << ' ';
            t = (nextAlph.size())? nextAlph.top() : firstRun;
            uint64_t phiInt = (intAtTop[t.psiInputInt]+1) % F.size();
            MoveStructureStartTable::IntervalPoint pPoint = {static_cast<uint64_t>(-1), phiInt, 0};
            uint64_t minL = static_cast<uint64_t>(-1), minLloc = static_cast<uint64_t>(-1);
            for (uint64_t i = 0; i < runLen; ++i) {
                pPoint = Phi.map(pPoint);
                uint64_t l = PLCPsamples[pPoint.interval] - pPoint.offset;
                //use <= for first index
                if (l < minL) {
                    minL = l;
                    minLloc = runLen - 1 - i;
                }
            }
            out << "( " << minLloc << ", " << minL << ")\n";
        }

        if (runs != F.size()) {
            std::cerr << "ERROR: runs found by recovering RLBWT not equal to runs in F!\n";
            exit(1);
        }
        if (v >= TIME) { Timer.stop(); } //Computing Min LCP per Run
    }

    void printRaw(const sdsl::int_vector<>& intAtTop) const {
        std::cout << "LCP\n";
        std::vector<uint64_t> lcp(totalLen);
        MoveStructureStartTable::IntervalPoint phiPoint{static_cast<uint64_t>(-1), intAtTop[0], 0};
        phiPoint.offset = Phi.data.get<2>(phiPoint.interval) - 1;
        phiPoint = Phi.map(phiPoint);
        for (uint64_t i = 0; i < totalLen; ++i) {
            lcp[totalLen - 1 - i] = PLCPsamples[phiPoint.interval] - phiPoint.offset;
            phiPoint = Phi.map(phiPoint);
        }
        for (uint64_t i = 0; i < totalLen; ++i) {
            std::cout << lcp[i] << '\n';
        }

        std::cout << "\n\n";
        std::cout << "PLCP";
        phiPoint = {static_cast<uint64_t>(-1), 0, 0};
        for (uint64_t i = 0; i < totalLen; ++i) {
            std::cout << '\t' << PLCPsamples[phiPoint.interval] - phiPoint.offset;
            ++phiPoint.offset;
            phiPoint.offset %= Phi.data.get<2>(phiPoint.interval);
            phiPoint.interval += (phiPoint.offset == 0);
        }
        std::cout << "\n";
    }

    void printPhiAndLCP(const sdsl::int_vector<>& PLCPsamples) const {
        std::cout << "runInd\tD_index\tD_offset\tintlen\tPLCPsamp\n";
        uint64_t numRuns = PLCPsamples.size();
        for (uint64_t i = 0; i < numRuns; ++i) {
            std::cout << i << '\t'
                << Phi.data.get<0>(i) << '\t'
                << Phi.data.get<1>(i) << '\t'
                << Phi.data.get<2>(i+1) - Phi.data.get<2>(i) << '\t'
                << PLCPsamples[i] << '\n';
        }
    }

    static bool validateRB3(const rb3_fmi_t* rb3);

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v=NULL, std::string name="") const {
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type bytes = 0;

        bytes += sdsl::serialize(totalLen, out, child, "totalLen");
        bytes += sdsl::serialize(F, out, child, "F");
        //bytes += sdsl::serialize(Flens, out, child, "Flens");
        bytes += sdsl::serialize(Psi, out, child, "Psi");
        bytes += sdsl::serialize(intAtTop, out, child, "intAtTop");
        //bytes += sdsl::serialize(PhiIntLen, out, child, "PhiIntLen");
        bytes += sdsl::serialize(Phi, out, child, "Phi");
        bytes += sdsl::serialize(PLCPsamples, out, child, "PLCPsamples");

        sdsl::structure_tree::add_size(child, bytes);
        return bytes;
    }

    void load(std::istream& in) {
        sdsl::load(totalLen, in);
        sdsl::load(F, in);
        //sdsl::load(Flens, in);
        sdsl::load(Psi, in);
        sdsl::load(intAtTop, in);
        //sdsl::load(PhiIntLen, in);
        sdsl::load(Phi, in);
        sdsl::load(PLCPsamples, in);
        //Phi.intLens = &PhiIntLen;
        //Psi.intLens = &Flens;
    }
};

bool LCPComputer::validateRB3(const rb3_fmi_t* rb3){
    if (!rb3->e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return false;
    }
    return true;
}

#endif //#ifndef R_SA_LCP_LCPCOMPUTER_H
