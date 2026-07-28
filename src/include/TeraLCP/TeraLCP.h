#ifndef R_SA_LCP_TeraLCP_H
#define R_SA_LCP_TeraLCP_H

#include"util/util.h"
#include"moveStructure/moveStructure.h"
#include"fm-index.h"
#include<queue>
#include<omp.h>
#include<atomic>
#include<mutex>
#include<stdexcept>
#include<string>
#include<fstream>
#include<cstring>
#include<cstdio>
#include<cerrno>
#include<numeric>
#include<algorithm>

static constexpr const char* lcp_index_extension = ".lcp_index";

class TeraLCP {
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

                //while (F[suffMatchEndIntPoint.interval] != 0 && F[suffMatchEndIntPoint.interval] != 5 &&
                //
                //removed N not matching special case for DNA/RNA.
                //TODO later: slightly optimize by starting at end of phi interval instead of beginning.
                while (F[suffMatchEndIntPoint.interval] != 0 &&
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

    // True once a complete lcp_index has been built. It is set false when the
    // constructor returns early after writing a -stop-after checkpoint, so the
    // driver knows not to emit a (partial) index. Not serialized; it is purely a
    // runtime status and survives the move-assignment used in main().
    bool complete = false;

    // Which slow construct phase to stop after, so a multi-week build can be
    // chunked across successive <=7-day jobs. NONE runs straight through.
    //   A: after endmarker repair (ComputeAuxAndRepairPsi)
    //   B: after Phi + samples (ConstructPhiAndSamples)
    // (Phase C -- LCP computation -- is the existing .lcp_index checkpoint.)
    enum class StopAfter { NONE, A, B };

    // ---- Checkpoint/resume plumbing ------------------------------------------
    // With -checkpoint <dir>, construct durably serializes each phase's outputs
    // into <dir> (one file per structure, via the same sdsl::serialize calls used
    // for the in-process spill and the final index). Each file is written to
    // <name>.tmp then renamed into place (atomic), and the manifest -- recording
    // the last completed phase plus the input identity (r, totalLen) -- is written
    // last, so a job killed mid-write leaves the previous checkpoint intact. On
    // startup construct reads the manifest and jumps to the first incomplete phase,
    // reloading that phase's inputs instead of recomputing the completed ones.

    struct CpManifest { int version = 0, phase = 0; uint64_t r = 0, totalLen = 0, numSequences = 0; bool found = false; };

    static void cpWriteU64Vec(std::ostream& out, const std::vector<uint64_t>& v) {
        uint64_t n = v.size();
        out.write(reinterpret_cast<const char*>(&n), sizeof(n));
        if (n) out.write(reinterpret_cast<const char*>(v.data()), n * sizeof(uint64_t));
    }
    static void cpReadU64Vec(std::istream& in, std::vector<uint64_t>& v) {
        uint64_t n = 0;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        v.resize(n);
        if (n) in.read(reinterpret_cast<char*>(v.data()), n * sizeof(uint64_t));
    }
    // Open a checkpoint component on its .tmp path; pair with cpCommit once written.
    static std::ofstream cpOpen(const std::string& dir, const std::string& name) {
        std::ofstream out(dir + "/" + name + ".tmp", std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            std::cerr << "ERROR: cannot open checkpoint file '" << dir << "/" << name << ".tmp' for writing: " << std::strerror(errno) << "\n";
            exit(1);
        }
        return out;
    }
    // Atomically publish a component written via cpOpen (rename .tmp -> name).
    static void cpCommit(const std::string& dir, const std::string& name) {
        const std::string tmp = dir + "/" + name + ".tmp", fin = dir + "/" + name;
        if (std::rename(tmp.c_str(), fin.c_str()) != 0) {
            std::cerr << "ERROR: cannot commit checkpoint '" << tmp << "' -> '" << fin << "': " << std::strerror(errno) << "\n";
            exit(1);
        }
    }
    static std::ifstream cpOpenIn(const std::string& dir, const std::string& name) {
        std::ifstream in(dir + "/" + name, std::ios::binary);
        if (!in.is_open()) {
            std::cerr << "ERROR: cannot open checkpoint file '" << dir << "/" << name << "' for reading: " << std::strerror(errno) << "\n";
            exit(1);
        }
        return in;
    }
    static CpManifest cpReadManifest(const std::string& dir) {
        CpManifest m;
        std::ifstream in(dir + "/manifest.txt");
        if (!in.is_open()) return m;
        m.found = true;
        std::string k;
        while (in >> k) {
            if (k == "version") in >> m.version;
            else if (k == "phase") in >> m.phase;
            else if (k == "r") in >> m.r;
            else if (k == "totalLen") in >> m.totalLen;
            else if (k == "numSequences") in >> m.numSequences;
        }
        return m;
    }
    // Written last; its rename commits the checkpoint at phase `phase` (1=A, 2=B).
    static void cpWriteManifest(const std::string& dir, int phase, uint64_t r, uint64_t totalLen, uint64_t numSequences) {
        std::ofstream out = cpOpen(dir, "manifest.txt");
        out << "version 1\n" << "phase " << phase << "\n" << "r " << r << "\n"
            << "totalLen " << totalLen << "\n" << "numSequences " << numSequences << "\n";
        out.close();
        cpCommit(dir, "manifest.txt");
    }

    // Boundary A: persist the endmarker-repair outputs (everything phases B and C
    // need that is not recomputed). Psi is const across B and C, so this single
    // snapshot serves both an A-resume and a B-resume.
    void cpWriteCheckpointA(const std::string& dir, uint64_t numSequences,
            const std::vector<uint64_t>& numTopRuns, const std::vector<uint64_t>& seqLens,
            uint64_t maxPhiIntLen, const sdsl::int_vector<>& PhiIntLen) const {
        { std::ofstream o = cpOpen(dir, "F.sdsl");        sdsl::serialize(F, o);        o.close(); cpCommit(dir, "F.sdsl"); }
        { std::ofstream o = cpOpen(dir, "Psi.sdsl");      sdsl::serialize(Psi, o);      o.close(); cpCommit(dir, "Psi.sdsl"); }
        { std::ofstream o = cpOpen(dir, "intAtTop.sdsl"); sdsl::serialize(intAtTop, o); o.close(); cpCommit(dir, "intAtTop.sdsl"); }
        { std::ofstream o = cpOpen(dir, "aux.bin");
          sdsl::serialize(totalLen, o);
          sdsl::serialize(numSequences, o);
          sdsl::serialize(maxPhiIntLen, o);
          cpWriteU64Vec(o, numTopRuns); //replace with sdsl::serialize?
          cpWriteU64Vec(o, seqLens); //replace with sdsl::serialize?
          sdsl::serialize(PhiIntLen, o);
          o.close(); cpCommit(dir, "aux.bin"); }
    }
    void cpLoadCheckpointA(const std::string& dir, uint64_t& numSequences,
            std::vector<uint64_t>& numTopRuns, std::vector<uint64_t>& seqLens,
            uint64_t& maxPhiIntLen, sdsl::int_vector<>& PhiIntLen) {
        { std::ifstream i = cpOpenIn(dir, "F.sdsl");        sdsl::load(F, i); }
        { std::ifstream i = cpOpenIn(dir, "Psi.sdsl");      sdsl::load(Psi, i); }
        { std::ifstream i = cpOpenIn(dir, "intAtTop.sdsl"); sdsl::load(intAtTop, i); }
        { std::ifstream i = cpOpenIn(dir, "aux.bin");
          sdsl::load(totalLen, i);
          sdsl::load(numSequences, i);
          sdsl::load(maxPhiIntLen, i);
          cpReadU64Vec(i, numTopRuns); //replace with sdsl::load?
          cpReadU64Vec(i, seqLens); //replace with sdsl::load?
          sdsl::load(PhiIntLen, i); }
    }
    // Boundary B: persist the Phi-and-samples outputs. F/Psi/intAtTop/intAtEnd/aux
    // are already durable from boundary A and the spill block.
    void cpWriteCheckpointB(const std::string& dir, uint64_t sampleInterval,
            const sdsl::int_vector<>& Psi_Index_Samples, const sdsl::int_vector<>& Psi_Offset_Samples) const {
        { std::ofstream o = cpOpen(dir, "Phi.sdsl"); sdsl::serialize(Phi, o); o.close(); cpCommit(dir, "Phi.sdsl"); }
        { std::ofstream o = cpOpen(dir, "samples.bin");
          sdsl::serialize(sampleInterval, o);
          sdsl::serialize(Psi_Index_Samples, o);
          sdsl::serialize(Psi_Offset_Samples, o);
          o.close(); cpCommit(dir, "samples.bin"); }
    }
    // Load everything phase C (ComputePLCPSamples) needs. intAtTop is intentionally
    // NOT loaded here -- like the straight-through path, it is reloaded only after
    // phase C (cpLoadIntAtTop) to keep peak memory the same.
    void cpLoadCheckpointB(const std::string& dir, uint64_t& numSequences,
            std::vector<uint64_t>& numTopRuns, std::vector<uint64_t>& seqLens,
            uint64_t& sampleInterval, sdsl::int_vector<>& Psi_Index_Samples,
            sdsl::int_vector<>& Psi_Offset_Samples, sdsl::int_vector<>& intAtEnd) {
        { std::ifstream i = cpOpenIn(dir, "F.sdsl");        sdsl::load(F, i); }
        { std::ifstream i = cpOpenIn(dir, "Psi.sdsl");      sdsl::load(Psi, i); }
        { std::ifstream i = cpOpenIn(dir, "intAtEnd.sdsl"); sdsl::load(intAtEnd, i); }
        { std::ifstream i = cpOpenIn(dir, "Phi.sdsl");      sdsl::load(Phi, i); }
        { std::ifstream i = cpOpenIn(dir, "samples.bin");
          sdsl::load(sampleInterval, i);
          sdsl::load(Psi_Index_Samples, i);
          sdsl::load(Psi_Offset_Samples, i); }
        { std::ifstream i = cpOpenIn(dir, "aux.bin");
          uint64_t maxPhiIntLen; // unused for phase C
          sdsl::load(totalLen, i);
          sdsl::load(numSequences, i);
          sdsl::load(maxPhiIntLen, i);
          cpReadU64Vec(i, numTopRuns);
          cpReadU64Vec(i, seqLens); }
    }
    void cpLoadIntAtTop(const std::string& dir) {
        std::ifstream i = cpOpenIn(dir, "intAtTop.sdsl");
        sdsl::load(intAtTop, i);
    }

    TeraLCP() = default;

	// Constructor to assist with matching statistic computation
	TeraLCP(const std::string& inFile
            #ifndef BENCHFASTONLY
            , verbosity v = QUIET
            #endif
            ) {
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.start("LCP index loading from file"); }
        #endif
		std::ifstream in = safeOpenFile<std::ifstream>(inFile);
        load(in);
        in.close();
        complete = true;
        #ifndef BENCHFASTONLY
        if (v >= TIME) { Timer.stop(); } //LCP index loading from file
        #endif
	}

    //input: a run length encoding of a multidollar BWT where all dollars are represented by 0
    //All characters between (and including) 0 and max_char are assumed to have more than 0 occurrences
    //in the text. max_char is ththe text
    TeraLCP(rb3_fmi_t * rb3, const std::string& safeTempName
            #ifndef BENCHFASTONLY
            , verbosity v = QUIET
            #endif
            , const std::string& checkpointDir = ""
            , StopAfter stopAfter = StopAfter::NONE
            ) {
        #ifndef BENCHFASTONLY
        if (v >= VERB) { std::cout << "Number of threads: " << omp_get_max_threads() << "\n"; }
        #endif

        // Checkpoint/resume: when a checkpoint dir already holds a manifest, jump to
        // the first incomplete phase and reload the completed phases' outputs.
        const bool useCp = !checkpointDir.empty();
        int resumeFrom = 0; // 0 = fresh, 1 = phase A done, 2 = phase B done
        if (useCp) {
            CpManifest m = cpReadManifest(checkpointDir);
            resumeFrom = m.found ? m.phase : 0;
            #ifndef BENCHFASTONLY
            if (v >= TIME && resumeFrom > 0)
                std::cout << "Resuming construct from checkpoint phase " << resumeFrom << " in '" << checkpointDir << "'\n";
            #endif
        }

        // State shared across phases (constructor locals in the straight-through
        // path; loaded from the checkpoint on a resume).
        uint64_t numSequences = 0;
        std::vector<uint64_t> numTopRuns, seqLens;
        uint64_t maxPhiIntLen = 0;
        sdsl::int_vector<> PhiIntLen;
        uint64_t sampleInterval = 0;
        sdsl::int_vector<> Psi_Index_Samples, Psi_Offset_Samples;
        sdsl::int_vector<> intAtEnd;
        std::ofstream tempOutFile;
        std::ifstream tempInFile;
        uint64_t rCount = 0; // run count (F.size()) captured before F is spilled

        #ifndef BENCHFASTONLY
        //sdsl::memory_monitor::granularity(std::chrono::milliseconds(1));
        //sdsl::memory_monitor::start();
        #endif

        // Open the temp spill file up front (fail fast on a bad path) whenever we
        // will run the spill block -- i.e. unless we are resuming past phase B.
        if (resumeFrom < 2) {
            tempOutFile.open(safeTempName);
            if (!tempOutFile.is_open()) {
                std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' failed to open for writing!" << std::endl;
                exit(1);
            }
        }

        //Psi.intLens = &Flens;
        if (resumeFrom < 1) {
            {
                #ifndef BENCHFASTONLY
                //auto event = sdsl::memory_monitor::event("Construct Psi");
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
            //for every suffix x at the top of a run in the BWT,
            //there is an input interval of psi, j, where
            //suffix x-1 is at the top of the input interval
            //intAtTop stores, for every input interval of psi
            //with suffix x-1 at the top of the input interval,
            //how many suffixes < x - 1 are at the top of a run in the BWT
            {
                #ifndef BENCHFASTONLY
                //auto event = sdsl::memory_monitor::event("Compute Auxiliary Data and Repair Endmarker Psis");
                #endif
                ComputeAuxAndRepairPsi(maxPhiIntLen, numTopRuns, seqLens, intAtTop, F, Psi, PhiIntLen, numSequences
                        #ifndef BENCHFASTONLY
                        , v
                        #endif
                        );
            }
            rCount = F.size();

            // ---- Checkpoint boundary A (after endmarker repair) ----
            if (useCp) {
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.start("Checkpoint write (phase A)"); }
                #endif
                cpWriteCheckpointA(checkpointDir, numSequences, numTopRuns, seqLens, maxPhiIntLen, PhiIntLen);
                cpWriteManifest(checkpointDir, 1, rCount, totalLen, numSequences);
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.stop(); } //Checkpoint write (phase A)
                #endif
                if (stopAfter == StopAfter::A) {
                    #ifndef BENCHFASTONLY
                    if (v >= TIME) { std::cout << "Stopping after phase A as requested; checkpoint written to '" << checkpointDir << "'\n"; }
                    #endif
                    complete = false;
                    return;
                }
            }
        }
        else {
            // Resume: phase A already complete, so the input FMD is not needed.
            rb3_fmi_free(rb3);
            // Only reload phase A's outputs when phase B has NOT yet run; on a
            // phase-B resume, cpLoadCheckpointB below loads everything phase C
            // needs (and intAtTop is reloaded after phase C), so loading A here
            // would be redundant I/O and would inflate peak memory.
            if (resumeFrom == 1) {
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.start("Checkpoint load (phase A)"); }
                #endif
                cpLoadCheckpointA(checkpointDir, numSequences, numTopRuns, seqLens, maxPhiIntLen, PhiIntLen);
                rCount = F.size();
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.stop(); } //Checkpoint load (phase A)
                #endif
            }
        }

        // ---- Spill block + phase B (Construct Phi and Samples) ----
        // Run when phase B has not yet completed (fresh run or an A-resume).
        if (resumeFrom < 2) {
            {
                #ifndef BENCHFASTONLY
                //auto event = sdsl::memory_monitor::event("Computing and storing intAtEnd and F");
                if (v >= TIME) { Timer.start("Computing and storing intAtEnd"); }
                #endif
                intAtEnd = sdsl::int_vector<>(F.size(), 0, sdsl::bits::hi(F.size() - 1) + 1);
                //intAtEnd[j] is the input interval of Psi where the suffix at the end of input interval j of Phi occurs
                //(at the top of, necessarily)
                for (uint64_t i = 0; i < F.size(); ++i)
                    intAtEnd[intAtTop[i]] = i;
                // Durably persist intAtEnd for a phase-B resume (F and intAtTop are
                // already durable from boundary A); the temp spill below is the
                // in-process copy used by the straight-through recovery.
                if (useCp) { std::ofstream o = cpOpen(checkpointDir, "intAtEnd.sdsl"); sdsl::serialize(intAtEnd, o); o.close(); cpCommit(checkpointDir, "intAtEnd.sdsl"); }
                sdsl::serialize(intAtEnd, tempOutFile);
                intAtEnd = sdsl::int_vector<>();
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

            {
                #ifndef BENCHFASTONLY
                //auto event = sdsl::memory_monitor::event("Construct Phi and Equidistant ISA Samples");
                #endif
                ConstructPhiAndSamples(Psi, PhiIntLen, Psi.data.c, numTopRuns, seqLens, intAtTop, numSequences, maxPhiIntLen, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples
                        #ifndef BENCHFASTONLY
                        , v
                        #endif
                        );
            }

            // ---- Checkpoint boundary B (after Phi + samples) ----
            if (useCp) {
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.start("Checkpoint write (phase B)"); }
                #endif
                cpWriteCheckpointB(checkpointDir, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples);
                cpWriteManifest(checkpointDir, 2, rCount, totalLen, numSequences);
                #ifndef BENCHFASTONLY
                if (v >= TIME) { Timer.stop(); } //Checkpoint write (phase B)
                #endif
                if (stopAfter == StopAfter::B) {
                    #ifndef BENCHFASTONLY
                    if (v >= TIME) { std::cout << "Stopping after phase B as requested; checkpoint written to '" << checkpointDir << "'\n"; }
                    #endif
                    complete = false;
                    return;
                }
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
            {
                #ifndef BENCHFASTONLY
                //auto event = sdsl::memory_monitor::event("Recover intAtEnd and F from disk");
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
        }
        else {
            // Resume: phases A and B already complete. Load everything phase C needs
            // (intAtTop is loaded after phase C, mirroring the straight-through path).
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.start("Checkpoint load (phase B)"); }
            #endif
            cpLoadCheckpointB(checkpointDir, numSequences, numTopRuns, seqLens, sampleInterval, Psi_Index_Samples, Psi_Offset_Samples, intAtEnd);
            rCount = F.size();
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Checkpoint load (phase B)
            #endif
        }


        {
            #ifndef BENCHFASTONLY
            //auto event = sdsl::memory_monitor::event("Compute PLCP Samples");
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

        //reload intAtTop
        {
            #ifndef BENCHFASTONLY
            //auto event = sdsl::memory_monitor::event("Recover intAtTop from disk");
            if (v >= TIME) { Timer.start("Recover intAtTop from disk"); }
            #endif
            if (resumeFrom < 2) {
                if (!tempInFile.is_open()) {
                    std::cerr << "ERROR: File provided for temporary writing/reading, '" << safeTempName << "' is no longer open for reading!" << std::endl;
                    exit(1);
                }
                sdsl::load(intAtTop, tempInFile);
                tempInFile.close();
            } else {
                // Phase-B resume: intAtTop comes from the checkpoint, not the temp file.
                cpLoadIntAtTop(checkpointDir);
            }
            #ifndef BENCHFASTONLY
            if (v >= TIME) { Timer.stop(); } //Recover intAtTop from disk
            #endif
        }

        complete = true;

        /*
        Timer.start("Computing minLCP per run");
        {
            auto event = sdsl::memory_monitor::event("Computing minLCP per run");
            ComputeMinLCPRun(intAtTop, F, Psi, lcpOut);
        }
        Timer.stop(); //Computing minLCP per run"
        */

        #ifndef BENCHFASTONLY
        //sdsl::memory_monitor::stop();
        //if (v >= TIME) {
            //std::cout << "peak usage = " << sdsl::memory_monitor::peak() << " bytes" << std::endl;
            //std::cout << "peak usage = " << static_cast<double>(sdsl::memory_monitor::peak())/1024 << " kibibytes" << std::endl;
        //}

        //std::ofstream cstofs("construction.html");
        //if (v >= VERB) { std::cout << "writing memory usage visualization to construction.html\n"; }
        //sdsl::memory_monitor::write_memory_log<sdsl::HTML_FORMAT>(cstofs);
        //cstofs.close();
        #endif

        /*
        needs to be rewritten, affects bench
        std::ofstream outpsi("TeraLCP.psi");
        //Flens.serialize(outpsi);
        Psi.serialize(outpsi);
        std::ofstream outphi("TeraLCP.phi");
        (*Phi.intLens).serialize(outphi);
        Phi.serialize(outphi);
        */
    }

    void ComputeMinLCPRunSequential(std::ofstream& out
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

    // RunInfo holds per-run lengths and run symbols (alphabet codes) in BWT order.
    struct RunInfo {
        std::vector<uint64_t> lengths;
        std::vector<uint64_t> symbols;
    };

    // Per-run LCP summary in BWT order with absolute BWT-row positions, used to
    // derive thresholds via an O(r) sweep. Filled by
    // ComputeMinLCPRunParallelDestructive when a non-null pointer is supplied, so
    // thresholds piggyback on the existing parallel per-run traversal.
    struct RunSummary {
        sdsl::int_vector<> topLCP;       // LCP at each run's first (top) row
        sdsl::int_vector<> minLCP;       // min LCP within each run
        sdsl::int_vector<> firstMinPos;  // smallest BWT row achieving the run min
        sdsl::int_vector<> lastMinPos;   // largest  BWT row achieving the run min
        sdsl::int_vector<> runStarts;    // BWT row of each run's first element (size runs+1)
    };

    // [Tier B] One fixed-size per-run record spilled to disk during the reverse
    // reordering pass, so the run-sized threshold summary never materializes in RAM.
    // Fields are absolute BWT rows / LCP values (same semantics as RunSummary).
    struct SpillRec { uint64_t runStart, topLCP, minLCP, firstMin, lastMin; };

    std::pair<sdsl::int_vector<>,sdsl::int_vector<>>
    ComputeMinLCPRunParallelDestructive(
            #ifndef BENCHFASTONLY
            const verbosity v,
            #endif
            RunSummary* summary = nullptr,
            std::ostream* spill = nullptr
            ) {
        // wantSummary drives the extra per-run captures (top LCP, first-min offset).
        // The captured summary is delivered either into RAM arrays (summary != null,
        // boundary path) or streamed as records to disk (spill != null, default path).
        const bool wantSummary = (summary != nullptr) || (spill != nullptr);
        if (v >= TIME) { Timer.start("Computing Min LCP per Run"); }
        uint64_t runs = F.size(), psilenwidth = Psi.data.c;
        uint64_t intAtTopBWT = intAtTop[0];
        //std::cout << intAtTopBWT << std::endl;
        assert(runs == Phi.num_intervals());

        F = sdsl::int_vector<>();
        Psi = MoveStructureTable();
        intAtTop = sdsl::int_vector<>();
        //bits allowed to play with: 
        //r log r + r log Psi_len_max >= r log n from ISA samples during LCP sampling
        //r log r + 2 r log Psi_len_max >= r log n from Psi removal
        //r log sigma from F removal

        //bits used here:
        //r log r : prevRunIntAtTop
        //r log lcp_max : thisRunMin
        //r log Psi_len_max : thisRunMinLoc
        //r log Psi_len_max : thisRunLength
        sdsl::int_vector<> prevRunIntAtTop, thisRunMin, thisRunMinLoc, thisRunLength;

        //each interval j of Phi is the suffix at the top of a run in the BWT, call it run i
        //then 
        //  prevrunIntAtTop[j] stores the interval j' where j' is the suffix at the top of run i - 1
        //  thisRunMin[j] stores the minLCP of run i - 1
        //  thisRunMinLoc[j] stores the location of the minLCP of run i - 1 (relative in the run)
        //  thisRunLength[j] stores the length of run i - 1
        prevRunIntAtTop = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(runs - 1) + 1);
        thisRunMin = sdsl::int_vector<>(runs, 0, PLCPsamples.width());
        thisRunMinLoc = sdsl::int_vector<>(runs, 0, psilenwidth);
        thisRunLength = sdsl::int_vector<>(runs, 0, psilenwidth);

        // Extra per-run captures for thresholds (only when a summary is wanted):
        //   thisRunTop[j]         : LCP at the run's top row (offset 0)
        //   thisRunMinLocFirst[j] : smallest in-run offset achieving the run min
        // (thisRunMinLoc already holds the largest such offset.)
        sdsl::int_vector<> thisRunTop, thisRunMinLocFirst;
        if (wantSummary) {
            thisRunTop = sdsl::int_vector<>(runs, 0, PLCPsamples.width());
            thisRunMinLocFirst = sdsl::int_vector<>(runs, 0, psilenwidth);
        }

        auto PLCP = [this] (const MoveStructureStartTable::IntervalPoint& p) {
            return PLCPsamples[p.interval] - p.offset;
        };

        if (v >= TIME) { Timer.start("Parallel per run computation"); }
        const uint64_t blockSize = std::min(std::max(((runs+3)/4)/omp_get_num_threads(), static_cast<uint64_t>(1)), static_cast<uint64_t>(1024));
        const uint64_t numBlocks = runs/blockSize + ((runs%blockSize) != 0);
        const uint64_t minBitWidth = std::min(std::min(psilenwidth, static_cast<uint64_t>(PLCPsamples.width())), static_cast<uint64_t>(prevRunIntAtTop.width()));
        const uint64_t dangerousInts = 64/minBitWidth + ((64%minBitWidth) != 0);
        if (v >= VERB) { std::cout << "Block size: " << blockSize << std::endl; }
        #pragma omp parallel for schedule(dynamic, 1)
        for (uint64_t block = 0; block < numBlocks; ++block) {
            const uint64_t start = block*blockSize;
            const uint64_t end = std::min(runs, start + blockSize);
            const uint64_t safeStart = start + dangerousInts,
                  safeEnd = end - dangerousInts;
            for (uint64_t run = start; run != end; ++run) {
                //compute values for run
                MoveStructureStartTable::IntervalPoint p = {Phi.data.get<2>(run), run, 0};
                uint64_t minLCP = static_cast<uint64_t> (-1);
                assert(PLCP(Phi.map(p)) < minLCP);
                uint64_t minLCPLoc = minLCP;
                uint64_t runLen = 0, l;
                //std::cout << "run " << run << std::endl;
                //std::cout << "suff " << p.position << std::endl;
                // For the threshold summary we also track the largest-i min (the
                // first/topmost offset in the run) and the run's top LCP.
                uint64_t minLCPLocLast = 0, topVal = 0;
                do {
                    p = Phi.map(p);
                    //std::cout << "suff " << p.position << std::endl;
                    //std::cout << "lcp " << PLCP(p) << std::endl;
                    ++runLen;
                    //currently outputs last location in run, use <= for first location in run
                    if ((l = PLCP(p)) < minLCP) {
                        //std::cout << "in if" << std::endl;
                        minLCPLoc = runLen;
                        minLCPLocLast = runLen;
                        minLCP = l;
                        //std::cout << "minLCPLoc " << minLCPLoc << " minLCP " << minLCP << std::endl;
                    } else if (l == minLCP) {
                        minLCPLocLast = runLen;
                    }
                    topVal = l; // last iteration is the run's top row (offset 0)
                } while (p.offset);
                assert(minLCPLoc != static_cast<uint64_t>(-1));
                minLCPLoc = runLen - minLCPLoc;
                //std::cout << "runLen " << runLen << " newmiNLCPLoc " << minLCPLoc << std::endl;
                const uint64_t minLocFirst = runLen - minLCPLocLast; // smallest (first) min offset

                //store values for run
                if (run >= safeStart && run < safeEnd) {
                    prevRunIntAtTop[run] = p.interval;
                    thisRunMin[run] = minLCP;
                    thisRunMinLoc[run] = minLCPLoc;
                    thisRunLength[run] = runLen;
                    if (wantSummary) {
                        thisRunTop[run] = topVal;
                        thisRunMinLocFirst[run] = minLocFirst;
                    }
                }
                else {
                    #pragma omp critical
                    {
                        prevRunIntAtTop[run] = p.interval;
                        thisRunMin[run] = minLCP;
                        thisRunMinLoc[run] = minLCPLoc;
                        thisRunLength[run] = runLen;
                        if (wantSummary) {
                            thisRunTop[run] = topVal;
                            thisRunMinLocFirst[run] = minLocFirst;
                        }
                    }
                }
            }
        }
        if (v >= TIME) { Timer.stop(); } //Parallel per run computation

        if (v >= TIME) { Timer.start("bit compress"); }
        sdsl::util::bit_compress(thisRunMin);
        sdsl::util::bit_compress(thisRunMinLoc);
        if (v >= TIME) { Timer.stop(); } //bit compress

        Phi = MoveStructureStartTable();
        PLCPsamples = sdsl::int_vector<>();

        // [Tier A] The base per-run minLCP/minLCPLoc are only consumed by the rlcp
        // return path (TeraLCP.cpp:319). writeThresholdsParallel discards the return,
        // and summary->minLCP already carries the same reordered per-run minima, so
        // when building a threshold summary we skip these two run-sized arrays
        // (~r*(lcp_width + log n) bits, ~30 GB at r~3.5e9) at the reordering peak.
        sdsl::int_vector<> minLCP, minLCPLoc;
        if (!wantSummary) {
            minLCP    = sdsl::int_vector<>(runs, 0, thisRunMin.width());
            minLCPLoc = sdsl::int_vector<>(runs, 0, sdsl::bits::hi(totalLen - 1) + 1);
        }

        const uint64_t posWidth = sdsl::bits::hi(totalLen) + 1;
        if (summary) {   // RAM summary (boundary-mode path); the spill path streams
                         // records to disk in the reordering loop below instead.
            summary->topLCP      = sdsl::int_vector<>(runs, 0, thisRunTop.width());
            summary->minLCP      = sdsl::int_vector<>(runs, 0, thisRunMin.width());
            summary->firstMinPos = sdsl::int_vector<>(runs, 0, posWidth);
            summary->lastMinPos  = sdsl::int_vector<>(runs, 0, posWidth);
            summary->runStarts   = sdsl::int_vector<>(runs + 1, 0, posWidth);
            summary->runStarts[runs] = totalLen;
        }

        if (v >= TIME) { Timer.start("Sequential reordering"); }
        //curr is the intAtTop of run currRun+1 % runs
        //startPos is the startPos of run currRun
        uint64_t curr = intAtTopBWT, currRun = runs - 1, startPos = totalLen, prev;
        uint64_t count = 0;
        do {
            startPos -= thisRunLength[curr];
            prev = prevRunIntAtTop[curr];
            // startPos is now runStarts[currRun]; convert in-run offsets to absolute
            // BWT rows for the threshold sweep. [Tier B] When spilling, append a
            // fixed-size record to disk instead of filling the run-sized summary
            // arrays; records emerge in reverse run order and the streaming sweep
            // reads them back to front.
            if (spill) {
                const SpillRec rec{ startPos, thisRunTop[curr], thisRunMin[curr],
                                    startPos + thisRunMinLocFirst[curr],
                                    startPos + thisRunMinLoc[curr] };
                spill->write(reinterpret_cast<const char*>(&rec), sizeof(rec));
            } else if (summary) {
                summary->runStarts[currRun]   = startPos;
                summary->topLCP[currRun]      = thisRunTop[curr];
                summary->minLCP[currRun]      = thisRunMin[curr];
                summary->firstMinPos[currRun] = startPos + thisRunMinLocFirst[curr];
                summary->lastMinPos[currRun]  = startPos + thisRunMinLoc[curr];
            } else {
                minLCP[currRun]    = thisRunMin[curr];
                //remove + startPos for relative positions
                //minLCPLoc[currRun] = thisRunMinLoc[currRun] + startPos;
                minLCPLoc[currRun] = thisRunMinLoc[curr];
            }
            ++count;
            //minLCPLoc[currRun] = startPos;
            curr = prev;
            assert(currRun != 0 || curr == intAtTopBWT);
            --currRun;
        } while (curr != intAtTopBWT); 
        assert(currRun == static_cast<uint64_t>(-1));
        assert(startPos == 0);
        if (v >= TIME) { Timer.stop(); } //Sequential reordering
        if (v >= TIME) { Timer.stop(); } //Computing Min LCP per Run
        return {minLCPLoc, minLCP};
    }

    /**
     * Performs the Phi walk and invokes callback(pos, val) for each BWT position
     * in order totalLen-1, totalLen-2, ..., 0. Streams LCP values without
     * materializing the full n-sized array, so callers can accumulate per-run
     * summaries in O(r) space.
     */
    template<typename Func>
    void phiWalkLCP(Func&& callback) const {
        const uint64_t runs = F.size();
        if (runs == 0) return;
        MoveStructureStartTable::IntervalPoint phiPoint{
            static_cast<uint64_t>(-1), intAtTop[0], 0};
        // Initialize to the predecessor of BWT position totalLen-1 w.r.t. Phi.
        // Phi advances in reverse SA order (lex largest first) and we want the
        // first callback to receive pos = totalLen-1. That predecessor is the last
        // position of the first Phi interval (intAtTop[0]); lenFirst is its length.
        uint64_t lenFirst = Phi.data.get<2>(phiPoint.interval + 1) - Phi.data.get<2>(phiPoint.interval);
        phiPoint.offset = (lenFirst > 0) ? (lenFirst - 1) : 0;
        phiPoint = Phi.map(phiPoint);
        for (uint64_t i = 0; i < totalLen; ++i) {
            uint64_t plcpSample = PLCPsamples[phiPoint.interval];
            uint64_t val = (phiPoint.offset <= plcpSample) ? (plcpSample - phiPoint.offset) : 0;
            if (val > totalLen) val = 0;
            callback(totalLen - 1 - i, val);
            phiPoint = Phi.map(phiPoint);
        }
    }

    /**
     * Extracts RunInfo (lengths and symbols in BWT order) from an open FMD.
     * Matches ConstructPsi / the lcp_index: raw rld_dec runs with symbol 0
     * (sentinel) are split into l runs of length 1 each; other symbols stay one
     * run per decode. Caller must ensure rb3 is valid and has e != nullptr.
     */
    static RunInfo runInfoFromFMD(const rb3_fmi_t* rb3) {
        if (!rb3 || !rb3->e) {
            throw std::runtime_error("runInfoFromFMD: invalid or multirope FMD (e is null)");
        }
        RunInfo info;
        rlditr_t itr;
        rld_itr_init(rb3->e, &itr, 0);
        int c = 0;
        int64_t l;
        while ((l = rld_dec(rb3->e, &itr, &c, 0)) > 0) {
            if (c == 0) {
                for (uint64_t i = 0; i < static_cast<uint64_t>(l); ++i) {
                    info.lengths.push_back(1);
                    info.symbols.push_back(0);
                }
            } else {
                info.lengths.push_back(static_cast<uint64_t>(l));
                info.symbols.push_back(static_cast<uint64_t>(c));
            }
        }
        return info;
    }

private:
    // Self-describing threshold-file format ("TLTHR v1"): a 32-byte header
    //   [0..7]  magic 93 'T' 'L' 'T' 'H' 'R' 00 01(version)
    //   [8]     field width (bytes per record), [9] flags (bit0 = little-endian),
    //   [10..11] reserved, [12..19] record count (uint64 LE), [20..31] reserved
    // followed by `count` little-endian records of `width` bytes each. The reader
    // autodetects this header; --thr-pfp emits the old headerless 5-byte
    // records instead. Format is specified in BenLangmead/thresh-tools (SPEC.md).
    struct ThrFileWriter {
        std::ofstream thrOut, thrPosOut;
        std::vector<char> buf1, buf2;
        unsigned width;
        ThrFileWriter(const std::string& basePath, uint64_t count, unsigned width_, bool legacy)
            : buf1(4 << 20), buf2(4 << 20), width(width_) {
            thrOut.open(basePath + ".thr", std::ios::binary);
            thrPosOut.open(basePath + ".thr_pos", std::ios::binary);
            if (!thrOut.is_open() || !thrPosOut.is_open())
                throw std::runtime_error("writeThresholds: failed to open output files");
            thrOut.rdbuf()->pubsetbuf(buf1.data(), buf1.size());
            thrPosOut.rdbuf()->pubsetbuf(buf2.data(), buf2.size());
            if (!legacy) { writeHeader(thrOut, width, count); writeHeader(thrPosOut, width, count); }
        }
        static void writeHeader(std::ofstream& o, unsigned width, uint64_t count) {
            char h[32] = {0};
            static const unsigned char MAGIC[8] = {0x93, 'T', 'L', 'T', 'H', 'R', 0x00, 0x01};
            std::memcpy(h, MAGIC, 8);
            h[8] = static_cast<char>(width);
            h[9] = 1;                          // flags: bit0 = little-endian
            std::memcpy(h + 12, &count, 8);    // record count (LE)
            o.write(h, 32);
        }
        void writeOne(uint64_t thr, uint64_t pos) {
            char b[8];
            std::memcpy(b, &thr, 8); thrOut.write(b, width);   // low `width` LE bytes
            std::memcpy(b, &pos, 8); thrPosOut.write(b, width);
        }
    };

    // Bytes per threshold field: 5 for --thr-pfp, an explicit override, or the
    // minimum width that holds the BWT length totalLen (positions/values are <= n).
    unsigned thrFieldWidth(bool legacy, unsigned widthOverride) const {
        if (legacy) return 5;
        if (widthOverride) return widthOverride;
        unsigned b = 1;
        while (b < 8 && totalLen >= (1ULL << (8 * b))) ++b;
        return b;
    }

    /**
     * Shared O(r)/O(sigma) sticky-register sweep that turns a per-run LCP summary
     * into .thr/.thr_pos. Generic over the array types so both the phi-walk
     * path (std::vector) and the fused parallel path (sdsl::int_vector) reuse it.
     * All position arrays use absolute BWT rows; runStarts has size runs+1.
     */
    template<typename SymArr, typename StartArr, typename TopArr,
             typename MinArr, typename PosArr>
    void writeThresholdsSweep(const std::string& basePath, bool boundaryMode,
                              const SymArr& symbols, const StartArr& runStarts,
                              const TopArr& topLCP, const MinArr& minLCP,
                              const PosArr& firstMinPos, const PosArr& lastMinPos,
                              unsigned width, bool legacy) const {
        const uint64_t runs = symbols.size();
        if (runs == 0) return;
        const uint64_t SENT = static_cast<uint64_t>(-1);

        // Per-character sticky registers over the currently open gap: running min
        // LCP and the earliest/latest BWT rows achieving it.
        uint64_t maxSym = 0;
        for (uint64_t i = 0; i < runs; ++i) maxSym = std::max<uint64_t>(maxSym, symbols[i]);
        std::vector<uint64_t> regMin(maxSym + 1, SENT),
                              regFirst(maxSym + 1, 0), regLast(maxSym + 1, 0);
        std::vector<char> seenChar(maxSym + 1, 0);

        ThrFileWriter w(basePath, runs, width, legacy);
        auto writeOne = [&](uint64_t thr, uint64_t pos) { w.writeOne(thr, pos); };

        for (uint64_t i = 0; i < runs; ++i) {
            uint64_t c = symbols[i];
            uint64_t thrVal, thrPosVal;
            if (!seenChar[c]) {
                thrVal = 0; thrPosVal = 0;          // first run of c: placeholder
                seenChar[c] = 1;
            } else {
                // Gap min = min over runs since the previous c-run (regMin[c]),
                // plus run i's first row (topLCP[i] at runStarts[i]).
                thrVal = regMin[c];
                uint64_t gFirst = regFirst[c], gLast = regLast[c];
                uint64_t topi = topLCP[i], starti = runStarts[i];
                if (topi < thrVal) {
                    thrVal = topi; gFirst = gLast = starti;
                } else if (topi == thrVal) {
                    if (starti > gLast) gLast = starti;
                }
                if (thrVal == SENT) {               // gap empty (adjacent c-runs)
                    thrVal = 0; thrPosVal = 0;
                } else if (boundaryMode) {
                    // Prefer the smallest run-boundary row within [gFirst, gLast];
                    // else fall back to the earliest min position.
                    thrPosVal = gFirst;
                    uint64_t lo = 0, hi = runs + 1;  // search runStarts[0..runs]
                    while (lo < hi) {
                        uint64_t mid = (lo + hi) / 2;
                        if (static_cast<uint64_t>(runStarts[mid]) < gFirst) lo = mid + 1;
                        else hi = mid;
                    }
                    if (lo <= runs) {
                        uint64_t cand = runStarts[lo];
                        if (cand >= gFirst && cand <= gLast) thrPosVal = cand;
                    }
                } else {
                    thrPosVal = gFirst;             // first-found (earliest)
                }
            }
            writeOne(thrVal, thrPosVal);

            // Run i now lies inside every other character's open gap; fold its
            // full min into their registers. Reset c's register: its next gap
            // begins after run i.
            uint64_t m = minLCP[i], mf = firstMinPos[i], ml = lastMinPos[i];
            for (uint64_t d = 0; d <= maxSym; ++d) {
                if (d == c) continue;
                if (m < regMin[d]) { regMin[d] = m; regFirst[d] = mf; regLast[d] = ml; }
                else if (m == regMin[d]) {
                    if (mf < regFirst[d]) regFirst[d] = mf;
                    if (ml > regLast[d]) regLast[d] = ml;
                }
            }
            regMin[c] = SENT; regFirst[c] = 0; regLast[c] = 0;
        }
    }

    /**
     * [Tier B] Streaming, O(sigma)-RAM variant of writeThresholdsSweep for the
     * default (non-boundary) mode. The per-run summary was spilled to spillPath
     * during the reverse reordering pass (records in reverse run order: run runs-1
     * first). We read it back in forward run order -- reading fixed-size records
     * from the end of the file toward the start -- and run the identical sticky-
     * register sweep, never holding the run-sized summary arrays in RAM. Boundary
     * mode is unsupported here (it needs random access to runStarts); the caller
     * routes boundary requests to the RAM path.
     */
    template<typename SymArr>
    void writeThresholdsSweepSpilled(const std::string& basePath, const SymArr& symbols,
                                     uint64_t runs, const std::string& spillPath,
                                     unsigned width, bool legacy) const {
        if (runs == 0) return;
        const uint64_t SENT = static_cast<uint64_t>(-1);

        uint64_t maxSym = 0;
        for (uint64_t i = 0; i < runs; ++i) maxSym = std::max<uint64_t>(maxSym, symbols[i]);
        std::vector<uint64_t> regMin(maxSym + 1, SENT),
                              regFirst(maxSym + 1, 0), regLast(maxSym + 1, 0);
        std::vector<char> seenChar(maxSym + 1, 0);

        ThrFileWriter w(basePath, runs, width, legacy);
        auto writeOne = [&](uint64_t thr, uint64_t pos) { w.writeOne(thr, pos); };

        // Backward block reader: records are in reverse run order, so we load blocks
        // from the end of the file and, within a block, consume highest file index
        // (== lowest run index) first to yield forward run order.
        std::ifstream in(spillPath, std::ios::binary);
        if (!in) throw std::runtime_error("writeThresholds: failed to reopen spill file");
        const uint64_t RECSZ = sizeof(SpillRec);
        const uint64_t BLK = 1u << 16;                 // records per block
        std::vector<SpillRec> blk(BLK);
        uint64_t fileRemaining = runs;                 // file records not yet loaded: [0, fileRemaining)
        uint64_t blkPos = 0;                           // consume blk[--blkPos] next
        auto nextRec = [&](SpillRec& out) {
            if (blkPos == 0) {
                uint64_t take = std::min<uint64_t>(BLK, fileRemaining);
                uint64_t startRec = fileRemaining - take;   // load file records [startRec, fileRemaining)
                in.seekg(static_cast<std::streamoff>(startRec * RECSZ), std::ios::beg);
                in.read(reinterpret_cast<char*>(blk.data()),
                        static_cast<std::streamsize>(take * RECSZ));
                fileRemaining = startRec;
                blkPos = take;
            }
            out = blk[--blkPos];
        };

        SpillRec rec;
        for (uint64_t i = 0; i < runs; ++i) {
            nextRec(rec);
            uint64_t c = symbols[i];
            uint64_t thrVal, thrPosVal;
            if (!seenChar[c]) {
                thrVal = 0; thrPosVal = 0;          // first run of c: placeholder
                seenChar[c] = 1;
            } else {
                thrVal = regMin[c];
                uint64_t gFirst = regFirst[c];
                uint64_t topi = rec.topLCP, starti = rec.runStart;
                if (topi < thrVal) { thrVal = topi; gFirst = starti; }
                if (thrVal == SENT) { thrVal = 0; thrPosVal = 0; }  // gap empty
                else                { thrPosVal = gFirst; }         // earliest min (non-boundary)
            }
            writeOne(thrVal, thrPosVal);

            uint64_t m = rec.minLCP, mf = rec.firstMin, ml = rec.lastMin;
            for (uint64_t d = 0; d <= maxSym; ++d) {
                if (d == c) continue;
                if (m < regMin[d]) { regMin[d] = m; regFirst[d] = mf; regLast[d] = ml; }
                else if (m == regMin[d]) {
                    if (mf < regFirst[d]) regFirst[d] = mf;
                    if (ml > regLast[d]) regLast[d] = ml;
                }
            }
            regMin[c] = SENT; regFirst[c] = 0; regLast[c] = 0;
        }
    }

    /**
     * Non-destructive path: compute thresholds from the LCP index via a single
     * phi walk that gathers a per-run LCP summary (top_lcp, min, first/last
     * argmin) in O(n) time / O(r) space, then emit .thr/.thr_pos with the O(r)
     * sticky-register sweep. writeThresholdsParallel is the faster default.
     *
     * A threshold for run i of character c is the minimum LCP over the
     * BWT-row gap [runStarts[prev_c+1], runStarts[i]] (rows after the previous run
     * of c, up to and including run i's first row); that range-minimum decomposes
     * into the per-run minima this gathers.
     */
    void writeThresholdsImpl(const std::string& basePath, bool boundaryMode, const RunInfo& runInfo,
                             unsigned width, bool legacy) const {
        const uint64_t runs = F.size();
        if (runs == 0) return;

        const std::vector<uint64_t>& runLen = runInfo.lengths;
        const std::vector<uint64_t>& symbols = runInfo.symbols;
        if (symbols.size() != runs || runLen.size() != runs)
            throw std::runtime_error("writeThresholds: run info size mismatch");

        // runStarts[r] = BWT row of run r's first element; runStarts[runs] == totalLen.
        std::vector<uint64_t> runStarts(runs + 1);
        runStarts[0] = 0;
        for (uint64_t r = 0; r < runs; ++r)
            runStarts[r + 1] = runStarts[r] + runLen[r];

        const uint64_t SENT = static_cast<uint64_t>(-1);

        // Per-run LCP summary, gathered below in one phi walk:
        //   topLCP[r]      : LCP at the first row of run r (row runStarts[r])
        //   minLCP[r]      : minimum LCP over all rows of run r
        //   firstMinPos[r] : smallest BWT row of run r achieving minLCP[r]
        //   lastMinPos[r]  : largest  BWT row of run r achieving minLCP[r]
        std::vector<uint64_t> topLCP(runs, 0), minLCP(runs, SENT),
                              firstMinPos(runs, 0), lastMinPos(runs, 0);

        // phiWalkLCP streams (pos, val) for pos = totalLen-1 down to 0, i.e. runs
        // in decreasing index order; we track the current run accordingly.
        {
            uint64_t curRun = runs - 1;
            phiWalkLCP([&](uint64_t pos, uint64_t val) {
                while (pos < runStarts[curRun]) --curRun;
                if (val < minLCP[curRun]) {
                    minLCP[curRun] = val;
                    firstMinPos[curRun] = lastMinPos[curRun] = pos;
                } else if (val == minLCP[curRun]) {
                    if (pos < firstMinPos[curRun]) firstMinPos[curRun] = pos;
                    if (pos > lastMinPos[curRun]) lastMinPos[curRun] = pos;
                }
                if (pos == runStarts[curRun]) topLCP[curRun] = val;  // first row of run
            });
        }

        writeThresholdsSweep(basePath, boundaryMode, symbols, runStarts,
                             topLCP, minLCP, firstMinPos, lastMinPos, width, legacy);
    }

public:
    /**
     * Writes .thr and .thr_pos binary files (format per BenLangmead/thresh-tools) at
     * basePath.thr and basePath.thr_pos. For each run i of character c, thr[i] is
     * the minimum LCP in the gap between the previous run of c and run i, and
     * thr_pos[i] is a BWT row achieving it; the first-ever run of c writes 0,0.
     *
     * boundaryMode: when true, prefer a thr_pos at a run boundary within the gap's
     * [first_min, last_min] range; otherwise use the earliest min position.
     *
     * runInfo must come from runInfoFromFMD() on the FMD matching this index. This
     * is the non-destructive phi-walk path; writeThresholdsParallel is faster.
     */
    void writeThresholds(const std::string& basePath, bool boundaryMode, const RunInfo& runInfo,
                         bool legacy, unsigned widthOverride) const {
        if (runInfo.lengths.size() != F.size()) {
            throw std::runtime_error("writeThresholds: FMD run count (" + std::to_string(runInfo.lengths.size())
                + ") != lcp_index run count (" + std::to_string(F.size())
                + "). FMD and lcp_index may be from different inputs.");
        }
        uint64_t fmdTotal = std::accumulate(runInfo.lengths.begin(), runInfo.lengths.end(), 0ULL);
        if (fmdTotal != totalLen) {
            throw std::runtime_error("writeThresholds: FMD total length (" + std::to_string(fmdTotal)
                + ") != lcp_index total length (" + std::to_string(totalLen)
                + "). FMD and lcp_index may be from different inputs.");
        }
        if (legacy && totalLen >= (1ULL << 40))
            throw std::runtime_error("pfp-thresholds-style 5-byte thresholds cannot represent BWT length >= 2^40 ("
                + std::to_string(totalLen) + "); omit --thr-pfp for the wide self-describing format");
        writeThresholdsImpl(basePath, boundaryMode, runInfo, thrFieldWidth(legacy, widthOverride), legacy);
    }

    /**
     * Fast default path: compute thresholds by piggybacking on the parallel
     * per-run traversal (ComputeMinLCPRunParallelDestructive), then the O(r)
     * sweep. DESTRUCTIVE (frees Phi/Psi/PLCPsamples), matching build-time use
     * where the index is consumed once. Produces identical .thr/.thr_pos to the
     * non-destructive writeThresholds().
     */
    void writeThresholdsParallel(const std::string& basePath, bool boundaryMode,
                                 const RunInfo& runInfo
                                 #ifndef BENCHFASTONLY
                                 , const verbosity v
                                 #endif
                                 , bool legacy, unsigned widthOverride) {
        const uint64_t runs = F.size();
        if (runs == 0) return;
        if (runInfo.symbols.size() != runs || runInfo.lengths.size() != runs)
            throw std::runtime_error("writeThresholds: run info size mismatch");
        if (legacy && totalLen >= (1ULL << 40))
            throw std::runtime_error("pfp-thresholds-style 5-byte thresholds cannot represent BWT length >= 2^40 ("
                + std::to_string(totalLen) + "); omit --thr-pfp for the wide self-describing format");
        const unsigned width = thrFieldWidth(legacy, widthOverride);
        if (boundaryMode) {
            // Boundary mode needs random access to runStarts (binary search in the
            // sweep); keep the proven RAM-summary path for it.
            RunSummary summary;
            ComputeMinLCPRunParallelDestructive(
                #ifndef BENCHFASTONLY
                v,
                #endif
                &summary);
            writeThresholdsSweep(basePath, true, runInfo.symbols, summary.runStarts,
                                 summary.topLCP, summary.minLCP,
                                 summary.firstMinPos, summary.lastMinPos, width, legacy);
            return;
        }
        // [Tier B] Default mode: spill the per-run summary to disk during the
        // reordering pass (so the 5 run-sized arrays never materialize), then stream
        // it back in forward run order for the sweep. Output is byte-identical to the
        // RAM path; peak RAM drops by the summary arrays (~80 GB at r~3.5e9).
        const std::string spillPath = basePath + ".thrspill.tmp";
        {
            std::ofstream spill(spillPath, std::ios::binary);
            if (!spill.is_open())
                throw std::runtime_error("writeThresholds: cannot open spill file " + spillPath);
            std::vector<char> sbuf(8 * 1024 * 1024);
            spill.rdbuf()->pubsetbuf(sbuf.data(), sbuf.size());
            ComputeMinLCPRunParallelDestructive(
                #ifndef BENCHFASTONLY
                v,
                #endif
                nullptr, &spill);
        }
        writeThresholdsSweepSpilled(basePath, runInfo.symbols, runs, spillPath, width, legacy);
        std::remove(spillPath.c_str());
    }

    /**
     * Writes the run-length BWT as its companion files basePath.bwt.heads (one
     * char per run) and basePath.bwt.len (5-byte little-endian length per run).
     * ropeBWT3 FMD symbol codes (0=$,1=A,2=C,3=G,4=T,5=N) are mapped to ASCII; the
     * sentinel is written as a 0 byte, matching pfp-thresholds.
     *
     * This transcoding is independent of LCP/thresholds, but it lives here on
     * purpose: the heads/len files and the .thr/.thr_pos files must describe the
     * *same* run-length BWT -- identical run partition (here, runInfoFromFMD's
     * sentinel handling), run order, and alphabet assumptions -- or a downstream
     * consumer's per-run arrays silently misalign. Producing them together from
     * one RunInfo, under a single set of assumptions about the RLBWT shape and
     * separators, removes that whole class of mismatch by construction. A
     * threshold-based matching index is the first consumer of these files, but
     * not necessarily the only one.
     */
    void writeBwtHeadsLen(const std::string& basePath, const RunInfo& runInfo) const {
        static const char ALPHA[6] = {'\0', 'A', 'C', 'G', 'T', 'N'};
        std::ofstream heads(basePath + ".bwt.heads", std::ios::binary);
        std::ofstream lens(basePath + ".bwt.len", std::ios::binary);
        if (!heads.is_open() || !lens.is_open())
            throw std::runtime_error("writeBwtHeadsLen: failed to open output files");
        const uint64_t runs = runInfo.symbols.size();
        if (runInfo.lengths.size() != runs)
            throw std::runtime_error("writeBwtHeadsLen: run info size mismatch");
        char b[5];
        for (uint64_t i = 0; i < runs; ++i) {
            uint64_t s = runInfo.symbols[i];
            heads.put(s < 6 ? ALPHA[s] : '?');
            uint64_t L = runInfo.lengths[i];
            std::memcpy(b, &L, 5);
            lens.write(b, 5);
        }
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

bool TeraLCP::validateRB3(const rb3_fmi_t* rb3){
    if (!rb3->e) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return false;
    }
    return true;
}

#endif //#ifndef R_SA_LCP_TeraLCP_H
