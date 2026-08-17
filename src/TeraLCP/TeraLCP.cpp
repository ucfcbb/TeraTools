#include"util/util.h"
#include"TeraLCP/TeraLCP.h"
#include<sys/stat.h>
#include<iterator>

static constexpr const char* rlcp_extension = ".rlcp";

void printUsage() {
    std::cout << 
        "TeraLCP computes the minimum LCP values in each run of the BWT of a text.\n"
        "It does so by computing the Psi and Phi move data structures in compressed space.\n"
        "It can optionally output the computed index to avoid later recomputation or for downstream analysis.\n"
        "\n"
        "Usage: TeraLCP <arguments>\n"
        "Options:\n"
        "  Input:\n"
        "    -f          [text,bwt,rlbwt,fmd,lcp_index]  REQUIRED       Format of input. 'text' is the original text, 'bwt' is the bwt of the text, 'rlbwt' is a run-length BWT given as <name>.bwt.heads (1 byte per run head) and <name>.bwt.len (5-byte little-endian run lengths), 'fmd' is the rlbwt in the ropebwt3 fmd format of the text, and 'lcp_index' is the index outputted by TeraLCP -oindex\n"
        "    -i          FILE                            REQUIRED       File name of input file\n"
        "    -t          FILE                            REQUIRED       Name of a file this program can read and write to temporarily\n"
        "\n"
        "  Output:\n"
        "    -oindex     FILE                            optional       Output constructed index to FILE" << lcp_index_extension << "\n"
        "    -orlcp      FILE                            optional       Output (position, minLCP) pairs per run to FILE" << rlcp_extension << "\n"
        "    -othresholds BASE                           optional       Output pfp-thresholds-style thresholds to BASE.thr and BASE.thr_pos,\n"
        "                                                              plus run-length BWT files BASE.bwt.heads and BASE.bwt.len\n"
        "    -threshbound                                optional       (-othresholds) prefer a run-boundary row for thr_pos within the min range\n"
        "    --thr-pfp                                   optional       (-othresholds) write pfp-thresholds-style headerless 5-byte threshold files (requires n < 2^40)\n"
        "    --thr-width N                               optional       (-othresholds) force threshold field width in bytes (default: auto from n)\n"
        "    --fmd       FILE                            optional       (-othresholds with -f lcp_index) FMD matching the index, for run metadata\n"
        "    --rlbwt-meta BASE                           optional       (-othresholds with -f lcp_index) grlBWT/rlbwt heads/len (BASE.bwt.heads/.bwt.len)\n"
        "                                                              for run metadata; use instead of --fmd to resume the separator pipeline\n"
        "\n"
        "  Checkpoint/resume (construct, -f fmd):\n"
        "    -checkpoint DIR                             optional       Read/write per-phase checkpoints in DIR. On startup, resume from the\n"
        "                                                                 first incomplete phase; on a fresh run, write a checkpoint after each\n"
        "                                                                 slow phase so a killed job can resume. DIR is created if missing.\n"
        "    -stop-after [A,B]                           optional       (requires -checkpoint) Stop after phase A (endmarker repair) or B (Phi+samples),\n"
        "                                                                 writing its checkpoint and exiting 0. Lets a >7-day build be chunked.\n"
        "\n"
        "  Behavior:\n"
        "    -p          INT                             optional       Limit the program to (nonnegative) INT threads. By default uses maximum available. Maximum on this hardware is " << omp_get_max_threads() << "\n"
        "    -mmap                                       optional       read input using memory mapping (only avaiable for fmd) default: no memory mapping\n"
        #ifndef BENCHFASTONLY
        "    -v          [quiet,time,verb]               optional       Verbosity, verb for most verbose output, time for timer info, and quiet for no output. time is default.\n"
        #else
        "    -bench                                      optional       Has no effect. required if no outputs specified.\n"
        #endif
        "    -h, --help                                  optional       Print this help message.\n"
        ;
    //add verification of psi and phi options
    //add sdsl::memory_monitor output option
    //add timer depth option? or remove many timer calls.
    //add move data structure output options?
    //TODO: either fix minor bug in TeraLCP when numthreads > max threads or limit num threads to at most max threads
    //TODO: add more BENCHFASTONLY ifs
    //TODO: add DNDEBUG ifs
}

struct options{
    enum inputFormat { text, bwt, rlbwt, fmd, lcp_index }inputFormat;
    std::string inputFile, tempFile, oindex="", orlcp="", othresholds="", fmdFile="", rlbwtMeta="";
    // Checkpoint/resume for construct: -checkpoint <dir> reads/writes per-phase
    // checkpoints; -stop-after A|B chunks a multi-week build across <=7-day jobs.
    std::string checkpointDir="";
    TeraLCP::StopAfter stopAfter = TeraLCP::StopAfter::NONE;

    bool threshbound = false;
    bool thrLegacy = false;      // --thr-pfp: emit old headerless 5-byte thresholds
    unsigned thrWidth = 0;       // --thr-width N: force field width (0 = auto from n)

    unsigned numThreads = omp_get_max_threads();
    bool mmap;
    #ifndef BENCHFASTONLY
    verbosity v = TIME;
    #endif
}o;

void processOptions(const int argc, const char* argv[]) {
    std::vector<bool> used(argc);
    used[0] = true;
    auto getArg = [&] (std::string arg, bool required, bool argument) -> std::string {
        return getArgument(argc, argv, used, arg, required, argument);
    };

    if (argc == 1 || getArg("-h", false, false) != "" || getArg("--help", false, false) != "") {
        printUsage();
        exit(0);
    }

    auto s = getArg("-f", true, true);
    if (s == "text") o.inputFormat=options::text;
    else if (s == "bwt") o.inputFormat=options::bwt;
    else if (s == "rlbwt") o.inputFormat=options::rlbwt;
    else if (s == "fmd") o.inputFormat=options::fmd;
    else if (s == "lcp_index") o.inputFormat=options::lcp_index;
    else {
        std::cout << "Invalid value passed to -f '" << s << "'\n";
        exit(1);
    }
    o.inputFile = getArg("-i", true, true);
    o.tempFile = getArg("-t", true, true);
    o.oindex = getArg("-oindex", false, true);
    if (o.oindex != "") 
        o.oindex += lcp_index_extension;
    o.orlcp = getArg("-orlcp", false, true);
    if (o.orlcp != "")
        o.orlcp += rlcp_extension;

    // Checkpoint/resume of construct across successive jobs.
    o.checkpointDir = getArg("-checkpoint", false, true);
    {
        std::string sa = getArg("-stop-after", false, true);
        if (sa == "") o.stopAfter = TeraLCP::StopAfter::NONE;
        else if (sa == "A" || sa == "a") o.stopAfter = TeraLCP::StopAfter::A;
        else if (sa == "B" || sa == "b") o.stopAfter = TeraLCP::StopAfter::B;
        else { std::cerr << "ERROR: -stop-after must be A or B\n"; exit(1); }
        if (o.stopAfter != TeraLCP::StopAfter::NONE && o.checkpointDir.empty()) {
            std::cerr << "ERROR: -stop-after requires -checkpoint <dir>\n"; exit(1);
        }
    }

    o.othresholds = getArg("-othresholds", false, true);
    o.threshbound = (getArg("-threshbound", false, false) != "");
    o.thrLegacy = (getArg("--thr-pfp", false, false) != "");
    { std::string tw = getArg("--thr-width", false, true); o.thrWidth = tw.empty() ? 0u : static_cast<unsigned>(std::stoul(tw)); }
    o.fmdFile = getArg("--fmd", false, true);
    // Run metadata for the -f lcp_index resume can instead come from grlBWT/rlbwt
    // companion files (BASE.bwt.heads/.bwt.len), so the rlbwt (separator) pipeline
    // can checkpoint via -oindex and resume thresholds without ever building an FMD.
    o.rlbwtMeta = getArg("--rlbwt-meta", false, true);

    s = getArg("-p", false, true);
    if (s != "")
        o.numThreads = std::stoul(s);
    o.mmap = ("-mmap" == getArg("-mmap", false, false));
#ifndef BENCHFASTONLY
    s = getArg("-v", false, true);
    if (s == "quiet")
        o.v = QUIET;
    else if (s == "time" || s == "")
        o.v = TIME;
    else if (s == "verb")
        o.v = VERB;
    else {
        std::cout << "Invalid value passed to -v '" << s << "'\n";
        exit(1);
    }
#else
    s = getArg("-bench", false, false);
    if (o.oindex == "" && o.orlcp == "" && s == "") {
        std::cout << "No output formats passed. If you want to construct the index but not output anything (for benchmarking purposes, typically), then '-bench' must be explictly passed.\n";
        exit(1);
    }
#endif
    for (int i = 0; i < argc; ++i)
        if (!used[i]) {
            std::cout << "Argument " << i << ", '" << argv[i] << "' not recognized or used as an argument for another option. It might have been passed more than once (invalid).\n";
            exit(1);
        }
    if (o.inputFormat == options::rlbwt) {
        // -i is a base path; the actual inputs are base.bwt.heads and base.bwt.len.
        testInFile(o.inputFile + ".bwt.heads");
        testInFile(o.inputFile + ".bwt.len");
    } else {
        testInFile(o.inputFile);
    }
    testOutFile(o.tempFile);
    testInFile(o.tempFile);
    testOutFile(o.oindex);
    testOutFile(o.orlcp);

    if (!o.checkpointDir.empty()) {
        // Create the checkpoint directory if needed (idempotent); ignore EEXIST.
        if (mkdir(o.checkpointDir.c_str(), 0777) != 0 && errno != EEXIST) {
            std::cerr << "ERROR: cannot create checkpoint directory '" << o.checkpointDir << "': " << strerror(errno) << "\n";
            exit(1);
        }
    }

    if (!o.fmdFile.empty())
        testInFile(o.fmdFile);
    if (!o.rlbwtMeta.empty()) {
        testInFile(o.rlbwtMeta + ".bwt.heads");
        testInFile(o.rlbwtMeta + ".bwt.len");
    }
}

// Builds an in-memory ropeBWT3 rld (FMD) from a run-length BWT given as its
// companion files base.bwt.heads (1 byte per run) and base.bwt.len (5-byte LE per
// run). Distinct head bytes are remapped to codes 0..k-1 in ascending byte order,
// so the smallest byte (e.g. a 0 sentinel) becomes code 0 and is treated as the
// endmarker by ConstructPsi/runInfoFromFMD -- matching the byte collation the BWT
// was built with (e.g. libsais/divsufsort on a separator-bearing text). This lets
// TeraLCP ingest BWTs over an alphabet ropeBWT3's DNA-only builder cannot produce
// (e.g. one carrying a '%' separator).
static rld_t* buildRldFromRlbwt(const std::string& base) {
    const std::string headsPath = base + ".bwt.heads", lenPath = base + ".bwt.len";
    std::ifstream hf(headsPath, std::ios::binary), lf(lenPath, std::ios::binary);
    if (!hf.is_open()) { std::cerr << "ERROR: cannot open " << headsPath << "\n"; std::exit(1); }
    if (!lf.is_open()) { std::cerr << "ERROR: cannot open " << lenPath << "\n"; std::exit(1); }
    std::vector<unsigned char> heads((std::istreambuf_iterator<char>(hf)),
                                      std::istreambuf_iterator<char>());
    const uint64_t runs = heads.size();
    if (runs == 0) { std::cerr << "ERROR: empty heads file " << headsPath << "\n"; std::exit(1); }
    std::vector<int> present(256, 0);
    for (unsigned char b : heads) present[b] = 1;
    std::vector<int> code(256, -1);
    int k = 0;
    for (int b = 0; b < 256; ++b) if (present[b]) code[b] = k++;
    if (k > RB3_ASIZE) {
        std::cerr << "ERROR: RLBWT alphabet size " << k << " exceeds RB3_ASIZE (" << RB3_ASIZE << ")\n";
        std::exit(1);
    }
    rld_t* e = rld_init(RB3_ASIZE, 3);
    rlditr_t ei;
    rld_itr_init(e, &ei, 0);
    for (uint64_t i = 0; i < runs; ++i) {
        uint64_t L = 0;
        lf.read(reinterpret_cast<char*>(&L), 5);
        if (!lf) { std::cerr << "ERROR: " << lenPath << " has fewer than " << runs << " run lengths\n"; std::exit(1); }
        rld_enc(e, &ei, static_cast<int64_t>(L), static_cast<uint8_t>(code[heads[i]]));
    }
    rld_enc_finish(e, &ei);
    return e;
}

int main(const int argc, const char*argv[]) {
    processOptions(argc, argv);
    omp_set_num_threads(o.numThreads);
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.start("TeraLCP"); }
    if (o.v >= TIME) { Timer.start("Program Initialization"); }
    #endif
    rb3_fmi_t fmi;
    {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("Reading Arguments"); }
        #endif
        if (o.inputFormat != options::fmd && o.inputFormat != options::rlbwt
                && o.inputFormat != options::lcp_index) {
            std::cerr << "Only fmd, rlbwt, and lcp_index input currently implemented!" << std::endl;
            exit(1);
        }
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //Reading Arguments 
        #endif
        if (o.inputFormat == options::fmd) {
            #ifndef BENCHFASTONLY
            if (o.v >= TIME) { Timer.start((o.mmap)? "Loading fmd with mmap" : "Loading fmd"); }
            #endif

            rb3_fmi_restore(&fmi, o.inputFile.c_str(), o.mmap);
            if (fmi.e == 0 && fmi.r == 0) {
                std::cerr << "ERROR: failed to load fmd from index file " << o.inputFile << std::endl;
                exit(1);
            }

            #ifndef BENCHFASTONLY
            if (o.v >= TIME) { Timer.stop(); } //(o.mmap)? "Loading fmd with mmap" : "Loading fmd" 
            #endif

            if (!TeraLCP::validateRB3(&fmi)) {
                std::cerr << "ERROR: invalid ropebwt3 inputted!" << std::endl;
                exit(1);
            }
        }
        else if (o.inputFormat == options::rlbwt) {
            #ifndef BENCHFASTONLY
            if (o.v >= TIME) { Timer.start("Building FMD from RLBWT heads/len"); }
            #endif
            rb3_fmi_init(&fmi, buildRldFromRlbwt(o.inputFile), 0);
            fmi.ssa = 0; fmi.sid = 0;
            #ifndef BENCHFASTONLY
            if (o.v >= TIME) { Timer.stop(); } //Building FMD from RLBWT heads/len
            #endif
            if (!TeraLCP::validateRB3(&fmi)) {
                std::cerr << "ERROR: failed to build FMD from RLBWT '" << o.inputFile << "'" << std::endl;
                exit(1);
            }
        }
    }
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.stop(); } //Program Initialization 
    #endif

    TeraLCP ourIndex;
    if (o.inputFormat == options::fmd || o.inputFormat == options::rlbwt) {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("LCP index construction"); }
        #endif
        ourIndex = TeraLCP(&fmi, o.tempFile
                #ifndef BENCHFASTONLY
                , o.v
                #endif
                , o.checkpointDir, o.stopAfter);
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //LCP index construction
        #endif
        // A -stop-after run wrote a checkpoint and returned without a full index;
        // there is nothing to emit, so finish cleanly (exit 0) for the driver.
        if (!ourIndex.complete) {
            #ifndef BENCHFASTONLY
            if (o.v >= TIME) { Timer.stop(); } //TeraLCP
            #endif
            return 0;
        }
    }
    else if (o.inputFormat == options::lcp_index) {
        ourIndex = TeraLCP(o.inputFile, o.v);
    }

    // pfp-thresholds-style thresholds (BASE.thr, BASE.thr_pos) plus the run-length BWT
    // companion files (BASE.bwt.heads, BASE.bwt.len). Run metadata is read from
    // the FMD (the index constructor frees the one used for construction, so we
    // reload it here). The fast path fuses threshold capture into the parallel
    // per-run pass and is destructive; fall back to the non-destructive phi-walk
    // path when -oindex or -orlcp still needs the in-memory index.
    if (o.othresholds != "") {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("thresholds output"); }
        #endif
        TeraLCP::RunInfo runInfo;
        // Run metadata (per-run symbols + lengths) comes from the rlbwt heads/len
        // when the input is rlbwt, OR when resuming from a saved lcp_index and the
        // caller supplied --rlbwt-meta BASE (the grlBWT/separator checkpoint path).
        const bool useRlbwtMeta = (o.inputFormat == options::rlbwt) || !o.rlbwtMeta.empty();
        if (useRlbwtMeta) {
            // Rebuild the rld from the heads/len files (the constructor freed the one
            // used for construction) to read run metadata.
            const std::string rlBase = (o.inputFormat == options::rlbwt) ? o.inputFile : o.rlbwtMeta;
            rb3_fmi_t fmiThr;
            rb3_fmi_init(&fmiThr, buildRldFromRlbwt(rlBase), 0);
            fmiThr.ssa = 0; fmiThr.sid = 0;
            runInfo = TeraLCP::runInfoFromFMD(&fmiThr);
            rb3_fmi_free(&fmiThr);
        } else {
            std::string fmdPath = (o.inputFormat == options::fmd) ? o.inputFile : o.fmdFile;
            if (fmdPath.empty()) {
                std::cerr << "ERROR: -othresholds with -f lcp_index requires --fmd FILE (FMD matching the index) or --rlbwt-meta BASE (grlBWT heads/len)\n";
                exit(1);
            }
            rb3_fmi_t fmiThr;
            rb3_fmi_restore(&fmiThr, fmdPath.c_str(), o.mmap);
            if (fmiThr.e == nullptr && fmiThr.r == nullptr) {
                std::cerr << "ERROR: failed to load FMD from '" << fmdPath << "' for thresholds\n";
                exit(1);
            }
            if (!TeraLCP::validateRB3(&fmiThr)) {
                std::cerr << "ERROR: invalid FMD (multirope or corrupted) for thresholds\n";
                rb3_fmi_free(&fmiThr);
                exit(1);
            }
            runInfo = TeraLCP::runInfoFromFMD(&fmiThr);
            rb3_fmi_free(&fmiThr);
        }
        const bool needIndexLater = (o.oindex != "" || o.orlcp != "");
        try {
            if (needIndexLater)
                ourIndex.writeThresholds(o.othresholds, o.threshbound, runInfo, o.thrLegacy, o.thrWidth);
            else
                ourIndex.writeThresholdsParallel(o.othresholds, o.threshbound, runInfo, o.v, o.thrLegacy, o.thrWidth);
            // Emit the run-length BWT companion files too, except when the run
            // metadata came from rlbwt heads/len (-f rlbwt or --rlbwt-meta): those
            // files already exist, and writeBwtHeadsLen's DNA code->ASCII map would
            // not match a remapped alphabet (e.g. one carrying '%').
            if (!useRlbwtMeta)
                ourIndex.writeBwtHeadsLen(o.othresholds, runInfo);
        } catch (const std::exception& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            exit(1);
        }
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //thresholds output
        #endif
    }


    if (o.orlcp != "") {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("min LCP per run computation"); }
        #endif
        std::ofstream lcpOut(o.orlcp);
        if (!lcpOut.is_open()) {
            std::cerr << "ERROR: File '" << o.orlcp << "' failed to open for writing!\n";
            exit(1);
        }
        auto l = ourIndex.ComputeMinLCPRunParallelDestructive(o.v);
        assert(l.first.size() == l.second.size());
        uint64_t runs = l.first.size();
        if (o.v >= TIME) { Timer.start("sequential output min LCP per run"); }
        for (uint64_t i = 0; i < runs; ++i) 
            lcpOut << "( " << l.first[i] << ", " << l.second[i] << ")\n";
        if (o.v >= TIME) { Timer.stop(); } //sequential output min LCP per run
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //min LCP per run computation
        #endif
    }

    /*
    if (o.v >= TIME) { Timer.start("Printing Raw"); }
    ourIndex.printRaw();
    if (o.v >= TIME) { Timer.stop(); } //Printing Raw 
    */

    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.start("Measure size"); }
    if (o.v >= VERB) { std::cout << "Size of our index: " << sdsl::size_in_bytes(ourIndex) << std::endl; }
    if (o.v >= TIME) { Timer.stop(); } //Measure size 
    #endif


    if (o.oindex != "") {
        std::ofstream indOut;
        indOut.open(o.oindex);
        if (!indOut.is_open()) {
            std::cerr << "ERROR: File '" << o.oindex << ".optbwtrl' failed to open for writing!\n";
            exit(1);
        }
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("Writing Index"); }
        #endif
        ourIndex.serialize(indOut);
        indOut.close();
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //Writing Index 
        #endif
    }

    /*
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.start("Writing Structure Tree"); }
    #endif
    sdsl::write_structure<sdsl::HTML_FORMAT>(ourIndex, treeOut);
    treeOut.close();
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.stop(); } //Writing Structure Tree 
    #endif
    */

    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.stop(); } //TeraLCP 
    #endif
}
