#include"lcpComputer/lcpComputer.h"

void printUsage() {
    std::cout << 
        "lcpComputer computes the minimum LCP values in each run of the BWT of a text.\n"
        "It does so by computing the Psi and Phi move data structures in compressed space.\n"
        "It can optionally output the computed index to avoid later recomputation or for downstream analysis.\n\n"
        "Usage: lcpComputer <arguments>\n"
        "Options:\n"
        "  Input:\n"
        "    -f       [text,bwt,rlbwt,fmd]       REQUIRED       Format of input. 'text' is the original text, 'bwt' is the bwt of the text, 'rlbwt' is the rlbwt of the text, and 'fmd' is the rlbwt in the ropebwt3 fmd format of the text\n"
        "    -i       FILE                       REQUIRED       File name of input file\n"
        "    -t       FILE                       REQUIRED       Name of a file this program can read and write to temporarily\n\n"
        "  Output:\n"
        "    -oindex  FILE                       optional       Output constructed index to FILE\n"
        "    -orlcp   FILE                       optional       Output (position, minLCP) pairs per run to FILE\n\n"
        "  Behavior:\n"
        "    -p       INT                        optional       Limit the program to (nonnegative) INT threads. By default uses maximum available. Maximum on this hardware is " << omp_get_max_threads() << "\n"
        "    -mmap                               optional       read input using memory mapping (only avaiable for fmd) default: no memory mapping\n"
        #ifndef BENCHFASTONLY
        "    -v       [quiet,time,verb]          optional       Verbosity, verb for most verbose output, time for timer info, and quiet for no output. time is default.\n"
        #else
        "    -bench                              optional       Has no effect. required if no outputs specified.\n"
        #endif
        "    -h, --help                          optional       Print this help message.\n"
        ;
    //add verification of psi and phi options
    //add sdsl::memory_monitor output option
    //add timer depth option? or remove many timer calls.
    //add move data structure output options?
    //TODO: either fix minor bug in lcpcomputer when numthreads > max threads or limit num threads to at most max threads
    //TODO: add more BENCHFASTONLY ifs
    //TODO: add DNDEBUG ifs
}

struct options{
    enum inputFormat { text, bwt, rlbwt, fmd }inputFormat;
    std::string inputFile, tempFile, oindex="", orlcp="";
    unsigned numThreads = omp_get_max_threads();
    bool mmap;
    #ifndef BENCHFASTONLY
    verbosity v = QUIET;
    #endif
}o;

void processOptions(const int argc, const char* argv[]) {
    std::vector<bool> used(argc);
    used[0] = true;
    auto getArg = [&] (std::string arg, bool required, bool argument) -> std::string {
        auto pos = std::find(argv, argv+argc, arg);
        if (!required && pos == argv+argc)
            return "";
        if (required && pos == argv+argc) {
            std::cout << arg << " was not passed as an argument, but it is required.\n";
            exit(1);
        }
        if (argument) {
            if(pos != argv+argc-1) {
                used[pos - argv] = true;
                used[pos - argv + 1] = true;
                return *(pos+1);
            }
            std::cout << arg << " was not passed an argument, but it requires one.\n";
            exit(1);
        }
        used[pos - argv] = true;
        return *pos;
    };

    auto testOutFile = [] (std::string name) {
        if (name == "") return;
        std::ofstream out(name, std::ios::app);
        if (!out.is_open()) {
            std::cout << "File '" << name << "' failed to open for writing.\n";
            exit(1);
        }
    };
    auto testInFile = [] (std::string name) {
        if (name == "") return;
        std::ofstream in(name, std::ios::app);
        if (!in.is_open()) {
            std::cout << "File '" << name << "' failed to open for writing.\n";
            exit(1);
        }
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
    else {
        std::cout << "Invalid value passed to -f '" << s << "'\n";
        exit(1);
    }
    o.inputFile = getArg("-i", true, true);
    o.tempFile = getArg("-t", true, true);
    o.oindex = getArg("-oindex", false, true);
    o.orlcp = getArg("-orlcp", false, true);
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
    testInFile(o.inputFile);
    testInFile(o.tempFile);
    testOutFile(o.tempFile);
    testOutFile(o.oindex);
    testOutFile(o.orlcp);
}

int main(const int argc, const char*argv[]) {
    processOptions(argc, argv);
    omp_set_num_threads(o.numThreads);
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.start("lcpComputer"); }
    if (o.v >= TIME) { Timer.start("Program Initialization"); }
    #endif
    rb3_fmi_t fmi;
    {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("Reading Arguments"); }
        #endif
        if (o.inputFormat != options::fmd) {
            std::cerr << "Only fmd input currently implemented!" << std::endl;
            exit(1);
        }
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.stop(); } //Reading Arguments 
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
        
        if (!LCPComputer::validateRB3(&fmi)) {
            std::cerr << "ERROR: invalid ropebwt3 inputted!" << std::endl;
            exit(1);
        }
    }
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.stop(); } //Program Initialization 

    if (o.v >= TIME) { Timer.start("LCP construction"); }
    #endif
    LCPComputer ourIndex(&fmi, o.tempFile
            #ifndef BENCHFASTONLY
            , o.v
            #endif
            );
    #ifndef BENCHFASTONLY
    if (o.v >= TIME) { Timer.stop(); } //LCP construction 
    #endif


    if (o.orlcp != "") {
        #ifndef BENCHFASTONLY
        if (o.v >= TIME) { Timer.start("min LCP per run computation"); }
        #endif
        std::ofstream lcpOut(o.orlcp);
        if (!lcpOut.is_open()) {
            std::cerr << "ERROR: File '" << o.orlcp << "' failed to open for writing!\n";
            exit(1);
        }
        ourIndex.ComputeMinLCPRun(lcpOut, o.v);
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
    if (o.v >= TIME) { Timer.stop(); } //lcpComputer 
    #endif
}
