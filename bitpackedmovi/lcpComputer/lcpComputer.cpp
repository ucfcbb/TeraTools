#include"optbwtrl/lcpComputer.h"

int main(int argc, char*argv[]) {
    Timer.start("lcpComputer");
    Timer.start("Program Initialization");
    rb3_fmi_t fmi;
    std::ofstream indOut, treeOut, lcpOut;
    std::string outputPref;
    {
        bool use_mmap;
        Timer.start("Reading Arguments");
        //const char characters[7] = "$ACGTN";
        if (argc != 3 && argc != 4) {
            std::cerr << "Usage: lcpComputer [-mmap] <input.fmd> <outputPrefix>\n\nIf -mmap is passed, the file is memory mapped before reading, otherwise traditional file io is used\n";
            std::cerr << argc-1 << " arguments passed instead of 1 or 2" << std::endl;
            exit(1);
        }
        if (argc == 4 && strcmp(argv[1], "-mmap") !=0) {
            std::cerr << "If three arguments are passed, the first must be \"-mmap\"\n";
            std::cerr << "Currently passed: \"" << argv[1] << "\"\n";
            exit(1);
        }
        std::string inputfmd = argv[argc-2];
        outputPref = argv[argc-1];
        use_mmap = (argc == 4);

        Timer.stop(); //Reading Arguments
        Timer.start((use_mmap)? "Loading fmd with mmap" : "Loading fmd");

        rb3_fmi_restore(&fmi, inputfmd.c_str(), use_mmap);
        if (fmi.e == 0 && fmi.r == 0) {
            std::cerr << "ERROR: failed to load fmd from index file " << inputfmd << std::endl;
            exit(1);
        }

        Timer.stop(); //(use_mmap)? "Loading fmd with mmap" : "Loading fmd"
        
        if (!LCPComputer::validateRB3(&fmi)) {
            std::cerr << "ERROR: invalid ropebwt3 inputted!" << std::endl;
            exit(1);
        }

        indOut.open(outputPref + ".optbwtrl");
        if (!indOut.is_open()) {
            std::cerr << "ERROR: File '" << outputPref << ".optbwtrl' failed to open for writing!\n";
            exit(1);
        }
        treeOut.open(outputPref + "_StructTree.html");
        if (!treeOut.is_open()) {
            std::cerr << "ERROR: File '" << outputPref << "_StructTree.html' failed to open for writing!\n";
            exit(1);
        }
        lcpOut.open(outputPref + ".rlcp");
        if (!lcpOut.is_open()) {
            std::cerr << "ERRPR: File '" << outputPref << ".rlcp' failed to open for writing!\n";
            exit(1);
        }
    }
    Timer.stop(); //Program Initialization

    Timer.start("LCP construction");
    LCPComputer ourIndex(&fmi, outputPref + ".TEMPFILE101Y90HDOSJFN", lcpOut);
    Timer.stop(); //LCP construction

    /*
    Timer.start("Printing Raw");
    ourIndex.printRaw();
    Timer.stop(); //Printing Raw
    */

    Timer.start("Measure size");
    std::cout << "Size of our index: " << sdsl::size_in_bytes(ourIndex) << std::endl;
    Timer.stop(); //Measure size

    Timer.start("Writing Index");
    ourIndex.serialize(indOut);
    indOut.close();
    Timer.stop(); //Writing Index

    Timer.start("Writing Structure Tree");
    sdsl::write_structure<sdsl::HTML_FORMAT>(ourIndex, treeOut);
    treeOut.close();
    Timer.stop(); //Writing Structure Tree

    Timer.stop(); //lcpComputer
}
