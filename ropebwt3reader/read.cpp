#include<iostream>
//#include<fstream>
//#include<string>
#include"fm-index.h"

int main(int argc, char *argv[]) {
    const char characters[7] = "$ACGTN";
    if (argc != 2 && argc != 3) {
        std::cerr << "Usage: read [-mmap] <input.fmd>\n\nIf -mmap is passed, the file is memory mapped before reading, otherwise traditional file io is used\n";
        std::cerr << argc-1 << " arguments passed instead of 1 or 2" << std::endl;
        return 1;
    }
    if (argc == 3 && strcmp(argv[1], "-mmap") !=0) {
        std::cerr << "If two arguments are passed, the first must be \"-mmap\"\n";
        std::cerr << "Currently passed: \"" << argv[1] << "\"\n";
        return 1;
    }
    char* inputfmd = argv[argc-1];

    const int use_mmap = (argc == 3);
    rb3_fmi_t fmi;
    rb3_fmi_restore(&fmi, inputfmd, use_mmap);
    if (fmi.e == 0 && fmi.r == 0) {
        std::cerr << "ERROR: failed to load fmd from index file " << inputfmd << std::endl;
        return 1;
    }

    /*
     * sends to cout: 
     * # runs
     * # bits needed for alphabet
     * run length encoded runs, each run is first the character, then the length
    */

    if (fmi.e) {
        std::cout << (uint64_t)rb3_fmi_get_r(&fmi) << ' ' << fmi.e->abits << ' ';
        int64_t l;
        int c = 0;
        rlditr_t itr;
        rld_itr_init(fmi.e, &itr, 0); //what does 0 mean in this function call? offset number of bits to start reading at?
        //std::string s;
        while ((l = rld_dec(fmi.e, &itr, &c, 0)) > 0) {
            //std::cerr << "run of length " << l << " of the character " << c << " equals " << characters[c] << std::endl;
            //s += std::string(l, characters[c]);

            std::cout << c << ' ' << l << ' ';
        }
        //std::cerr << "Final BWT:\n" << s << std::endl;
    } else if (fmi.r) {
        std::cerr << "ERROR: fmd is a multirope (mrope)? I don't know what that is." << std::endl;
        return 1;
    }

    return 0;
}
