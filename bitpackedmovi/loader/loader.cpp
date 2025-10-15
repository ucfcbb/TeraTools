#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: loader [BWT,MINLCP] <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 2" << std::endl;
        exit(1);
    }
    if (strcmp(argv[1], "BWT") && strcmp(argv[1], "PHI") && strcmp(argv[1], "MINLCP")) {
        std::cerr << "Second parameter must be either 'BWT' or 'MINLCP', '" << argv[1] << "' passed." << std::endl;
        exit(1);
    }
    OptBWTRL index(argv[2]);

    if (strcmp(argv[1], "BWT") == 0)
        index.printRaw();
    else if (strcmp(argv[1], "PHI") == 0)
        index.printPhiAndLCP();
    else if (strcmp(argv[1], "MINLCP") == 0)
        index.ComputeMinLCPRun(std::cout);

    /*
    if (index.validateAllExceptRLBWT())
        std::cout << "Passed" << std::endl;
    else {
        std::cout << "Failed!" << std::endl;
        std::cerr << "ERROR: Failed!" << std::endl;
        exit(1);
    }
    */
    return 0;
}
