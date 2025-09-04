#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: loader <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 1" << std::endl;
        exit(1);
    }
    OptBWTRL index(argv[1]);

    if (index.validateAllExceptRLBWT())
        std::cout << "Passed" << std::endl;
    else {
        std::cout << "Failed!" << std::endl;
        std::cerr << "ERROR: Failed!" << std::endl;
        exit(1);
    }

    return 0;
}
