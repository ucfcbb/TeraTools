#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: repeats {SM,LM} [Length Threshold] <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 1 or 2" << std::endl;
        exit(1);
    }
    if (strcmp(argv[1], "SM") && strcmp(argv[1], "LM")) {
        std::cerr << "First argument must be either \"SM\" or \"LM\"."
           << " If passed SM, program will output supermaximal repeats."
           << " If passed LM, it will output regular (a.k.a. maximal a.k.a. locally maximal) repeats." << std::endl;
        exit(1);
    }
    uint64_t len = 1;
    if (argc == 4)
        len = atoi(argv[2]);
    std::cout << "Outputting" << ((strcmp(argv[1], "LM"))? " supermaximal repeats" : " repeats")
        << " of length " << len << std::endl;

    


    Timer.start("Loading");
    OptBWTRL index(argv[argc-1]);
    Timer.stop(); //Loading

    if (strcmp(argv[1], "SM") == 0) {
        Timer.start("Computing Supermaximal Repeats");
        index.superMaximalRepeats(std::cout, len);
        Timer.stop();
    }

    if (strcmp(argv[1], "LM") == 0) {
        Timer.start("Computing repeats");
        index.repeats(std::cout, len);
        Timer.stop();
    }
    return 0;
}
