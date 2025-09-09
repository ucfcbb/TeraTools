#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        std::cerr << "Usage: repeats <input.optbwtrl> [Length Threshold]\n"
            << argc-1 << " arguments passed instead of 1 or 2" << std::endl;
        exit(1);
    }
    
    uint64_t len = 1;
    if (argc == 3)
        len = atoi(argv[2]);
    std::cout << "Outputting supermaximal repeats of length " << len << std::endl;


    Timer.start("Loading");
    OptBWTRL index(argv[1]);
    Timer.stop(); //Loading

    Timer.start("Computing Supermaximal Repeats");
    index.superMaximalRepeats(std::cout, len);
    Timer.stop();
    return 0;
}
