#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: stats <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 1" << std::endl;
        exit(1);
    }
    Timer.start("Loading");
    OptBWTRL index(argv[1]);
    Timer.stop(); //Loading

    Timer.start("Computing Supermaximal Repeats");
    index.superMaximalRepeats(std::cout);
    Timer.stop();
    return 0;
}
