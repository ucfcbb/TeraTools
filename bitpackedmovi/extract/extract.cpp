#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 4 && argc != 5 && argc != 6) {
        std::cerr << "Usage: extract <seq> [start position] [length] <input.len> <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 3, 4, or 5" << std::endl;
        exit(1);
    }

    uint64_t seq = atoi(argv[1]), 
             pos = (argc > 4)? atoi(argv[2]) : static_cast<uint64_t>(-1), 
             len = (argc > 5)? atoi(argv[3]) : static_cast<uint64_t>(-1);

    std::ifstream in(argv[argc-2]);
    if (!in.is_open()) {
        std::cerr << "Contig name file: \"" << argv[argc-2] << "\" failed to open!" << std::endl;
        exit(1);
    }

    uint64_t contigNum = seq/2;

    auto maxLen = std::numeric_limits<std::streamsize>::max();
    for (uint64_t i = 0; i < contigNum; ++i) 
        in.ignore(maxLen, '\n');
    std::string contigName; 
    getline(in, contigName);

    //Timer.start("Loading");
    OptBWTRL index(argv[argc-1]);
    //Timer.stop(); //Loading
    
    std::cout << ">" << contigName << "\n" << index.extract(seq, pos, len) << '\n';

    return 0;
}
