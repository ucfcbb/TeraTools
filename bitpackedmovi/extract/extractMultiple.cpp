#include<optbwtrl/optbwtrl.h>



int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: extractMultiple <inputFile> <input.len> <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 3" << std::endl;
        exit(1);
    }


    std::ifstream in(argv[1]), inLen(argv[2]), inDex(argv[3]);

    if (!in.is_open()) {
        std::cerr << "ERROR: Input file for strings to extract: \"" << argv[1] << "\" failed to open!" << std::endl;
        exit(1);
    }
    if (!inLen.is_open()) {
        std::cerr << "ERROR: Input file for contig names: \"" << argv[2] << "\" failed to open!" << std::endl;
        exit(1);
    }

    timer errTimer('|', "Timer:", std::cerr, true);

    //get all queries
    errTimer.start("Getting queries");
    struct Query {
        uint64_t seq, pos, len;
    };
    std::vector<Query> queries;
    Query curr;
    while (in >> curr.seq >> curr.pos >> curr.len) 
        queries.push_back(curr);
    errTimer.stop();

    //get all contig names
    errTimer.start("Getting contig names");
    std::map<uint64_t,std::string> contigNames;
    for (auto a : queries)
        contigNames.emplace(a.seq/2, "");

    uint64_t line = 0;
    auto maxLen = std::numeric_limits<std::streamsize>::max();
    for (auto it = contigNames.begin(); it != contigNames.end(); ++it) {
        for (; line < it->first; ++line)
            in.ignore(maxLen, '\n');
        getline(in, it->second);
        ++line;
    }
    errTimer.stop();

    errTimer.start("Loading");
    OptBWTRL index(argv[argc-1]);
    errTimer.stop(); //Loading

    //do all queries
    //there is a much more efficient way to do this but I haven't implemented it
    errTimer.start("Computing Queries");
    for (auto query: queries) {
        std::cout << ">" << contigNames[query.seq/2] << "seq[start,end]: " << query.seq << '[' << query.pos << ',' << query.pos + query.len - 1 << ']' << "\n"
            << index.extract(query.seq, query.pos, query.len) << "\n\n";
    }
    errTimer.stop();
}
