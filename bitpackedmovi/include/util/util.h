#ifndef R_SA_LCP_UTIL_H
#define R_SA_LCP_UTIL_H
#include<iostream>
#include<stack>
#include<chrono>
#include<string>
#include<vector>
#include<fstream>
#include<algorithm>

#ifndef BENCHFASTONLY
enum verbosity { QUIET, TIME, VERB };
#endif

std::string getArgument(const int argc, const char* argv[], std::vector<bool>& used, const std::string& arg, bool required, bool argument) {
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
}

void testOutFile(const std::string& name) {
    if (name == "") return;
    std::ofstream out(name, std::ios::app);
    if (!out.is_open()) {
        std::cout << "File '" << name << "' failed to open for writing.\n";
        exit(1);
    }
};

void testInFile(const std::string& name) {
    if (name == "") return;
    std::ofstream in(name, std::ios::app);
    if (!in.is_open()) {
        std::cout << "File '" << name << "' failed to open for writing.\n";
        exit(1);
    }
};


class timer {
    std::stack<std::pair<std::chrono::time_point<std::chrono::high_resolution_clock>, std::string>> processes; 
    const char scopeChar;
    const std::string prefix;
    std::ostream& os;
    const bool flush;
    public:

    void start(std::string processName) {
        os << prefix << std::string(1+processes.size(), scopeChar) << "Starting '" << processName << "'\n";
        if (flush)
            os.flush();
        processes.emplace(std::chrono::high_resolution_clock::now(), processName);
    }
    double stop() {
        auto end = std::chrono::high_resolution_clock::now();
        auto beginNamePair = processes.top();
        double sec = std::chrono::duration<double>(end - beginNamePair.first).count();
        processes.pop();
        os << prefix << std::string(1+processes.size(), scopeChar) << "Ending '" << beginNamePair.second 
            << "'. It took " << sec << " seconds\n";
        if (flush)
            os.flush();
        return sec;
    }
    void stopAllProcesses() {
        while (processes.size())
            stop();
    }

    timer() = delete;
    timer(const char c, const std::string& pref, std::ostream& s, bool f): scopeChar(c), prefix(pref), os(s), flush(f) {}
    ~timer() {
        uint64_t notDone = processes.size();
        if (notDone)
            os << prefix << std::string(1+processes.size(), scopeChar) << "I AM BEING DESTRUCTED BEFORE THE FOLLOWING " << notDone << " PROCESSES ARE COMPLETED!" << '\n'; 
        if (flush)
            os.flush();
        stopAllProcesses();
        if (notDone)
            os << prefix << std::string(1+processes.size(), scopeChar) << "I WAS DESTRUCTED BEFORE THE LAST " << notDone << " PROCESSES WERE COMPLETED!" << '\n'; 
        if (flush)
            os.flush();
    }
}Timer('|', "Timer:", std::cout, true);

struct uRange {
    uint64_t min, max;
};

std::ostream& operator<<(std::ostream& os, uRange range) {
    os << '[' << range.min << ',' << range.max << ']';
    return os;
}
#endif //#ifndef R_SA_LCP_UTIL_H
