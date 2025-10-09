#include<iostream>
#include<stack>
#include<chrono>

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

