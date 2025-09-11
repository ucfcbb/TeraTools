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
    void stop() {
        auto end = std::chrono::high_resolution_clock::now();
        auto beginNamePair = processes.top();
        processes.pop();
        os << prefix << std::string(1+processes.size(), scopeChar) << "Ending '" << beginNamePair.second 
            << "'. It took " << std::chrono::duration<double>(end - beginNamePair.first).count() << " seconds\n";
        if (flush)
            os.flush();
    }
    void stopAllProcesses() {
        while (processes.size())
            stop();
    }

    timer() = delete;
    timer(const char c, const std::string& pref, std::ostream& s, bool f): scopeChar(c), prefix(pref), os(s), flush(f) {}
    ~timer() {
        stopAllProcesses();
    }
}Timer('|', "Timer:", std::cout, true);

struct uRange {
    uint64_t min, max;
};

std::ostream& operator<<(std::ostream& os, uRange range) {
    os << '[' << range.min << ',' << range.max << ']';
    return os;
}

