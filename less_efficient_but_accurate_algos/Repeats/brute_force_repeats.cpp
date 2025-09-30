#include <bits/stdc++.h>
using namespace std;

using vi = vector<int>;
#define all(x) begin(x), end(x)

struct SuffixArray {
    vi sa, lcp;
    SuffixArray(vi &s, int lim = 256) {
        int n = s.size() + 1, k = 0, a, b;
        vi x(all(s) + 1), y(n), ws(max(n, lim)), rank(n);
        x.back() = 0;
        sa = lcp = y, iota(all(sa), 0);
        for (int j = 0, p = 0; p < n; j = max(1, j * 2), lim = p) {
            p = j, iota(all(y), n - j);
            for (int i = 0; i < n; i++)
                if (sa[i] >= j) y[p++] = sa[i] - j;
            fill(all(ws), 0);
            for (int i = 0; i < n; i++) ws[x[i]]++;
            for (int i = 1; i < lim; i++) ws[i] += ws[i - 1];
            for (int i = n; i--;) sa[--ws[x[y[i]]]] = y[i];
            swap(x, y), p = 1, x[sa[0]] = 0;
            for (int i = 1; i < n; i++) a = sa[i - 1], b = sa[i], x[b] = (y[a] == y[b] && y[a + j] == y[b + j]) ? p - 1 : p++;
        }
        for (int i = 1; i < n; i++) rank[sa[i]] = i;
        for (int i = 0, j; i < n - 1; lcp[rank[i++]] = k)
            for (k &&k--, j = sa[rank[i] - 1]; s[i + k] == s[j + k]; k++);
    }
};

// Helper to build character id maps
void build_char_map(const string &combined, int numStrs, map<char, int> &id) {
    auto sorted = combined;
    sort(all(sorted));
    for (int i = numStrs; i < sorted.size(); i++)
        if (!id.count(sorted[i]))
            id[sorted[i]] = id.size() + numStrs + 1;
}

// Helper to convert combined string to integer vector
vi build_num_vec(const string &combined, const map<char, int> &id, int numStrs) {
    vi num;
    int c = 1;
    for (char x : combined) {
        if (x == '$')
            num.push_back(c++);
        else
            num.push_back(id.at(x));
    }
    return num;
}

// Helper to print suffix array and LCP
void print_sa_lcp(const SuffixArray &sa, int n) {
    cout << "Suffix Array (SA):\n  ";
    for (int i = 1; i <= n; i++) cout << sa.sa[i] + 1 << (i < n ? " " : "\n");
    cout << "Longest Common Prefix (LCP):\n  ";
    for (int i = 1; i <= n; i++) cout << sa.lcp[i] << (i < n ? " " : "\n");
    cout << endl;
}

// Helper to print BWT
void print_bwt(const SuffixArray &sa, const string &combined, int n) {
    cout << "Burrows-Wheeler Transform (BWT):\n  ";
    for (int i = 1; i <= n; i++) {
        int idx = sa.sa[i] - 1;
        char bwt_char = (idx < 0) ? combined.back() : combined[idx];
        cout << bwt_char;
    }
    cout << "\n\n";
}

// Helper to compute and print PHI/IPHI
void print_phi_iphi(const SuffixArray &sa, int n) {
    cout << "PHI array:\n  ";
    vector<int> phi(n);
    for (int k = 1; k <= n; ++k) phi[sa.sa[k]] = sa.sa[((k - 2 + n) % n) + 1];
    for (int i = 0; i < n; i++) cout << phi[i] + 1 << (i < n - 1 ? " " : "\n");
    cout << "IPHI array:\n  ";
    vector<int> iphi(n);
    for (int k = 1; k <= n; ++k) iphi[sa.sa[k]] = sa.sa[(k % n) + 1];
    for (int i = 0; i < n; i++) cout << iphi[i] + 1 << (i < n - 1 ? " " : "\n");
    cout << endl;
}

// Helper to compute and print LF
void print_lf(const SuffixArray &sa, int n) {
    vector<int> rank(n);
    for (int i = 1; i <= n; i++) rank[sa.sa[i]] = i - 1;
    cout << "LF mapping:\n  ";
    vector<int> LF(n);
    for (int i = 0; i < n; i++) LF[i] = rank[(sa.sa[i + 1] - 1 + n) % n];
    for (int i = 0; i < n; i++) cout << LF[i] + 1 << (i < n - 1 ? " " : "\n");
    cout << endl;
}

int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        std::cerr << "Usage: ./a.out {SM,LM} [Length Threshold] < input.in\n"
                  << argc - 1 << " arguments passed instead of 1 or 2" << std::endl;
        exit(1);
    }
    if (strcmp(argv[1], "SM") && strcmp(argv[1], "LM")) {
        std::cerr << "First argument must be either \"SM\" or \"LM\"."
                  << " If passed SM, program will output supermaximal repeats."
                  << " If passed LM, it will output regular (a.k.a. maximal a.k.a. locally maximal) repeats." << std::endl;
        exit(1);
    }
    int numStrs = 0;
    vector<int> str_starts;
    string combined;
    {
        int idx = 0;
        string s;
        while (cin >> s) {
            ++numStrs;
            str_starts.push_back(idx);
            combined += s + (char)('$');
            idx += s.size() + 1;
        }
    }

    auto getSOpair = [&str_starts](int sa_val) -> pair<int, int> {
        auto it = upper_bound(str_starts.begin(), str_starts.end(), sa_val) - 1;
        int sa_a = it - str_starts.begin();
        int sa_o = sa_val - *it;
        return {sa_a, sa_o};
    };

    map<char, int> id;
    build_char_map(combined, numStrs, id);
    int n = combined.size();
    vi num = build_num_vec(combined, id, numStrs);
    SuffixArray sa(num, id.size() + numStrs + 1);

    auto getMaximalRepeats_BRUTE = [&](const string &text, int threshold) -> void {
        cout << "seq1\tpos1\tseq2\tpos2\tlen\n";
        int n = text.size();
        for (int len = threshold; len < n; ++len) {
            unordered_map<string, vector<int>> occs;
            for (int i = 0; i + len <= n; ++i) {
                string sub = text.substr(i, len);
                if (sub.find('$') != string::npos) continue;
                occs[sub].push_back(i);
            }
            for (auto &kv : occs) {
                const string &sub = kv.first;
                const vector<int> &positions = kv.second;
                if (positions.size() == 1) continue;

                vi valids;
                for (int posIdx = 0; posIdx < positions.size(); ++posIdx) {
                    for (int otherIdx = posIdx + 1; otherIdx < positions.size(); ++otherIdx) {
                        int pos = positions[posIdx];
                        int other_pos = positions[otherIdx];

                        bool left = true, right = true;
                        if (pos > 0) {
                            if (text[pos - 1] != '$' && text[other_pos - 1] == text[pos - 1]) {
                                left = false;
                            }
                        }

                        char r = text[pos + len];
                        if (other_pos + len < n && r != '$' && text[other_pos + len] == r) {
                            right = false;
                        }

                        if (left && right) {
                            auto [sa_a1, sa_o1] = getSOpair(pos);
                            auto [sa_a2, sa_o2] = getSOpair(other_pos);
                            cout << sa_a1 << "\t" << sa_o1 << "\t" << sa_a2 << "\t" << sa_o2 << "\t" << len << "\n";
                        }
                    }
                }
            }
        }
    };

    auto getSupermaximalRepeats_BRUTE = [&](const string &text, int threshold) -> void {
        cout << "seq\tpos\tlen\tocc\n";

        int n = text.size();
        vector<int> ISA(n);
        for (int i = 0; i < n; i++) ISA[sa.sa[i + 1]] = i;

        unordered_map<string, vector<int>> next_occs;

        for (int len = n; len >= threshold; --len) {
            unordered_map<string, vector<int>> occs;
            for (int i = 0; i + len <= n; ++i) {
                string sub = text.substr(i, len);
                if (sub.find('$') != string::npos) continue;
                occs[sub].push_back(i);
            }

            for (auto &kv : occs) {
                const string &sub = kv.first;
                const vector<int> &positions = kv.second;
                if (positions.size() == 1) continue;
                for (int pos : positions) {
                    bool is_supermaximal = true;
                    if (pos > 0 && text[pos - 1] != '$') {
                        string extended_left = text[pos - 1] + sub;
                        int count = next_occs.count(extended_left) ? next_occs[extended_left].size() : 0;
                        if (count != 1) {
                            is_supermaximal = false;
                        }
                    }
                    if (is_supermaximal && pos + len < n && text[pos + len] != '$') {
                        string extended_right = sub + text[pos + len];
                        int count = next_occs.count(extended_right) ? next_occs[extended_right].size() : 0;
                        if (count != 1) {
                            is_supermaximal = false;
                        }
                    }
                    if (is_supermaximal) {
                        cout << getSOpair(pos).first << "\t" << getSOpair(pos).second << "\t" << len << "\t";
                        cout << positions.size();
                        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
                        for (auto p : positions) {
                            if (p == pos) continue;
                            pq.push({-ISA[p], p});
                        }
                        while (!pq.empty()) {
                            auto [_, p] = pq.top();
                            pq.pop();
                            cout << "\t" << getSOpair(p).first << "\t" << getSOpair(p).second;
                        }
                        cout << "\n";
                    }
                }
            }
            next_occs = move(occs);
        }
    };

    uint64_t len = 1;
    if (argc == 4)
        len = atoi(argv[2]);

    if (strcmp(argv[1], "SM") == 0) {
        getSupermaximalRepeats_BRUTE(combined, len);
    }
    if (strcmp(argv[1], "LM") == 0) {
        getMaximalRepeats_BRUTE(combined, len);
    }
}