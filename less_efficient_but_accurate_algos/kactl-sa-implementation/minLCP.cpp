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

int main() {
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

    unordered_map<char, int> BWTInt;
    BWTInt['$'] = 0;
    BWTInt['A'] = 1;
    BWTInt['C'] = 2;
    BWTInt['G'] = 3;
    BWTInt['T'] = 4;
    BWTInt['N'] = 5;

    map<char, int> id;
    build_char_map(combined, numStrs, id);

    int n = combined.size();
    vi num = build_num_vec(combined, id, numStrs);

    SuffixArray sa(num, id.size() + numStrs + 1);

    auto BWT = [&] (int i) {
        int s = sa.sa[i + 1];
        int ind = (s)? s - 1 : n - 1;
        return combined[ind];
    };

    for (int i = 0; i < n; ++i) {
        if (i != 0 && BWT(i) == BWT(i-1)) 
            continue;
        int minLCP = sa.lcp[i + 1];
        int minLCPLoc = 0;
        for (int j = 1; i+j < n && BWT(i+j) == BWT(i); ++j) {
            if (sa.lcp[i + j + 1] <= minLCP) {
                minLCP = sa.lcp[i + j + 1];
                minLCPLoc = j;
            }
        }
        std::cout << "( " << minLCPLoc << ", " << minLCP << ")\n";
    }

    /*
    // Compute LF mapping
    vector<int> rank(n);
    for (int i = 1; i <= n; i++) rank[sa.sa[i]] = i - 1;
    vector<int> LF(n);
    for (int i = 0; i < n; i++) LF[i] = rank[(sa.sa[i + 1] - 1 + n) % n];

    // Print header
    cout << "i\tSA_S\tSA_O\tLCP\tLF\tBWT\n";
    for (int i = 0; i < n; i++) {
        int sa_val = sa.sa[i + 1];
        int lcp_val = sa.lcp[i + 1];
        int lf_val = LF[i];
        int idx = sa_val - 1;
        char bwt_char = (idx < 0) ? combined.back() : combined[idx];
        // Find sa_a and sa_o
        int sa_a = -1, sa_o = -1;
        if (sa_val < n) {
            auto it = upper_bound(str_starts.begin(), str_starts.end(), sa_val) - 1;
            sa_a = it - str_starts.begin();
            sa_o = sa_val - *it;
        }
        cout << i << "\t" << sa_a << "\t" << sa_o << "\t" << lcp_val << "\t" << lf_val << "\t" << BWTInt[bwt_char] << "\n";
    }

    cout << "\n\ni\t";
    for (int i = 0; i < n; i++) cout << i << (i < n - 1 ? "\t" : "\n");

    cout << "TEXT\t";
    for (int i = 0; i < n; i++) {
        cout << BWTInt[combined[i]] << (i < n - 1 ? "\t" : "\n");
    }

    // ISA: inverse suffix array
    cout << "ISA\t";
    vector<int> ISA(n);
    for (int i = 0; i < n; i++) ISA[sa.sa[i + 1]] = i;
    for (int i = 0; i < n; i++) cout << ISA[i] << (i < n - 1 ? "\t" : "\n");

    // PLCP: permuted LCP
    cout << "PLCP\t";
    vector<int> PLCP(n);
    for (int i = 0; i < n; i++)
        PLCP[i] = sa.lcp[ISA[i] + 1];
    for (int i = 0; i < n; i++) cout << PLCP[i] << (i < n - 1 ? "\t" : "\n");

    auto getSOpair = [&str_starts](int sa_val) -> pair<int, int> {
        auto it = upper_bound(str_starts.begin(), str_starts.end(), sa_val) - 1;
        int sa_a = it - str_starts.begin();
        int sa_o = sa_val - *it;
        return {sa_a, sa_o};
    };

    vector<int> PHI_S(n), PHI_O(n);
    int s = sa.sa[n];
    auto a = getSOpair(s);

    PHI_S[sa.sa[1]] = a.first;
    PHI_O[sa.sa[1]] = a.second;
    for (int i = 2; i <= n; i++) {
        auto a = getSOpair(sa.sa[i - 1]);
        PHI_S[sa.sa[i]] = a.first;
        PHI_O[sa.sa[i]] = a.second;
    }
    cout << "PHI_S\t";
    for (int i = 0; i < n; i++) cout << PHI_S[i] << (i < n - 1 ? "\t" : "\n");
    cout << "PHI_O\t";
    for (int i = 0; i < n; i++) cout << PHI_O[i] << (i < n - 1 ? "\t" : "\n");

    vector<int> IPHI_S(n), IPHI_O(n);
    s = sa.sa[1];
    a = getSOpair(s);
    IPHI_S[sa.sa[n]] = a.first;
    IPHI_O[sa.sa[n]] = a.second;
    for (int i = 1; i < n; i++) {
        auto a = getSOpair(sa.sa[i + 1]);
        IPHI_S[sa.sa[i]] = a.first;
        IPHI_O[sa.sa[i]] = a.second;
    }
    cout << "IPHI_S\t";
    for (int i = 0; i < n; i++) cout << IPHI_S[i] << (i < n - 1 ? "\t" : "\n");
    cout << "IPHI_O\t";
    for (int i = 0; i < n; i++) cout << IPHI_O[i] << (i < n - 1 ? "\t" : "\n");
    */
}
