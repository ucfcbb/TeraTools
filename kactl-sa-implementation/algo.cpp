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

// Helper to read input strings and combine with separator
string read_and_combine(int &numStrs) {
    cin >> numStrs;
    string combined;
    for (int i = 0; i < numStrs; i++) {
        string s;
        cin >> s;
        combined += s + (char)('$');
    }
    return combined;
}

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
    int numStrs;
    string combined = read_and_combine(numStrs);

    map<char, int> id;
    build_char_map(combined, numStrs, id);

    int n = combined.size();
    vi num = build_num_vec(combined, id, numStrs);

    // cout << "--- Suffix Array Construction ---" << endl;
    SuffixArray sa(num, id.size() + numStrs + 1);
    // cout << "Suffix array built successfully.\n"
    // << endl;

    print_bwt(sa, combined, n);
    print_sa_lcp(sa, n);
    print_phi_iphi(sa, n);
    print_lf(sa, n);
}