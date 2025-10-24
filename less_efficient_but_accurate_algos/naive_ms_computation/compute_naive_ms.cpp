#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include "build_bwt.h"


template<typename T>
void printVectors(vector<T>& v,const string label="vector") {
	cout << label << " contains: " << endl;
	for(auto val: v){
		cout << val << " ";
	}
	cout << endl << endl;
}

// Prints the BWT column, i.e. the last column (L) of
// the sorted suffixes of the input text
void printBWT(SuffixArray& sa, const string& combined) {
	cout << "BWT (L): " << endl;
	for(size_t i = 1; i < sa.sa.size(); ++i) {
		cout << (sa.sa[i] != 0 ? combined[sa.sa[i]-1]: combined[combined.size()-1]);
	}
	cout << endl << endl;

}

void readPatternFile(const string& patternFile, string& pattern){
	if (!filesystem::exists(patternFile)){
		cout << patternFile << " does not exist!!!" << endl;
		exit(EXIT_FAILURE);
	}
	ifstream in(patternFile);
	getline(in, pattern);
	return;
}

int main(int argc, char* argv[]) {
	// input will be text
	// input will be pattern/query
	//
	if (argc != 3) {
		std::cout << "./compute_naive_ms <path/to/text/file> <path/to/pattern/file>" << std::endl;
		exit(1);
	}
	unordered_map<char, int> BWTInt;
	BWTInt['$'] = 0;
	BWTInt['A'] = 1;
	BWTInt['C'] = 2;
	BWTInt['G'] = 3;
	BWTInt['T'] = 4;
	BWTInt['N'] = 5;

	std::string textFile = argv[1];
	std::string patternFile = argv[2];

	std::cout << "Text file: " << textFile << std::endl;
	std::cout << "Pattern file: " << patternFile << std::endl;

    int numStrs = 0;
    vector<int> str_starts;
    string combined;
	readTextFile(textFile, str_starts, combined, numStrs);

	string pattern;
	readPatternFile(patternFile, pattern);
	cout << "Lenth of pattern: " << pattern.size() << endl;

	std::cout << "numStrs = " << numStrs << std::endl;
	std::cout << "combined = " << combined << std::endl;

    map<char, int> id;
    build_char_map(combined, numStrs, id);
	for(auto kv: id) {
		std::cout << kv.first << " : " << kv.second << std::endl;
	}
	std::cout << std::endl;

    int n = combined.size();
	cout << "n = " << n << endl;
    vi num = build_num_vec(combined, id, numStrs);
	// cout << "num.length = " << num.size() << endl;
	// printVectors<int>(num, "num");

	// Compute SA and LCP
    SuffixArray sa(num, id.size() + numStrs + 1);
	// printBWT(sa, combined);
	vector<int> bwt(sa.sa.size(), 0);
	compute_bwt(sa, combined, n, bwt);
	print_bwt(sa, combined, n);

	for(int i = 1; i <= bwt.size(); ++i){
		cout << bwt[i] << " ";
	}
	cout << endl;
	
	// Compute C array 
	// the maximum value in the values of the key,value pair
	vector<int> count(id.size() + numStrs + 2);
	compute_c(num, count);
	printVectors<int>(count, "count");

	cout << "SuffixArray.length = " << sa.sa.size() << endl;
	printVectors<int>(sa.sa, "Suffix Array");
	cout << "LCP.length = " << sa.lcp.size() << endl;
	printVectors<int>(sa.lcp, "LCP");

    // Compute LF mapping
    vector<int> rank(n);
    for (int i = 1; i <= n; i++) rank[sa.sa[i]] = i - 1;
    vector<int> LF(n);
    for (int i = 0; i < n; i++) LF[i] = rank[(sa.sa[i + 1] - 1 + n) % n];
	// cout << "Rank.length = " << rank.size() << endl;
	// printVectors<int>(rank, "rank");
	cout << "LF.length = " << LF.size() << endl;
	printVectors<int>(LF, "LF");

	// compute ms
	vector<int> ms_row(pattern.size(), 0);
	vector<int> ms_length(pattern.size(), 0);
	// Compute row
	int pos = 1; 
	int len = 0; 
	for (int i = pattern.size()-1; i >= 0; --i) {
		if (pos == 0 || pattern[i] != bwt[pos]) {
			// scan up until we find a character 
			// that matches pattern[i]
			int up = pos;
			while (up > 0 && pattern[i] != bwt[up]) { 
				len = min(sa.lcp[up], len);
				--up;
			}
			// whichever suffix has the maximal match
			// that will be the new pos
			pos = up;
		}

		if (pos == 0){
			// c[i] is the number of occurences of characters smaller than i that in the Text
			pos = count[id[pattern[i]]];
			assert(len == 0);
			len = 0;
		} else {
			pos = LF[pos-1] + 1;
			++len;
		}
		ms_row[i] = pos;
		ms_length[i] = len;
	}

	cout << "MS.row: " << endl;
	for(auto val: ms_row){
		cout << val << " ";
	}
	cout << endl;

	cout << "MS.suff: " << endl;
	for(auto val: ms_row){
		cout << sa.sa[val] << " ";
	}
	cout << endl;
	return 0;
}

