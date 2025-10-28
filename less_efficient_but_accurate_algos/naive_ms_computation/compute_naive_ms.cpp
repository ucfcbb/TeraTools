#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cassert>

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

// Here, the definition of MS.row[i] is the suffix of the text that preceeds the suffix that has the longest prefix match with i-th suffix of the pattern
pair<vector<int>, vector<int>> computePredOnlyMS(const string& pattern, vector<int>& bwt, vector<int>& LF, SuffixArray& sa, vector<int>& count, map<char, int>& id) {
	// Compute row and length
	// Here, the definition of row is the predecessor suffix
	vector<int> ms_row(pattern.size(), 0);
	vector<int> ms_length(pattern.size(), 0);
	int pos = 1;  // starting at the top
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
	return {ms_row, ms_length};
}

// Here, the definition of MS.row[i] is the suffix of the text whose prefix has the longest match with suffix i-th suffix of the pattern
pair<vector<int>, vector<int>> computeTrueMS(const string& pattern, vector<int>& bwt, vector<int>& LF, SuffixArray& sa, vector<int>& count, map<char, int>& id) {
    // Compute Matching Statistics
	vector<int> ms_row(pattern.size(), 0);
	vector<int> ms_length(pattern.size(), 0);
	int n = sa.sa.size() - 1;
	int pos = n;  // start at bottom so that it's easy to check with Moni's output as per Nate
	int len = 0; 

	for (int i = pattern.size()-1; i >= 0; --i) {
		if (pos == 0 || pos == n + 1 || pattern[i] != bwt[pos]) {
			// scan up until we find a character 
			// that matches pattern[i]
			int up = pos;
			int upLen = len;
			while (up > 0 && pattern[i] != bwt[up]) { 
				upLen = min(sa.lcp[up], upLen);
				--up;
			}

			// scan down until we find a character that matches pattern[i]
			int down = pos+1; // check with the next suffix's lcp
			int downLen = len;
			while (down < n + 1 && pattern[i] != bwt[down]) { 
				downLen = min(sa.lcp[down], downLen);
				++down;
			}

			// whichever suffix has the maximal match
			// that will be the new pos to be LF'd
			cout << "### Longest common extension ###" << endl;
			cout << "At i = " << i << endl;
			cout << "up = " << up << ", upLen = " << upLen << endl;
			cout << "down = " << down << ", downLen = " << downLen << endl;

			if (up == 0) {
				assert(down < n+ 1);
				pos = down;
				len = downLen;
			} else if (down == n + 1) {
				assert(up > 0);
				pos = up;
				len = upLen;
			} else {
				assert(up > 0);
				assert(down < n+1);
				pos = (upLen >= downLen) ? up : down;
				len = max(upLen, downLen);
			}
		}

		if (pos == 0 || pos == n + 1) {
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
	return {ms_row, ms_length};
}

void printUsage() {
	std::cout << "Computes the matching statistics naively" << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << "./compute_naive_ms <path/to/text/file> <path/to/pattern/file> [optional]\n" << std::endl;
	std::cout << "If <pattern file> is not provided, just computes the BWT (and auxiliary data structures) of the text." << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	// input will be text
	// input will be pattern/query
	string textFile;
	string patternFile;
	bool computeMS = false;
	switch (argc) {
		case 2:
			textFile = argv[1];
			break;
		case 3:
			textFile = argv[1];
			patternFile = argv[2];
			computeMS = true;
			break;
		default:
			printUsage();
			exit(1);
			break;
	}

	unordered_map<char, int> BWTInt;
	BWTInt['$'] = 0;
	BWTInt['A'] = 1;
	BWTInt['C'] = 2;
	BWTInt['G'] = 3;
	BWTInt['T'] = 4;
	BWTInt['N'] = 5;

	// Read both text input files
	cout << "Text file: " << textFile << endl;
	int numStrs = 0;
	vector<int> str_starts;
	string combined;

	readTextFile(textFile, str_starts, combined, numStrs);
	cout << "combined/concatenated text = " << combined << endl;
	cout << "numStrs = " << numStrs << endl;

    map<char, int> id;
    build_char_map(combined, numStrs, id);
    int n = combined.size();
	cout << "n (length of concatenated text) = " << n << endl;
    vi num = build_num_vec(combined, id, numStrs);

	// Compute SA and LCP
    SuffixArray sa(num, id.size() + numStrs + 1);

	vector<int> bwt(sa.sa.size(), 0);
	compute_bwt(sa, combined, n, bwt);
	print_bwt(sa, combined, n);

	// Compute C (count) array 
	// the maximum value in the values of the key,value pair
	vector<int> count(id.size() + numStrs + 2);
	compute_c(num, count);
	// printVectors<int>(count, "count");

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

	// compute Inverse LF
	vector<int> InvLF(n);
    for (int i = 0; i < n; i++) InvLF[LF[i]] = i;
	printVectors<int>(InvLF, "Inverse LF");

	if (computeMS) {
		// Read pattern input files
		std::cout << "Pattern file: " << patternFile << std::endl;
		string pattern;
		readPatternFile(patternFile, pattern);
		cout << "pattern = " << pattern << std::endl;
		cout << "Lenth of pattern: " << pattern.size() << endl;

		// auto matchingStatistics = computePredOnlyMS(pattern, bwt, LF, sa, count, id);
		auto matchingStatistics = computeTrueMS(pattern, bwt, LF, sa, count, id);

		// Compute Matching Statistics
		cout << "MS.row: " << endl;
		for(auto val: matchingStatistics.first) {
			cout << val << " ";
		}
		cout << endl;

		cout << "MS.suff: " << endl;
		for(auto val: matchingStatistics.first) {
			cout << sa.sa[val] << " ";
		}
		cout << endl;

		cout << "MS.len: " << endl;
		for(auto val: matchingStatistics.second) {
			cout << val << " ";
		}
		cout << endl;
	}

	return 0;
}

