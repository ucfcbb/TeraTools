#include "TeraLCP/TeraLCP.h"
#include "matchingStatistics/msIndex.h"
#include "moveStructure/moveStructure.h"
#include <chrono>

// void computeMSRow(const std::string& inpMSIndex, const std::string& query = ""){
// 	/* Assume, we'll have the following data-structures from msIndex
// 	 * i - position - [1, n]
//      * x is the interval that contains i
//      * LF(i, x) = (LF[i], x')
//      * Phi(i, x) = (Phi[i], x')
//      * RLBWT[x]
// 	 */
    
//     // TeraLCP lcp_comp(inputOptbwtrl);
	
// 	MSIndex msIndex;
// 	std::ifstream in(inpMSIndex);
// 	std::cout << "Loading " << inpMSIndex << "..." << std::endl;
// 	msIndex.load(in);
// 	std::cout << "Successfully loaded" << inpMSIndex << "." << std::endl;
// 	in.close();

//     std::vector row(query.size(), 0);
//     // Encode the query into integers 
// 	std::unordered_map<char, int> BWTInt;
// 	BWTInt['$'] = 0;
// 	BWTInt['A'] = 1;
// 	BWTInt['C'] = 2;
// 	BWTInt['G'] = 3;
// 	BWTInt['T'] = 4;
// 	BWTInt['N'] = 5;

// 	// Initial interval
// 	MoveStructureTable::IntervalPoint inp;
//     inp.interval = 0;
//     inp.offset = 0;
//     for (int i = query.size() - 1; i >= 0; --i) {
// 		row[i] = inp.interval;
//         if (BWTInt[query[i]] != msIndex.rlbwt[inp.interval]) {
//             // find the next run that has the same character as query[i]
//             // we do this by "repositioning"
//             // In Movi, they search up and down, here we don't need to do that?
//         }
//         // perform LF(i, x) to find (LF[i], x')
// 		MoveStructureTable::IntervalPoint newinp = msIndex.LF.map(inp);
//     }
// }

int main(int argc, char*argv[]) {
    if (argc != 3) {
		std::cout << "msComputer <path/to/msIndex/file> <path/to/pattern/file>" << std::endl;
        exit(1);
    }

	std::string inpMSIndex = argv[1];
	MSIndex msIndex;
	std::ifstream in(inpMSIndex);
	std::cout << "Loading " << inpMSIndex << "..." << std::endl;
	msIndex.load(in);
	std::cout << "Successfully loaded" << inpMSIndex << "." << std::endl;
	in.close();

	std::string patternFile = argv[2];
	std::ifstream patternIn(patternFile);
	std::string pattern;
	getline(patternIn, pattern);
	patternIn.close();
	std::cout << "Successfully loaded pattern from " << patternFile << "." << std::endl;
	const char* pattern_chars = pattern.c_str();
	uint64_t m = pattern.size();
	std::cout << "Length of pattern: " << m << std::endl;

	std::ofstream outPhi(patternFile + ".ms_phi.txt");
	std::ofstream outPsi(patternFile + ".ms_psi.txt");
	std::ofstream outDual(patternFile + ".ms_dual.txt");

    std::cout << "ms_phi: " << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	auto ms_result_phi = msIndex.ms_phi(pattern_chars, m);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	std::cout << "\t Time taken: " << duration.count() << " milliseconds" << std::endl;
	std::vector<uint64_t> ms_len_phi = std::get<0>(ms_result_phi);
	std::vector<uint64_t> ms_pos_phi = std::get<1>(ms_result_phi);
	outPhi << "\tms_len: ";
	for (auto len : ms_len_phi) {
		outPhi << len << " ";
	}
	outPhi << std::endl;
	outPhi << "\tms_pos: ";
	for (auto pos : ms_pos_phi) {
		outPhi << pos << " ";
	}
	outPhi << std::endl;
	#ifdef STATS
	msIndex.print_ms_stats();
	#endif

    std::cout << "ms_psi: " << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    auto ms_result_psi = msIndex.ms_psi(pattern_chars, m);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\t Time taken: " << duration.count() << " milliseconds" << std::endl;
	std::vector<uint64_t> ms_len_psi = std::get<0>(ms_result_psi);
	std::vector<uint64_t> ms_pos_psi = std::get<1>(ms_result_psi);
	outPsi << "\tms_len: ";
	for (auto len : ms_len_psi) {
		outPsi << len << " ";
	}
	outPsi << std::endl;
	outPsi << "\tms_pos: ";
	for (auto pos : ms_pos_psi) {
		outPsi << pos << " ";
	}
	outPsi << std::endl;
	#ifdef STATS
	msIndex.print_ms_stats();
	#endif

    std::cout << "ms_dual: " << std::endl;
    start_time = std::chrono::high_resolution_clock::now();
    auto ms_result_dual = msIndex.ms_dual(pattern_chars, m);
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\t Time taken: " << duration.count() << " milliseconds" << std::endl;
	std::vector<uint64_t> ms_len_dual = std::get<0>(ms_result_dual);
	std::vector<uint64_t> ms_pos_dual = std::get<1>(ms_result_dual);
	outDual << "\tms_len: ";
	for (auto len : ms_len_dual) {
		outDual << len << " ";
	}
	outDual << std::endl;
	outDual << "\tms_pos: ";
	for (auto pos : ms_pos_dual) {
		outDual << pos << " ";
	}
	outDual << std::endl;
	#ifdef STATS
	msIndex.print_ms_stats();
	#endif

    return 0;
}
