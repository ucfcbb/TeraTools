#include "lcpComputer/lcpComputer.h"
#include "matchingStatistics/msIndex.h"
#include "moveStructure/moveStructure.h"
#include <unordered_map>

void testLastWeek(const std::string& inputOptbwtrl){ 
	LCPComputer lcp_comp(inputOptbwtrl);
    // std::ifstream in(inputOptbwtrl);

    // lcp_comp.load(in);
    // in.close();
    
    // lcp_comp.printPhiAndLCP();
    // lcp_comp.printRaw();
    
    // size_t steps = 0;

    // MoveStructure::IntervalPoint phiPoint;
    // phiPoint.interval = 0;
    // phiPoint.offset = 0;

    // do {
    //     phiPoint = lcp_comp.Phi.map(phiPoint);
    //     std::cout << phiPoint.interval << '\t' << phiPoint.offset << "\t-->\t";
    //     std::cout << lcp_comp.PLCPsamples[phiPoint.interval] - phiPoint.offset << std::endl;
    //     steps++;
    // } while (phiPoint.interval != 0 || phiPoint.offset != 0);
    // std::cout << steps << std::endl;

	/*
		MoveStructure::IntervalPoint psiPoint;
		psiPoint.interval = 0;
		psiPoint.offset = 0;
		steps = 0;
		do {
			psiPoint = lcp_comp.Psi.map(psiPoint);
			std::cout << psiPoint.interval << '\t' << psiPoint.offset << "\t-->\t";
			std::cout << lcp_comp.F[psiPoint.interval] << std::endl;
			steps++;
		} while (psiPoint.interval != 0 || psiPoint.offset != 0);
		std::cout << steps << std::endl;
	*/

    // int_vector<16> bwt_to_psi_int_at_top;
    // int_vector<16> bwt_to_psi_offset_at_top;
    // int_vector<16> bwt_to_phi_int_at_bottom;
    // int_vector<16> bwt_to_phi_offset_at_bottom;
}


void computeMSRow(const std::string& inpMSIndex, const std::string& query = ""){
	/* Assume, we'll have the following data-structures from msIndex
	 * i - position - [1, n]
     * x is the interval that contains i
     * LF(i, x) = (LF[i], x')
     * Phi(i, x) = (Phi[i], x')
     * RLBWT[x]
	 */
    
    // LCPComputer lcp_comp(inputOptbwtrl);
	
	MSIndex ms;
	std::ifstream in(inpMSIndex);
	std::cout << "Loading msIndex..." << std::endl;
	ms.load(in);
	std::cout << "Loaded msIndex successfully." << std::endl;
	in.close();
	// ms.constructFromLCPIndexFileWriteAndClear(in, out);
	// out.close();

    std::vector row(query.size(), 0);
    // Encode the query into integers 
	std::unordered_map<char, int> BWTInt;
	BWTInt['$'] = 0;
	BWTInt['A'] = 1;
	BWTInt['C'] = 2;
	BWTInt['G'] = 3;
	BWTInt['T'] = 4;
	BWTInt['N'] = 5;

	// Initial interval
	/*
	MoveStructureTable::IntervalPoint inp;
    inp.interval = 0;
    inp.offset = 0;
    for (int i = query.size() - 1; i >= 0; --i) {
		row[i] = inp.interval;
        // if (BWTInt[query[i]] != msIndex::RLBWT[inp.interval]) {
            // find the next run that has the same character as query[i]
            // we do this by "repositioning"
            // In Movi, they search up and down, here we don't need to do that?
        // }
        // perform LF(i, x) to find (LF[i], x')
		// inp = msIndex::LF(inp.offset, inp.interval) 
    }*/
}


int main(int argc, char*argv[]) {
    if (argc != 2) {
		std::cout << "msComputer <path/to/msIndex/file>" << std::endl;
        exit(1);
    }

	std::string inpMSIndex = argv[1];
	computeMSRow(inpMSIndex);
    return 0;
}
