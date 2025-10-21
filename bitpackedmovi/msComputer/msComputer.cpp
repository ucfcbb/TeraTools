#include "lcpComputer/lcpComputer.h"
#include "matchingStatistics/msIndex.h"
#include "moveStructure/moveStructure.h"

void testLastWeek(){
    // std::string inputOptbwtrl = argv[1];
    // LCPComputer lcp_comp(inputOptbwtrl);
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


void computeMSRow(const std::string& query){
	/* Assume, we'll have the following data-structures
	 * i - position - [1, n]
     * x is the interval that contains i
     * LF(i, x) = (LF[i], x')
     * Phi(i, x) = (Phi[i], x')
     * RLBWT[x]
	 */

	MoveStructureTable::IntervalPoint inp;
    inp.interval = 0;
    inp.offset = 0;
    
    std::vector row(query.size(), 0);
    // Encode the query into integers 
	//
    for (int i = query.size() - 1; i >= 0; ++i) {
        if (query[i] != RLBWT[inp.interval]){
            // find the next run that has the same character as query[i]
            // we do this by "repositioning"
            // In Movi, they search up and down, here we don't need to do that?
        }
        // perform LF(i, x)
		// .map()
    }
}


int main(int argc, char*argv[]) {
    if (argc != 2) {
		std::cout << "msComputer <path/to/optbwtrl/file>" << std::endl;
        exit(1);
    }

    return 0;
}
