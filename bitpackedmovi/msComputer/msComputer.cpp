#include"optbwtrl/lcpComputer.h"

int main(int argc, char*argv[]) {
    if (argc != 2) {
        exit(1);
    }

    LCPComputer lcp_comp;
    std::string inputOptbwtrl = argv[1];
    std::ifstream in(inputOptbwtrl);

    lcp_comp.load(in);
    in.close();
    
    // lcp_comp.printPhiAndLCP();
    // lcp_comp.printRaw();
    
    size_t steps = 0;

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

    // int_vector<16> bwt_to_psi_int_at_top;
    // int_vector<16> bwt_to_psi_offset_at_top;
    // int_vector<16> bwt_to_phi_int_at_bottom;
    // int_vector<16> bwt_to_phi_offset_at_bottom;

    return 0;
}