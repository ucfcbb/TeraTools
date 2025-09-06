#include<optbwtrl/optbwtrl.h>


int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: stats <input.optbwtrl>\n"
            << argc-1 << " arguments passed instead of 1" << std::endl;
        exit(1);
    }
    Timer.start("Loading");
    OptBWTRL index(argv[1]);
    Timer.stop(); //Loading

    Timer.start("Summing");
    std::cout << "The sum of LCPs at the top of runs in the BWT is: " << index.sumLCPTopRun() << std::endl;
    Timer.stop(); //Summing

    Timer.start("Interval Traversals");
    OptBWTRL::Histograms a = index.IntervalTraversals();
    Timer.stop();

    uint64_t maxSize = std::max(a.LF.size(), a.INVMOVE.Phi.size());
    maxSize = std::max(maxSize, a.INVMOVE.InvPhi.size());
    
    std::cout << "i\tLF\tINVMOVEPHI\tINVMOVEINVPHI\n";
    for (uint64_t i = 0; i < maxSize; ++i) {
        std::cout << i << '\t' 
            << ((i < a.LF.size())? a.LF[i] : 0) << '\t'
            << ((i < a.INVMOVE.Phi.size())? a.INVMOVE.Phi[i] : 0) << '\t'
            << ((i < a.INVMOVE.InvPhi.size())? a.INVMOVE.InvPhi[i] : 0) << '\n';
    }
    return 0;
}
