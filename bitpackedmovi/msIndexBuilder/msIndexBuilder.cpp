#include<optbwtrl/optbwtrl.h>
#include<optbwtrl/lcpComputer.h>

static constexpr char* ms_index_extension = ".ms_index";

int main(int argc, char*argv[]) {
    // Check builder.cpp to parse rb3 arguments
    // Optbwtr object

    // Build RLBWT from ropebwt3, take totalLen, rlbwt, runlens from Optbwtr
    uRange alphRange; 
    std::vector<uint64_t> alphCounts;
    RLBWTconstruction(rb3, alphRange, alphCounts);
    // rb3_fmi_free(rb3);

    // Take LF from Optbwtr
    std::vector<MoveStructure::IntervalPoint> alphStarts;
    LFconstruction(alphRange, alphCounts, alphStarts);
    // Set intLens pointer to data for intLens

    // serialize rlbwt, runlens, LF

    // lcpComputer
    LCPComputer lcp_comp(rb3);
    // call serialize

}

uint64_t totalLen = 0;
sdsl::int_vector<> rlbwt, runlens;
MoveStructure LF;

// sdsl::int_vector<> F;
// sdsl::int_vector<> Flens;
// MoveStructure Psi;

// sdsl::int_vector<> intAtTop;

// sdsl::int_vector<> PhiIntLen;
// MoveStructure Phi;
// sdsl::int_vector<> PLCPsamples;