#ifndef CCProtonPi0_PDG_h 
#define CCProtonPi0_PDG_h

namespace PDG
{
    // leptons
    const int electron = 11;
    const int nu_e     = 12;
    const int muon     = 13;
    const int nu_mu    = 14;
    const int tau      = 15;
    const int nu_tau   = 16;

    // gamma
    const int gamma    = 22;

    // mesons
    const int pi       = 211;
    const int pi0      = 111;
    const int K        = 321;
    const int K0       = 311;

    // baryons
    const int proton   = 2212;
    const int neutron  = 2112;

    // other
    const int bindino  = 2000000101; // nuclear binding energy
    const int nucleus  = 1000000000; // excited nucleus (10LZZZAAAI)
}

#endif
