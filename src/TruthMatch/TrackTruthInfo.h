#ifndef TrackTruthInfo_h
#define TrackTruthInfo_h

#include <map>

#include <GaudiKernel/SmartRef.h>

namespace Minerva {
    class Track;
    class TG4Trajectory;
}

class TrackTruthInfo {
  public:
    TrackTruthInfo();
    ~TrackTruthInfo() {}

    void ParseTruth(const SmartRef<Minerva::Track>& track,
                    const std::map<int,Minerva::TG4Trajectory*>& trajectories);

    double GetTotalNormEnergy() const;
    double GetDataNormEnergy() const;
    double GetMCNormEnergy() const;
    double GetMCXtalkNormEnergy() const;
    double GetTotalTruthEnergy() const;
    double GetSharedTruthEnergy() const;

    double GetDataFraction() const;
    double GetXtalkFraction() const;
    double GetSharedFraction() const;

        /// return the number of primary particles contributing to
        /// this track. Note that if there are two particles with the
        /// same PDG code, then they are counted as one
    unsigned int GetTruthParticleCount() const;
        /// index = 0 return the particle with most energy deposit
        /// index = 1 second to most, and so on...
    int GetTruthPdg(unsigned int index = 0) const;

        /// Right after the GetTruthPdg() call, this returns the
        /// energy fraction by the particle
    double GetEdepFraction() const;
    
  private:
    double fTotalNormEnergy;
    double fDataNormEnergy;
    double fMCNormEnergy;
    double fMCXtalkNormEnergy;
    double fTotalTruthEnergy;
    double fSharedTruthEnergy;
    mutable double fEdepFraction;

    std::map<double,int> fEdepPdgMap;
    std::map<int,double> fPdgEdepMap;
    
};

#endif
