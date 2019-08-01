#ifndef DigitVectorTruthInfo_h
#define DigitVectorTruthInfo_h

#include <map>

#include <GaudiKernel/SmartRef.h>

namespace Minerva {
    class TG4Trajectory;
    class IDDigit;
    class IDCluster;
}

class DigitVectorTruthInfo {
  public:
    DigitVectorTruthInfo();
    ~DigitVectorTruthInfo() {}

    void ParseTruth(const SmartRefVector<Minerva::IDCluster>& clusters,
                    const std::map<int,Minerva::TG4Trajectory*>& trajectories, int goalcode = 0);
    
    double GetTotalNormEnergy() const;
    double GetDataNormEnergy() const;
    double GetMCNormEnergy() const;
    double GetMCXtalkNormEnergy() const;
    double GetTotalTruthEnergy() const;
    double GetSharedTruthEnergy() const;

    double GetDataFraction() const;
    double GetXtalkFraction() const;
    double GetSharedFraction() const;

    double GetEdepByPdg(int pdg) const;
    double GetEdepByTrackId(int trackId) const;
    int    GetMostEvisPdg() const;
    int    GetMostEvisId() const;
    int    GetParentEvisPdg() const;
    int    GetParentEvisId() const;
    
    
    const std::map<int,double>& GetPdgEdepMap() const;
    const std::map<int,double>& GetNextToPrimaryIdEdepMap() const;
    
    void PrintMap() const;
    
  private:
    double fTotalNormEnergy;
    double fDataNormEnergy;
    double fMCNormEnergy;
    double fMCXtalkNormEnergy;
    double fTotalTruthEnergy;
    double fSharedTruthEnergy;

    std::map<int,double> fPdgEdepMap;
    std::map<int,double> fIdEdepMap;
    std::map<int,double> fParentPdgEdepMap;
    std::map<int,double> fParentIdEdepMap;
    std::map<int,double> fNextToPrimaryIdEdepMap; 
};

#endif
