#ifndef HitVectorTruthInfo_h
#define HitVectorTruthInfo_h

#include <map>

#include <GaudiKernel/SmartRefVector.h>

namespace Minerva {
    class TG4Trajectory;
    class MCHit;
    class IDDigit;
}

class HitVectorTruthInfo {
  public:
    HitVectorTruthInfo();
    ~HitVectorTruthInfo() {}

    void ParseTruth(const SmartRefVector<Minerva::MCHit>& hits,
                    const std::map<int,Minerva::TG4Trajectory*>& trajectories);

    void ParseTruth(const SmartRef<Minerva::IDDigit>& digit,
                    const std::map<int,Minerva::TG4Trajectory*>& trajectories);
    
    double GetTotalTruthEnergy() const;
    double GetEdepByPdg(int pdg) const;
    double GetEdepByTrackId(int trackId) const;
    
    const std::map<int,double>& GetPdgEdepMap() const;
    const std::map<int,double>& GetNextToPrimaryIdEdepMap() const;

    void PrintMap() const;
    
  private:
    double fTotalTruthEnergy;

    std::map<int,double> fPdgEdepMap;
    std::map<int,double> fNextToPrimaryIdEdepMap;
    
};

#endif
