#include <set>
#include <iterator>
#include <iostream>

#include <Event/MCHit.h>
#include <Event/IDDigit.h>
#include <Event/MCIDDigit.h>
#include <Event/TG4Trajectory.h>


#include "HitVectorTruthInfo.h"
#include "TraverseHistory.h"

HitVectorTruthInfo::HitVectorTruthInfo()
    : fTotalTruthEnergy(0.0)
{}

double HitVectorTruthInfo::GetTotalTruthEnergy() const {
    return fTotalTruthEnergy;
}

double HitVectorTruthInfo::GetEdepByPdg(int pdg) const {
    std::map<int,double>::const_iterator target = fPdgEdepMap.find(pdg);
    if (target != fPdgEdepMap.end()) return target->second;

    return 0.0; /* If the particle pdg does not exist */
}

double HitVectorTruthInfo::GetEdepByTrackId(int trackId) const {
    std::map<int,double>::const_iterator target = fNextToPrimaryIdEdepMap.find(trackId);
    if (target != fNextToPrimaryIdEdepMap.end()) return target->second;

    return 0.0; /* If trackId does not exist */
}

const std::map<int,double>& HitVectorTruthInfo::GetPdgEdepMap() const {
    return fPdgEdepMap;
}

const std::map<int,double>& HitVectorTruthInfo::GetNextToPrimaryIdEdepMap() const {
    return fNextToPrimaryIdEdepMap;
}

void HitVectorTruthInfo::ParseTruth(const SmartRef<Minerva::IDDigit>& digit_ref,
                                    const std::map<int, TG4Trajectory*>& trajectories)
{
    const Minerva::IDDigit* digit = digit_ref;
    const Minerva::MCIDDigit* mcdigit
        = dynamic_cast<const Minerva::MCIDDigit*>(digit);

    ParseTruth(mcdigit->hits(),trajectories);
}

void HitVectorTruthInfo::ParseTruth(const SmartRefVector<Minerva::MCHit>& hits,
                                    const std::map<int, TG4Trajectory*>& trajectories)
{

    for (SmartRefVector<Minerva::MCHit>::const_iterator h = hits.begin();
         h != hits.end(); ++h) {
        const double energy = (*h)->energy();
        fTotalTruthEnergy += energy;
        
        TraverseHistory history((*h)->GetTrackId());
        history.DoTraverse(trajectories);
        int primaryPdg      = history.GetPrimaryPdg();
        int nextToPrimaryId = history.GetNextToPrimaryId(); 
        
        fPdgEdepMap[primaryPdg] += energy;
        fNextToPrimaryIdEdepMap[nextToPrimaryId] += energy;
    }
 
}

void HitVectorTruthInfo::PrintMap() const
{
    std::cout << "Primary particle PDG | energy deposit " << std::endl;
    for (std::map<int,double>::const_iterator p = fPdgEdepMap.begin();
         p != fPdgEdepMap.end(); ++p) {
        std::cout << "\t--> " << p->first << " " << p->second << std::endl;
    }
}
