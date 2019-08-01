#include <set>
#include <iterator>
#include <iostream>
#include <limits>

#include <Event/Track.h>
#include <Event/Node.h>
#include <Event/IDCluster.h>
#include <Event/IDDigit.h>
#include <Event/MCIDDigit.h>
#include <Event/MCHit.h>
#include <Event/TG4Trajectory.h>


#include "TrackTruthInfo.h"
#include "DigitVectorTruthInfo.h"


namespace {
    const double epsilon = std::numeric_limits<double>::epsilon();
}

TrackTruthInfo::TrackTruthInfo()
    : fTotalNormEnergy(0.0),
      fDataNormEnergy(0.0),
      fMCNormEnergy(0.0),
      fMCXtalkNormEnergy(0.0),
      fTotalTruthEnergy(0.0),
      fSharedTruthEnergy(0.0),
      fEdepFraction(1.e9)
{}

double TrackTruthInfo::GetTotalNormEnergy() const {
    return fTotalNormEnergy;
}
double TrackTruthInfo::GetDataNormEnergy() const {
    return fDataNormEnergy;
}
double TrackTruthInfo::GetMCNormEnergy() const {
    return fMCNormEnergy;
}
double TrackTruthInfo::GetMCXtalkNormEnergy() const {
    return fMCXtalkNormEnergy;
}
double TrackTruthInfo::GetTotalTruthEnergy() const {
    return fTotalTruthEnergy;
}
double TrackTruthInfo::GetSharedTruthEnergy() const {
    return fSharedTruthEnergy;
}

double TrackTruthInfo::GetDataFraction() const {
    return fDataNormEnergy/fTotalNormEnergy;
}
double TrackTruthInfo::GetXtalkFraction() const {
    return fMCXtalkNormEnergy/fMCNormEnergy;
}
double TrackTruthInfo::GetSharedFraction() const {
    if (fTotalTruthEnergy < epsilon) return 0.0;

    return fSharedTruthEnergy/fTotalTruthEnergy;
    
}
unsigned int TrackTruthInfo::GetTruthParticleCount() const {
    return fEdepPdgMap.size();
}
int TrackTruthInfo::GetTruthPdg(unsigned int index) const {

        /* This should never happen for MC */
    if (fEdepPdgMap.empty()) {
        fEdepFraction = 0.0;
        return -1;
    }

        /* Handle invalid index */
    if (index > fEdepPdgMap.size()-1) {
        fEdepFraction = 0.0;
        return -1;
    }

    std::map<double,int>::const_reverse_iterator p = fEdepPdgMap.rbegin();
    std::advance(p,index);
        
    fEdepFraction = p->first/fTotalTruthEnergy;

    return p->second;
    
}
double TrackTruthInfo::GetEdepFraction() const {
    return fEdepFraction;
}


void TrackTruthInfo::ParseTruth(const SmartRef<Minerva::Track>& track,
                                const std::map<int, Minerva::TG4Trajectory*>& trajectories)
{

    if (trajectories.empty()) return;


    SmartRefVector<Minerva::IDCluster> trackClusters;
    
    const Minerva::Track::NodeContainer& nodes = track->nodes();
    for (Minerva::Track::NodeContainer::const_iterator node = nodes.begin();
             node != nodes.end(); ++node) {
        SmartRef<Minerva::IDCluster> cluster((*node)->idcluster());
        trackClusters.push_back(cluster);
    }
    
    DigitVectorTruthInfo digitVectorInfo;
    digitVectorInfo.ParseTruth(trackClusters,trajectories);
    
    fTotalNormEnergy   += digitVectorInfo.GetTotalNormEnergy();
    fDataNormEnergy    += digitVectorInfo.GetDataNormEnergy();
    fMCNormEnergy      += digitVectorInfo.GetMCNormEnergy();
    fMCXtalkNormEnergy += digitVectorInfo.GetMCXtalkNormEnergy(); 
    fTotalTruthEnergy  += digitVectorInfo.GetTotalTruthEnergy();
    fSharedTruthEnergy += digitVectorInfo.GetSharedTruthEnergy();
    
    const std::map<int, double>& pdgEdepMap = digitVectorInfo.GetPdgEdepMap();
    for (std::map<int, double>::const_iterator p = pdgEdepMap.begin();
         p != pdgEdepMap.end(); ++p) {
        const int primaryPdg = p->first;
        const double amount  = p->second;
        std::map<int,double>::iterator target = fPdgEdepMap.find(primaryPdg);
        if (target == fPdgEdepMap.end()) {
            fPdgEdepMap.insert(std::make_pair(primaryPdg,amount));
        } else {
            target->second += amount;
        }
    }
            

    
    std::cout << "\t\tTotal truth energy: " << fTotalTruthEnergy << std::endl;
    for (std::map<int,double>::iterator p = fPdgEdepMap.begin();
         p != fPdgEdepMap.end(); ++p) {
        const int pdg       = p->first;
        const double amount = p->second;
        std::cout << "\t\t\tTrackTruthInfo: " << pdg << " " << amount << std::endl;
        fEdepPdgMap.insert(std::make_pair(amount,pdg));
    }


    if (fEdepPdgMap.empty()) {
        std::cerr << "Warning: no MC truth particle produce this track!" << std::endl;
    }

}
