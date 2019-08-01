
#include <GaudiKernel/SmartRef.h>

#include <Event/IDCluster.h>

#include "ClusterVectorInfo.h"


ClusterVectorInfo::ClusterVectorInfo(const SmartRefVector<Minerva::IDCluster>& clusters,
                                     bool unusedOnly,
                                     bool noLowActivity,
                                     bool noLowLight)
    : fAllClusters(clusters),
      fDigitCountX(0),
      fDigitCountU(0),
      fDigitCountV(0),
      fUnusedOnly(unusedOnly),
      fNoLowActivity(noLowActivity),
      fNoLowLight(noLowLight)
      
{
    
    for (SmartRefVector<Minerva::IDCluster>::const_iterator cluster = fAllClusters.begin();
         cluster != fAllClusters.end(); ++cluster ){
        if ( fUnusedOnly &&    (*cluster)->history() != Minerva::IDCluster::Unused ) continue;
        if ( fNoLowActivity && (*cluster)->type() == Minerva::IDCluster::LowActivity  ) continue;
        if ( fNoLowLight &&    (*cluster)->pe()/(*cluster)->iddigs() <= 3.0 ) continue;
        if ( (*cluster)->view() == Minerva::IDCluster::X ) {
            fDigitCountX += (*cluster)->iddigs();
            fXClusters.push_back(*cluster);
        }
        if ( (*cluster)->view() == Minerva::IDCluster::U ) {
            fDigitCountU += (*cluster)->iddigs();
            fUClusters.push_back(*cluster);
        }
        if ( (*cluster)->view() == Minerva::IDCluster::V ) {
            fDigitCountV += (*cluster)->iddigs();
            fVClusters.push_back(*cluster);
        }
    }
    

}

unsigned int ClusterVectorInfo::GetN() const {
    return GetNx() + GetNu() + GetNv();
}

unsigned int ClusterVectorInfo::GetNx() const {
    return fXClusters.size();
}

unsigned int ClusterVectorInfo::GetNu() const {
    return fUClusters.size();
}

unsigned int ClusterVectorInfo::GetNv() const {
    return fVClusters.size();
}

unsigned int ClusterVectorInfo::GetDigitCount() const {
    return fDigitCountX + fDigitCountU + fDigitCountV;
}

unsigned int ClusterVectorInfo::GetDigitCountX() const {
    return fDigitCountX;
}

unsigned int ClusterVectorInfo::GetDigitCountU() const {
    return fDigitCountU;
}

unsigned int ClusterVectorInfo::GetDigitCountV() const {
    return fDigitCountV;
}

const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetXClusters() const {
    return fXClusters;
}

const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetUClusters() const {
    return fUClusters;
}

const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetVClusters() const {
    return fVClusters;
}
