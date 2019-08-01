#ifndef ClusterVectorInfo_h
#define ClusterVectorInfo_h

#include <Event/MinervaEventFwd.h>

class ClusterVectorInfo {
  public:
    explicit ClusterVectorInfo(const SmartRefVector<Minerva::IDCluster>& clusters,
                               bool unusedOnly = false,
                               bool noLowActivity = false,
                               bool noLowLight = false);

    unsigned int GetN() const;
    unsigned int GetNx() const;
    unsigned int GetNu() const;
    unsigned int GetNv() const;

    unsigned int GetDigitCount() const;
    unsigned int GetDigitCountX() const;
    unsigned int GetDigitCountU() const;
    unsigned int GetDigitCountV() const;

    const SmartRefVector<Minerva::IDCluster>& GetXClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetUClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetVClusters() const;
    
  private:
    SmartRefVector<Minerva::IDCluster> fAllClusters;
    SmartRefVector<Minerva::IDCluster> fXClusters;
    SmartRefVector<Minerva::IDCluster> fUClusters;
    SmartRefVector<Minerva::IDCluster> fVClusters;

    unsigned int fDigitCountX;
    unsigned int fDigitCountU;
    unsigned int fDigitCountV;
    
    bool fUnusedOnly;
    bool fNoLowActivity;
    bool fNoLowLight;
};

#endif
