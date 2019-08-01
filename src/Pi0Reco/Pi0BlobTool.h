#ifndef CCProtonPi0_Pi0BlobTool_h
#define CCProtonPi0_Pi0BlobTool_h

#include <iostream>
#include <vector>

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class Pi0BlobTool
{
    public:
        Pi0BlobTool();
        bool isBlobGood( const Minerva::IDBlob* pi0_blob);

    private:
        std::vector<Minerva::IDCluster*> All_clusters;
        std::vector<Minerva::IDCluster*> X_clusters;
        std::vector<Minerva::IDCluster*> U_clusters;
        std::vector<Minerva::IDCluster*> V_clusters;
        std::vector<Minerva::Node*> nodes;
        SmartRef<Minerva::Track> Pi0Blob_Track;
       
        bool IsBlobStartTrackLike(std::vector<Minerva::IDCluster*> &clusters);
        void CheckClusterVectors();
        void FillClusterVectors();
        void FillNodesVector_and_SetTrack(const Minerva::IDBlob* pi0_blob);
        void CheckClusterVector( std::vector<Minerva::IDCluster*> &clusters);
        bool isClusterOrderedinZ( std::vector<Minerva::IDCluster*> &clusters);
        int getNPlanes(std::vector<Minerva::IDCluster*> &clusters);
        void FillClusterPosVariance(std::vector<Minerva::IDCluster*> &clusters);
        void getPosVarianceStats(std::vector<double> &pos_vars);

};

#endif

