#ifndef CCProtonPi0_Pi0BlobTool_cpp
#define CCProtonPi0_Pi0BlobTool_cpp

#include "MinervaDet/DeDetector.h"
#include "MinervaDet/DePlane.h"
#include "Event/IDCluster.h"
#include "Event/IDBlob.h"
#include "Event/Node.h"
#include "Event/Track.h"

#include "Pi0BlobTool.h"

using namespace std;

Pi0BlobTool::Pi0BlobTool()
{
    // Do Nothing!
}

bool Pi0BlobTool::isBlobGood( const Minerva::IDBlob* pi0_blob)
{
    cout<<"Enter Pi0BlobTool::isBlobGood()"<<endl;
    
    bool blobGood = true;
    // Sanity Check
    if (pi0_blob == NULL) return false;
   
    FillNodesVector_and_SetTrack(pi0_blob);
    
    FillClusterVectors();

    CheckClusterVectors();
    
    bool isTrackLike = IsBlobStartTrackLike(X_clusters);
    
    if (isTrackLike){ 
        cout<<"Track Like"<<endl;
        blobGood = false;
    } else{
       cout<<"NOT Track Like"<<endl;
    }
    cout<<"Exit Pi0BlobTool::isBlobGood()"<<endl;
    
    return blobGood;
}

void Pi0BlobTool::CheckClusterVectors()
{
    //CheckClusterVector(All_clusters);
    CheckClusterVector(X_clusters);
    //CheckClusterVector(U_clusters);
    //CheckClusterVector(V_clusters);
}

void Pi0BlobTool::CheckClusterVector(std::vector<Minerva::IDCluster*> &clusters)
{
    if(clusters.size() == 0){
        cout<<"Empty Cluster -- returning"<<endl;
        return;
    }
    bool isZOrdered = isClusterOrderedinZ(clusters);
    if( isZOrdered) cout<<" Z is Ordered!"<<endl;
    else cout<<" WARNING: Z is NOT Ordered!"<<endl;
   
    int nPlanes = getNPlanes(clusters);
    cout<<"nPlanes = "<<nPlanes<<endl;

    // Debugging 
    for (unsigned int i = 0; i < clusters.size(); i++){
        cout<<"View = "<<clusters[i]->view();
        cout<<"\tZ = "<<clusters[i]->z();
        cout<<"\tPos = "<<clusters[i]->position()<<endl;
    }
    cout<<"----------------------------"<<endl;

//    if (nPlanes > 4){
//        cout<<"nPlanes = "<<nPlanes<<" entering FillClusterPosVariance()"<<endl;
//        FillClusterPosVariance(clusters);
//    }
    //cout<<"============================"<<endl;

}

void Pi0BlobTool::FillClusterPosVariance(std::vector<Minerva::IDCluster*> &clusters)
{
    double prevZ = clusters[0]->z();
    double prevPos = clusters[0]->position();
    
    double currentZ;
    double currentPos;
    std::vector<double> pos_vars;
    double pos_variance;

    for( unsigned int i = 1; i < clusters.size(); i++){
        currentZ = clusters[i]->z();
        currentPos = clusters[i]->position();
        cout<<"prevZ = "<<prevZ<<" currentZ = "<<currentZ<<endl;
        
        if (currentZ != prevZ){
            pos_variance = currentPos - prevPos; 
            cout<<"prevPos = "<<prevPos<<" currentPos = "<<currentPos<<" pos_variance = "<<pos_variance<<endl;
            pos_vars.push_back(pos_variance);
        }
        // Update prevZ
        prevZ = currentZ;
        prevPos = currentPos;
    }


    // Get Pos Variance Stats
    getPosVarianceStats(pos_vars);
    
}

void Pi0BlobTool::getPosVarianceStats(std::vector<double> &pos_vars)
{
    double total = 0;
    double avg;

    for (unsigned int i = 0; i < pos_vars.size(); i++){
        total = total + pos_vars[i];            
    }

    avg = total / (double)pos_vars.size();
    cout<<"Average Pos Variance = "<<avg<<endl;

    // Get Distance to Average
    double dist_to_avg;
    for (unsigned int i = 0; i < pos_vars.size(); i ++){
        dist_to_avg = pos_vars[i]-avg;
        cout<<"dist_to_avg = "<<dist_to_avg<<" std::abs(dist_to_avg) = "<<std::abs(dist_to_avg)<<endl;
    }
}

int Pi0BlobTool::getNPlanes(std::vector<Minerva::IDCluster*> &clusters)
{
    int nPlanes = 1;
    double prevZ = clusters[0]->z();
    double currentZ;

    for (unsigned int i = 1; i < clusters.size(); i++){
        currentZ = clusters[i]->z();
        if (currentZ != prevZ) nPlanes++;
        prevZ = currentZ;
    }

    return nPlanes;
}


bool Pi0BlobTool::isClusterOrderedinZ( std::vector<Minerva::IDCluster*> &clusters)
{
    double currentZ;
    double prevZ;
    for (unsigned int i = 1; i < clusters.size(); i++){
        prevZ = clusters[i-1]->z();
        currentZ = clusters[i]->z();
        if (currentZ >= prevZ) continue;
        else{ 
            cout<<"prevZ > currentZ"<<prevZ<<" > "<<currentZ<<endl;
            return false;
        }
    }

    return true;
}

void Pi0BlobTool::FillClusterVectors()
{
    std::vector<Minerva::IDCluster*> clusters = Pi0Blob_Track->idclusters();
    Minerva::IDCluster* tempCluster = NULL;
    
    for( unsigned int i = 0; i < clusters.size(); i ++){
        tempCluster = new Minerva::IDCluster(*(clusters[i]));
        
        All_clusters.push_back(tempCluster);
        
        // Fill Each View
        if ( tempCluster->view() == Minerva::IDCluster::X) X_clusters.push_back(tempCluster);
        else if ( tempCluster->view() == Minerva::IDCluster::U) U_clusters.push_back(tempCluster);
        else if ( tempCluster->view() == Minerva::IDCluster::V) V_clusters.push_back(tempCluster);
    }
}

void Pi0BlobTool::FillNodesVector_and_SetTrack(const Minerva::IDBlob* pi0_blob)
{
    Minerva::Node* tempNode = NULL;
    SmartRefVector<Minerva::IDCluster> clusters = pi0_blob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator iter_c;
    
    for (iter_c = clusters.begin(); iter_c != clusters.end(); ++iter_c){
        tempNode = new Minerva::Node(*iter_c);
        nodes.push_back(tempNode);
    }

    // Set Track
    Pi0Blob_Track = new Minerva::Track(nodes);
}

bool Pi0BlobTool::IsBlobStartTrackLike(std::vector<Minerva::IDCluster*> &clusters)
{
    // We need at least 3 clusters to work with
    if (clusters.size() < 3 ) return false;

    std::vector<double> plane_z;
    std::vector<double> cluster_pos;
    
    // Fill Vectors 
    for (unsigned int i = 0; i < clusters.size(); i++){
        plane_z.push_back( clusters[i]->z() );
        cluster_pos.push_back( clusters[i]->position() );
    }

    const int nPlanes_to_check = 4;
    double prevZ = plane_z[0];
    double prevPos = cluster_pos[0];
    double currentZ;
    double currentPos;
    std::vector<double> z_change;
    std::vector<double> pos_change;
    for (int i = 1; i < nPlanes_to_check; i++){
        currentZ = plane_z[i];
        currentPos = cluster_pos[i];
        
        z_change.push_back(currentZ - prevZ);
        pos_change.push_back(currentPos - prevPos);

        // Reset Prev Values
        prevZ = currentZ;
        prevPos = currentPos;
    }


    // Check Z-Change
    for(unsigned int i = 0; i < z_change.size(); i++){
        if (z_change[i] == 0) return false;
    }


    // Check Pos-Change
    double total_pos_change = 0;
    double avg_pos_change;
    double min_change = 99999;
    double max_change = -99999;
    for(unsigned int i = 0; i < pos_change.size(); i++){
        total_pos_change = total_pos_change + pos_change[i];
        if (pos_change[i] < min_change) min_change = pos_change[i];
        if (pos_change[i] > max_change) max_change = pos_change[i];
        cout<<"i = "<<i<<" pos_change = "<<pos_change[i]<<endl;
    }

    avg_pos_change = total_pos_change / (double)pos_change.size();

    cout<<"min_change = "<<min_change<<" max_change = "<<max_change<<endl;
    cout<<"avg_pos_change = "<<avg_pos_change<<endl;

    if (std::abs(max_change - min_change) > 100 ) return false;  
    else return true;
}



#endif

