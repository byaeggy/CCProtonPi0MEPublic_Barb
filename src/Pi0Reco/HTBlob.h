/*
    HTBlob.h Duplicated from CCPi0 Package on 2014-05-24
        Purpose: Make CCProtonPi0 Package independent of CCPi0 Package
        Future: Common Tools and Functions will be combined under AnaUtils or
                PionUtils
                
    Original Author:    Trung Le
    Author:             Ozgur Altinok  - ozgur.altinok@tufts.edu
    Date:               2014_05_24
    Last Revision:      2014_05_24
*/
#ifndef HTBLOB_H
#define HTBLOB_H 1

// inheritance
#include "MinervaUtils/MinervaHistoTool.h"
#include "CCProtonPi0/IHoughBlob.h"

// forwards
class IMinervaMathTool;
class TH1D;

class HTBlob: public MinervaHistoTool, virtual public IHoughBlob {

  public:

    // constructor
    HTBlob( const std::string& type, const std::string& name, const IInterface* parent );

    // destructor
    virtual ~HTBlob();

    virtual StatusCode initialize();    //< Algorithm initialization
    virtual StatusCode finalize();    //< Algorithm finalization

    StatusCode GetViewClustersAboveThreshold( SmartRefVector<Minerva::IDCluster> &idClusterVec,
                                              SmartRefVector<Minerva::IDCluster> &idClusterView, 
                                              Minerva::IDCluster::View view, double pe_cut ) const;

    StatusCode PseudoCone(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorX,
                          Gaudi::XYZVector direction, Gaudi::XYZPoint vert ) const;
    
    StatusCode XUVMatch(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorU,
                        SmartRefVector<Minerva::IDCluster> &ClusVectorV, double match ) const;
    
    StatusCode XUVMatch(SmartRef<Minerva::IDCluster> Cluster, SmartRefVector<Minerva::IDCluster> &Seed, 
                        SmartRefVector<Minerva::IDCluster> &ClusVectorU,	SmartRefVector<Minerva::IDCluster> &ClusVectorV,
                        double zmin, double zmax, double match) const;
    
    StatusCode AddClusterInsideCone(SmartRef<Minerva::IDCluster> UnuCluster, std::vector<Minerva::IDBlob*> &idBlobs, 
                                    Gaudi::XYZPoint vert) const;
    
    StatusCode Angle( SmartRef<Minerva::IDCluster> Cluster, Gaudi::XYZVector direction, Gaudi::XYZPoint vert, double &radius ) const;
    
    StatusCode Create2dHTSeed( SmartRefVector<Minerva::IDCluster> &idClusterView, SmartRefVector<Minerva::IDCluster> &HT2dClusters, 
                               double r, double theta, Gaudi::XYZPoint ref, double &spX, double &spZ) const;
    
    StatusCode GetDirection( Minerva::IDBlob *idBlob ) const;
    
    bool GetDirection( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vertex ) const;
    
    bool GetStartPosition( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vert, bool is_vertex ) const;

    StatusCode idBlobdEdx( Minerva::IDBlob *idblob, double &dEdx ) const;
    
    StatusCode isPhoton( SmartRefVector<Minerva::IDCluster> Seed, Gaudi::XYZPoint vtX ) const;
    
    void getBlobEnergyTime_Old( Minerva::IDBlob *idblob, double &energy, double& tracker_evis,double& ecal_evis, double& hcal_evis, double& scal_evis) const;
    double getBlobEnergyTime_New( Minerva::IDBlob *idblob, std::vector<double>& evis_v, std::vector<double>& energy_v) const;
    
    StatusCode invariantMass( Minerva::IDBlob* idblob1, Minerva::IDBlob* idblob2, double &mass, Gaudi::XYZPoint vert ) const;

    StatusCode energyCorrection( double &energy, double z, double radius ) const;

  private:

    IMinervaMathTool* m_mathTool;
	
    int m_planesdEdx;
    
    //vertex blob to event vertex
    double m_minDistanceStripPhoton;
    double m_minDistanceModulePhoton;
    
    // calibration constants
    double m_kT;
    double m_kE;
    double m_kS_X;
    double m_kS_UV;
    double m_kH;
    
    StatusCode FinddEdxPlanes(TH1D *h, int &index, double &dEdx) const;
    double get_kT(double evis) const;
    double getShowerEnergy(std::vector<double>& evis_v, std::vector<double>& energy_v) const;
   
};

#endif // HBLOB_H
