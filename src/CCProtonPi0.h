/*
   ================================================================================
   CCProtonPi0

   Reconstruction Package:
   Exclusive Channel for muon,proton,pi0 on the final state
   Uses Nightly Build

   Main Package:
   AnalysisFramework/Ana/CCProtonPi0/

Setup:
> getpack -u Ana/CCProtonPi0
> cmt config
> cmt make

Usage:
There is a corresponding Options file under Tools/SystemTests
Use Tools/ProductionScripts/ana_scripts/ProcessAna.py


Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCPROTONPI0_H 
#define CCPROTONPI0_H 1

// C++ libraries
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <utility>
#include <fstream>
#include <functional>

// ROOT Libraries
#include "TString.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TGraphErrors.h"
// Inheritance
#include "AnaUtils/MinervaAnalysisTool.h"

// MLVFTool
#include "AnaUtils/IMLVFTool.h"
#include "RecInterfaces/IRecoObjectTimeTool.h"

// Libriaries not included in MinervaAnalysisTool.h
#include "AnaUtils/IMuonUtils.h"
#include "AnaUtils/IAnalysisToolUtils.h"
#include "TruthMatcher/ITruthMatcher.h"
#include "EnergyRecTools/ICalorimetryUtils.h"
#include "AnaUtils/IDeadTimeTool.h"
#include "AnaUtils/AnaFilterTags.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

// Local
#include "Helper/PDG.h"
#include "Pi0Reco/AngleScan.h"
#include "Pi0Reco/AngleScan_U.h"
#include "Pi0Reco/AngleScan_V.h"
#include "Pi0Reco/ClusterVectorInfo.h"
#include "Pi0Reco/OneParLineFit.h"
#include "Pi0Reco/TwoParLineFit.h"
#include "Pi0Reco/Pi0BlobTool.h"
#include "TruthMatch/DigitVectorTruthInfo.h"
#include "TruthMatch/TraverseHistory.h"
#include "TruthMatch/TrackTruthInfo.h"
#include "TruthMatch/HitVectorTruthInfo.h"
#include "CCProtonPi0/IHoughBlob.h"
#include "CCProtonPi0/IHoughTool.h"

// Gaudi
#include "GaudiKernel/PhysicalConstants.h"

// Minerva Analysis Framework
#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/MCTrack.h"
#include "BadChannels/IGetDeadTime.h"
#include "BlobFormation/IBlobCreatorUtils.h"
#include "BlobFormation/IIDAnchoredBlobCreator.h"
#include "BlobFormation/IIDBlobSeedingTool.h"
#include "CalTools/IGetCalAttenuation.h"
#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "EnergyRecTools/IExtraEnergyTool.h"
#include "Event/GenMinHeader.h"
#include "Event/TimeSlice.h"
#include "Event/MCIDDigit.h"
#include "Event/MCHit.h"
#include "Event/TG4Trajectory.h"
#include "Event/TG4TrajectoryPointExtraInfo.h"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "GiGaCnv/IGiGaGeomCnvSvc.h"
#include "MinervaDet/DeDetector.h"
#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaUtils/IHitTaggerTool.h"
#include "MinervaUtils/IMinervaMathTool.h"
#include "MinervaUtils/IMinervaObjectAssociator.h"
#include "ODDet/DeOuterDetector.h"
#include "ParticleMaker/IParticleMakerTool.h"
#include "ParticleMaker/IParticleTool.h"
#include "ParticleMaker/ILikelihoodParticleTool.h"
#include "ProngMaker/IMichelTool.h"
#include "ProngMaker/IImprovedMichelTool.h"
#include "ProngMaker/IODProngClassificationTool.h"
#include "RecInterfaces/IAnchoredTrackFormation.h"
#include "RecInterfaces/IFiducialPointTool.h"
#include "RecInterfaces/IRecoObjectTimeTool.h"
#include "RecInterfaces/ITrackLinearPropagator.h"
#include "RecUtils/Cone.h"
#include "RecUtils/IConeUtilsTool.h"
#include "RecoStudies/IVertexEnergyStudyTool.h"
#include "VertexCreation/IVertexFitter.h"


class IAnchoredTrackFormation;
class IBlobCreatorUtils;
class ICalorimetryUtils;
class ITrackLinearPropagator;
class IConeUtilsTool;
class IEnergyCorrectionTool;
class IExtraEnergyTool;
class IGeomUtilSvc;
class IGetCalAttenuation;
class IGetDeadTime;
class IGiGaGeomCnvSvc;
class IHitTaggerTool;
class IHoughBlob;
class IHoughTool;
class IIDAnchoredBlobCreator;
class IIDBlobSeedingTool;
class IMichelTool;
class IImprovedMichelTool;
class IMinervaCoordSysTool;
class IMinervaMathTool;
class IMinervaObjectAssociator;
class IODProngClassificationTool;
class IParticleMakerTool;
class IParticleTool;
class IProtonUtils;
class IRecoObjectTimeTool;
class IVertexEnergyStudyTool;
class IVertexFitter;
class TRandom3;
class IRecoObjectTimeTool;

namespace Minerva {
    class DeDetector;
    class DeOuterDetector;
}

//! This class is for Reconstruct Pi0 using muon match vertex
class CCProtonPi0 : public MinervaAnalysisTool
{
private:
    typedef std::vector<Minerva::NeutrinoInt*> NeutrinoVect;

public:

    //! Standard constructor
    CCProtonPi0(const std::string& type, 
            const std::string& name, 
            const IInterface* parent );

    //! Destructor (mandatory for inheritance)
    virtual ~CCProtonPi0(){};

    StatusCode initialize();
    StatusCode finalize();

    //! Reconstruct the event (mandatory for inheritance)
    StatusCode reconstructEvent(Minerva::PhysicsEvent* event, 
            Minerva::GenMinInteraction* truthEvent = NULL ) const;

    //! Attach an interpretations to the event (mandatory for inheritance)
    StatusCode interpretEvent(  const Minerva::PhysicsEvent* event, 
            const Minerva::GenMinInteraction* truthEvent, 
            NeutrinoVect& interaction_hyp ) const;

    StatusCode tagTruth( Minerva::GenMinInteraction* truthEvent ) const;

    bool truthIsPlausible( const Minerva::PhysicsEvent*) const;

private:
    // mutable Member Variables which can be modified by const functions
    mutable Gaudi::LorentzVector m_muon_4P;
    mutable Gaudi::LorentzVector m_proton_4P;
    mutable Gaudi::LorentzVector m_pion_4P;
    mutable Gaudi::LorentzVector m_pi0_4P;
    mutable SmartRef<Minerva::Prong>    m_MuonProng;
    mutable SmartRef<Minerva::Particle> m_MuonParticle;
    mutable SmartRef<Minerva::IDBlob>   m_Pi0Blob1;
    mutable SmartRef<Minerva::IDBlob>   m_Pi0Blob2;
    mutable Minerva::ProngVect    m_ProtonProngs;
    mutable Minerva::ParticleVect m_ProtonParticles;
    mutable Minerva::ProngVect    m_PionProngs;
    mutable Minerva::ParticleVect m_PionParticles; 
    mutable Minerva::ProngVect    m_HadronProngs;
    mutable Minerva::ParticleVect m_HadronParticles;
    
// Counters for Functions - Debugging Purposes
    mutable double N_tagTruth;
    mutable double N_reconstructEvent;

    // Truth Match
    mutable std::map<int, Minerva::TG4Trajectory*> fTrajectoryMap;
    mutable const Minerva::TG4Trajectory* fPizero;
    mutable const Minerva::TG4Trajectory* fGamma1;
    mutable const Minerva::TG4Trajectory* fGamma2;

    std::vector<std::string>   m_dedx_uncertainties;

    // Fiducial Volume
    double m_fidHexApothem;
    double m_fidUpStreamZ;
    double m_fidDownStreamZ;

    // Analyzable Volume
    double m_recoHexApothem;
    double m_recoUpStreamZ;
    double m_recoDownStreamZ;

    // Signal Definition Requirements
    double m_W_max;
    double m_Enu_min;
    double m_Enu_max;
    double m_muon_P_min;
    double m_muon_theta_max;

    // Analysis Parameters
    double m_beamAngleBias;

    // Keep Events After Each Cut Set
    bool m_keepAfter_vertex_cuts;
    bool m_keepAfter_muon_cuts;
    bool m_keepAfter_michel_cuts;
    bool m_keepAfter_proton_cuts;
    bool m_keepAfter_pi0_cuts;

    bool m_removeEvents_withMichel;

    // Algorihm Flow
    bool m_writeFSParticle_Table;
    bool m_DoPlausibilityCuts;
    bool m_DoTruthMatch;

    // Optional Studies
    bool m_study_michel_electron;
    bool m_study_shower_energy;
    bool m_study_shower_recovery;
    bool m_study_unused_energy;

    // Prong and Cluster Colors
    int m_Color_muonProng;
    int m_Color_protonProng;
    int m_Color_primaryVertex;
    int m_Color_secondaryVertex;
    int m_Color_endPointVertex;
    int m_Color_unattachedProng;
    int m_Color_Gamma1Prong;
    int m_Color_Gamma2Prong;
    int m_Color_GammaOtherProng;
    int m_Color_clusterUnused;
    int m_Color_clusterUsed;
    int m_Color_VertexBlob;
    int m_Color_RejectedBlob;

    // VertexBlob
    double m_vertex_blob_radius;
    bool fSkipLowEnergyClusterVtxEnergy;
    bool fThresholdVertexEnergy;

    // ConeBlobs
    double m_rejectedClustersTime; ///< window time  to allow clusters

    bool m_ApplyAttenuationCorrection;
    bool   m_AllowUVMatchWithMoreTolerance;
    double m_UVMatchTolerance;
    double m_UVMatchMoreTolerance;

        bool m_TrytoRecover_1Shower;
        bool m_TrytoRecover_3Shower;
	bool m_TrytoRecover_4Shower;
        bool m_recoverSingleShower_SmallAngle;
        bool m_recoverSingleShower_SearchView;
        bool m_recoverShower_invMass;
        bool m_recoverShower_Direction;
    

    TRandom3*                 m_randomGen;
    unsigned long int         m_randomSeed;

    Minerva::DeDetector*        m_InnerDetector;
    Minerva::DeOuterDetector*   m_OuterDetector;

    // Duplicated variables -- will fix in future
    const Minerva::DeDetector*      m_idDet;                ///< Inner detector
    const Minerva::DeOuterDetector* m_odDet;                ///< Outer detector

    IAnchoredTrackFormation* m_anchoredTracker;
    IBlobCreatorUtils* m_blobUtils;
    ITrackLinearPropagator* m_trackPropagator;
    ICalorimetryUtils* m_caloUtils;
    IConeUtilsTool* m_coneUtilsTool;
    IEnergyCorrectionTool* m_energyCorrectionTool;
    IExtraEnergyTool* m_extraEnergyTool;
    IGeomUtilSvc* m_GeomUtilSvc;
    IGetCalAttenuation* m_AttenuationCorrectionTool;
    IGetDeadTime* m_getDeadTimeTool;
    IGiGaGeomCnvSvc* m_gigaCnvSvc;
    IHitTaggerTool* m_hitTagger;
    IIDBlobSeedingTool* m_blobSeedingTool;
    IHoughBlob* m_idHoughBlob;
    IHoughTool* m_idHoughTool;
    IIDAnchoredBlobCreator* m_stopPointBlobTool;
    IMichelTool* m_michelTool_Default;
    IMichelTool* m_michelTool_Large;
    IImprovedMichelTool* m_michelTool_Improved;
    IMinervaCoordSysTool* m_coordSysTool;
    IMinervaMathTool* m_mathTool;
    IMinervaObjectAssociator* m_objectAssociator;
    IODProngClassificationTool* m_odMatchTool;
    IParticleMakerTool* m_particleMaker;
    ILikelihoodParticleTool* m_LikelihoodPIDTool;
    std::string                 m_LPIDToolAlias;
    IProtonUtils* m_protonUtils;
    IRecoObjectTimeTool* m_recoTimeTool;
    IVertexEnergyStudyTool* m_vertexEnergyStudyTool;
    IVertexFitter* m_vertexFitter;

IRecoObjectTimeTool* m_timeTool;
	std::string m_timeToolAlias;

        IParticleTool* m_particleTool;
		std::string           m_particleToolName;
		std::string           m_particleToolAlias;

    //! Private Functions
    Minerva::IDClusterVect getClusters( Minerva::PhysicsEvent* event ) const;
    StatusCode ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const;
    StatusCode getNearestPlane( double z, int & module_return, int & plane_return) const;
    StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
    bool AreBlobsDirectionGood(Minerva::PhysicsEvent *event) const;
    bool ConeBlobs( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const;
    bool PreFilterPi0( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const;
    bool ShouldReconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent) const;
    bool TrackEndPointHasMichels(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent) const;
    bool ImprovedMichel_EventHasMichel(Minerva::PhysicsEvent *event) const;
    bool ImprovedMichel_HasMichelAtPoint(SmartRef<Minerva::Vertex> point) const;
    bool ImprovedMichel_HasMichelAtPoint(double x, double y, double z) const;
    bool VertexHasMichels(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent) const;
    bool checkMichel(Minerva::GenMinInteraction* truthEvent) const;
    bool createTrackedParticles(Minerva::PhysicsEvent *event ) const;
    bool createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex, bool make_primary_short_tracks ) const;
    bool getProtonProng(Minerva::PhysicsEvent *event ) const;
    bool getPionProngs(Minerva::PhysicsEvent *event ) const; 
    bool hasEventMinosMatchedMuon(Minerva::PhysicsEvent *event) const;
    bool hasEventVertex(Minerva::PhysicsEvent *event) const;
    bool isInteractionNC(Minerva::GenMinInteraction* truthEvent) const;
    bool isMichelProngGood(Minerva::Prong &michelProng) const;
    bool isMotherPrimary(std::vector<int>& motherList, int mother ) const;
    bool isMuonChargeNegative(Minerva::PhysicsEvent *event) const;
    bool isMuonPlausible(Minerva::PhysicsEvent *event) const;
    bool isTrueVertexFiducial(Minerva::GenMinInteraction* truthEvent) const;
    bool setMuonData (Minerva::PhysicsEvent *event ) const;
    bool setPi0Data( Minerva::PhysicsEvent *event) const;
    bool setProtonData( Minerva::PhysicsEvent *event ) const;
    bool setPionData( Minerva::PhysicsEvent *event) const;
    bool tagSignal( Minerva::GenMinInteraction* truthEvent) const;
    bool tagSignal_BeforeFSI( Minerva::GenMinInteraction* truthEvent) const;
    bool vertexInFiducialVolume(Minerva::PhysicsEvent *event) const;
    bool vertexInRecoVolume(Minerva::PhysicsEvent *event) const;
    double CalculatePhotonScore(double dEdX, double Energy) const;
    double CalculateConeAngle(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    double MeanEnergyClusters_1stHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    double MeanEnergyClusters_2ndHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    double MedianEnergyClusters_1stHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    double MedianEnergyClusters_2ndHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    std::vector<std::pair<double, double> >EnergyDistanceClusters(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;

    double CalProbFitLine(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const;
    double GetVertexEnergyCalConst(double evis) const;
    double GetVertexEnergy(const Minerva::PhysicsEvent *event) const;
    void SaveExtraEnergyLeftover(Minerva::PhysicsEvent *event) const;

     //! Private Functions
        Gaudi::LorentzVector Get_Neutrino_4P(double Enu) const;
        StatusCode HoughBlob( Minerva::PhysicsEvent *event, SmartRefVector<Minerva::IDCluster> idClusters, std::vector<Minerva::IDBlob*>& outBlobs) const;
        SmartRefVector<Minerva::IDCluster> FilterInSphereClusters( Minerva::PhysicsEvent *event, const SmartRefVector<Minerva::IDCluster>& clusters, const double sphereRadius) const;

    
    bool isSignal_CCNeutrino(Minerva::GenMinInteraction* truthEvent) const;
    bool isSignal_AcceptanceGood(Minerva::GenMinInteraction* truthEvent, double muon_theta_beam) const;
    bool isSignal_KinematicsGood(Minerva::GenMinInteraction* truthEvent) const;

    // Function to improve the blob direction
    double GetTotalpe(SmartRefVector<Minerva::IDCluster> clusters) const;
    void GetBlobDirection(Minerva::PhysicsEvent *event, Minerva::IDBlob *blob) const;
    double Calc_Enu_1Track(double vertex_energy, double extra_energy) const;
    double Calc_Enu_1Track_Alt() const;
    double Calc_Enu_2Track(double vertex_energy, double extra_energy) const;
    void Calc_TrueEventKinematics(Minerva::GenMinInteraction* truthEvent) const;
    double CalcMinBlobSeparation( const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const;
    double Calc_Longitudinal_Momentum(Gaudi::LorentzVector particle_4P) const;
    double Calc_Px_wrt_Beam(Gaudi::LorentzVector particle_4P) const;
    double Calc_Py_wrt_Beam(Gaudi::LorentzVector particle_4P) const;
    double Calc_QSq(double Enu, double muon_E, double muon_P, double muon_angle_beam) const;
    double Calc_WSq(double Enu, double QSq, double muon_E) const;
    double FindShowerMaxZ(SmartRef<Minerva::IDBlob> blob) const; 
    double TwoParLineFitBlobVtxDistance(const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const;
    double calcDistance( double x1, double y1, double z1,double x2, double y2, double z2) const;
    int getMichelPion(std::vector<int>& piList, int ID ) const;
    std::pair<int,double> OneParLineFitBlob(const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const;
    void PrintZPositionMaterial() const;
    void ApplyAttenuationCorrection(Minerva::IDBlob* blob) const;
    int Calculate_BlobdEdX( Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob, unsigned int blob_number) const;
    void Calculate_dEdx( Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob, unsigned int blob_number) const;
    void ColorUnusedIDClusters(Minerva::PhysicsEvent *event) const;
    void DiscardFarTracks(Minerva::PhysicsEvent *event) const;
    void DispersedBlob( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent ) const;
    void FillTrajectoryMap() const;
    void GetMuonExtraEnergy(Minerva::PhysicsEvent *event) const;
    void MarkFoundPi0Blobs(std::vector<Minerva::IDBlob*> &foundBlobs) const;
    void Check_Quality_of_Blobs(std::vector<Minerva::IDBlob*> &foundBlobs, std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent) const;
    void RefitVertex_Using_AnchoredShortTracks(Minerva::PhysicsEvent *event) const;
    void SaveEventTime(Minerva::PhysicsEvent *event) const;
    void SaveTruthUnusedClusterEnergyInsideDetector(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> ecalClusters, SmartRefVector<Minerva::IDCluster> hcalClusters, SmartRefVector<Minerva::IDCluster> otherClusters) const;
    void SaveTruthClusterEnergy_FoundBlobs(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent) const;
    void SaveTruthUnusedClusterEnergy_NearVertex(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> vertexClusters) const;
    void SaveTruthUnusedClusterEnergy_Dispersed(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> clusters) const;
    void SaveTruthUnusedClusterEnergy_Rejected(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> clusters) const;
    void setSignalKinematics(Minerva::GenMinInteraction* truthEvent) const;
    void SetVertexCount(Minerva::PhysicsEvent *event) const;
    void SetExtraEnergyTotal(Minerva::PhysicsEvent* event) const;
    void VertexBlob( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent ) const;
    void correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, double& p_calCorrection, double& p_visEnergyCorrection ) const;
    void correctPionProngEnergy(  SmartRef<Minerva::Prong>& pionProng, double& p_calCorrection, double& p_visEnergyCorrection ) const;
    void processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const;
    void saveMichelElectron(Minerva::GenMinInteraction* truthEvent, int muon_ID) const;
    void saveMichelProngToNTuple(Minerva::PhysicsEvent* event, Minerva::Prong &michelProng, bool isVertex, bool isLargeSearch) const;
    void saveMichelProngTruth(Minerva::GenMinInteraction* truthEvent, Minerva::Prong &michelProng, bool isVertex, bool isLargeSearch) const;
    void setBlobData(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction *truthEvent) const;
    void setEventKinematics(const Minerva::PhysicsEvent *event, Minerva::NeutrinoInt* nuInt) const;

    void setPi0GenieRecord(Minerva::GenMinInteraction* truthEvent) const;
    void setTargetMaterial(Minerva::GenMinInteraction* truthEvent) const;
    void setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const;
    void setTrackProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const;
    void setVertexData(Minerva::PhysicsEvent* event ) const;
    void tagBackground(Minerva::GenMinInteraction* truthEvent) const;
    void tagBackground_Compact(Minerva::GenMinInteraction* truthEvent) const;
    void tagBackgroundWithPi0(Minerva::GenMinInteraction* truthEvent) const;
    void tagPrimaryMuon(Minerva::PhysicsEvent *event) const;
    void writeBackgroundType(Minerva::GenMinInteraction* truthEvent) const;
    void writeEventRecord(Minerva::GenMinInteraction* truthEvent, bool isSignal) const;
    void writeFSParticleTable(bool isSignal) const;
    void setSignal_SecondaryTrajectoryKinematics(Minerva::GenMinInteraction *truthEvent, int &Pi0_ID) const;
    void setSignal_PrimaryTrajectoryKinematics(Minerva::GenMinInteraction *truthEvent, int &Pi0_ID) const;
    int GetMCHitPDG(const SmartRef<Minerva::MCHit> mc_hit) const;
    int GetDigitPDG(const Minerva::MCIDDigit* mcdigit) const;
    int GetPDGfromVector(std::vector<int> &all_pdg) const;
    void FillUsableClusters(SmartRefVector<Minerva::IDCluster> &usableClusters, Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const;
    void ProcessRejectedClusters(SmartRefVector<Minerva::IDCluster> &rejectedClusters,Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const;
    bool RecoverShowers_InvMass( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const;
    bool RecoverShowers_Direction( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const;
    bool Save_1ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const;
    bool Save_3ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const;
    void Save_1ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const;
    void Save_3ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const;
    bool Save_AllShowersInfo( std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::PhysicsEvent *event) const;
    void Save_AllShowersTruthMatch( std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::GenMinInteraction* truthEvent) const;
void Save_4ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const;
bool Save_4ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const;
    bool isShowerGood(Minerva::IDBlob* shower, Minerva::PhysicsEvent* event) const;
    bool RecoverSingleShower_SmallAngle(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const;
    bool RecoverSingleShower_View_U(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const;
    bool RecoverSingleShower_View_V(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const;
    double GetTotalExtraEnergy(const Minerva::PhysicsEvent* event) const;
    double getClusterEnergy( SmartRefVector<Minerva::IDCluster> clusters) const;
    std::vector<double> EstimateShowerEndPoint(const Gaudi::LorentzVector& shower_4P, const std::vector<double> init_pos, double max_z) const;
    void SaveTrajectories(Minerva::PhysicsEvent* event) const;

    // --------------------------------------------------------------------
    // Study: Shower Energy Functions
    //      See Studies/CCProtonPi0_Study_ShowerEnergy.cpp
    // --------------------------------------------------------------------
    void SaveBlobDigitEnergy(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const;
    bool Use_SCAL_minZ(Minerva::IDBlob* blob, double minZ) const;
    bool isHitInsideSCAL(double x, double y, double z) const;
    double Find_SCAL_MinZ(Minerva::IDBlob* blob) const;
    void SaveBlobTrueEvisBySubDetector(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const;
    void SaveSCALHits(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const;
    void SaveSCALHits_Improved(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const;
    void SaveSCAL_minZ_Info(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const;
};

#endif // CCPROTONPI0_H 

