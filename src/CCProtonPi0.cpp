/*
   See CCProtonPi0.h header for Class Information
   */
#ifndef CCProtonPi0_cpp 
#define CCProtonPi0_cpp 1

#include "CCProtonPi0.h"
#include <AnaUtils/IMuonUtils.h>
#include <AnaUtils/IAnalysisToolUtils.h>
#include <AnaUtils/IDeadTimeTool.h>
#include <AnaUtils/IRecoilUtils.h>
#include <AnaUtils/AnaFilterTags.h>
#include <EnergyRecTools/ICalorimetryUtils.h>
#include <TruthMatcher/ITruthMatcher.h>

TF1 *f1_photon    = new TF1("f1_photon",   "[0]*TMath::Power(x,6)+[1]*TMath::Power(x,5)+[2]*TMath::Power(x,4)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,2)+[5]*x+[6]", 0.0, 3000.0);
TF1 *std1_photon  = new TF1("std1_photon", "[0]*TMath::Power(x,6)+[1]*TMath::Power(x,5)+[2]*TMath::Power(x,4)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,2)+[5]*x+[6]", 0.0, 3000.0);
TF1 *f1_photon1   = new TF1("f1_photon1",  "[0]*x+[1]", 0.0, 3000.0);
TF1 *std1_photon1 = new TF1("std1_photon1","[0]*x+[1]", 0.0, 3000.0);

TF1 *lineX = new TF1("lineX", "pol1", 5000, 10000);
TF1 *lineU = new TF1("lineU", "pol1", 5000, 10000);
TF1 *lineV = new TF1("lineV", "pol1", 5000, 10000);
TF1 *fit_lineX;
TF1 *fit_lineU;
TF1 *fit_lineV;

TGraphErrors *graphX;
TGraphErrors *graphU;
TGraphErrors *graphV;

namespace
{
    class increasingDistanceFromVertex : public std::binary_function
                                         <
                                         SmartRef<Minerva::IDCluster>,
                                         SmartRef<Minerva::IDCluster>,
                                         bool
                                         > {
                                             public:
                                                 explicit increasingDistanceFromVertex(Gaudi::XYZPoint pos) : fVertex(pos) {}

                                                 bool operator() (const SmartRef<Minerva::IDCluster>& lhs,
                                                         const SmartRef<Minerva::IDCluster>& rhs) const {


                                                     const double z0 = fVertex.Z();

                                                     return std::abs(lhs->z()-z0) < std::abs(rhs->z()-z0);
                                                 }

                                             private:
                                                 Gaudi::XYZPoint fVertex;
                                         };

}

DECLARE_TOOL_FACTORY( CCProtonPi0 );

using namespace Minerva;

//==============================================================================
// Standard constructor 
//==============================================================================
CCProtonPi0::CCProtonPi0(const std::string& type, const std::string& name, const IInterface* parent ) :
    MinervaAnalysisTool( type, name, parent ) 
{
    declareInterface<IInteractionHypothesis>(this);
    // Mandatory declaration of analysis signature: CCProtonPi0
    m_anaSignature = "CCProtonPi0";

    // Private Properties
    declareProperty("WriteFSParticleTable", m_writeFSParticle_Table = false);

    declareProperty("KeepAfterVertexCuts", m_keepAfter_vertex_cuts = true);
    declareProperty("KeepAfterMuonCuts", m_keepAfter_muon_cuts = true);
    declareProperty("KeepAfterMichelCuts", m_keepAfter_michel_cuts = true);
    declareProperty("KeepAfterProtonCuts", m_keepAfter_proton_cuts = true);
    declareProperty("KeepAfterPi0Cuts", m_keepAfter_pi0_cuts = true);
//declareProperty("KeepAfterPionCuts",   m_keepAfter_pion_cuts   = true);        // I have to change this in the opts file
    declareProperty("RemoveEventsWithMichel", m_removeEvents_withMichel = false);

    declareProperty("DoPlausibilityCuts",   m_DoPlausibilityCuts    =   true);
    declareProperty("DoTruthMatch",         m_DoTruthMatch          =   true);

    declareProperty("BeamAngleBias",       m_beamAngleBias = 0.006*CLHEP::radian);

    // Signal Definition Requirements
    declareProperty("Signal_Min_NeutrinoE", m_Enu_min = 1500.0);
    declareProperty("Signal_Max_NeutrinoE", m_Enu_max = 20000.0);
    declareProperty("Signal_Max_W", m_W_max = 20000.0);
    declareProperty("Signal_Max_MuonAngle", m_muon_theta_max = 25.0);

    // Optional Studies
    declareProperty("StudyMichelElectron", m_study_michel_electron = true);
    declareProperty("StudyShowerEnergy", m_study_shower_energy = false);
    declareProperty("StudyShowerRecovery", m_study_shower_recovery = false);
    declareProperty("StudyUnusedEnergy", m_study_unused_energy = false);

  // Pion Tools
    declareProperty( "ParticleToolName", m_particleToolName  = "dEdXTool" );
    declareProperty( "ParticleToolAlias",m_particleToolAlias = "dEdXTool" );
    //--------------------------------------------------------------------------
    // Arachne Colors
    //--------------------------------------------------------------------------
    // Prongs and Tracks
    declareProperty("MuonProngColor",           m_Color_muonProng       = 0x228B22); //-- green
    declareProperty("ProtonProngColor",         m_Color_protonProng     = 0x0000FF); //-- blue
    declareProperty("TrackEndProngColor",       m_Color_endPointVertex  = 0xFFA500); //-- orange
    declareProperty("UnattachedProngColor",     m_Color_unattachedProng = 0xFFFF00); //-- yellow

    // Blobs and Clusters
    declareProperty("BlobRejectedColor",        m_Color_RejectedBlob    = 0xFF0000); //-- red
    declareProperty("ClusterUnusedColor",       m_Color_clusterUnused   = 0xFF00FF); //-- magenta
    declareProperty("ClusterUsedColor",         m_Color_clusterUsed     = 0x008080); //-- cyan
    declareProperty("Gamma1ProngColor",         m_Color_Gamma1Prong     = 0x00FF00); //-- light green (lime)
    declareProperty("Gamma2ProngColor",         m_Color_Gamma2Prong     = 0x00BFFF); //-- light blue 
    declareProperty("GammaOtherProngColor",     m_Color_GammaOtherProng = 0xFFA500); //-- orange

    // Vertex
    declareProperty("PrimaryVertexProngColor",  m_Color_primaryVertex   = 0x000000); //-- black
    declareProperty("SecondaryVertexProngColor",m_Color_secondaryVertex = 0x000000); //-- black
    declareProperty("VertexSphereColor",        m_Color_VertexBlob    = 0xFA8072); //-- salmon 

    // VertexBlob
    declareProperty( "SkipLowEnergyClusterVtxEnergy", fSkipLowEnergyClusterVtxEnergy = true);
    declareProperty( "ThresholdVertexEnergy",         fThresholdVertexEnergy = 1.0*CLHEP::MeV);
    declareProperty( "VertexBlobRadius", m_vertex_blob_radius = 90.0 * CLHEP::millimeter ); 

    // ConeBlobs
    declareProperty( "RejectedClustersTime",  m_rejectedClustersTime  = 25 * CLHEP::ns );

    declareProperty( "ApplyAttenuationCorrection", m_ApplyAttenuationCorrection = true);
    declareProperty( "UVMatchTolerance", m_UVMatchTolerance = 10.0*CLHEP::mm);
    declareProperty( "UVMatchMoreTolerance", m_UVMatchMoreTolerance = 100.0*CLHEP::mm);
    declareProperty( "AllowUVMatchWithMoreTolerance", m_AllowUVMatchWithMoreTolerance = true);

    declareProperty( "TrytoRecover_1Shower",           m_TrytoRecover_1Shower = true);
    declareProperty( "TrytoRecover_3Shower",           m_TrytoRecover_3Shower = true);
    declareProperty( "TrytoRecover_4Shower",           m_TrytoRecover_4Shower = true);
    declareProperty( "RecoverShower_InvMass",          m_recoverShower_invMass = false);
    declareProperty( "RecoverShower_Direction",        m_recoverShower_Direction = true);
    declareProperty( "RecoverSingleShower_SmallAngle", m_recoverSingleShower_SmallAngle = true);
    declareProperty( "RecoverSingleShower_SearchView", m_recoverSingleShower_SearchView = true);

    // dedx uncertainties
    m_dedx_uncertainties.push_back("Nominal");
    m_dedx_uncertainties.push_back("Mass_Up");
    m_dedx_uncertainties.push_back("Mass_Down");
    m_dedx_uncertainties.push_back("MEU_Up");
    m_dedx_uncertainties.push_back("MEU_Down");
    m_dedx_uncertainties.push_back("BetheBloch_Up");
    m_dedx_uncertainties.push_back("BetheBloch_Down");
    m_dedx_uncertainties.push_back("Birks");

    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    info() << " CCProtonPi0 Hypothesis added " << endmsg;

}

//==============================================================================
// Initialize
//==============================================================================
StatusCode CCProtonPi0::initialize()
{
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) return Error( "Failed to initialize!", sc ); 

    info()<<"\tInitialized MinervaAnalysisTool"<<endmsg;

    //Seed the RNG
    m_randomGen = new TRandom3( m_randomSeed );   
    if ( m_randomGen == NULL ) return StatusCode::FAILURE;

    // Reset Counters for Functions - Debugging Purposes
    // They can be used to track events in DST and Ana
    N_tagTruth = 0;
    N_reconstructEvent = 0;

    // Fiducial Volume
    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5991.37*CLHEP::mm;   // ~module 27, plane 1
    m_fidDownStreamZ = 8363.92*CLHEP::mm;   // ~module 79, plane 2

    // Reconstructable Volume
    m_recoHexApothem  = 1000.0*CLHEP::mm; 
    m_recoUpStreamZ   = 5810.0*CLHEP::mm;
    m_recoDownStreamZ = 8600.0*CLHEP::mm;

    //--------------------------------------------------------------------------
    // Initialize Analysis Tools and Services
    //--------------------------------------------------------------------------
    try {
        m_InnerDetector = getDet<Minerva::DeDetector>("/dd/Structure/Minerva/Detector/InnerDetector");
        m_OuterDetector = getDet<Minerva::DeOuterDetector>("/dd/Structure/Minerva/Detector/OuterDetector");
        m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool");
        m_protonUtils = tool<IProtonUtils>("ProtonUtils"); 
        m_particleMaker = tool<IParticleMakerTool>("ParticleMakerTool"); 
        m_LikelihoodPIDTool = tool<ILikelihoodParticleTool>("LikelihoodPIDTool", m_LPIDToolAlias);
        m_vertexFitter = tool<IVertexFitter>("VertexFitterKalman"); 
        m_vertexEnergyStudyTool = tool<IVertexEnergyStudyTool>("VertexEnergyStudyTool"); 
        m_anchoredTracker = tool<IAnchoredTrackFormation>("AnchoredShortTracker"); 
        m_coneUtilsTool = tool<IConeUtilsTool>("ConeUtilsTool"); 
        m_odMatchTool = tool<IODProngClassificationTool>("ODTrackMatchTool"); 
        m_energyCorrectionTool = tool<IEnergyCorrectionTool>("EnergyCorrectionTool"); 
        m_hitTagger = tool<IHitTaggerTool>( "HitTaggerTool" ); 
        m_michelTool_Improved = tool<IImprovedMichelTool>("ImprovedMichelTool");
        m_michelTool_Default = tool<IMichelTool>("MichelTool","MichelTool_Default");
        m_michelTool_Large = tool<IMichelTool>("MichelTool","MichelTool_Large");
        m_objectAssociator = tool<IMinervaObjectAssociator>("MinervaObjectAssociator");
        m_caloUtils = tool<ICalorimetryUtils>("CalorimetryUtils");  
        m_stopPointBlobTool = tool<IIDAnchoredBlobCreator>("VertexBlobCreator");
        m_extraEnergyTool = tool<IExtraEnergyTool>("ExtraEnergyTool");
        m_getDeadTimeTool = tool<IGetDeadTime>( "GetDeadTime" );
        m_trackPropagator = tool<ITrackLinearPropagator>("TrackLinearPropagator");
        m_recoTimeTool = tool<IRecoObjectTimeTool>( "RecoObjectTimeTool" );
        m_mathTool = tool<IMinervaMathTool>("MinervaMathTool");
        m_blobUtils = tool<IBlobCreatorUtils>("BlobCreatorUtils");
        m_idHoughBlob = tool<IHoughBlob>("HTBlob");
        m_idHoughTool = tool<IHoughTool>("HTtool");
        m_AttenuationCorrectionTool = tool<IGetCalAttenuation>("GetCalAttenuation");
        m_gigaCnvSvc = svc<IGiGaGeomCnvSvc>("GiGaGeo", true);
        m_blobSeedingTool = tool<IIDBlobSeedingTool>("BlobSeedingTool");
	m_timeTool = tool<IRecoObjectTimeTool>( "RecoObjectTimeTool", m_timeToolAlias );
    } catch(GaudiException& e) {
        error()<<"Could not obtained all Tools and Services!"<<endmsg;
        return StatusCode::FAILURE;
    }

    service("GeomUtilSvc", m_GeomUtilSvc, true);
    m_idDet = m_GeomUtilSvc->getIDDet();
    m_odDet = m_GeomUtilSvc->getODDet();

    // declare common branches
    declareMinosMuonBranches();
    declareGenieWeightBranches();

    //--------------------------------------------------------------------------
    // Select the branches you want in your AnaTuple
    //--------------------------------------------------------------------------
    const double SENTINEL = -9.9; 

    // Truth Branches
    declareDoubleTruthBranch("eventID",0.0);

    declareIntTruthBranch("N_FSParticles", -1 );
    declareIntTruthBranch("N_proton", -1 );
    declareIntTruthBranch("N_pi0", -1 );
    declareIntTruthBranch("N_other", -1 );

    declareIntTruthBranch("InNucleus_N_pi0_initial", -1 );
    declareIntTruthBranch("InNucleus_N_piplus_initial", -1 );
    declareIntTruthBranch("InNucleus_N_piminus_initial", -1 );
    declareIntTruthBranch("InNucleus_N_pi0_final", -1 );
    declareIntTruthBranch("InNucleus_N_piplus_final", -1 );
    declareIntTruthBranch("InNucleus_N_piminus_final", -1 );
   
    declareIntTruthBranch( "vertex_module", 500);
    declareIntTruthBranch( "vertex_plane", 0);
    declareIntTruthBranch( "target_material", -1); 

    // Event Kinematics
    declareDoubleTruthBranch("QSq_exp", SENTINEL);
    declareDoubleTruthBranch("WSq_exp", SENTINEL);
    declareDoubleTruthBranch("W_exp", SENTINEL);

    // Signal Kinematics
    declareIntTruthBranch("gamma1_TrackID", SENTINEL);
    declareIntTruthBranch("gamma2_TrackID", SENTINEL);
    declareContainerDoubleTruthBranch("muon_4P",          4, SENTINEL );
    declareContainerDoubleTruthBranch("proton_4P",        4, SENTINEL );
    declareContainerDoubleTruthBranch("pi0_4P",           4, SENTINEL );
    declareContainerDoubleTruthBranch("gamma1_4P",        4, SENTINEL );
    declareContainerDoubleTruthBranch("gamma2_4P",        4, SENTINEL );
    declareContainerDoubleTruthBranch("gamma1_init_pos",  3, SENTINEL );
    declareContainerDoubleTruthBranch("gamma2_init_pos",  3, SENTINEL );
    declareContainerDoubleTruthBranch("gamma1_final_pos", 3, SENTINEL );
    declareContainerDoubleTruthBranch("gamma2_final_pos", 3, SENTINEL );
    declareContainerDoubleTruthBranch("gamma1_final_pos_estimated", 3, SENTINEL );
    declareContainerDoubleTruthBranch("gamma2_final_pos_estimated", 3, SENTINEL );
    declareBoolTruthBranch("isGamma1_conv_inside");
    declareBoolTruthBranch("isGamma2_conv_inside");
    declareDoubleTruthBranch("muon_P",SENTINEL);
    declareDoubleTruthBranch("pi0_P",SENTINEL);
    declareDoubleTruthBranch("pi0_KE",SENTINEL);
    declareDoubleTruthBranch("proton_P",SENTINEL);
    declareDoubleTruthBranch("muon_theta",SENTINEL);
    declareDoubleTruthBranch("muon_theta_beam",SENTINEL);
    declareDoubleTruthBranch("muon_thetaX_beam",SENTINEL);
    declareDoubleTruthBranch("muon_thetaY_beam",SENTINEL);
    declareDoubleTruthBranch("pi0_theta",SENTINEL);
    declareDoubleTruthBranch("pi0_theta_beam",SENTINEL);
    declareDoubleTruthBranch("proton_theta",SENTINEL);
    declareDoubleTruthBranch("proton_theta_beam",SENTINEL);

    declareIntTruthBranch("pi0_status", -9 );
    declareIntTruthBranch("pi0_Mother", -9 );
    declareIntTruthBranch("pi0_MotherStatus", -9 );   
    declareIntTruthBranch("pi0_GrandMother", -9 );
    declareIntTruthBranch("pi0_GrandMotherStatus", -9 );

    // Truth Michel Information
    declareIntTruthBranch("N_trueMichelElectrons", -1 );
    declareDoubleTruthBranch("michelElectron_P", SENTINEL );
    declareDoubleTruthBranch("michelElectron_E", SENTINEL );
    declareContainerDoubleTruthBranch("michelMuon_endPoint", 3, SENTINEL ); 
    declareDoubleTruthBranch("michelMuon_P", SENTINEL );
    declareDoubleTruthBranch("michelMuon_length", SENTINEL );
    declareDoubleTruthBranch("michelMuon_end_dist_vtx", SENTINEL );
    declareDoubleTruthBranch("michelPion_P", SENTINEL );
    declareDoubleTruthBranch("michelPion_length", SENTINEL );
    declareDoubleTruthBranch("michelPion_begin_dist_vtx", SENTINEL );

    // ------------------------------------------------------------------------
    // Signal and Background
    // ------------------------------------------------------------------------
    declareBoolTruthBranch("isSignal");
    declareBoolTruthBranch("isSignal_BeforeFSI");
    declareBoolTruthBranch("isSignalOut_Acceptance");
    declareBoolTruthBranch("isSignalOut_Kinematics");
    declareBoolTruthBranch("isSignal_EventRecord");
    declareBoolTruthBranch("isFidVol");
    declareBoolTruthBranch("isNC");
    declareBoolTruthBranch("ReconstructEvent");

    // Cross Section Variables Before FSI
    declareDoubleTruthBranch("muon_P_BeforeFSI", -1); 
    declareDoubleTruthBranch("muon_theta_beam_BeforeFSI", -1); 
    declareDoubleTruthBranch("pi0_P_BeforeFSI", -1); 
    declareDoubleTruthBranch("pi0_KE_BeforeFSI", -1); 
    declareDoubleTruthBranch("pi0_theta_beam_BeforeFSI", -1); 
    declareDoubleTruthBranch("Enu_BeforeFSI", -1); 
    declareDoubleTruthBranch("QSq_exp_BeforeFSI", -1); 
    declareDoubleTruthBranch("W_exp_BeforeFSI", -1); 

    // BackgroundWithPi0
    declareIntTruthBranch("Bckg_nPi0_Primary", -1);
    declareIntTruthBranch("Bckg_nPi0_Secondary", -1);
    declareIntTruthBranch("Bckg_nPi0_Total", -1);
    declareBoolTruthBranch("isBckg_NoPi0");
    declareBoolTruthBranch("isBckg_SinglePi0");
    declareBoolTruthBranch("isBckg_MultiPi0");

    // Background Types
    declareBoolTruthBranch("isBckg_Compact_WithPi0");
    declareBoolTruthBranch("isBckg_Compact_QELike");
    declareBoolTruthBranch("isBckg_Compact_SinglePiPlus");
    declareBoolTruthBranch("isBckg_Compact_Other");

    // Background Types All
    declareIntTruthBranch("Bckg_nOther", -1);
    declareIntTruthBranch("Bckg_nPiCharged", -1);
    declareIntTruthBranch("Bckg_nPiCharged_ChargeExchanged", -1);
    declareBoolTruthBranch("isBckg_NC");
    declareBoolTruthBranch("isBckg_AntiNeutrino");
    declareBoolTruthBranch("isBckg_QELike");
    declareBoolTruthBranch("isBckg_SingleChargedPion");
    declareBoolTruthBranch("isBckg_SingleChargedPion_ChargeExchanged");
    declareBoolTruthBranch("isBckg_DoublePionWithPi0");
    declareBoolTruthBranch("isBckg_DoublePionWithoutPi0");
    declareBoolTruthBranch("isBckg_MultiPionWithPi0");
    declareBoolTruthBranch("isBckg_MultiPionWithoutPi0");
    declareBoolTruthBranch("isBckg_Other");

    
    // Background Branching - For Each Background Type
    declareBoolTruthBranch("isBckg_withMichel");

     // ------------------------------------------------------------------------
    // Unused Clusters Visible Energy
    // ------------------------------------------------------------------------
    // ECAL -- Filled in PreFilterPi0()
    declareDoubleTruthBranch("ecal_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_total_truth", -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_pizero",      -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_piplus",      -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_piminus",     -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_muon",        -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_proton",      -1.0);
    declareDoubleTruthBranch("ecal_unused_evis_neutron",     -1.0);
    // HCAL -- Filled in PreFilterPi0()
    declareDoubleTruthBranch("hcal_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_total_truth", -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_pizero",      -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_piplus",      -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_piminus",     -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_muon",        -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_proton",      -1.0);
    declareDoubleTruthBranch("hcal_unused_evis_neutron",     -1.0);
    // Tracker + ECAL + HCAL -- Filled in PreFilterPi0()
    declareDoubleTruthBranch("other_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("other_unused_evis_total_truth", -1.0);
    declareDoubleTruthBranch("other_unused_evis_pizero",      -1.0);
    declareDoubleTruthBranch("other_unused_evis_piplus",      -1.0);
    declareDoubleTruthBranch("other_unused_evis_piminus",     -1.0);
    declareDoubleTruthBranch("other_unused_evis_muon",        -1.0);
    declareDoubleTruthBranch("other_unused_evis_proton",      -1.0);
    declareDoubleTruthBranch("other_unused_evis_neutron",     -1.0);

    // Near Vertex -- Filled in VertexBlob()
    declareIntTruthBranch("vertex_unused_evis_most_pdg",       -1);
    declareDoubleTruthBranch("vertex_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_total_truth", -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_pizero",      -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_piplus",      -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_piminus",     -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_muon",        -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_proton",      -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_neutron",     -1.0);
    declareDoubleTruthBranch("vertex_unused_evis_gamma",       -1.0);

    // Rejected -- Filled in ConeBlobs()
    declareIntTruthBranch("Rejected_unused_evis_most_pdg",       -1);
    declareDoubleTruthBranch("Rejected_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("Rejected_unused_evis_total_truth", -1.0);

    // Dispersed Clusters -- Filled in DispersedBlob(), After Pi0 Reconstruction
    declareIntTruthBranch("dispersed_unused_evis_most_pdg",       -1);
    declareDoubleTruthBranch("dispersed_unused_evis_total_norm",  -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_total_truth", -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_pizero",      -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_piplus",      -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_piminus",     -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_muon",        -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_proton",      -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_neutron",     -1.0);
    declareDoubleTruthBranch("dispersed_unused_evis_gamma",       -1.0);

    // Truth Match for Found Pi0 Blobs
    declareIntTruthBranch("blob1_evis_most_pdg", -1);
    declareIntTruthBranch("blob1_evis_most_TrackID",   -1);
    declareIntTruthBranch("blob1_evis_parent_pdg",     -1);
    declareIntTruthBranch("blob1_evis_parent_TrackID", -1);
    declareDoubleTruthBranch("blob1_evis_total_norm", -1); 
    declareDoubleTruthBranch("blob1_evis_total_truth", -1);
    declareDoubleTruthBranch("blob1_evis_pizero", -1);
    declareDoubleTruthBranch("blob1_evis_piplus", -1);
    declareDoubleTruthBranch("blob1_evis_piminus", -1);
    declareDoubleTruthBranch("blob1_evis_muon", -1);
    declareDoubleTruthBranch("blob1_evis_proton", -1);
    declareDoubleTruthBranch("blob1_evis_neutron", -1);

    declareIntTruthBranch("blob2_evis_most_pdg", -1);
    declareIntTruthBranch("blob2_evis_most_TrackID",   -1);
    declareIntTruthBranch("blob2_evis_parent_pdg",     -1);
    declareIntTruthBranch("blob2_evis_parent_TrackID", -1);
    declareDoubleTruthBranch("blob2_evis_total_norm", -1); 
    declareDoubleTruthBranch("blob2_evis_total_truth", -1);
    declareDoubleTruthBranch("blob2_evis_pizero", -1);
    declareDoubleTruthBranch("blob2_evis_piplus", -1);
    declareDoubleTruthBranch("blob2_evis_piminus", -1);
    declareDoubleTruthBranch("blob2_evis_muon", -1);
    declareDoubleTruthBranch("blob2_evis_proton", -1);
    declareDoubleTruthBranch("blob2_evis_neutron", -1);

    declareDoubleTruthBranch("total_captured_evis_total_norm", -1);
    declareDoubleTruthBranch("total_captured_evis_total_truth", -1);
    declareDoubleTruthBranch("total_captured_evis_pizero", -1);
    declareDoubleTruthBranch("allClusters_evis_pizero", -1);

    // ------------------------------------------------------------------------       
    // Event Branches
    // ------------------------------------------------------------------------       
    // MC Trajectory Information
    declareIntEventBranch("detmc_ntrajectory", -1);
    declareIntEventBranch("detmc_ntrajectory2", -1);

    declareContainerIntEventBranch("detmc_traj_pdg");
    declareContainerIntEventBranch("detmc_traj_id");
    declareContainerIntEventBranch("detmc_traj_mother");
    declareContainerIntEventBranch("detmc_traj_proc");
    declareContainerIntEventBranch("detmc_traj_status");

    declareContainerDoubleEventBranch("detmc_traj_x0");
    declareContainerDoubleEventBranch("detmc_traj_y0");
    declareContainerDoubleEventBranch("detmc_traj_z0");
    declareContainerDoubleEventBranch("detmc_traj_t0");

    declareContainerDoubleEventBranch("detmc_traj_px0");
    declareContainerDoubleEventBranch("detmc_traj_py0");
    declareContainerDoubleEventBranch("detmc_traj_pz0");
    declareContainerDoubleEventBranch("detmc_traj_E0");

    declareContainerDoubleEventBranch("detmc_traj_xf");
    declareContainerDoubleEventBranch("detmc_traj_yf");
    declareContainerDoubleEventBranch("detmc_traj_zf");
    declareContainerDoubleEventBranch("detmc_traj_tf");

    declareContainerDoubleEventBranch("detmc_traj_pxf");
    declareContainerDoubleEventBranch("detmc_traj_pyf");
    declareContainerDoubleEventBranch("detmc_traj_pzf");
    declareContainerDoubleEventBranch("detmc_traj_Ef");

    declareContainerDoubleEventBranch("detmc_traj_prepxf");
    declareContainerDoubleEventBranch("detmc_traj_prepyf");
    declareContainerDoubleEventBranch("detmc_traj_prepzf");
    declareContainerDoubleEventBranch("detmc_traj_preEf");

    // Cut Results
    declareIntEventBranch( "Cut_Vertex_None", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Reconstructable", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Fiducial", -1 );
    declareIntEventBranch( "Cut_Muon_None",-1);
    declareIntEventBranch( "Cut_Muon_Charge",-1);
    declareIntEventBranch( "Cut_Vertex_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_Vertex_Large_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_EndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_secEndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_Particle_None", -1 );
    declareIntEventBranch( "Cut_Proton_None", -1 );
    declareIntEventBranch( "Cut_Proton_Bad", -1 );
    declareIntEventBranch( "Cut_PreFilter_Pi0", -1 );
    declareIntEventBranch( "Cut_ConeBlobs", -1 );
    declareIntEventBranch( "Cut_BlobDirectionBad", -1 );
    declareIntEventBranch( "Cut_Pi0_Bad", -1 );

    // Michel Veto
    declareBoolEventBranch( "ImprovedMichel_EventHasMichel");
    declareBoolEventBranch( "ImprovedMichel_EventHasMatchedMichel");
    declareBoolEventBranch( "ImprovedMichel_VertexHasMichel");
    declareBoolEventBranch( "ImprovedMichel_EndPointHasMichel");
    declareBoolEventBranch( "ImprovedMichel_secEndPointHasMichel");
    declareBoolEventBranch( "ImprovedMichel_Gamma1HasMichel");
    declareBoolEventBranch( "ImprovedMichel_Gamma2HasMichel");
    
    //Charge Pion cuts
    // Filled in CreateTrackedParticles()
    /*declareContainerDoubleEventBranch("all_prongs_LLRScore",    10, -999.9);
    declareIntEventBranch("N_Primary_Hadrons",              -1);
    declareDoubleEventBranch("hadron_visible_energy",       -9.9);

    // GetPionProngs
    declareDoubleEventBranch("prongs_pionScore",   -9.9);
    declareDoubleEventBranch("prongs_protonScore", -9.9);
     */
    // General Reco
    declareDoubleEventBranch( "time", -1.0 );
    declareDoubleEventBranch( "reco_eventID", -1.0 );

    // Vertex Information 
    declareIntEventBranch("vtx_total_count",-1);
    declareIntEventBranch("vtx_secondary_count",-1);
    declareIntEventBranch("vtx_primary_index",-1);
    declareIntEventBranch("vtx_primary_multiplicity",-1);
    declareIntEventBranch("vtx_fit_converged",-1);
    declareDoubleEventBranch("vtx_fit_chi2", -1);
    declareContainerDoubleEventBranch("fit_vtx", 3, -1);

    // Michel Prong
    declareDoubleEventBranch( "michelProng_distance",  -1.0);
    declareDoubleEventBranch( "michelProng_energy",    -1.0);
    declareDoubleEventBranch( "michelProng_time_diff", -1.0);
    declareDoubleEventBranch( "michelProng_end_Z",     -1.0);
    declareDoubleEventBranch( "michelProng_begin_Z",   -1.0);

    // Number of Proton Candidates

    // Discard Far Tracks  - Update nTracks
    declareIntEventBranch("nTracks", -1);
    declareIntEventBranch("nTracks_Close",-1);
    declareIntEventBranch("nTracks_Far",-1);
    declareIntEventBranch("nTracks_Discarded",-1);
    declareContainerIntEventBranch("nTracks_Secondary_Vtx");

    // PreFilterPi0()
    declareIntEventBranch("preFilter_Result", -1);
    declareDoubleEventBranch("preFilter_rejectedEnergy", -1.0);
    declareDoubleEventBranch("preFilter_evis_nearvtx", -1.0);
    declareDoubleEventBranch("preFilter_evis_total", -1.0);
    declareDoubleEventBranch("preFilter_evis_NuclearTarget",  -1.0);
    declareDoubleEventBranch("preFilter_evis_Tracker",  -1.0);
    declareDoubleEventBranch("preFilter_evis_ECAL",  -1.0);
    declareDoubleEventBranch("preFilter_evis_HCAL",  -1.0);
    declareDoubleEventBranch("preFilter_evis_TotalExceptNuclearTarget", -1.0);

    // VertexBlob()
    declareDoubleEventBranch("vertex_blob_evis", SENTINEL );
    declareDoubleEventBranch("vertex_blob_energy", SENTINEL );

    // ConeBlobs()
    declareContainerIntEventBranch("shower_ndigits",         "anglescan_nfoundBlobs");
    declareContainerIntEventBranch("shower_nclusters",       "anglescan_nfoundBlobs");
    declareContainerIntEventBranch("shower_isGood",          "anglescan_nfoundBlobs");
    declareContainerIntEventBranch("shower_pdg",             "anglescan_nfoundBlobs");
    declareContainerDoubleEventBranch("shower_energy",       "anglescan_nfoundBlobs");
    declareContainerDoubleEventBranch("shower_dist_fitline", "anglescan_nfoundBlobs");
    declareContainerDoubleEventBranch("shower_score",        "anglescan_nfoundBlobs");
    declareIntEventBranch("Good_Blobs", -1);

    declareContainerIntEventBranch("AllShower_nclusters");
    declareContainerDoubleEventBranch("AllShower_energy");
    declareContainerDoubleEventBranch("AllShower_theta");
    declareContainerDoubleEventBranch("AllShower_dist_vtx");
    declareContainerDoubleEventBranch("AllShower_dEdX");

    
    declareDoubleEventBranch("ConeBlobs_usable_evis_Tracker", SENTINEL );
    declareDoubleEventBranch("Coneblobs_usable_evis_ECAL", SENTINEL );
    declareDoubleEventBranch("Coneblobs_usable_evis_HCAL", SENTINEL );
    declareIntEventBranch("anglescan_ncandx", -1);
    declareIntEventBranch("anglescan_ncand", -1);
    declareIntEventBranch("anglescan_nfoundBlobs", -1);
    declareBoolEventBranch("is_blobs_recovered");
    declareBoolEventBranch("is_blobs_recovered_direction");
    declareBoolEventBranch("is_blobs_recovered_invMass");
    declareBoolEventBranch("is_blobs_recovered_small_angle");
    declareBoolEventBranch("is_blobs_recovered_search_view_U");
    declareBoolEventBranch("is_blobs_recovered_search_view_V");

    // ConeBlobs() -- Testing Showers GoodPosition and GoodDirection
    declareIntEventBranch("shower_has_goodPosition", -1);
    declareIntEventBranch("shower_has_goodDirection", -1);
    
    // ConeBlobs() -- Recovered Showers
    
    declareContainerIntTruthBranch("AllShowers_evis_most_pdg");
    declareContainerDoubleTruthBranch("AllShowers_evis_total_norm");
    declareContainerDoubleTruthBranch("AllShowers_evis_total_truth");
    declareContainerDoubleTruthBranch("AllShowers_evis_pizero");
    declareContainerDoubleTruthBranch("AllShowers_evis_piplus");
    declareContainerDoubleTruthBranch("AllShowers_evis_piminus");
    declareContainerDoubleTruthBranch("AllShowers_evis_minus");
    declareContainerDoubleTruthBranch("AllShowers_evis_proton");
    declareContainerDoubleTruthBranch("AllShowers_evis_neutron");

    declareIntEventBranch("OneShower_nClusters", -1);
    declareDoubleEventBranch("OneShower_energy",   SENTINEL);
    declareDoubleEventBranch("OneShower_theta",    SENTINEL);
    declareDoubleEventBranch("OneShower_dist_vtx", SENTINEL);
    declareIntTruthBranch("OneShower_evis_most_pdg", -1);
    declareDoubleTruthBranch("OneShower_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("OneShower_evis_neutron",     SENTINEL);

    declareIntEventBranch("ThreeShower_s1_nClusters", -1);
    declareDoubleEventBranch("ThreeShower_s1_energy",   SENTINEL);
    declareDoubleEventBranch("ThreeShower_s1_theta",    SENTINEL);
    declareDoubleEventBranch("ThreeShower_s1_dist_vtx", SENTINEL);
    declareIntTruthBranch("ThreeShower_s1_evis_most_pdg", -1);
    declareDoubleTruthBranch("ThreeShower_s1_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s1_evis_neutron",     SENTINEL);

    declareIntEventBranch("ThreeShower_s2_nClusters", -1);
    declareDoubleEventBranch("ThreeShower_s2_energy",   SENTINEL);
    declareDoubleEventBranch("ThreeShower_s2_theta",    SENTINEL);
    declareDoubleEventBranch("ThreeShower_s2_dist_vtx", SENTINEL);
    declareIntTruthBranch("ThreeShower_s2_evis_most_pdg", -1);
    declareDoubleTruthBranch("ThreeShower_s2_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s2_evis_neutron",     SENTINEL);

    declareIntEventBranch("ThreeShower_s3_nClusters", -1);
    declareDoubleEventBranch("ThreeShower_s3_energy",   SENTINEL);
    declareDoubleEventBranch("ThreeShower_s3_theta",    SENTINEL);
    declareDoubleEventBranch("ThreeShower_s3_dist_vtx", SENTINEL);
    declareIntTruthBranch("ThreeShower_s3_evis_most_pdg", -1);
    declareDoubleTruthBranch("ThreeShower_s3_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("ThreeShower_s3_evis_neutron",     SENTINEL);

    declareIntEventBranch("FourShower_s1_nClusters", -1);
    declareDoubleEventBranch("FourShower_s1_energy",   SENTINEL);
    declareDoubleEventBranch("FourShower_s1_theta",    SENTINEL);
    declareDoubleEventBranch("FourShower_s1_dist_vtx", SENTINEL);
    declareIntTruthBranch("FourShower_s1_evis_most_pdg", -1);
    declareDoubleTruthBranch("FourShower_s1_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s1_evis_neutron",     SENTINEL);

    declareIntEventBranch("FourShower_s2_nClusters", -1);
    declareDoubleEventBranch("FourShower_s2_energy",   SENTINEL);
    declareDoubleEventBranch("FourShower_s2_theta",    SENTINEL);
    declareDoubleEventBranch("FourShower_s2_dist_vtx", SENTINEL);
    declareIntTruthBranch("FourShower_s2_evis_most_pdg", -1);
    declareDoubleTruthBranch("FourShower_s2_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s2_evis_neutron",     SENTINEL);

    declareIntEventBranch("FourShower_s3_nClusters", -1);
    declareDoubleEventBranch("FourShower_s3_energy",   SENTINEL);
    declareDoubleEventBranch("FourShower_s3_theta",    SENTINEL);
    declareDoubleEventBranch("FourShower_s3_dist_vtx", SENTINEL);
    declareIntTruthBranch("FourShower_s3_evis_most_pdg", -1);
    declareDoubleTruthBranch("FourShower_s3_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s3_evis_neutron",     SENTINEL);

    declareIntEventBranch("FourShower_s4_nClusters", -1);
    declareDoubleEventBranch("FourShower_s4_energy",   SENTINEL);
    declareDoubleEventBranch("FourShower_s4_theta",    SENTINEL);
    declareDoubleEventBranch("FourShower_s4_dist_vtx", SENTINEL);
    declareIntTruthBranch("FourShower_s4_evis_most_pdg", -1);
    declareDoubleTruthBranch("FourShower_s4_evis_total_norm",  SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_total_truth", SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_pizero",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_piplus",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_piminus",     SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_muon",        SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_proton",      SENTINEL);
    declareDoubleTruthBranch("FourShower_s4_evis_neutron",     SENTINEL);

    // processBlobs() -- Called from ConeBlobs()
    declareIntEventBranch("g1blob_1ParFit_ndof", -1);
    declareIntEventBranch("g2blob_1ParFit_ndof", -1);
    declareDoubleEventBranch("g1blob_1ParFit_fval", -1.0);
    declareDoubleEventBranch("g2blob_1ParFit_fval", -1.0);    
    declareDoubleEventBranch("g1blob_2ParFit_vtx_distance",-1.0);
    declareDoubleEventBranch("g2blob_2ParFit_vtx_distance",-1.0);

    // ODActivity() -- Called from ConeBlobs()
    declareIntEventBranch( "od_energeticTower", SENTINEL );
    declareDoubleEventBranch( "od_upstreamFrame", SENTINEL );
    declareDoubleEventBranch( "od_downstreamFrame", SENTINEL );
    declareDoubleEventBranch( "od_upstreamFrame_z", SENTINEL );
    declareDoubleEventBranch( "od_downstreamFrame_z", SENTINEL );
    declareDoubleEventBranch( "od_highStory", SENTINEL );
    declareDoubleEventBranch( "od_lowStory", SENTINEL );
    declareDoubleEventBranch( "od_highStory_t", SENTINEL );
    declareDoubleEventBranch( "od_lowStory_t", SENTINEL );
    declareDoubleEventBranch( "od_maxEnergy", SENTINEL );
    declareContainerDoubleEventBranch( "od_distanceBlobTower");
    declareContainerDoubleEventBranch( "od_towerEnergy" );
    declareContainerDoubleEventBranch( "od_towerNClusters" );
    declareContainerDoubleEventBranch( "od_towerTime" );
    declareContainerDoubleEventBranch( "od_idBlobTime" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobOD" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobMuon" );

    // Blob Data -- Filled in setBlobData()
    declareIntEventBranch("gamma1_blob_ndigits",-1);
    declareIntEventBranch("gamma1_blob_nclusters",-1);
    declareDoubleEventBranch("gamma1_blob_theta", -1);
    declareDoubleEventBranch("gamma1_blob_energy",-1);
    declareDoubleEventBranch("gamma1_blob_minsep", -1.0);

    declareIntEventBranch("gamma2_blob_ndigits",-1);
    declareIntEventBranch("gamma2_blob_nclusters",-1);
    declareDoubleEventBranch("gamma2_blob_theta", -1);
    declareDoubleEventBranch("gamma2_blob_energy",-1);
    declareDoubleEventBranch("gamma2_blob_minsep", -1.0);

    // Calculate_dEdx() -- Called from setBlobData()
    declareIntEventBranch("g1dedx_nplane",      -1);
    declareIntEventBranch("g1dedx_doublet",     -1);
    declareIntEventBranch("g1dedx_empty_plane", -1);
    declareDoubleEventBranch("g1dedx_cos_theta", -2.0);
    declareDoubleEventBranch("g1dedx_total",   -1.0);
    declareDoubleEventBranch("g1dedx",         -1.0);
    declareDoubleEventBranch("g1dedx_total1",  -1.0);
    declareDoubleEventBranch("g1dedx1",        -1.0);
    declareContainerIntEventBranch("g1dedx_cluster_occupancy");
    declareContainerDoubleEventBranch("g1dedx_cluster_energy");
    declareContainerDoubleEventBranch("g1dedx_rev_cluster_energy");

    declareIntEventBranch("g2dedx_nplane",      -1);
    declareIntEventBranch("g2dedx_doublet",     -1);
    declareIntEventBranch("g2dedx_empty_plane", -1);
    declareDoubleEventBranch("g2dedx_cos_theta", -2.0);
    declareDoubleEventBranch("g2dedx_total",    -1.0);
    declareDoubleEventBranch("g2dedx",          -1.0);
    declareDoubleEventBranch("g2dedx_total1",   -1.0);
    declareDoubleEventBranch("g2dedx1",         -1.0);
    declareContainerIntEventBranch("g2dedx_cluster_occupancy");
    declareContainerDoubleEventBranch("g2dedx_cluster_energy");
    declareContainerDoubleEventBranch("g2dedx_rev_cluster_energy"); 

    // Extra Energy 
    declareDoubleEventBranch("Extra_Energy_Total", SENTINEL );
    declareDoubleEventBranch("Extra_Energy_Muon", SENTINEL );
    declareDoubleEventBranch("Extra_Energy_Rejected", SENTINEL );
    declareDoubleEventBranch("Extra_Energy_Leftover", SENTINEL);

    // Primary Vertex
    declareIntEventBranch("vtx_module", -99);
    declareIntEventBranch("vtx_plane",-1);
    declareDoubleEventBranch("vtx_x",0.0);
    declareDoubleEventBranch("vtx_y",0.0);
    declareDoubleEventBranch("vtx_z",0.0);

    // Muon Kinematics  -- Filled in setMuonData()
    declareIntEventBranch("muon_hasMinosMatchTrack", -1);
    declareIntEventBranch("muon_hasMinosMatchStub", -1);
    declareIntEventBranch("muon_minervaTrack_types", -1);
    declareIntEventBranch("muon_N_minosTracks", -1);
    declareIntEventBranch("muon_minosTrackQuality", -1);
    declareIntEventBranch("muon_roadUpstreamPlanes", -1);
    declareIntEventBranch("muon_charge", -99);
    declareDoubleEventBranch("muon_roadUpstreamEnergy", 0.0);
    declareDoubleEventBranch("muon_E", 0.0);
    declareDoubleEventBranch("muon_P", 0.0);
    declareDoubleEventBranch("muon_KE", 0.0);
    declareDoubleEventBranch("muon_px", 0.0);
    declareDoubleEventBranch("muon_py", 0.0);
    declareDoubleEventBranch("muon_pz", 0.0);  
    declareDoubleEventBranch("muon_phi", 0.0);
    declareDoubleEventBranch("muon_phi_beam", 0.0);
    declareDoubleEventBranch("muon_theta", 0.0);
    declareDoubleEventBranch("muon_theta_beam", 0.0);
    declareDoubleEventBranch("muon_thetaX_beam", 0.0);
    declareDoubleEventBranch("muon_thetaY_beam", 0.0);
    declareDoubleEventBranch("muon_theta_beam_biasUp", 0.0);
    declareDoubleEventBranch("muon_theta_beam_biasDown", 0.0); 
    declareDoubleEventBranch("muon_muScore", -1.0);
    declareDoubleEventBranch("muon_qp", 99.0);
    declareDoubleEventBranch("muon_qpqpe", 99.0);
    declareDoubleEventBranch("muon_E_shift", 0.0);
    declareDoubleEventBranch("muon_TrackVertexTime", -9999.9);
    declareContainerDoubleEventBranch( "muon_theta_allNodes" );
    declareContainerDoubleEventBranch( "muon_thetaX_allNodes" );
    declareContainerDoubleEventBranch( "muon_thetaY_allNodes" );

    
    // Proton Kinematics -- Filled in setProtonData()
    declareIntEventBranch("nProtonCandidates", -1);
    declareContainerIntEventBranch(   "all_protons_kinked", 10, -1);
    declareContainerIntEventBranch(   "all_protons_odMatch", 10, -1);
    declareContainerDoubleEventBranch("all_protons_length",10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_startPointX",10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_startPointY",10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_startPointZ", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_endPointX", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_endPointY", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_endPointZ", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_protonScore", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_pionScore", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_LLRScore", 10,  SENTINEL);
    declareContainerDoubleEventBranch("all_protons_chi2_ndf", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_theta", 10,  SENTINEL);
    declareContainerDoubleEventBranch("all_protons_theta_beam", 10,  SENTINEL);
    declareContainerDoubleEventBranch("all_protons_phi", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_phi_beam", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_KE", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_E", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_P", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_px", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_py", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_pz", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_p_calCorrection", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_protons_p_visEnergy", 10,SENTINEL);
    declareContainerDoubleEventBranch("all_protons_p_dEdXTool", 10,SENTINEL);

    for(unsigned int i = 0; i < m_dedx_uncertainties.size(); i++) {
        std::string name = m_dedx_uncertainties[i];
        declareContainerDoubleEventBranch("all_protons_energy_shift_" + name, 10, SENTINEL);
        declareContainerDoubleEventBranch("all_protons_score1_shift_"   + name, 10, SENTINEL);
    }

    // Leading(Interaction Proton) Kinematics 
    declareDoubleEventBranch("proton_px",SENTINEL);
    declareDoubleEventBranch("proton_py",SENTINEL);
    declareDoubleEventBranch("proton_pz",SENTINEL);
    declareDoubleEventBranch("proton_E",SENTINEL);
    declareDoubleEventBranch("proton_P",SENTINEL);
    declareDoubleEventBranch("proton_KE", SENTINEL);
    declareDoubleEventBranch("proton_theta", SENTINEL);
    declareDoubleEventBranch("proton_theta_beam", SENTINEL);
    declareDoubleEventBranch("proton_phi", SENTINEL);
    declareDoubleEventBranch("proton_phi_beam", SENTINEL);
    declareDoubleEventBranch("proton_length", SENTINEL);
    declareDoubleEventBranch("proton_protonScore", SENTINEL);
    declareDoubleEventBranch("proton_pionScore", SENTINEL);
    declareDoubleEventBranch("proton_LLRScore", SENTINEL);
    declareIntEventBranch("proton_kinked",-1);
    declareIntEventBranch("proton_leadingIndice",-1);

    for(unsigned int i = 0; i < m_dedx_uncertainties.size(); i++) {
        std::string name = m_dedx_uncertainties[i];
        declareDoubleEventBranch("proton_energy_shift_" + name, SENTINEL);
        declareDoubleEventBranch("proton_score1_shift_"   + name, SENTINEL);
    }

    // Pion Kinematics -- Filled in setPionData()
    /*declareIntEventBranch("nPionCandidates", -1);
    declareContainerIntEventBranch(   "all_pions_kinked",  10, -1);
    declareContainerIntEventBranch(   "all_pions_odMatch", 10, -1);
    declareContainerDoubleEventBranch("all_pions_length",      10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_startPointX", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_startPointY", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_startPointZ", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_endPointX",   10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_endPointY",   10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_endPointZ",   10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_protonScore", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_pionScore",   10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_LLRScore",    10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_chi2_ndf",    10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_theta",       10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_thetaX",      10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_thetaY",      10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_phi",         10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_KE",          10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_E",           10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_P",           10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_px",          10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_py",          10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_pz",          10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_p_calCorrection", 10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_p_visEnergy",     10, SENTINEL);
    declareContainerDoubleEventBranch("all_pions_p_dEdXTool",      10, SENTINEL);

    // Leading(Interaction Pion) Kinematics
    declareDoubleEventBranch("pion_px",          SENTINEL);
    declareDoubleEventBranch("pion_py",          SENTINEL);
    declareDoubleEventBranch("pion_pz",          SENTINEL);
    declareDoubleEventBranch("pion_E",           SENTINEL);
    declareDoubleEventBranch("pion_P",           SENTINEL);
    declareDoubleEventBranch("pion_theta",       SENTINEL);
    declareDoubleEventBranch("pion_KE",          SENTINEL);
    declareDoubleEventBranch("pion_phi",         SENTINEL);
    declareDoubleEventBranch("pion_thetaX",      SENTINEL);
    declareDoubleEventBranch("pion_thetaY",      SENTINEL);
    declareDoubleEventBranch("pion_length",      SENTINEL);
    declareDoubleEventBranch("pion_protonScore", SENTINEL);
    declareDoubleEventBranch("pion_pionScore",   SENTINEL);
    declareDoubleEventBranch("pion_LLRScore",    SENTINEL);
    declareIntEventBranch("pion_kinked",        -1);
    declareIntEventBranch("pion_leadingIndice", -1);
*/

    // Pi0 & Gamma1,2 Kinematics -- Filled in setPi0Data()
    declareDoubleEventBranch("pi0_px",               SENTINEL);
    declareDoubleEventBranch("pi0_py",               SENTINEL);
    declareDoubleEventBranch("pi0_pz",               SENTINEL);
    declareDoubleEventBranch("pi0_E",                SENTINEL);
    declareDoubleEventBranch("pi0_E_Cal",            SENTINEL);
    declareDoubleEventBranch("pi0_P",                SENTINEL);
    declareDoubleEventBranch("pi0_KE",               SENTINEL);
    declareDoubleEventBranch("pi0_invMass",          SENTINEL);
    declareDoubleEventBranch("pi0_invMass_Old",      SENTINEL);
    declareDoubleEventBranch("pi0_theta",            SENTINEL);
    declareDoubleEventBranch("pi0_phi",              SENTINEL);
    declareDoubleEventBranch("pi0_thetaX",           SENTINEL);
    declareDoubleEventBranch("pi0_thetaY",           SENTINEL);
    declareDoubleEventBranch("pi0_openingAngle",     SENTINEL);
    declareDoubleEventBranch("pi0_cos_openingAngle", SENTINEL);
     
    declareBoolEventBranch("gamma1_isMichel_begin");
    declareBoolEventBranch("gamma1_isMichel_end");
    declareDoubleEventBranch("gamma1_px",                  SENTINEL);
    declareDoubleEventBranch("gamma1_py",                  SENTINEL);
    declareDoubleEventBranch("gamma1_pz",                  SENTINEL);
    declareDoubleEventBranch("gamma1_E",                   SENTINEL);
    declareDoubleEventBranch("gamma1_E_Old",               SENTINEL);
    declareDoubleEventBranch("gamma1_P",                   SENTINEL);
    declareDoubleEventBranch("gamma1_theta",               SENTINEL);
    declareDoubleEventBranch("gamma1_theta_beam",SENTINEL);    
    declareDoubleEventBranch("gamma1_phi",                 SENTINEL);
    declareDoubleEventBranch("gamma1_phi_beam",SENTINEL);    
    declareDoubleEventBranch("gamma1_dEdx",                SENTINEL);
    declareDoubleEventBranch("gamma1_time",                SENTINEL);
    declareDoubleEventBranch("gamma1_dist_vtx",            SENTINEL);
    declareDoubleEventBranch("gamma1_maxZ",                SENTINEL);
    declareDoubleEventBranch("gamma1_evis_trkr",           SENTINEL);
    declareDoubleEventBranch("gamma1_evis_ecal",           SENTINEL);
    declareDoubleEventBranch("gamma1_evis_scal_X",         SENTINEL);
    declareDoubleEventBranch("gamma1_evis_scal_UV",        SENTINEL);
    declareDoubleEventBranch("gamma1_evis_hcal",           SENTINEL);
    declareDoubleEventBranch("gamma1_energy_trkr",         SENTINEL);
    declareDoubleEventBranch("gamma1_energy_ecal",         SENTINEL);
    declareDoubleEventBranch("gamma1_energy_scal_X",       SENTINEL);
    declareDoubleEventBranch("gamma1_energy_scal_UV",      SENTINEL);
    declareDoubleEventBranch("gamma1_energy_hcal",         SENTINEL);
    declareContainerDoubleEventBranch("gamma1_direction",3,SENTINEL);
    declareContainerDoubleEventBranch("gamma1_vertex",3,SENTINEL);
    declareContainerDoubleEventBranch("gamma1_end_vertex",3,SENTINEL);
    //just added
    declareDoubleEventBranch("gamma1_cone_angle",          SENTINEL);
    declareDoubleEventBranch("gamma1_pvalue",              SENTINEL);
    declareDoubleEventBranch("gamma1_MeanE_1stHalf",       SENTINEL);
    declareDoubleEventBranch("gamma1_MeanE_2ndHalf",       SENTINEL);
    declareDoubleEventBranch("gamma1_MedianE_1stHalf",     SENTINEL);
    declareDoubleEventBranch("gamma1_MedianE_2ndHalf",     SENTINEL);
    declareContainerDoubleEventBranch("gamma1_clusters_distance");
    declareContainerDoubleEventBranch("gamma1_clusters_energy");

    declareBoolEventBranch("gamma2_isMichel_begin");
    declareBoolEventBranch("gamma2_isMichel_end");
    declareDoubleEventBranch("gamma2_px",SENTINEL);
    declareDoubleEventBranch("gamma2_py",SENTINEL);
    declareDoubleEventBranch("gamma2_pz",SENTINEL);
    declareDoubleEventBranch("gamma2_E",SENTINEL);
    declareDoubleEventBranch("gamma2_E_Old",SENTINEL);
    declareDoubleEventBranch("gamma2_P",SENTINEL);
    declareDoubleEventBranch("gamma2_theta",SENTINEL);
    declareDoubleEventBranch("gamma2_theta_beam",SENTINEL);
    declareDoubleEventBranch("gamma2_phi",SENTINEL);
    declareDoubleEventBranch("gamma2_phi_beam",SENTINEL);
    declareDoubleEventBranch("gamma2_dEdx",SENTINEL);
    declareDoubleEventBranch("gamma2_time",SENTINEL);
    declareDoubleEventBranch("gamma2_dist_vtx",SENTINEL);
    declareDoubleEventBranch("gamma2_evis_trkr", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_ecal", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_scal_X", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_scal_UV", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_hcal", SENTINEL);
    declareDoubleEventBranch("gamma2_energy_trkr", SENTINEL);
    declareDoubleEventBranch("gamma2_energy_ecal", SENTINEL);
    declareDoubleEventBranch("gamma2_energy_scal_X", SENTINEL);
    declareDoubleEventBranch("gamma2_energy_scal_UV", SENTINEL);
    declareDoubleEventBranch("gamma2_energy_hcal", SENTINEL);
    declareDoubleEventBranch("gamma2_maxZ", SENTINEL);
    declareContainerDoubleEventBranch("gamma2_direction",3,SENTINEL);
    declareContainerDoubleEventBranch("gamma2_vertex",3,SENTINEL); 
    declareContainerDoubleEventBranch("gamma2_end_vertex",3,SENTINEL); 
     // just added 
    declareDoubleEventBranch("gamma2_cone_angle",          SENTINEL);
    declareDoubleEventBranch("gamma2_pvalue",              SENTINEL);
    declareDoubleEventBranch("gamma2_MeanE_1stHalf",       SENTINEL);
    declareDoubleEventBranch("gamma2_MeanE_2ndHalf",       SENTINEL);
    declareDoubleEventBranch("gamma2_MedianE_1stHalf",     SENTINEL);
    declareDoubleEventBranch("gamma2_MedianE_2ndHalf",     SENTINEL);
    declareContainerDoubleEventBranch("gamma2_clusters_distance");
    declareContainerDoubleEventBranch("gamma2_clusters_energy");

    //-------------------------------------------------------------------------
    // NeutrinoInt Branches 
    //-------------------------------------------------------------------------
    // Hadron Recoil
    declareDoubleBranch( m_hypMeths, "hadron_recoil_CCInc", 0.0);
    declareDoubleBranch( m_hypMeths, "hadron_recoil", 0.0);

    // Event Kinematics -- Filled in setEventKinematics()
   /* declareDoubleBranch( m_hypMeths, "vertex_energy",         SENTINEL);
    declareDoubleBranch( m_hypMeths, "Extra_Energy_Total",    SENTINEL);
    declareDoubleBranch( m_hypMeths, "neutrino_E",            SENTINEL);
    declareDoubleBranch( m_hypMeths, "QSq",                   SENTINEL);
    declareDoubleBranch( m_hypMeths, "WSq",                   SENTINEL);
    declareDoubleBranch( m_hypMeths, "neutrino_E_1Track_Alt", SENTINEL);
    declareDoubleBranch( m_hypMeths, "QSq_1Track_Alt",        SENTINEL);
    declareDoubleBranch( m_hypMeths, "WSq_1Track_Alt",        SENTINEL);
    */
    // Truth Match for Prongs
    declareIntBranch(m_hypMeths,    "isMuonInsideOD",        -1);
    declareIntBranch(m_hypMeths,    "ntrajMuonProng",        -1);
    declareIntBranch(m_hypMeths,    "trajMuonProngPrimary",  -1);
    declareIntBranch(m_hypMeths,    "trajMuonProngPDG",      -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPx", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPy", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPz", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngEnergy", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngMomentum", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPSelf", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonTheta",         -1);
    declareDoubleBranch(m_hypMeths, "trajMuonPhi",           -1);
    declareDoubleBranch(m_hypMeths, "endMuonTrajMomentum",   SENTINEL);
    declareDoubleBranch(m_hypMeths, "endMuonTrajXPosition",  SENTINEL);
    declareDoubleBranch(m_hypMeths, "endMuonTrajYPosition",  SENTINEL);
    declareDoubleBranch(m_hypMeths, "endMuonTrajZPosition",  SENTINEL);

    declareContainerIntBranch(m_hypMeths,    "isProtonInsideOD",        10, -1);
    declareContainerIntBranch(m_hypMeths,    "ntrajProtonProng",         10, -1);
    declareContainerIntBranch(m_hypMeths,    "trajProtonProngPrimary",  10, -1);
    declareContainerIntBranch(m_hypMeths,    "trajProtonProngPDG",      10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPx", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPy", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPz", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngEnergy", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngMomentum", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPSelf", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonTheta",         10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonPhi",           10, -1);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajMomentum",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajXPosition",  10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajYPosition",  10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajZPosition",  10, SENTINEL);

    // ------------------------------------------------------------------------
    // OPTIONAL STUDY BRANCHES
    // ------------------------------------------------------------------------
    if (m_study_michel_electron){ 
        // Reco Values 
        declareDoubleEventBranch( "vtx_michelProng_Large_distance", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_Large_energy", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_Large_time_diff", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_Large_end_Z", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_Large_begin_Z", -1.0);

        declareDoubleEventBranch( "vtx_michelProng_distance", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_energy", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_time_diff", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_end_Z", -1.0);
        declareDoubleEventBranch( "vtx_michelProng_begin_Z", -1.0);

        declareDoubleEventBranch( "track_michelProng_distance", -1.0);
        declareDoubleEventBranch( "track_michelProng_energy", -1.0);
        declareDoubleEventBranch( "track_michelProng_time_diff", -1.0);
        declareDoubleEventBranch( "track_michelProng_end_Z", -1.0);
        declareDoubleEventBranch( "track_michelProng_begin_Z", -1.0);

        // Truth Values
        declareIntTruthBranch("vtx_michel_large_evis_most_pdg", -1);
        declareDoubleTruthBranch("vtx_michel_large_evis_total_truth", -1.0);

        declareIntTruthBranch("vtx_michel_evis_most_pdg", -1);
        declareDoubleTruthBranch("vtx_michel_evis_total_truth", -1.0);

        declareIntTruthBranch("track_michel_evis_most_pdg", -1);
        declareDoubleTruthBranch("track_michel_evis_total_truth", -1.0);
    }


    if (m_study_unused_energy){
        // ------------------------------------------------------------------------       
        // Unused Clusters Visible Energy
        // ------------------------------------------------------------------------       
        // ECAL -- Filled in PreFilterPi0()
        declareDoubleTruthBranch("ecal_unused_evis_total_norm", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_total_truth", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_pizero", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_piplus", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_piminus", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_muon", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_proton", -1.0);
        declareDoubleTruthBranch("ecal_unused_evis_neutron", -1.0);
        // HCAL -- Filled in PreFilterPi0()
        declareDoubleTruthBranch("hcal_unused_evis_total_norm", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_total_truth", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_pizero", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_piplus", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_piminus", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_muon", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_proton", -1.0);
        declareDoubleTruthBranch("hcal_unused_evis_neutron", -1.0);
        // Tracker + ECAL + HCAL -- Filled in PreFilterPi0()
        declareDoubleTruthBranch("other_unused_evis_total_norm", -1.0);
        declareDoubleTruthBranch("other_unused_evis_total_truth", -1.0);
        declareDoubleTruthBranch("other_unused_evis_pizero", -1.0);
        declareDoubleTruthBranch("other_unused_evis_piplus", -1.0);
        declareDoubleTruthBranch("other_unused_evis_piminus", -1.0);
        declareDoubleTruthBranch("other_unused_evis_muon", -1.0);
        declareDoubleTruthBranch("other_unused_evis_proton", -1.0);
        declareDoubleTruthBranch("other_unused_evis_neutron", -1.0);

        // Near Vertex -- Filled in VertexBlob()
        declareIntTruthBranch("vertex_unused_evis_most_pdg", -1);
        declareDoubleTruthBranch("vertex_unused_evis_total_norm", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_total_truth", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_pizero", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_piplus", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_piminus", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_muon", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_proton", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_neutron", -1.0);
        declareDoubleTruthBranch("vertex_unused_evis_gamma", -1.0);

        // Rejected -- Filled in ConeBlobs()
        declareIntTruthBranch("Rejected_unused_evis_most_pdg", -1);
        declareDoubleTruthBranch("Rejected_unused_evis_total_norm", -1.0);
        declareDoubleTruthBranch("Rejected_unused_evis_total_truth", -1.0);
    }

    if (m_study_shower_recovery){
        declareIntEventBranch("OneShower_nClusters", -1);
        declareDoubleEventBranch("OneShower_energy", SENTINEL);
        declareDoubleEventBranch("OneShower_theta", SENTINEL);
        declareDoubleEventBranch("OneShower_dist_vtx", SENTINEL);
        declareIntTruthBranch("OneShower_evis_most_pdg", -1);
        declareDoubleTruthBranch("OneShower_evis_total_norm", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_total_truth", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_pizero", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_piplus", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_piminus", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_muon", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_proton", SENTINEL);
        declareDoubleTruthBranch("OneShower_evis_neutron", SENTINEL);

        declareIntEventBranch("ThreeShower_s1_nClusters", -1);
        declareDoubleEventBranch("ThreeShower_s1_energy", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s1_theta", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s1_dist_vtx", SENTINEL);
        declareIntTruthBranch("ThreeShower_s1_evis_most_pdg", -1);
        declareDoubleTruthBranch("ThreeShower_s1_evis_total_norm", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_total_truth", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_pizero", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_piplus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_piminus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_muon", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_proton", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s1_evis_neutron", SENTINEL);

        declareIntEventBranch("ThreeShower_s2_nClusters", -1);
        declareDoubleEventBranch("ThreeShower_s2_energy", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s2_theta", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s2_dist_vtx", SENTINEL);
        declareIntTruthBranch("ThreeShower_s2_evis_most_pdg", -1);
        declareDoubleTruthBranch("ThreeShower_s2_evis_total_norm", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_total_truth", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_pizero", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_piplus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_piminus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_muon", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_proton", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s2_evis_neutron", SENTINEL);

        declareIntEventBranch("ThreeShower_s3_nClusters", -1);
        declareDoubleEventBranch("ThreeShower_s3_energy", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s3_theta", SENTINEL);
        declareDoubleEventBranch("ThreeShower_s3_dist_vtx", SENTINEL);
        declareIntTruthBranch("ThreeShower_s3_evis_most_pdg", -1);
        declareDoubleTruthBranch("ThreeShower_s3_evis_total_norm", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_total_truth", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_pizero", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_piplus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_piminus", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_muon", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_proton", SENTINEL);
        declareDoubleTruthBranch("ThreeShower_s3_evis_neutron", SENTINEL);
    }

    if (m_study_shower_energy){
        // Blob Digit Energy -- Filled in SaveBlobDigitEnergy 
        declareContainerDoubleEventBranch("gamma1_blob_all_digit_E");
        declareContainerDoubleEventBranch("gamma1_blob_pi0_digit_E");
        declareContainerDoubleEventBranch("gamma1_blob_pi_digit_E");
        declareContainerDoubleEventBranch("gamma1_blob_proton_digit_E");
        declareContainerDoubleEventBranch("gamma1_blob_neutron_digit_E");
        declareContainerDoubleEventBranch("gamma1_blob_muon_digit_E");

        declareContainerDoubleEventBranch("gamma2_blob_all_digit_E");
        declareContainerDoubleEventBranch("gamma2_blob_pi0_digit_E");
        declareContainerDoubleEventBranch("gamma2_blob_pi_digit_E");
        declareContainerDoubleEventBranch("gamma2_blob_proton_digit_E");
        declareContainerDoubleEventBranch("gamma2_blob_neutron_digit_E");
        declareContainerDoubleEventBranch("gamma2_blob_muon_digit_E");

        // Blob True Evis by SubDetector -- Filled in SaveBlobTrueEvisBySubDetector()
        declareDoubleEventBranch("gamma1_trkr_true_evis",-1);
        declareDoubleEventBranch("gamma1_scal_true_evis",-1);
        declareDoubleEventBranch("gamma1_ecal_true_evis",-1);
        declareDoubleEventBranch("gamma1_hcal_true_evis",-1);

        declareDoubleEventBranch("gamma2_trkr_true_evis",-1);
        declareDoubleEventBranch("gamma2_scal_true_evis",-1);
        declareDoubleEventBranch("gamma2_ecal_true_evis",-1);
        declareDoubleEventBranch("gamma2_hcal_true_evis",-1);

        // Blob nHits for Tracker and SCAL -- Filled in SaveBlobTrueEvisBySubdetector()
        declareDoubleEventBranch("gamma1_center_nHits_all",-1);
        declareDoubleEventBranch("gamma1_center_nHits_scal",-1);
        declareDoubleEventBranch("gamma1_center_nHits_trkr",-1);
        declareDoubleEventBranch("gamma1_side_nHits_all",-1);
        declareDoubleEventBranch("gamma1_side_nHits_scal",-1);
        declareDoubleEventBranch("gamma1_side_nHits_trkr",-1);

        declareDoubleEventBranch("gamma2_center_nHits_all",-1);
        declareDoubleEventBranch("gamma2_center_nHits_scal",-1);
        declareDoubleEventBranch("gamma2_center_nHits_trkr",-1);
        declareDoubleEventBranch("gamma2_side_nHits_all",-1);
        declareDoubleEventBranch("gamma2_side_nHits_scal",-1);
        declareDoubleEventBranch("gamma2_side_nHits_trkr",-1);

        // Blob nHits for Tracker and SCAL -- Filled in SaveSCALHits
        declareDoubleEventBranch("gamma1_trkr_nHits_true",-1);
        declareDoubleEventBranch("gamma1_trkr_nHits_reco",-1);
        declareDoubleEventBranch("gamma1_scal_nHits_true",-1);
        declareDoubleEventBranch("gamma1_scal_nHits_reco",-1);

        declareDoubleEventBranch("gamma2_trkr_nHits_true",-1);
        declareDoubleEventBranch("gamma2_trkr_nHits_reco",-1);
        declareDoubleEventBranch("gamma2_scal_nHits_true",-1);
        declareDoubleEventBranch("gamma2_scal_nHits_reco",-1);

        // Blob nHits for Tracker and SCAL -- Filled in SaveSCALHits_Improved
        declareDoubleEventBranch("gamma1_improved_trkr_nHits_true",-1);
        declareDoubleEventBranch("gamma1_improved_trkr_nHits_reco",-1);
        declareDoubleEventBranch("gamma1_improved_scal_nHits_true",-1);
        declareDoubleEventBranch("gamma1_improved_scal_nHits_reco",-1);

        declareDoubleEventBranch("gamma2_improved_trkr_nHits_true",-1);
        declareDoubleEventBranch("gamma2_improved_trkr_nHits_reco",-1);
        declareDoubleEventBranch("gamma2_improved_scal_nHits_true",-1);
        declareDoubleEventBranch("gamma2_improved_scal_nHits_reco",-1);

        // Blob minZ for new method -- Filled in  SaveSCAL_minZ_Info
        declareDoubleEventBranch("gamma1_scal_minZ_evis",-1);
        declareDoubleEventBranch("gamma1_scal_minZ_nDigits",-1);
        declareDoubleEventBranch("gamma2_scal_minZ_evis",-1);
        declareDoubleEventBranch("gamma2_scal_minZ_nDigits",-1);
    }

    return sc;
}

//==============================================================================
//
// reconstructEvent() --
//
//==============================================================================
StatusCode CCProtonPi0::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    if(truthEvent) info()<<"This is MC Event!"<<endmsg;

    Minerva::GenMinHeader* header(NULL);
    if (exist<Minerva::GenMinHeader>(Minerva::GenMinHeaderLocation::Default)) {
        header = get<Minerva::GenMinHeader>(Minerva::GenMinHeaderLocation::Default); 
    }

    if (header) {
        int run    = header->RunNumber();
        int subrun = header->SubRunNumber();
        int gate   = header->SpillNumber();
        info()<<"Reconstruct (run, subrun, gate): (" <<run<<", "<<subrun<<", "<< gate<<")"<<endmsg;
    }

    //--------------------------------------------------------------------------
    // Count Events Enters reconstructEvent()
    //      N_reconstructEvent is used as the reco_eventID
    //--------------------------------------------------------------------------
    debug()<<"reconstructEvent() Count = "<<N_reconstructEvent<<" tagTruth() Count = "<<N_tagTruth-1<<endmsg;
    event->setDoubleData("reco_eventID", N_reconstructEvent);
    N_reconstructEvent++;

    bool shouldRecoEvent = ShouldReconstructEvent(event,truthEvent);
    if(!shouldRecoEvent){
        debug()<<"Skipping Event!"<<endmsg;
        return StatusCode::SUCCESS;
    }

    //--------------------------------------------------------------------------
    // Reset Mutable Member Variables
    //--------------------------------------------------------------------------
    // SmartRefs
    m_MuonProng = NULL;
    m_MuonParticle = NULL;
    m_Pi0Blob1 = NULL;
    m_Pi0Blob2 = NULL;

    // Vectors
    //m_PionProngs.clear();
    //m_PionParticles.clear();
    m_ProtonProngs.clear();
    m_ProtonParticles.clear();

    //m_HadronProngs.clear();
    //m_HadronParticles.clear();

    fTrajectoryMap.clear();

    FillTrajectoryMap();

    SaveTrajectories(event);
    //==========================================================================
    // Vertex Reconstruction
    //==========================================================================
    debug() << "START: Vertex Reconstruction..." << endmsg;

    if (!hasEventVertex(event) ){
        if( m_keepAfter_vertex_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    if ( !vertexInRecoVolume(event) ){
        if( m_keepAfter_vertex_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    RefitVertex_Using_AnchoredShortTracks(event);

    if ( !vertexInFiducialVolume(event) ){
        if( m_keepAfter_vertex_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    SetVertexCount(event);

    setVertexData(event);

    debug() << "FINISH: Vertex Reconstruction!" << endmsg;

    //==========================================================================
    // Muon Reconstruction
    //==========================================================================
    debug() << "START: Muon Reconstruction..." << endmsg;

    if( !hasEventMinosMatchedMuon(event) ){
        if( m_keepAfter_muon_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    if ( !isMuonChargeNegative(event)){
        if( m_keepAfter_muon_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    tagPrimaryMuon(event);

    GetMuonExtraEnergy(event);

    bool muonFilled = setMuonData(event);
    if( !muonFilled ){ 
        error()<<"Muon NTuple Branches did not filled!"<<endmsg;
        return StatusCode::SUCCESS;
    }

    debug() << "FINISH: Muon Reconstruction" << endmsg;

    //==========================================================================
    // Michel Electrons
    //==========================================================================
    debug()<<"START: Michel Electron Search"<<endmsg;

    bool has_michel_event = ImprovedMichel_EventHasMichel(event);
    bool has_michel_vertex = VertexHasMichels(event, truthEvent);
    bool has_michel_trackend = TrackEndPointHasMichels(event, truthEvent);

    if (m_removeEvents_withMichel && has_michel_event) {
        if( m_keepAfter_michel_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    if (m_removeEvents_withMichel && has_michel_vertex ){
        if( m_keepAfter_michel_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    if ( m_removeEvents_withMichel && has_michel_trackend ){
        if( m_keepAfter_michel_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    debug()<<"FINISH: Michel Electron Search"<<endmsg;

    //==========================================================================
    // Discard Far Tracks 
    //      They might be tracked photons
    //      Updates nTracks
    //==========================================================================

    ProngVect before_DiscardFarTracks = event->primaryProngs();
    debug()<<"before_DiscardFarTracks Size = "<<before_DiscardFarTracks.size()<<endmsg;
    DiscardFarTracks(event);
    ProngVect after_DiscardFarTracks = event->primaryProngs();
    debug()<<"after_DiscardFarTracks Size = "<<after_DiscardFarTracks.size()<<endmsg;

    //==========================================================================
    // Proton Reconstruction
    //==========================================================================
    debug() << "START: Proton Reconstruction" << endmsg;

    ProngVect primaryProngs = event->primaryProngs();
    int nPrimaryProngs = primaryProngs.size();

    if (nPrimaryProngs > 1){     

        bool makeParticles = createTrackedParticles(event);
        if (!makeParticles){
            debug() << "Creation of Particles are FAILED!"<< endmsg;
            event->setIntData("Cut_Particle_None",1);
            if( m_keepAfter_proton_cuts ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }

        bool foundProton = getProtonProng(event);
        if( !foundProton ) {
            debug() << "Didn't find any contained in the tracker bit-positive prong with a proton particle!" << endmsg;
            event->setIntData("Cut_Proton_None",1);
            if( m_keepAfter_proton_cuts ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
    }

    // Set Proton Kinematics
    if( m_ProtonParticles.size() > 0){
        bool protonFilled = setProtonData( event );
        if( !protonFilled ){ 
            debug()<<"Proton Momentum is NaN, rejecting event!"<<endmsg;
            event->setIntData("Cut_Proton_Bad",1);
            if( m_keepAfter_proton_cuts ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
    }else{
        debug()<<"No Proton Particle, Setting -9.9 to m_proton_4P"<<endmsg;
        m_proton_4P.SetPxPyPzE(-9.9,-9.9,-9.9,-9.9);
    }

    debug() <<"Found "<<m_ProtonParticles.size()<<" Good Proton Candidates!"<<endmsg;
    event->setIntData("nProtonCandidates", (int)m_ProtonParticles.size());

    //--------------------------------------------------------------------------
    // Debugging: Check values
    debug()<<"m_proton_4P = ( "<<m_proton_4P.px()<<", "<<m_proton_4P.py()<<", "<<m_proton_4P.pz()<<", "<<m_proton_4P.E()<<" )"<<endmsg;
    //--------------------------------------------------------------------------

    debug()<<"FINISH: Proton Reconstruction"<<endmsg;

    //==========================================================================
    // Pi0 Reconstruction
    //==========================================================================
    debug()<<"START: Pi0 Reconstruction"<<endmsg;

    if ( !PreFilterPi0(event, truthEvent) ){
        event->setIntData("Cut_PreFilter_Pi0",1);
        if( m_keepAfter_pi0_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }

    SmartRefVector<Minerva::IDCluster> beforeVertexBlob; 
    beforeVertexBlob = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    debug()<<"N(Unused) Before VertexBlob = "<<beforeVertexBlob.size()<<endmsg; 

    VertexBlob(event,truthEvent);

    SmartRefVector<Minerva::IDCluster> beforeConeBlobs; 
    beforeConeBlobs = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    debug()<<"N(Unused) Before ConeBlobs = "<<beforeConeBlobs.size()<<endmsg; 


    // MAKE CUT - If ConeBlobs Can NOT Find Two Blobs
    bool FoundTwoBlobs = ConeBlobs(event, truthEvent);
    if ( !FoundTwoBlobs ){
        event->setIntData("Cut_ConeBlobs",1);
        if( m_keepAfter_pi0_cuts ) return interpretFailEvent(event);
        else return StatusCode::SUCCESS;  
    }

    SmartRefVector<Minerva::IDCluster> afterConeBlobs; 
    afterConeBlobs = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    debug()<<"N(Unused) After ConeBlobs = "<<afterConeBlobs.size()<<endmsg; 

    // Set Data for Blobs found in ConeBlobs()
    setBlobData(event, truthEvent);

    if ( !AreBlobsDirectionGood(event) ){
        event->setIntData("Cut_BlobDirectionBad",1);
        if( m_keepAfter_pi0_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }

    // Get Calorimetric Unused Energy after Pi0 Reconstruction
    DispersedBlob(event, truthEvent);

    // Set Pi0 Kinematics
    bool pi0Filled = setPi0Data( event );
    if( !pi0Filled ){
        debug()<<"Pi0 Momentum is NaN, rejecting event!"<<endmsg;
        event->setIntData("Cut_Pi0_Bad",1);
        if( m_keepAfter_pi0_cuts ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }

    debug()<<"FINISH: Pi0 Reconstruction"<<endmsg;

    //--------------------------------------------------------------------------
    //
    // Particle Reconstruction Finished 
    //
    //--------------------------------------------------------------------------

    // Get Extra Energy after Particle Reconstructions
    SaveExtraEnergyLeftover(event);
    SetExtraEnergyTotal(event);

    // Write FS Particle Table and Event Record for Reconstructed Events
    if(truthEvent){
        bool isSignal = truthEvent->filtertaglist()->isFilterTagTrue("isSignal");
        if (m_writeFSParticle_Table){
            info()<<"FS Particle Table and Event Record for Reconstructed Event!"<<endmsg;
            writeBackgroundType(truthEvent);
            writeFSParticleTable(isSignal);
            writeEventRecord(truthEvent,isSignal);
        }
    }

    ColorUnusedIDClusters(event);
    SaveEventTime(event);
    //fillCommonPhysicsAnaBranches( event );

    //--------------------------------------------------------------------------
    // Call the interpretEvent function.
    //--------------------------------------------------------------------------
    NeutrinoVect interactions;
    StatusCode interpret = this->interpretEvent( event, truthEvent, interactions );

    // If there were any neutrino interactions reconstructed, mark the event.
    if( interactions.size() == 1 ){
        markEvent( event );
    }
    else {
        warning()<<"interpretEvent() returns "<<interactions.size()<<" interactions!"<<endmsg;
        return StatusCode::SUCCESS; // We didn't crash.
    }

    // Add the interactions to the event. Use MinervaHistoTool::addInteractionHyp.
    StatusCode sc = addInteractionHyp( event, interactions );

    return sc;    
}

//==============================================================================
//
// interpretEvent()
//
//==============================================================================
StatusCode CCProtonPi0::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truthEvent, NeutrinoVect& interaction_hyp ) const
{
    if( truthEvent ){
        debug() << "This is a MC event." << endmsg;
    }

    if( !event ){
        warning() << "NULL Event" << endmsg;
        return StatusCode::FAILURE; // we crashed!
    }

    //--------------------------------------------------------------------------
    // Create interaction hypothesis
    //--------------------------------------------------------------------------
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    interaction_hyp.push_back( nuInt );

    //--------------------------------------------------------------------------
    // Identify and Store - Muon and Proton Prongs
    //--------------------------------------------------------------------------
    ProngVect primaryProngs = event->primaryProngs();

    SmartRef<Prong> muonProng = (Prong*)NULL;
    ProngVect protonProngs;
    ParticleVect protonParticles;

    for (ProngVect::iterator itProng = primaryProngs.begin(); itProng != primaryProngs.end(); ++itProng) {

        bool isPrimaryMuon = false;
        bool isPrimaryProton = false;

        if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryMuon") ) {
            (*itProng)->filtertaglist()->checkFilterTag( "PrimaryMuon", isPrimaryMuon );
        }

        if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryProton") ) {
            (*itProng)->filtertaglist()->checkFilterTag( "PrimaryProton", isPrimaryProton );
        }    

        if (isPrimaryMuon && !isPrimaryProton){
            muonProng = *itProng;
        }

        if (isPrimaryProton && !isPrimaryMuon){
            protonProngs.push_back(*itProng);
            protonParticles.push_back((*itProng)->bestParticle());
        }

        if (isPrimaryMuon && isPrimaryProton ) {
            warning()<<"Prong is two primary particles!"<<endmsg;
        }
    }

    // Neutrino Interaction Vertex 
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

    Gaudi::XYZTVector vtx_position( vertex->position().x(), 
            vertex->position().y(),
            vertex->position().z(),
            event->time() );     
    nuInt->setVertex( vtx_position );
    nuInt->setScore( 1.0 );

    // Primary Lepton
    nuInt->setLeptonEnergy( m_muon_4P );
    fillMinosMuonBranches(nuInt, m_MuonProng);

    // Interaction Parameters
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );

    //--------------------------------------------------------------
    // Hadronic Recoil
    //--------------------------------------------------------------
    double hadron_recoilE_CCInc = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "CCInclusive", -1); 
    double hadron_recoilE = RecoilUtils->calcRecoilEFromClusters( event, muonProng, "Default");
 
    debug()<<"hadron_recoilE_CCInc = "<<hadron_recoilE_CCInc<<" at event time "<<event->time()<<" used muon: "<<event->filtertaglist()->isFilterTagTrue("used_muon_time")<<endmsg;
    debug()<<"hadron_recoilE_Default = "<<hadron_recoilE<<" at event time "<<event->time()<<" used muon: "<<event->filtertaglist()->isFilterTagTrue("used_muon_time")<<endmsg;

    nuInt->setDoubleData("hadron_recoil_CCInc", hadron_recoilE_CCInc);
    nuInt->setDoubleData("hadron_recoil", hadron_recoilE);
    
    debug()<<"Finished hadronic recoil"<<endmsg;

    //--------------------------------------------------------------------------
    // Truth Matching for Muon and Proton Prongs
    //--------------------------------------------------------------------------
    if( haveNeutrinoMC() && m_DoTruthMatch ) {      
        // Muon Truth Matching
        Minerva::ProngVect muonProngs;
        muonProngs.push_back( muonProng );
        setTrackProngTruth(nuInt,muonProngs);
        debug() << "       set the muon prong track truth information" << endmsg;

        // Proton Turth Matching
        if (protonProngs.size() == 0 ){
            debug() << "       No Proton Candidate! No Truth Matching" << endmsg;
        }else{
            setTrackProngTruth(nuInt,protonProngs);
            debug() << "       set the proton prong track truth information" << endmsg;
        }
    }

    return StatusCode::SUCCESS;
}


//==============================================================================
// tagTruth()
//==============================================================================
StatusCode CCProtonPi0::tagTruth( Minerva::GenMinInteraction* truthEvent ) const 
{
    //--------------------------------------------------------------------------
    // Count Events Enters tagTruth()
    //      N_tagTruth is used as the eventID to locate event in DST
    //      Required for Arachne Scannings
    //--------------------------------------------------------------------------
    debug()<<"tagTruth() Count = "<<N_tagTruth<<endmsg;
    truthEvent->setDoubleData("eventID", N_tagTruth);
    N_tagTruth++;

    if (!truthEvent) {
        warning()<<"Passed a null truthEvent to tagTruth()!"<<endmsg;
        return StatusCode::SUCCESS;
    }

    isTrueVertexFiducial(truthEvent);
    Calc_TrueEventKinematics(truthEvent);

    //--------------------------------------------------------------------------
    // Decide To Reconstruct Event or NOT
    //--------------------------------------------------------------------------

    // Return immediately - If True Vertex is NOT inside Fid Volume
    if ( !isTrueVertexFiducial(truthEvent) ){
        debug() <<"Do NOT Reconstruct Event - True vtx is not in Fiducial Volume!" << endmsg;
        truthEvent->filtertaglist()->setOrAddFilterTag( "ReconstructEvent", false );
	      debug() << "Event Tagged as no Reconstructable" << endmsg;
        return StatusCode::SUCCESS;
    }

    //--------------------------------------------------------------------------
    // Fill GENIE Weight Branches
    //--------------------------------------------------------------------------
    StatusCode sc = fillGenieWeightBranches( truthEvent );
    if (sc.isFailure() ) {
        warning()<<"Genie weight branch filling failed!"<<endmsg;
        return sc;
    }

    //--------------------------------------------------------------------------
    // Tag Event as Signal or Background
    //--------------------------------------------------------------------------
    bool isSignal = tagSignal(truthEvent);
    if (isSignal){
        setPi0GenieRecord(truthEvent);  
        setSignalKinematics(truthEvent);
    }else{
        // Two Different Backround Sets
        tagBackground(truthEvent);
        tagBackground_Compact(truthEvent);
        tagBackgroundWithPi0(truthEvent);
    }

    //--------------------------------------------------------------------------
    // Create Final State Particle Table
    //--------------------------------------------------------------------------
    if (m_writeFSParticle_Table){ 
        writeFSParticleTable(isSignal);
        writeEventRecord(truthEvent,isSignal);
    }

    //------------------------------------------------------------
    // Set Target Material 
    //------------------------------------------------------------
    setTargetMaterial(truthEvent);

    //--------------------------------------------------------------
    // Set Vertex module and plane
    //--------------------------------------------------------------   
    Gaudi::LorentzVector t_vtx = truthEvent->Vtx();    
    int vertex_module;
    int vertex_plane;

    getNearestPlane(t_vtx.z(), vertex_module, vertex_plane);  

    truthEvent->setIntData("vertex_module", vertex_module);
    truthEvent->setIntData("vertex_plane", vertex_plane);

    //--------------------------------------------------------------   
    // Other Studies
    //--------------------------------------------------------------   
    tagSignal_BeforeFSI(truthEvent);
    //PrintZPositionMaterial();

    //--------------------------------------------------------------------------
    // Decide whether to Reconstruct the Event or NOT
    //    Set ReconstructEvent true to Reconstruct All Events
    //--------------------------------------------------------------------------
    truthEvent->filtertaglist()->setOrAddFilterTag( "ReconstructEvent", true);

    return StatusCode::SUCCESS;
}

void CCProtonPi0::PrintZPositionMaterial() const
{
    G4Navigator* m_geantNav = new G4Navigator;
    m_geantNav->SetWorldVolume(m_gigaCnvSvc->world());
    
    for (double z = 4293; z < 9987; ++z){
        G4ThreeVector p(0.0, 0.0, z); // The point you want to look up
        G4ThreeVector v(1,0,0); // Need a direction to pass to LocateGlobalPointAndSetup
        G4VPhysicalVolume* physVol = m_geantNav->LocateGlobalPointAndSetup(p, &v);
        G4Material* material = physVol->GetLogicalVolume()->GetMaterial();

        std::string material_name = material->GetName();
        //int material_nElements = material->GetNumberOfElements();
        //double material_Z = material_nElements == 1 ? material->GetZ() : -1;
        
        debug()<<"pos_z = "<<z<<" Name = "<<material_name<<endmsg;
    }

    delete m_geantNav;
}

//==============================================================================
//  truthIsPlausible()
//      truthIsPlausible() is now a pure virtual function of MinervaAnalysisTool, so you must implement it.
//      It is called automatically by PhysicsEventAnalysisAlg AFTER reconstructEvent() and interpretEvent() are run.
//      If it returns false, the event is not included in the analysis ntuple
//==============================================================================
bool CCProtonPi0::truthIsPlausible( const Minerva::PhysicsEvent*  ) const
{
    // I am checking only Muon Prong right now
    if (m_MuonProng == NULL){
        return true;
    }else{
        double mc_frac = -1.0;
        return muonIsPlausible( m_MuonProng, mc_frac);
    }
}

//==============================================================================
//  Finalize
//==============================================================================
StatusCode CCProtonPi0::finalize()
{
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );

    return sc;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  Private Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//------------------------------------------------------------------------------
// interpret Events which fails the reconstructor cuts
//------------------------------------------------------------------------------
StatusCode CCProtonPi0::interpretFailEvent( Minerva::PhysicsEvent* event ) const
{  
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    NeutrinoVect nuInts;
    nuInts.push_back( nuInt );
    markEvent(event);
    addInteractionHyp(event,nuInts);
    //fillCommonPhysicsAnaBranches(event);
    //fillNuMIBranches(event);

    return StatusCode::SUCCESS;
}

//==============================================================================
// Set Vertex Data 
//==============================================================================
void CCProtonPi0::setVertexData( Minerva::PhysicsEvent* event ) const 
{
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

    Gaudi::XYZTVector vtx_position( vertex->position().x(), 
            vertex->position().y(),
            vertex->position().z(),
            event->time() );     
    //event->setVertex( vtx_position );
    //event->setScore( 1.0 );
    event->setDoubleData("vtx_x", vtx_position.x() );
    event->setDoubleData("vtx_y", vtx_position.y() );
    event->setDoubleData("vtx_z", vtx_position.z() );

    // Get Vertex Module and Planes
    int vtx_module, vtx_plane;
    debug()<<"Calling getNearestPlane, vtx is "<<vtx_position.z()<<endmsg;
    getNearestPlane(vtx_position.z(), vtx_module, vtx_plane); 
    event->setIntData("vtx_module", vtx_module);
    event->setIntData("vtx_plane", vtx_plane);
}

//==============================================================================
// Set Muon particle data
//==============================================================================
bool CCProtonPi0::setMuonData( Minerva::PhysicsEvent *event ) const 
{
    const double SENTINEL = -9.9;
    // Sanity Check
    if (!m_MuonProng || !m_MuonParticle) {
        fatal()<< "NO Muon!!" <<endmsg;
        return false;
    }

    //--------------------------------------------------------------------------
    // Get Muon Track Information
    //--------------------------------------------------------------------------
    SmartRefVector<Track>::iterator iterTrk;
    SmartRefVector<Track> muonTracks = m_MuonProng->minervaTracks();

    int muon_minervaTrack_types = 0;
    int muon_N_minosTracks = 0;
    int nlong = 0;
    int nshort = 0;

    // Number of Short and Long Tracks
    for (iterTrk = muonTracks.begin(); iterTrk != muonTracks.end(); ++iterTrk) {
        if ((*iterTrk)->type() == Track::Long) nlong++;
        if ((*iterTrk)->type() == Track::Short) nshort++;
    }

    if (nlong > 0 && nshort == 0) muon_minervaTrack_types = 1;
    else if (nlong == 0 && nshort > 0) muon_minervaTrack_types = 2;
    else if (nlong >0 && nshort >0) muon_minervaTrack_types = 3;

    // Get number of MINOS Tracks
    muon_N_minosTracks = m_MuonProng->minosTracks().size();   

    //--------------------------------------------------------------------------
    // Get Muon Kinematics
    //--------------------------------------------------------------------------

    // 4-Momentum and Angle Information
    Gaudi::LorentzVector muon_4p    = m_MuonParticle->momentumVec();
    double muon_px    = muon_4p.px();
    double muon_py    = muon_4p.py();
    double muon_pz    = muon_4p.pz();    
    double muon_E     = muon_4p.E();
    double muon_p     = muon_4p.P();
    double muon_KE    = muon_4p.E() - MinervaUnits::M_mu; 
    double muon_phi = muon_4p.phi();
    double muon_theta = muon_4p.theta();
    double muon_phi_beam = m_coordSysTool->phiWRTBeam(muon_4p);
    double muon_theta_beam = m_coordSysTool->thetaWRTBeam(muon_4p);
    double muon_thetaX_beam = m_coordSysTool->thetaXWRTBeam(muon_4p);
    double muon_thetaY_beam = m_coordSysTool->thetaYWRTBeam(muon_4p);
    double muon_theta_beam_biasUp = m_coordSysTool->thetaWRTBeam(muon_4p,m_beamAngleBias) - muon_theta_beam;
    double muon_theta_beam_biasDown = m_coordSysTool->thetaWRTBeam(muon_4p, -1.0*m_beamAngleBias) - muon_theta_beam;

    // Get Theta in All Track Nodes
    Minerva::Track* muonTrack = m_MuonProng->minervaTracks()[0];
    std::vector<double> thetaNodes;
    std::vector<double> thetaXNodes;
    std::vector<double> thetaYNodes;

    for(uint i = 0; i < muonTrack->nNodes(); i++){
        const Minerva::Node* node = muonTrack->nodes()[i];
        double ax = node->state().ax();
        double ay = node->state().ay();
        double theta = m_coordSysTool->thetaWRTBeam( ax, ay, 1 );
        double thetaX = m_coordSysTool->thetaXWRTBeam( ax, ay, 1 );
        double thetaY = m_coordSysTool->thetaYWRTBeam( ay, 1 );
        thetaNodes.push_back(theta);
        thetaXNodes.push_back(thetaX);
        thetaYNodes.push_back(thetaY);
    }

    // Muon Score
    double muon_muScore = m_MuonParticle->score();

    // Muon Q/P / Q/P-Error
    int is_minos_track = -1;
    int is_minos_stub = -1;
    int muon_minosTrackQuality = -99;
    double muon_qp = SENTINEL; 
    double muon_qpqpe = MuonUtils->minosQPQPE(m_MuonProng);
    if ( m_MuonProng->MinosTrack() ) {
        is_minos_track = 1;
        SmartRef<MinosRecoTrack> minosTrack = m_MuonProng->minosTracks()[0];
        muon_qp = minosTrack->qp();
        muon_minosTrackQuality = minosTrack->quality();      
    }
    else if ( m_MuonProng->MinosStub() ) {
        is_minos_stub = 1;
    }
    else{
        warning()<<"No MINOS track or stub!  How is this a muon?"<<endmsg;
    }

    // Muon energy in road upstream info
    int muon_roadUpstreamPlanes = -99;
    double muon_roadUpstreamEnergy = AnaToolUtils->energyInTheRoadUpstream(m_MuonProng, muon_roadUpstreamPlanes);  

    // Muon Energy Shift
    //      the shifts are so small compared to MINOS-matched momentum that we can 
    //      approximate the momentum shift as the energy shift (with ~0.2% at 1.5 GeV/c)
    double muon_E_shift = MuonUtils->calculateMomentumCorrection(m_MuonProng);

    //Get Muon Charge
    int muon_charge = -99;
    MuonUtils->muonCharge(m_MuonProng,muon_charge); 

    // Write Muon 4-Momentum to Global Variable
    m_muon_4P = muon_4p;

    //--------------------------------------------------------------------------
    // Debugging: Check values
    debug()<<"m_muon_4P = ( "<<m_muon_4P.px()<<", "<<m_muon_4P.py()<<", "<<m_muon_4P.pz()<<", "<<m_muon_4P.E()<<" )"<<endmsg;
    debug()<<"P4(Muon) = ( "<<muon_4p.px()<<", "<<muon_4p.py()<<", "<<muon_4p.pz()<<", "<<muon_4p.E()<<" )"<<endmsg;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Fill Muon Branches
    //--------------------------------------------------------------------------
    event->setIntData("muon_hasMinosMatchTrack", is_minos_track );
    event->setIntData("muon_hasMinosMatchStub", is_minos_stub );
    event->setIntData("muon_minervaTrack_types", muon_minervaTrack_types);
    event->setIntData("muon_N_minosTracks", muon_N_minosTracks);
    event->setIntData("muon_minosTrackQuality", muon_minosTrackQuality);
    event->setIntData("muon_roadUpstreamPlanes", muon_roadUpstreamPlanes);
    event->setIntData("muon_charge",muon_charge);
    event->setDoubleData("muon_roadUpstreamEnergy", muon_roadUpstreamEnergy);
    event->setDoubleData("muon_px",muon_px);
    event->setDoubleData("muon_py",muon_py);
    event->setDoubleData("muon_pz",muon_pz);
    event->setDoubleData("muon_E",muon_E);
    event->setDoubleData("muon_P",muon_p);
    event->setDoubleData("muon_KE",muon_KE);
    event->setDoubleData("muon_phi",muon_phi);
    event->setDoubleData("muon_phi_beam",muon_phi_beam);
    event->setDoubleData("muon_theta", muon_theta);
    event->setDoubleData("muon_theta_beam", muon_theta_beam);
    event->setDoubleData("muon_thetaX_beam", muon_thetaX_beam);
    event->setDoubleData("muon_thetaY_beam", muon_thetaY_beam);
    event->setDoubleData("muon_theta_beam_biasUp", muon_theta_beam_biasUp);
    event->setDoubleData("muon_theta_beam_biasDown", muon_theta_beam_biasDown);
    event->setDoubleData("muon_muScore", muon_muScore);
    event->setDoubleData("muon_qp",muon_qp );
    event->setDoubleData("muon_qpqpe",muon_qpqpe);
    event->setDoubleData("muon_E_shift",muon_E_shift);
    event->setContainerDoubleData("muon_theta_allNodes", thetaNodes);
    event->setContainerDoubleData("muon_thetaX_allNodes", thetaXNodes);
    event->setContainerDoubleData("muon_thetaY_allNodes", thetaYNodes);

    return true;
}

double CCProtonPi0::Calc_QSq(double Enu, double muon_E, double muon_P, double muon_angle_beam) const
{
    double muon_mass = MinervaUnits::M_muon;
    double QSq = 2*Enu*(muon_E - muon_P*cos(muon_angle_beam))-(muon_mass*muon_mass);

    return QSq;
}

double CCProtonPi0::Calc_WSq(double Enu, double QSq, double muon_E) const
{
    double neutron_mass = MinervaUnits::M_neutron;
    double WSq = neutron_mass*neutron_mass + 2*neutron_mass*(Enu - muon_E) - QSq; 

    return WSq;
}

/*
double CCProtonPi0::Calc_Enu_1Track(double vertex_energy, double extra_energy) const
{
    debug()<<"Calculating Enu_1Track Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tP4(Pi0) = "<<m_pi0_4P<<endmsg;
    debug()<<"\tVertex Energy = "<<vertex_energy<<endmsg;
    debug()<<"\tExtra Energy = "<<extra_energy<<endmsg;

    double Enu = m_muon_4P.E() + m_pi0_4P.E() + vertex_energy + extra_energy;

    debug()<<"\tNeutrino Energy = "<<Enu<<endmsg;

    return Enu;
}


// Using Equation from DocDB: 9749 v2
// Back-of-napkin derivations for E_nu of CC(nu, meson) reactions
double CCProtonPi0::Calc_Enu_1Track_Alt() const
{
    debug()<<"Calculating Enu_1Track_Alt Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tP4(Pi0) = "<<m_pi0_4P<<endmsg;
    debug()<<"\tT(proton) = from Transverse Momentum Balance"<<endmsg;

    // Main Components
    const double Mn = MinervaUnits::M_n;    // Neutron Rest Mass [MeV]
    double Enu;     // Neutrino Energy
    double TN;      // Hadron Kinetic Energy
    double Emu = m_muon_4P.E();
    double Epi0 = m_pi0_4P.E();
    double Pmu_long = Calc_Longitudinal_Momentum(m_muon_4P);
    double Ppi0_long = Calc_Longitudinal_Momentum(m_pi0_4P);
    double Px_mu_wrt_Beam = Calc_Px_wrt_Beam(m_muon_4P);
    double Py_mu_wrt_Beam = Calc_Py_wrt_Beam(m_muon_4P);
    double Px_pi0_wrt_Beam = Calc_Px_wrt_Beam(m_pi0_4P);
    double Py_pi0_wrt_Beam = Calc_Py_wrt_Beam(m_pi0_4P);

    // TN Calculation in Parts -- TN is very long to calculate in a single Line
    double Emu_Pmu_long = Emu - Pmu_long;
    double Epi0_Ppi0_long = Epi0 - Ppi0_long;
    double P_long = std::pow((Emu_Pmu_long + Epi0_Ppi0_long),2);
    // Transverse Component (Vectoral Calculation Page 4)
    double Pmu_trans = std::pow(Px_mu_wrt_Beam,2) + std::pow(Py_mu_wrt_Beam,2);
    double Ppi0_trans = std::pow(Px_pi0_wrt_Beam,2) + std::pow(Py_pi0_wrt_Beam,2);
    double Pmu_trans_Ppi0_trans = 2 * ( (Px_mu_wrt_Beam * Px_pi0_wrt_Beam) +
            (Py_mu_wrt_Beam * Py_pi0_wrt_Beam) );
    double P_trans = Pmu_trans + Ppi0_trans + Pmu_trans_Ppi0_trans;

    // Form Hadron Kinetic Energy
    TN = 0.5 * (P_long + P_trans) / (Mn - Emu_Pmu_long - Epi0_Ppi0_long);

    Enu = Emu + Epi0 + TN;

    return Enu;
}

double CCProtonPi0::Calc_Enu_2Track(double vertex_energy, double extra_energy) const
{
    debug()<<"Calculating Enu_2Track Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tP4(Proton) = "<<m_proton_4P<<endmsg;
    debug()<<"\tP4(Pi0) = "<<m_pi0_4P<<endmsg;
    debug()<<"\tVertex Energy = "<<vertex_energy<<endmsg;
    debug()<<"\tExtra Visible Energy = "<<extra_energy<<endmsg;

    double Mp = MinervaUnits::M_p; // Proton Rest Mass [MeV]
    double Emu = m_muon_4P.E();
    double Epi0 = m_pi0_4P.E();
    double KEproton = m_proton_4P.E() - Mp;

    // Calculate Enu -- Use eq. in Research Log Book page 29
    double Enu = Emu + Epi0 + KEproton + vertex_energy + extra_energy;

    debug()<<"\tNeutrino Energy = "<<Enu<<endmsg;

    return Enu;
}

double CCProtonPi0::Calc_Longitudinal_Momentum(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double P_long = particle_4P.P() * std::cos(Theta_wrt_Beam);
    return P_long;
}

double CCProtonPi0::Calc_Px_wrt_Beam(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double Phi_wrt_Beam = m_coordSysTool->phiWRTBeam(particle_4P);
    double Px_wrt_Beam = particle_4P.P() * std::sin(Theta_wrt_Beam) * std::cos(Phi_wrt_Beam);
    return Px_wrt_Beam;
}

double CCProtonPi0::Calc_Py_wrt_Beam(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double Phi_wrt_Beam = m_coordSysTool->phiWRTBeam(particle_4P);
    double Py_wrt_Beam = particle_4P.P() * std::sin(Theta_wrt_Beam) * std::sin(Phi_wrt_Beam);
    return Py_wrt_Beam;
}

//==============================================================================
// Calculates and Set Event Kinematics
//      Uses 4-Momentum of Final State Particles
//==============================================================================
void CCProtonPi0::setEventKinematics(const Minerva::PhysicsEvent* event, Minerva::NeutrinoInt* nuInt) const
{
    double vertex_energy = GetVertexEnergy(event);
    double extra_energy = GetTotalExtraEnergy(event);
    double Enu;
    double QSq;
    double WSq;

    // Enu Calculation depends on Proton Reconstruction
    if (m_ProtonParticles.size() == 0) Enu = Calc_Enu_1Track(vertex_energy, extra_energy);
    else Enu = Calc_Enu_2Track(vertex_energy, extra_energy);

    QSq = Calc_QSq(Enu);
    WSq = Calc_WSq(Enu, QSq);

    // Fill NTuples
    nuInt->setDoubleData("Vertex_Energy",vertex_energy);
    nuInt->setDoubleData("Extra_Energy_Total",extra_energy);
    nuInt->setDoubleData("neutrino_E",Enu);
    nuInt->setDoubleData("QSq",QSq);
    nuInt->setDoubleData("WSq",WSq);


    // Calculate Enu Alternatively for 1 Track
    // Estimate Proton KE from Transverse Momentum Balance
    if (m_ProtonParticles.size() == 0){
        double Enu_Alt = Calc_Enu_1Track_Alt();
        double QSq_Alt = Calc_QSq(Enu_Alt);
        double WSq_Alt = Calc_WSq(Enu_Alt, QSq_Alt);

        nuInt->setDoubleData("neutrino_E_1Track_Alt",Enu_Alt);
        nuInt->setDoubleData("QSq_1Track_Alt",QSq_Alt);
        nuInt->setDoubleData("WSq_1Track_Alt",WSq_Alt);
    }

}
*/
//==============================================================================
// Find the plane nearest to a point
//==============================================================================
StatusCode CCProtonPi0::getNearestPlane(double z, int & module_return, int & plane_return) const
{
    // Got From Brandon's CCNuPionInc - 2014_04_14
    // testing new MinervaDet routine
    // works as advertised, but a couple problems:
    // 1) Does not go to more downstream plane (upstream if backwards track) if point is in passive material

    DePlane const * pPlane = m_InnerDetector->getClosestDePlane(z);
    if (!pPlane) {
        error()<<"getClosestDePlane failed on z "<<z<<endmsg;
        return StatusCode::FAILURE;
    }

    PlaneID planeid = pPlane->getPlaneID();
    module_return = planeid.module();
    plane_return = planeid.plane();

    return StatusCode::SUCCESS;

}


//==============================================================================
// Created particles for negative bit Minos prongs
//==============================================================================
bool CCProtonPi0::createTrackedParticles(Minerva::PhysicsEvent *event ) const
{
    ProngVect prongs =  event->primaryProngs();

    bool makeParticles = false;

    // Check if the prongs are odMatch
    Minerva::EventMgr* mgr = getEventMgr(prongs[0]);
    Minerva::ODClusterVect odClusters = mgr->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
    m_odMatchTool->classifyProngs(prongs,odClusters);

    // Loop over prongs
    for(unsigned int p = 0; p < prongs.size(); p++) {
        debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() << endmsg;

        // Make sure the prong is not a bad object
        bool pass = true; std::string tag = "BadObject";
        if( prongs[p]->filtertaglist()->filterTagExists(tag) ) {
            if( prongs[p]->filtertaglist()->checkFilterTag(tag,pass) ) {
                error() << "This prong = " << prongs[p]->typeBitsToString() << " has been flagged as a \"BadObject\", skipping!" << endmsg;
                continue;
            }
        }

        // Skip prongs with particles
        debug()<< "The prong has n particles = " << prongs[p]->particles().size() << endmsg;
        if( !prongs[p]->particles().empty() ){
            debug()<<"Skipping Prong: Prong has attached particles!"<<endmsg;
            continue;
        }

        // Select the particle hypotheses candidates and particle tools
        std::vector< Minerva::Particle::ID > hypotheses;
        IParticleMakerTool::NameAliasListType toolsToUse;

        // Do NOT use prongs which have OD Match
        if( !prongs[p]->OdMatch() ) {
            hypotheses.push_back( Minerva::Particle::Proton );
            hypotheses.push_back( Minerva::Particle::Pion );
            toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
        }


        // Make particles
        Minerva::Prong* prong = prongs[p];
        bool isParticleCreated = m_particleMaker->makeParticles(prong,hypotheses,toolsToUse);

        if( isParticleCreated ) {
            debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() 
                << " has " << prong->particles().size() << " number of particles attached." << endmsg; 
            makeParticles = true;

            if (prong->particles().size() == 2){
                Gaudi::LorentzVector P4_0 = prong->particles()[0]->momentumVec();
                Gaudi::LorentzVector P4_1 = prong->particles()[1]->momentumVec();
            }
        } else{
            debug() << "Did not make particles for the prong type = " << prongs[p]->typeBitsToString() << endmsg;
        }
        hypotheses.clear();
    }

    return makeParticles;
}

/*bool CCProtonPi0::createTrackedParticles(Minerva::PhysicsEvent *event ) const
{
  ProngVect prongs = event->primaryProngs();
  ProngVect HadronProngs;
  bool makeParticles = false; 
  double SENTINEL = -999.9;

  //--------------------------------------------
  // Likelihood PID Tool
  //--------------------------------------------
  Minerva::ParticleVect pionLikelihood;
  std::vector<Minerva::Particle::ID> pionHypotheses;
  pionHypotheses.push_back(Minerva::Particle::Pion);

std::vector<double> HadronScoreLLR(10,SENTINEL);

double Hadron_Visible_Energy = 0.0;

    // Check if the prongs are odMatch
    Minerva::EventMgr* mgr = getEventMgr(prongs[0]);
    Minerva::ODClusterVect odClusters = mgr->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
    m_odMatchTool->classifyProngs(prongs,odClusters);

    // Loop over prongs
    for(unsigned int p = 0; p < prongs.size(); p++) {
        debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() << endmsg;

        // Make sure the prong is not a bad object
        bool pass = true; std::string tag = "BadObject";
        if( prongs[p]->filtertaglist()->filterTagExists(tag) ) {
            if( prongs[p]->filtertaglist()->checkFilterTag(tag,pass) ) {
                error() << "This prong = " << prongs[p]->typeBitsToString() << " has been flagged as a \"BadObject\", skipping!" << endmsg;
                continue;
            }
        }

        // Skip prongs with particles
        debug()<< "The prong has n particles = " << prongs[p]->particles().size() << endmsg;
        if( !prongs[p]->particles().empty() ){
            debug()<<"Skipping Prong: Prong has attached particles!"<<endmsg;
            continue;
        }

        // Select the particle hypotheses candidates and particle tools
        //std::vector< Minerva::Particle::ID > hypotheses;
        //IParticleMakerTool::NameAliasListType toolsToUse;

        // Do NOT use prongs which have OD Match
          //hypotheses.push_back( Minerva::Particle::Proton );
            //hypotheses.push_back( Minerva::Particle::Pion );
            //toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
         debug() << "The prong reach the OD, skipping" << endmsg;
      continue;
        }
       
  // Skip the Muon Prong
    if ( prongs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
      debug() << "The prong is a PrimaryMuon, skipping" << endmsg;
      continue;
    }

 // Check if the prong has associated a track, otherwise it could return an Unknown particle
    Minerva::TrackVect tracks = prongs[p]->minervaTracks();
    debug() << " The Prong has " << tracks.size() << " tracks" << endmsg;
    if ( tracks.empty() ){
    	debug() << " The Prong hasn't a track, skipping the prong" << endmsg;
    	continue;
    }

    // Check if the prong is not a muon (MINOS track or stub)
    if ( prongs[p]->MinosTrack() || prongs[p]->MinosStub() ){
    	debug() << " This is a MINOS matched prong, skipping it" << endmsg;
    	continue;
    }

    // Check if the prong is forked
    if ( prongs[p]->Forked() ){
      debug() << "  This prong is forked, skipping it" << endmsg;
      continue;
    }

    // Count nodes
    int nodes = 0;
    Minerva::TrackVect::iterator it_tracks;
    for (it_tracks = tracks.begin(); it_tracks != tracks.end(); ++it_tracks){
      nodes += (*it_tracks)->nNodes();
    }

    if ( nodes < 3){
      debug() << " This prong has less than 3 nodes, skipping it" << endmsg;
      continue;
    }

    // Check for the direction of the prongs.
    Minerva::Track::Direction direction = prongs[p]->minervaTracks().back()->direction();
    bool isODMatchAtTrackEnd = false;
    if ( (direction == Minerva::Track::Forward && prongs[p]->OdMatchDS() ) || (direction == Minerva::Track::Backward && prongs[p]->OdMatchUS() ) ){
      isODMatchAtTrackEnd = true;
    }

    if (prongs[p]->ExitID() && !isODMatchAtTrackEnd){
      debug() << "The prong is exiting the ID but don't match with the OD, skipping it" << endmsg;
      continue;
    }

    if ( isODMatchAtTrackEnd ){
      debug() << "The prong is matching the OD at the end of the track, skipping it" << endmsg;
      continue;
    }

    HadronProngs.push_back(prongs[p]);
    Hadron_Visible_Energy += prongs[p]->minervaVisibleEnergySum();
       

// Make particles
        Minerva::Prong* prong = prongs[p];
        bool isParticleCreated = m_particleMaker->makeParticles(prong,hypotheses,toolsToUse);

        if( isParticleCreated ) {
            debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() 
                << " has " << prong->particles().size() << " number of particles attached." << endmsg; 
            makeParticles = true;

            if (prong->particles().size() == 2){
                Gaudi::LorentzVector P4_0 = prong->particles()[0]->momentumVec();
                Gaudi::LorentzVector P4_1 = prong->particles()[1]->momentumVec();
            }
        } else{
            debug() << "Did not make particles for the prong type = " << prongs[p]->typeBitsToString() << endmsg;
        }
        hypotheses.clear();
    }
 
debug() << "Making particles" << endmsg;
    prongs[p]->filtertaglist()->setOrAddFilterTag("PrimaryHadron", true);
    Minerva::Prong* prong = prongs[p];

    bool isLikelihoodPion = m_LikelihoodPIDTool->makeParticles(prong, pionLikelihood, pionHypotheses);

    if ( isLikelihoodPion){
      makeParticles = true;
      m_HadronProngs.push_back(prong);
      //m_HadronParticles.push_back(pionLikelihood[0]);   // I'm only attaching the first particle.
      debug() << "The prong has associated " << pionLikelihood.size() << " Likelihood particles" << endmsg;

      int index = 0;
      for(unsigned int i = 0; i < pionLikelihood.size(); i++ ){
	debug() << "Pion Likelihood score = " << pionLikelihood[i]->score() << endmsg;
	//m_HadronParticles.push_back(pionLikelihood[i]);
	HadronScoreLLR[i] = pionLikelihood[i]->score();    // this is a vector for all scores in the prong
	// check for the highest score and to be bigger than 0, to be more likely to be a pion
	double max_score = 0.0;
	//int index = 0;
	if (HadronScoreLLR[i] >= max_score){
	  max_score = HadronScoreLLR[i];
	  index = i;
          // NOT SURE IF IT IS GOING TO WORK DAMED IT!!!!
	}
	//m_HadronParticles.push_back(pionLikelihood[index]);

      } // End of loop over particles associated to a prong
      m_HadronParticles.push_back(pionLikelihood[index]);
    }

  } // End of the loop over prongs

  event->setContainerDoubleData("all_prongs_LLRScore", HadronScoreLLR);
  event->setIntData("N_Primary_Hadrons", HadronProngs.size());
  event->setDoubleData("hadron_visible_energy", Hadron_Visible_Energy);

    return makeParticles;
}

//=======================================================================================================
// Get Pion Prong
//     It selects the prongs from m_HadronProngs vector, and select thoose that has
//     Likelihood score >= 0.
//=======================================================================================================
bool CCProtonPi0::getPionProngs(Minerva::PhysicsEvent *event) const
{
  debug() << " Entering to getPionProngs " << endmsg;

  //-- Check for the size of the prongs
  if (m_HadronProngs.size() == 0 ){
    debug() << "The Hadron Prong size is 0" << endmsg;
    return false;
  }

  bool isPionParticle = false;
  //-- Loop over the Hadron Prongs
  for (unsigned int i = 0; i < m_HadronProngs.size(); i++){

    SmartRef<Minerva::Prong> prong = m_HadronProngs[i];
    SmartRef<Minerva::Particle> particle = m_HadronParticles[i];

    if (particle->score() >= 0){
      debug() << "A GOOD pion candidate, the likelihood score is = " << particle->score() << endmsg;
      isPionParticle = true;
      m_PionProngs.push_back(prong);
      m_PionParticles.push_back(particle);
    }

  } // End of the loop over Hadron Prongs

  return isPionParticle;
}
*/
//==============================================================================
// Return the momentum analyzable contained ( proton candidate ) prong/particle
// Uses Global Variables: m_ProtonProngs and m_ProtonParticles
//==============================================================================
bool CCProtonPi0::getProtonProng(  Minerva::PhysicsEvent *event) const
{
    ProngVect primaryProngs = event->primaryProngs();

    debug() <<"N(primaryProngs) =  " << primaryProngs.size() << endmsg;

    // Initialize
    Minerva::ProngVect tmpProngs;
    bool isProtonExist = false;
    bool pass = true;
    std::string tag = "PrimaryProton";

    // Get All Proton Candidates
    for(unsigned int p = 0; p < primaryProngs.size(); p++) {

        debug() <<"Checking prong "<<p<<endmsg;

        // Skip Muon Prong
        if( primaryProngs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
            debug() <<"Muon Prong, skipping this prong! "<< endmsg;
            continue;
        }

        // Temp Storage for ProngVect, Single Prong
        SmartRef<Minerva::Prong> prong       = (Minerva::Prong*)NULL;
        SmartRef<Minerva::Particle> particle = (Minerva::Particle*)NULL;

        // Push current Prong to temp ProngVect
        tmpProngs.push_back( primaryProngs[p] );

        // Find Proton using m_protonUtils
        bool isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle); 

        tmpProngs.clear();

        if( isProton ) {
            debug() <<"Found a proton prong!"<< endmsg;

            Gaudi::LorentzVector P4 = particle->momentumVec(); 
            debug()<<"Proton Particle Score: " << particle->score() << endmsg;
            debug()<<"Proton P4 = ( "<<P4.px()<<", "<<P4.py()<<", "<<P4.pz()<<", "<<P4.E()<<" )"<<endmsg;

            prong->filtertaglist()->setOrAddFilterTag(tag,pass);
            particle->filtertaglist()->setOrAddFilterTag(tag,pass);
            prong->updateBestParticle(particle);
            m_hitTagger->applyColorTag(prong, m_Color_protonProng);
            m_ProtonProngs.push_back( prong );
            m_ProtonParticles.push_back( particle );
            isProtonExist = true;
        }

    }

    return isProtonExist;
}


//==============================================================================
// Set proton particle data
//==============================================================================
bool CCProtonPi0::setProtonData( Minerva::PhysicsEvent *event ) const 
{
    if ( m_ProtonParticles.size() == 0 ) {
        warning()<< "m_ProtonParticles is empty! Exiting..." <<endmsg;
        return false;
    }

    // Sanity Check
    if ( m_ProtonProngs.size() != m_ProtonParticles.size()){
        warning()<< "m_ProtonProngs and m_ProtonParticles NOT SAME SIZE" <<endmsg;
        return false;
    }

    const double SENTINEL = -9.9;
    // PID Score from Likelihood Parameters
    Minerva::ParticleVect protonLikelihood;
    Minerva::ParticleVect pionLikelihood;
    std::vector<Minerva::Particle::ID> protonHypotheses;
    std::vector<Minerva::Particle::ID> pionHypotheses;
    protonHypotheses.push_back( Minerva::Particle::Proton );
    pionHypotheses.push_back( Minerva::Particle::Pion );

    // Declare Variables
    int leadingProtonIndice = 0;
    double tempProtonMomentum = 0.0;
    std::vector<double> p_calCorrection(10,SENTINEL);
    std::vector<double> p_visEnergyCorrection(10,SENTINEL);
    std::vector<double> p_dedx(10,SENTINEL);

    std::vector<double> proton_theta(10,SENTINEL);
    std::vector<double> proton_theta_beam(10,SENTINEL);
    std::vector<double> proton_phi(10,SENTINEL);
    std::vector<double> proton_phi_beam(10,SENTINEL);

    std::vector<double> E(10,SENTINEL);
    std::vector<double> px(10,SENTINEL);
    std::vector<double> py(10,SENTINEL);
    std::vector<double> pz(10,SENTINEL);
    std::vector<double> p(10,SENTINEL);

    std::vector<double> ekin(10,SENTINEL);

    std::vector<double> proton_end_x(10,SENTINEL);
    std::vector<double> proton_end_y(10,SENTINEL);
    std::vector<double> proton_end_z(10,SENTINEL);

    std::vector<double> proton_start_x(10,SENTINEL);
    std::vector<double> proton_start_y(10,SENTINEL);
    std::vector<double> proton_start_z(10,SENTINEL);
    std::vector<double> length(10,SENTINEL);

    std::vector<int> kinked(10,-1);
    std::vector<int> odMatch(10,-1);

    std::vector<double> pionScore(10,SENTINEL);
    std::vector<double> protonScore(10,SENTINEL);
    std::vector<double> protonScoreLLR(10,SENTINEL);
    std::vector<double> chi2(10,SENTINEL);

    // Get Vertex Z Position
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   
    double vertexZ = vertex->position().z();

    std::map< std::string, std::vector<double> > dedxMapEnergy;
    std::map< std::string, std::vector<double> > dedxMapScore;

    for(unsigned int i = 0; i < m_dedx_uncertainties.size(); i++) {
        std::string name = m_dedx_uncertainties[i];
        std::vector<double> tmp(10,-1);
        dedxMapEnergy.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
        dedxMapScore.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
    }

    // Loop over all Proton Candidates
    for(unsigned int i = 0; i < m_ProtonProngs.size(); i++) {
        SmartRef<Minerva::Prong> prong = m_ProtonProngs[i];
        SmartRef<Minerva::Particle> particle = m_ProtonParticles[i];

        double theta = prong->minervaTracks().front()->theta();
        double phi   = prong->minervaTracks().front()->phi();
        double E_original = particle->momentumVec().E();
        
        proton_theta[i] = theta;
        proton_phi[i] = phi;
        proton_theta_beam[i]  = m_coordSysTool->thetaWRTBeam(particle->momentumVec());
        proton_phi_beam[i]    = m_coordSysTool->phiWRTBeam(particle->momentumVec());

        Gaudi::LorentzVector dEdXprotonfourVec;
        StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(prong,particle,vertexZ,dEdXprotonfourVec);
        if( sc ) particle->setMomentumVec(dEdXprotonfourVec);
        else dEdXprotonfourVec = particle->momentumVec();

        p_dedx[i] = sqrt( std::pow(dEdXprotonfourVec.E(),2) - std::pow(MinervaUnits::M_p,2) );
        correctProtonProngEnergy(prong,p_calCorrection[i],p_visEnergyCorrection[i]);

        E[i]  = sqrt( std::pow(p_calCorrection[i],2) + std::pow(MinervaUnits::M_p,2) );
        px[i] = p_calCorrection[i]*sin(theta)*cos(phi);
        py[i] = p_calCorrection[i]*sin(theta)*sin(phi);
        pz[i] = p_calCorrection[i]*cos(theta);
        p[i]  = sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );

        if( prong->minervaTracks().front()->direction() == Minerva::Track::Backward && pz[i] > 0. ) pz[i] *= -1.0;
        Gaudi::LorentzVector protonfourVec(px[i],py[i],pz[i],E[i]);

        ekin[i] = protonfourVec.E() - MinervaUnits::M_p;

        proton_end_x[i]   = (prong->minervaTracks().back())->lastState().x();
        proton_end_y[i]   = (prong->minervaTracks().back())->lastState().y();
        proton_end_z[i]   = (prong->minervaTracks().back())->lastState().z();

        proton_start_x[i] = (prong->minervaTracks().back())->firstState().x();
        proton_start_y[i] = (prong->minervaTracks().back())->firstState().y();
        proton_start_z[i] = (prong->minervaTracks().back())->firstState().z();

        length[i] = calcDistance(   proton_end_x[i],proton_end_y[i], proton_end_z[i],
                proton_start_x[i],proton_start_y[i],proton_start_z[i]);
        kinked[i]  = (int)prong->Kinked();
        odMatch[i] = (int)prong->OdMatch(); 
        chi2[i]    = particle->getDoubleData("chi2_ndf");
        m_protonUtils->getParticleScores(prong,protonScore[i],pionScore[i]);

        // --------------------------------------------------------------------
        // dEdX Uncertainties
        // --------------------------------------------------------------------
        for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
            std::string name = m_dedx_uncertainties[j];

            double score_varied = particle->hasDoubleData("score1_" + name) ? particle->getDoubleData("score1_" + name) : 0.0;
            double E_varied     = particle->hasDoubleData("FitE_"   + name) ? particle->getDoubleData("FitE_"   + name) : 0.0;

            double score_shift  = particle->getDoubleData("score1") - score_varied;
            double E_shift      = E_original - E_varied;

            (dedxMapEnergy.find(name))->second.at(i) = E_shift;
            (dedxMapScore.find(name))->second.at(i) = score_shift;
        }

        m_LikelihoodPIDTool->makeParticles( prong, protonLikelihood, protonHypotheses );
        m_LikelihoodPIDTool->makeParticles( prong, pionLikelihood, pionHypotheses );

        // Phil's Likelihood Method for PID
        if(!protonLikelihood.empty()){
            protonScoreLLR[i] = protonLikelihood[i]->score();
        }else{
            debug()<<"particlesLikelihood is empty!"<<endmsg;
        }

        // Find Interaction Proton -- Highest Momentum Proton
        if (p[i]  > tempProtonMomentum){
            tempProtonMomentum = p[i];
            leadingProtonIndice = i;
        }

        //----------------------------------------------------------------------
        // Debugging: Check values for All Protons
        debug()<<"Proton Candidate = "<<i<<endmsg;
        debug()<<"Original Energy = "<<E_original<<endmsg;
        debug()<<"Corrected Momentum = "<<p[i]<<endmsg;
        debug()<<"P4(Proton) = ( "<<px[i]<<", "<<py[i]<<", "<<pz[i]<<", "<<E[i]<<" )"<<endmsg;
        debug()<<"  protonScore = "<<protonScore[i]<<" pionScore = "<<pionScore[i]<<endmsg;
        debug()<<"  protonScoreLLR = "<<protonScoreLLR[i]<<endmsg;
        debug()<<"  length[i] = "<<length[i]<<endmsg;
        debug()<<"  ekin[i] = "<<ekin[i]<<endmsg;
        debug()<<"  kinked[i] = "<<kinked[i]<<endmsg;
        debug()<<"  odMatch[i] = "<<odMatch[i]<<endmsg;
        for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
            std::string name = m_dedx_uncertainties[j];

            debug()<<" proton_energy_shift_"<<name<<"[i] = "<<(dedxMapEnergy.find(name))->second.at(i)<<endmsg;
            debug()<<" proton_score1_shift_"<<name<<"[i] = "<<(dedxMapScore.find(name))->second.at(i)<<endmsg;
        }
        //----------------------------------------------------------------------

        // Check Momentum 
        if ( isnan(p[i]) ) return false;
    }


    //----------------------------------------------------------------------
    // Debugging: Check values for Leading Proton
    debug()<<"Proton Candidate = "<<leadingProtonIndice<<endmsg;
    debug()<<"Corrected Momentum = "<<p[leadingProtonIndice]<<endmsg;
    debug()<<"P4(Proton) = ( "<<px[leadingProtonIndice]<<", "<<py[leadingProtonIndice]<<", "<<pz[leadingProtonIndice]<<", "<<E[leadingProtonIndice]<<" )"<<endmsg;
    debug()<<"  protonScore = "<<protonScore[leadingProtonIndice]<<" pionScore = "<<pionScore[leadingProtonIndice]<<endmsg;
    debug()<<"  protonScoreLLR = "<<protonScoreLLR[leadingProtonIndice]<<endmsg;
    debug()<<"  length = "<<length[leadingProtonIndice]<<endmsg;
    debug()<<"  ekin = "<<ekin[leadingProtonIndice]<<endmsg;
    debug()<<"  kinked = "<<kinked[leadingProtonIndice]<<endmsg;
    debug()<<"  odMatch = "<<odMatch[leadingProtonIndice]<<endmsg;
    for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
        std::string name = m_dedx_uncertainties[j];

        debug()<<" proton_energy_shift_"<<name<<" = "<<(dedxMapEnergy.find(name))->second.at(leadingProtonIndice)<<endmsg;
        debug()<<" proton_score1_shift_"<<name<<" = "<<(dedxMapScore.find(name))->second.at(leadingProtonIndice)<<endmsg;
    }
    //----------------------------------------------------------------------


    // Set Leading Proton 4-Momentum
    m_proton_4P.SetPxPyPzE( px[leadingProtonIndice], py[leadingProtonIndice], pz[leadingProtonIndice],E[leadingProtonIndice]);

    event->setContainerDoubleData("all_protons_p_calCorrection",p_calCorrection);
    event->setContainerDoubleData("all_protons_p_visEnergy",p_visEnergyCorrection);
    event->setContainerDoubleData("all_protons_p_dEdXTool",p_dedx);
    event->setContainerDoubleData("all_protons_endPointX",proton_end_x);
    event->setContainerDoubleData("all_protons_endPointY",proton_end_y);
    event->setContainerDoubleData("all_protons_endPointZ",proton_end_z);
    event->setContainerDoubleData("all_protons_startPointX",proton_start_x);
    event->setContainerDoubleData("all_protons_startPointY",proton_start_y);
    event->setContainerDoubleData("all_protons_startPointZ",proton_start_z);
    event->setContainerDoubleData("all_protons_length",length);
    event->setContainerDoubleData("all_protons_px",px);
    event->setContainerDoubleData("all_protons_py",py);
    event->setContainerDoubleData("all_protons_pz",pz);
    event->setContainerDoubleData("all_protons_E",E);
    event->setContainerDoubleData("all_protons_P",p);
    event->setContainerDoubleData("all_protons_KE",ekin);
    event->setContainerIntData("all_protons_kinked",kinked);
    event->setContainerIntData("all_protons_odMatch",odMatch);
    event->setContainerDoubleData("all_protons_protonScore",protonScore);
    event->setContainerDoubleData("all_protons_pionScore",pionScore);
    event->setContainerDoubleData("all_protons_LLRScore",protonScoreLLR);
    event->setContainerDoubleData("all_protons_chi2_ndf",chi2);
    event->setContainerDoubleData("all_protons_theta",proton_theta);
    event->setContainerDoubleData("all_protons_theta_beam",proton_theta_beam);
    event->setContainerDoubleData("all_protons_phi",proton_phi);
    event->setContainerDoubleData("all_protons_phi_beam",proton_phi_beam);
    for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
        std::string name = m_dedx_uncertainties[j];

        std::vector<double> E_shift_vect     = dedxMapEnergy.find(name)->second;
        std::vector<double> score_shift_vect = dedxMapScore.find(name)->second;

        event->setContainerDoubleData("all_protons_energy_shift_" + name, E_shift_vect);
        event->setContainerDoubleData("all_protons_score1_shift_" + name, score_shift_vect);
    }

    event->setDoubleData("proton_px",px[leadingProtonIndice]);
    event->setDoubleData("proton_py",py[leadingProtonIndice]);
    event->setDoubleData("proton_pz",pz[leadingProtonIndice]);
    event->setDoubleData("proton_E",E[leadingProtonIndice]);
    event->setDoubleData("proton_P",p[leadingProtonIndice]);
    event->setDoubleData("proton_KE",ekin[leadingProtonIndice]);
    event->setDoubleData("proton_theta",proton_theta[leadingProtonIndice]);
    event->setDoubleData("proton_theta_beam",proton_theta_beam[leadingProtonIndice]);
    event->setDoubleData("proton_phi",proton_phi[leadingProtonIndice]);
    event->setDoubleData("proton_phi_beam",proton_phi_beam[leadingProtonIndice]);
    event->setDoubleData("proton_length", length[leadingProtonIndice]);
    event->setDoubleData("proton_protonScore", protonScore[leadingProtonIndice]);
    event->setDoubleData("proton_pionScore", pionScore[leadingProtonIndice]);
    event->setDoubleData("proton_LLRScore", protonScoreLLR[leadingProtonIndice]);
    event->setIntData("proton_kinked",kinked[leadingProtonIndice]);
    event->setIntData("proton_leadingIndice",leadingProtonIndice);
    for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
        std::string name = m_dedx_uncertainties[j];

        std::vector<double> E_shift_vect     = dedxMapEnergy.find(name)->second;
        std::vector<double> score_shift_vect = dedxMapScore.find(name)->second;

        event->setDoubleData("proton_energy_shift_" + name, E_shift_vect[leadingProtonIndice]);
        event->setDoubleData("proton_score1_shift_" + name, score_shift_vect[leadingProtonIndice]);
    }

    return true;
}
//==============================================================================
// Correct Pion Energy
//==============================================================================
/*void CCProtonPi0::correctPionProngEnergy(SmartRef<Minerva::Prong>& pionProng, double& p_calCorrection, double& p_visEnergyCorrection) const
{
    //-- initialize
    double calEnergy = 0.0;
    double visEnergy = 0.0;

    //-- get the pion prong subprongs
    Minerva::ProngVect subProngs = pionProng->subProngs();
    debug() << "Pion nSubProngs" << subProngs.size() << endmsg;
    if ( !subProngs.empty() ){
        debug() << "The pion subProngs are = " << subProngs.size() << endmsg;

        //-- loop over sub prongs
        for (unsigned int p = 0; p < subProngs.size(); p++){
            //-- get the subprong link vertex
            Minerva::Vertex* vertex      = (Minerva::Vertex*)NULL;
            Minerva::VertexVect vertices = subProngs[p]->linkVertices();

            //--skip subprongs attached to secondar vertex -- add back in energy at the dEdX stage
            for (unsigned int v = 0; v < vertices.size(); ++v){
                if (vertices[v]->type() == Minerva::Vertex::Kinked) vertex = vertices[v];
            }
            if (vertex) continue;

            //-- apply polyline calorimetric energy correction to subProngs
            calEnergy += CaloUtils->applyCalConsts(subProngs[p], "Default", true, true);
            debug() << " The calorimetric energy = " << calEnergy << " MeV" << endmsg;

            //-- get the total visible energy correction to subProngs;
            Minerva::IDClusterVect clusters = subProngs[p]->getAllIDClusters();
            for (unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
            debug() << " The pion visible energy = " << visEnergy << " MeV" << endmsg;
        }// end of subprongs
    } // end if !subProngs.empty()

    //-- get the pion prong fuzz energy
    Minerva::IDBlobVect fuzzBlobs = pionProng->getAllIDBlobs();
    for(unsigned int blob = 0; blob < fuzzBlobs.size(); blob++) {
        if( !fuzzBlobs[blob]->hasIntData("fuzzAttachedTrack") ) continue;

        //-- apply polyline calorimetric energy correction to the fuzz energy
        calEnergy += CaloUtils->applyCalConsts(fuzzBlobs[blob], "Default", true);

        //-- get the total visible energy of the fuzz blobs
        Minerva::IDClusterVect clusters = fuzzBlobs[blob]->clusters();
        for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
    }// end of loop over blobs in the pion prong

    //-- get the prong first track
    SmartRef<Minerva::Track> track = (*pionProng->minervaTracks().begin());

    //-- update the particles four momentum
    Minerva::ParticleVect particles = pionProng->particles();
    for(unsigned int part = 0; part < particles.size(); part++) {
        if( particles[part]->idcode() != Minerva::Particle::Pion ) continue;

        Gaudi::XYZTVector fourMomentum = particles[part]->momentumVec();
        debug() << "   particle idcode = " << particles[part]->idcode() << ", and four momomentum = " << fourMomentum << endmsg;

        double E_calCorrection = fourMomentum.E() + calEnergy;
        p_calCorrection = sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << "  update energy using calorimetric correction = " << E_calCorrection << endmsg;

        double E_visEnergyCorrection = fourMomentum.E() + visEnergy;
        p_visEnergyCorrection = sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << " update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
    }

    return;
}
*/

//==============================================================================
// Correct Proton Energy
//==============================================================================
void CCProtonPi0::correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
        double& p_calCorrection, 
        double& p_visEnergyCorrection ) const
{
    //-- initialize
    double calEnergy = 0.0;
    double visEnergy = 0.0;

    //-- get the proton prong subprongs
    Minerva::ProngVect subProngs = protonProng->subProngs();
    debug()<<"nSubProngs = "<<subProngs.size()<<endmsg; 
    if( !subProngs.empty() ) {
        debug() << "  The number of subProngs = " << subProngs.size() << endmsg;

        //-- loop over sub prongs
        for(unsigned int p = 0; p < subProngs.size(); p++) {

            //-- get the subprong link vertex
            Minerva::Vertex* vertex      = (Minerva::Vertex*)NULL;
            Minerva::VertexVect vertices = subProngs[p]->linkVertices();

            //-- skip subprongs attached to secondary vertex -- add back in energy at the dEdX stage 
            for(unsigned int v = 0; v < vertices.size(); ++v) {
                if( vertices[v]->type() == Minerva::Vertex::Kinked ) vertex = vertices[v];
            }
            if( vertex ) continue;

            //-- apply polyline calorimetric energy correction to the subProngs
            calEnergy += CaloUtils->applyCalConsts(subProngs[p]);
            debug() << "  The calorimetric energy = " << calEnergy << " MeV" << endmsg;

            //-- get the total visible energy of the sub prongs
            Minerva::IDClusterVect clusters = subProngs[p]->getAllIDClusters();
            for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
            debug() << "  The visible energy = " << visEnergy << " MeV" << endmsg;
        }  
    }    

    //-- get the proton prong fuzz energy
    Minerva::IDBlobVect fuzzBlobs = protonProng->getAllIDBlobs();
    for(unsigned int blob = 0; blob < fuzzBlobs.size(); blob++) {
        if( !fuzzBlobs[blob]->hasIntData("fuzzAttachedTrack") ) continue;

        //-- apply polyline calorimetric energy correction to the fuzz energy  
        calEnergy += CaloUtils->applyCalConsts(fuzzBlobs[blob]);

        //-- get the total visible energy of the fuzz blobs
        Minerva::IDClusterVect clusters = fuzzBlobs[blob]->clusters();
        for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
    }

    //-- get the prong first track
    SmartRef<Minerva::Track> track = (*protonProng->minervaTracks().begin());

    //-- update the particles four momentum
    Minerva::ParticleVect particles = protonProng->particles();
    for(unsigned int part = 0; part < particles.size(); part++) {
        if( particles[part]->idcode() != Minerva::Particle::Proton ) continue;

        Gaudi::XYZTVector fourMomentum = particles[part]->momentumVec();
        debug() << "   particle idcode = " << particles[part]->idcode() << ", and four momomentum = " << fourMomentum << endmsg;

        double E_calCorrection = fourMomentum.E() + calEnergy;
        p_calCorrection = sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << "  update energy using calorimetric correction = " << E_calCorrection << endmsg;

        double E_visEnergyCorrection = fourMomentum.E() + visEnergy;
        p_visEnergyCorrection = sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );     
        debug() << " update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
    }     

    return;
}

//==============================================================================
// setPi0Data
//==============================================================================
bool CCProtonPi0::setPi0Data( Minerva::PhysicsEvent *event ) const 
{    
    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<<"Passed NULL IDBlob"<<endmsg;
        return false;
    }

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    if (m_ApplyAttenuationCorrection) {
        ApplyAttenuationCorrection(m_Pi0Blob1);
        ApplyAttenuationCorrection(m_Pi0Blob2);
    }

    // ------------------------------------------------------------------------
    // Get Shower Energy Using New Method
    const double SENTINEL = -9.9;
    std::vector<double> gamma1_evis_v(5,SENTINEL); 
    std::vector<double> gamma1_energy_v(5,SENTINEL); 
    std::vector<double> gamma2_evis_v(5,SENTINEL); 
    std::vector<double> gamma2_energy_v(5,SENTINEL); 
    double g1energy = m_idHoughBlob->getBlobEnergyTime_New(m_Pi0Blob1, gamma1_evis_v, gamma1_energy_v);
    double g2energy = m_idHoughBlob->getBlobEnergyTime_New(m_Pi0Blob2, gamma2_evis_v, gamma2_energy_v);
    // ------------------------------------------------------------------------

    //********* Here I should multiply the by the factor of the energy resolution y*g1energy *****************//


    // ------------------------------------------------------------------------
    // Get Energy Using Old Method
    double g1energy_Old = 0.0;
    double g2energy_Old = 0.0;
    double dummy_trkr = 0.0;
    double dummy_ecal = 0.0;
    double dummy_hcal = 0.0;
    double dummy_scal = 0.0;
    m_idHoughBlob->getBlobEnergyTime_Old( m_Pi0Blob1, g1energy_Old, dummy_trkr, dummy_ecal, dummy_hcal, dummy_scal );
    m_idHoughBlob->getBlobEnergyTime_Old( m_Pi0Blob2, g2energy_Old, dummy_trkr, dummy_ecal, dummy_hcal, dummy_scal );
    // ------------------------------------------------------------------------


    // Make sure Gamma1 is the more energetic one
    if (g2energy > g1energy) {
        warning()<<" Gamma2 Energy is Higher than Gamma1"<<endmsg;
    }

    debug()<<" ------------------------------------ "<<endmsg;
    debug()<<"Gamma 1 Energy:"<<endmsg;
    debug()<<"  New = "<<g1energy<<endmsg;
    debug()<<"  Old = "<<g1energy_Old<<endmsg;

    debug()<<"Gamma 2 Energy:"<<endmsg;
    debug()<<"  New = "<<g2energy<<endmsg;
    debug()<<"  Old = "<<g2energy_Old<<endmsg;
    debug()<<" ------------------------------------ "<<endmsg;

    // Get Blob Time
    double time1 = m_Pi0Blob1->time();
    double time2 = m_Pi0Blob2->time();

    // Get Blob dEdx
    double dEdx1 = 0.0;
    double dEdx2 = 0.0;
    m_idHoughBlob->idBlobdEdx( m_Pi0Blob1, dEdx1 );
    m_idHoughBlob->idBlobdEdx( m_Pi0Blob2, dEdx2 );

    // Get Blob Directions
    Gaudi::XYZVector direction1 = m_Pi0Blob1->direction();
    Gaudi::XYZVector direction2	= m_Pi0Blob2->direction();

    double Blob1_cone_angle = CalculateConeAngle(m_Pi0Blob1, vtx_position); // angle in degrees
    double Blob2_cone_angle = CalculateConeAngle(m_Pi0Blob2, vtx_position);
    double Blob1_MeanE_1stHalf = MeanEnergyClusters_1stHalf(m_Pi0Blob1, vtx_position);
    double Blob2_MeanE_1stHalf = MeanEnergyClusters_1stHalf(m_Pi0Blob2, vtx_position);
    double Blob1_MeanE_2ndHalf = MeanEnergyClusters_2ndHalf(m_Pi0Blob1, vtx_position);
    double Blob2_MeanE_2ndHalf = MeanEnergyClusters_2ndHalf(m_Pi0Blob2, vtx_position);

    /*
    double Blob1_MedianE_1stHalf = MedianEnergyClusters_1stHalf(m_Pi0Blob1, vtx_position);
    double Blob2_MedianE_1stHalf = MedianEnergyClusters_1stHalf(m_Pi0Blob2, vtx_position);
    double Blob1_MedianE_2ndHalf = MedianEnergyClusters_2ndHalf(m_Pi0Blob1, vtx_position);
    double Blob2_MedianE_2ndHalf = MedianEnergyClusters_2ndHalf(m_Pi0Blob2, vtx_position);
    */

    //-------- Obtain de the energy and distace of the clusters in the blobs
    std::vector<std::pair<double, double> > blob1_clusters_distance_energy;
    std::vector<std::pair<double, double> > blob2_clusters_distance_energy;
    blob1_clusters_distance_energy = EnergyDistanceClusters(m_Pi0Blob1, vtx_position);
    blob2_clusters_distance_energy = EnergyDistanceClusters(m_Pi0Blob2, vtx_position);

    std::vector<double> blob1_clusters_distance;
    std::vector<double> blob1_clusters_energy;
    for (int i = 0; i < blob1_clusters_distance_energy.size(); i++){
      double distance = blob1_clusters_distance_energy[i].first;
      double energy = blob1_clusters_distance_energy[i].second;
      blob1_clusters_distance.push_back(distance);
      blob1_clusters_energy.push_back(energy);
    }

    std::vector<double> blob2_clusters_distance;
    std::vector<double> blob2_clusters_energy;
    for (int i = 0; i < blob2_clusters_distance_energy.size(); i++){
      double distance = blob2_clusters_distance_energy[i].first;
      double energy = blob2_clusters_distance_energy[i].second;
      blob2_clusters_distance.push_back(distance);
      blob2_clusters_energy.push_back(energy);
    }


    //----- Fit the Blob to a straight line between the centroid and the interaction vertex in the X-view
    double Blob1_pvalue = CalProbFitLine(m_Pi0Blob1, vtx_position);
    double Blob2_pvalue = CalProbFitLine(m_Pi0Blob2, vtx_position);
    info()<< "\tFit to a Line "<< "gamma1_pvalue = "<<Blob1_pvalue <<endmsg;
    info()<< "\tFit to a Line" << "gamma2_pvalue = "<<Blob2_pvalue <<endmsg;

    // Save Gamma1 and Gamma2 Starting Position (blob vertex)
    std::vector<double> vertex1 ;
    std::vector<double> vertex2;
    vertex1.push_back(m_Pi0Blob1->startPoint().x());
    vertex1.push_back(m_Pi0Blob1->startPoint().y());
    vertex1.push_back(m_Pi0Blob1->startPoint().z());

    vertex2.push_back(m_Pi0Blob2->startPoint().x());
    vertex2.push_back(m_Pi0Blob2->startPoint().y());
    vertex2.push_back(m_Pi0Blob2->startPoint().z());

    // Calculate Distance between Gamma Vertex and Event Vertex
    double gamma1_dist_vtx = calcDistance(  vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
            vertex1[0], vertex1[1], vertex1[2]);

    double gamma2_dist_vtx = calcDistance(  vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
            vertex2[0], vertex2[1], vertex2[2]);

    // Save Gamma1 and Gamma2 Direction
    std::vector<double> direc_1;
    std::vector<double> direc_2;
    direc_1.push_back(direction1.x());
    direc_1.push_back(direction1.y());
    direc_1.push_back(direction1.z());

    direc_2.push_back(direction2.x());
    direc_2.push_back(direction2.y());
    direc_2.push_back(direction2.z());

    // Calculate Gamma1 and Gamma2 3-Momentums
    TVector3 g1mom(direction1.x(),direction1.y(),direction1.z());
    TVector3 g2mom(direction2.x(),direction2.y(),direction2.z());
    g1mom *= g1energy;
    g2mom *= g2energy;

    // Get Pi0 Kinematics
    double M_pi0 = MinervaUnits::M_pi0;
    double pi0_E_Cal = g1energy + g2energy;
    TVector3 pimom = g1mom + g2mom;
    double pi0_P = pimom.Mag();
    double pi0_E = sqrt(pi0_P*pi0_P + M_pi0*M_pi0); 
    double pi0_KE = pi0_E - M_pi0;
    Gaudi::LorentzVector pi0_4P(pimom.Px(), pimom.Py(), pimom.Pz(), pi0_E);
    Gaudi::LorentzVector gamma1_4P(g1mom.Px(), g1mom.Py(), g1mom.Pz(), g1energy);
    Gaudi::LorentzVector gamma2_4P(g2mom.Px(), g2mom.Py(), g2mom.Pz(), g2energy);

    debug()<<"pi0_4P = "<<pi0_4P<<endmsg;
    debug()<<"gamma1_4P = "<<gamma1_4P<<endmsg;
    debug()<<"gamma2_4P = "<<gamma2_4P<<endmsg;

    // Estimate Shower End Point
    double gamma1_maxZ = FindShowerMaxZ(m_Pi0Blob1);
    double gamma2_maxZ = FindShowerMaxZ(m_Pi0Blob2);
    std::vector<double> end_vertex1 = EstimateShowerEndPoint(gamma1_4P, vertex1, gamma1_maxZ);
    std::vector<double> end_vertex2 = EstimateShowerEndPoint(gamma2_4P, vertex2, gamma2_maxZ);

    // Check for Michels in Shower Begin & End Points
    Minerva::Prong dummy_prong_1;
    Minerva::Prong dummy_prong_2;
    Minerva::Prong dummy_prong_3;
    Minerva::Prong dummy_prong_4;
    bool gamma1_isMichel_begin = m_michelTool_Default->findMichel( vertex1[0], vertex1[1], vertex1[2], time1, dummy_prong_1);
    bool gamma1_isMichel_end = m_michelTool_Default->findMichel( end_vertex1[0], end_vertex1[1], end_vertex1[2], time1, dummy_prong_2);

    bool gamma2_isMichel_begin = m_michelTool_Default->findMichel( vertex2[0], vertex2[1], vertex2[2], time2, dummy_prong_3);
    bool gamma2_isMichel_end = m_michelTool_Default->findMichel( end_vertex2[0], end_vertex2[1], end_vertex2[2], time2, dummy_prong_4);

    bool ImprovedMichel_Gamma1_Michel_Begin = ImprovedMichel_HasMichelAtPoint(vertex1[0], vertex1[1], vertex1[2]);
    bool ImprovedMichel_Gamma1_Michel_End = ImprovedMichel_HasMichelAtPoint(end_vertex1[0], end_vertex1[1], end_vertex1[2]);

    bool ImprovedMichel_Gamma2_Michel_Begin = ImprovedMichel_HasMichelAtPoint(vertex2[0], vertex2[1], vertex2[2]);
    bool ImprovedMichel_Gamma2_Michel_End = ImprovedMichel_HasMichelAtPoint(end_vertex2[0], end_vertex2[1], end_vertex2[2]);
    
    if (ImprovedMichel_Gamma1_Michel_Begin || ImprovedMichel_Gamma1_Michel_End){
        event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMatchedMichel", true);
    }
 
    if (ImprovedMichel_Gamma2_Michel_Begin || ImprovedMichel_Gamma2_Michel_End){
        event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMatchedMichel", true);
    }

    // Check Momentum
    if ( isnan(pi0_4P.P()) ) return false;

    // Get Angles wrt Beam Coordinates
    double pi0_theta_beam = m_coordSysTool->thetaWRTBeam(pi0_4P);
    double pi0_theta_beam_biasUp = m_coordSysTool->thetaWRTBeam(pi0_4P,m_beamAngleBias) - pi0_theta_beam;
    double pi0_theta_beam_biasDown = m_coordSysTool->thetaWRTBeam(pi0_4P,-1.0*m_beamAngleBias) - pi0_theta_beam;
    double pi0_phi_beam = m_coordSysTool->phiWRTBeam(pi0_4P);
    double gamma1_theta_beam = m_coordSysTool->thetaWRTBeam(gamma1_4P);
    double gamma1_phi_beam = m_coordSysTool->phiWRTBeam(gamma1_4P);
    double gamma2_theta_beam = m_coordSysTool->thetaWRTBeam(gamma2_4P);
    double gamma2_phi_beam = m_coordSysTool->phiWRTBeam(gamma2_4P);

    // Get Angles wrt Lab Coordinates
    double pi0_theta = pi0_4P.theta();
    double pi0_phi = pi0_4P.phi(); 
    double gamma1_theta = gamma1_4P.theta();
    double gamma1_phi = gamma1_4P.phi();
    double gamma2_theta = gamma2_4P.theta();
    double gamma2_phi = gamma2_4P.phi();

    // Calculate Opening Angle
    const double openingAngle       = (g1mom.Angle(g2mom))*TMath::RadToDeg();
    const double cos_openingAngle   = direction1.Dot(direction2);

    // Calculate invariant Mass of Pi0
    const double invMass = std::sqrt(2*g1energy*g2energy*(1-cos_openingAngle));
    const double invMass_Old = std::sqrt(2*g1energy_Old*g2energy_Old*(1-cos_openingAngle));

    // Set Pi0 4 - Momentum
    m_pi0_4P.SetPxPyPzE(pimom.x(),pimom.y(),pimom.z(),pi0_E);

    //--------------------------------------------------------------------------
    // Debugging: Check values
    debug()<<"m_pi0_4P = "<<m_pi0_4P<<endmsg;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //    Fill Branches
    //--------------------------------------------------------------------------
    // Pi0 Information
    event->setDoubleData("pi0_openingAngle", openingAngle);
    event->setDoubleData("pi0_cos_openingAngle", cos_openingAngle );
    event->setDoubleData("pi0_invMass", invMass);
    event->setDoubleData("pi0_invMass_Old", invMass_Old);
    event->setDoubleData("pi0_px" ,pimom.Px());
    event->setDoubleData("pi0_py", pimom.Py());
    event->setDoubleData("pi0_pz", pimom.Pz());
    event->setDoubleData("pi0_E", pi0_E);
    event->setDoubleData("pi0_E_Cal", pi0_E_Cal);
    event->setDoubleData("pi0_P", pi0_P);
    event->setDoubleData("pi0_KE", pi0_KE);
    event->setDoubleData("pi0_theta_beam", pi0_theta_beam);
    event->setDoubleData("pi0_theta_beam_biasUp", pi0_theta_beam_biasUp);
    event->setDoubleData("pi0_theta_beam_biasDown", pi0_theta_beam_biasDown);
    event->setDoubleData("pi0_phi_beam", pi0_phi_beam);
    event->setDoubleData("pi0_theta", pi0_theta);
    event->setDoubleData("pi0_phi", pi0_phi);
    event->setDoubleData("pi0_thetaX", std::atan2(pimom.X(),pimom.Z())*TMath::RadToDeg());
    event->setDoubleData("pi0_thetaY", std::atan2(pimom.Y(),pimom.Z())*TMath::RadToDeg());

    // Gamma1 Information
    event->filtertaglist()->setOrAddFilterTag( "gamma1_isMichel_begin", gamma1_isMichel_begin );
    event->filtertaglist()->setOrAddFilterTag( "gamma1_isMichel_end", gamma1_isMichel_end );
    event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_Gamma1HasMichel", ImprovedMichel_Gamma1_Michel_Begin || ImprovedMichel_Gamma1_Michel_End);
    event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_Gamma2HasMichel", ImprovedMichel_Gamma2_Michel_Begin || ImprovedMichel_Gamma2_Michel_End);
    event->setDoubleData("gamma1_px",g1mom.Px());
    event->setDoubleData("gamma1_py", g1mom.Py());
    event->setDoubleData("gamma1_pz", g1mom.Pz());
    event->setDoubleData("gamma1_E", g1energy);
    event->setDoubleData("gamma1_E_Old", g1energy_Old);
    event->setDoubleData("gamma1_P", g1mom.Mag());
    event->setDoubleData("gamma1_theta_beam", gamma1_theta_beam);
    event->setDoubleData("gamma1_phi_beam",  gamma1_phi_beam);
    event->setDoubleData("gamma1_theta", gamma1_theta);
    event->setDoubleData("gamma1_phi",  gamma1_phi);
    event->setDoubleData("gamma1_dEdx", dEdx1 );
    event->setDoubleData("gamma1_time", time1 );
    event->setDoubleData("gamma1_dist_vtx",gamma1_dist_vtx);
    event->setDoubleData("gamma1_maxZ",gamma1_maxZ);
    event->setContainerDoubleData("gamma1_direction", direc_1 );
    event->setContainerDoubleData("gamma1_vertex", vertex1 );
    event->setContainerDoubleData("gamma1_end_vertex", end_vertex1);
    event->setDoubleData("gamma1_evis_trkr", gamma1_evis_v[0]);
    event->setDoubleData("gamma1_evis_ecal", gamma1_evis_v[1]);
    event->setDoubleData("gamma1_evis_scal_X", gamma1_evis_v[2]);
    event->setDoubleData("gamma1_evis_scal_UV", gamma1_evis_v[3]);
    event->setDoubleData("gamma1_evis_hcal", gamma1_evis_v[4]);
    event->setDoubleData("gamma1_energy_trkr", gamma1_energy_v[0]);
    event->setDoubleData("gamma1_energy_ecal", gamma1_energy_v[1]);
    event->setDoubleData("gamma1_energy_scal_X", gamma1_energy_v[2]);
    event->setDoubleData("gamma1_energy_scal_UV", gamma1_energy_v[3]);
    event->setDoubleData("gamma1_energy_hcal", gamma1_energy_v[4]);
    event->setDoubleData("gamma1_cone_angle", Blob1_cone_angle);
    event->setDoubleData("gamma1_pvalue", Blob1_pvalue);
    event->setDoubleData("gamma1_MeanE_1stHalf", Blob1_MeanE_1stHalf);
    event->setDoubleData("gamma1_MeanE_2ndHalf", Blob1_MeanE_2ndHalf);
    //event->setDoubleData("gamma1_MedianE_1stHalf", Blob1_MedianE_1stHalf);
    //event->setDoubleData("gamma1_MedianE_2ndHalf", Blob1_MedianE_2ndHalf);
    event->setContainerDoubleData("gamma1_clusters_distance", blob1_clusters_distance);
    event->setContainerDoubleData("gamma1_clusters_energy",   blob1_clusters_energy);

    // Gamma2 Information
    event->filtertaglist()->setOrAddFilterTag( "gamma2_isMichel_begin", gamma2_isMichel_begin );
    event->filtertaglist()->setOrAddFilterTag( "gamma2_isMichel_end", gamma2_isMichel_end );
    event->setDoubleData("gamma2_px", g2mom.Px());
    event->setDoubleData("gamma2_py", g2mom.Py());
    event->setDoubleData("gamma2_pz", g2mom.Pz());
    event->setDoubleData("gamma2_E", g2energy);
    event->setDoubleData("gamma2_E_Old", g2energy_Old);
    event->setDoubleData("gamma2_P", g2mom.Mag());
    event->setDoubleData("gamma2_theta_beam", gamma2_theta_beam);
    event->setDoubleData("gamma2_phi_beam",  gamma2_phi_beam);
    event->setDoubleData("gamma2_theta", gamma2_theta);
    event->setDoubleData("gamma2_phi",  gamma2_phi);
    event->setDoubleData("gamma2_dEdx", dEdx2 );
    event->setDoubleData("gamma2_time", time2 );
    event->setDoubleData("gamma2_dist_vtx", gamma2_dist_vtx);
    event->setDoubleData("gamma2_maxZ",gamma2_maxZ);
    event->setContainerDoubleData("gamma2_direction", direc_2 );
    event->setContainerDoubleData("gamma2_vertex", vertex2 );
    event->setContainerDoubleData("gamma2_end_vertex", end_vertex2);
    event->setDoubleData("gamma2_evis_trkr", gamma2_evis_v[0]);
    event->setDoubleData("gamma2_evis_ecal", gamma2_evis_v[1]);
    event->setDoubleData("gamma2_evis_scal_X", gamma2_evis_v[2]);
    event->setDoubleData("gamma2_evis_scal_UV", gamma2_evis_v[3]);
    event->setDoubleData("gamma2_evis_hcal", gamma2_evis_v[4]);
    event->setDoubleData("gamma2_energy_trkr", gamma2_energy_v[0]);
    event->setDoubleData("gamma2_energy_ecal", gamma2_energy_v[1]);
    event->setDoubleData("gamma2_energy_scal_X", gamma2_energy_v[2]);
    event->setDoubleData("gamma2_energy_scal_UV", gamma2_energy_v[3]);
    event->setDoubleData("gamma2_energy_hcal", gamma2_energy_v[4]);
    event->setDoubleData("gamma2_cone_angle", Blob2_cone_angle);
    event->setDoubleData("gamma2_pvalue", Blob2_pvalue);
    event->setDoubleData("gamma2_MeanE_1stHalf", Blob2_MeanE_1stHalf);
    event->setDoubleData("gamma2_MeanE_2ndHalf", Blob2_MeanE_2ndHalf);
    //event->setDoubleData("gamma1_MedianE_1stHalf", Blob1_MedianE_1stHalf);
    //event->setDoubleData("gamma1_MedianE_2ndHalf", Blob1_MedianE_2ndHalf);
    event->setContainerDoubleData("gamma2_clusters_distance", blob2_clusters_distance);
    event->setContainerDoubleData("gamma2_clusters_energy",   blob2_clusters_energy);

    return true;
}

//---------------------------------------------------------------------------------------------
// Calculte the best straight line that describe the clusters in the blob, in each view,
// then calculate the distance for the line to the muon vertex at z0
// select the minimum distance, if the distance is less than 200 mm, it should be a good blob
// almost pointing to the interaction vertex
//--------------------------------------------------------------------------------------------

double CCProtonPi0::CalProbFitLine(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const
{
  //Gaudi::XYZVector direction = blob->direction();
  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  //-- Fill a vector of points with the (z,x,0, energy) with each cluster in the blob
  //std::vector<LineFit::Point> data_X;
  int n = clusters.size();
  double blob_energy = blob->energy();
  double Xc_X[n], ErrXc_X[n], Zc_X[n], ErrZc_X[n];
  double Xc_U[n], ErrXc_U[n], Zc_U[n], ErrZc_U[n];
  double Xc_V[n], ErrXc_V[n], Zc_V[n], ErrZc_V[n];
  int counterX = 0, counterU = 0, counterV = 0;
  for (c = clusters.begin(); c != clusters.end(); ++c){
    //double x, z, energy;
    if ((*c)->view() == Minerva::IDCluster::X ){   // Only X-view
      double x = (*c)->position();
      double x_err = (*c)->errPosition();
      double z = (*c)->z();
      double z_err = 17.0;   //-- in mm, it is approx the width of a plane
      double energy = (*c)->energy();

      Xc_X[counterX] = x; ErrXc_X[counterX] = x_err*(blob_energy/energy); Zc_X[counterX] = z; ErrZc_X[counterX] = z_err;
      counterX++;
      //data_X.push_back(LineFit::Point(z,x,0.0, energy));
    }

    if ((*c)->view() == Minerva::IDCluster::U ){   // Only X-view
      double x = (*c)->position();
      double x_err = (*c)->errPosition();
      double z = (*c)->z();
      double z_err = 17.0;   //-- in mm, it is approx the width of a plane
      double energy = (*c)->energy();

      Xc_U[counterU] = x; ErrXc_U[counterU] = x_err*(blob_energy/energy); Zc_U[counterU] = z; ErrZc_U[counterU] = z_err;
      counterU++;
      //data_X.push_back(LineFit::Point(z,x,0.0, energy));
    }

    if ((*c)->view() == Minerva::IDCluster::V ){   // Only X-view
      double x = (*c)->position();
      double x_err = (*c)->errPosition();
      double z = (*c)->z();
      double z_err = 17.0;   //-- in mm, it is approx the width of a plane
      double energy = (*c)->energy();

      Xc_V[counterV] = x; ErrXc_V[counterV] = x_err*(blob_energy/energy); Zc_V[counterV] = z; ErrZc_V[counterV] = z_err;
      counterV++;
      //data_X.push_back(LineFit::Point(z,x,0.0, energy));
    }
  }

  double x0 = vert.x(), z0 = vert.z(), u0 = m_mathTool->calcUfromXY(vert.x(), vert.y()), v0 = m_mathTool->calcVfromXY(vert.x(), vert.y());
  std::vector<double> delta;

  if (counterX != 0){
      graphX = new TGraphErrors(counterX, Zc_X, Xc_X, ErrZc_X, ErrXc_X);
      graphX->Fit("lineX", "R");
      fit_lineX = graphX->GetFunction("lineX");

      double deltaX = x0 - fit_lineX->Eval(z0);
      double chi2X = fit_lineX->GetChisquare();
      double ndofX = fit_lineX->GetNDF();

      for (int i = 0; i < counterX; i++){
          info() << "z = " << Zc_X[i] << ", Err_z = " << ErrZc_X[i] << ", x = " << Xc_X[i] << ", Err_x = " << ErrXc_X[i] << endmsg;
      }

      delta.push_back(TMath::Abs(deltaX));
      info() << "Blob fit line: chi2X = " << chi2X << " , ndof = " << ndofX << " , deltaX = " << deltaX << " z0 = " << z0 << " x0 = " << x0 << endmsg;
  }

  if (counterU != 0){
      graphU = new TGraphErrors(counterU, Zc_U, Xc_U, ErrZc_U, ErrXc_U);
      graphU->Fit("lineU", "R");
      fit_lineU = graphU->GetFunction("lineU");

      double deltaU = u0 - fit_lineU->Eval(z0);
      double chi2U = fit_lineU->GetChisquare();
      double ndofU = fit_lineU->GetNDF();

      for (int i = 0; i < counterU; i++){
          info() << "z = " << Zc_U[i] << ", Err_z = " << ErrZc_U[i] << ", u = " << Xc_U[i] << ", Err_u = " << ErrXc_U[i] << endmsg;
      }
      delta.push_back(TMath::Abs(deltaU));
      info() << "Blob fit line: chi2U = " << chi2U << " , ndof = " << ndofU << " , deltaU = " << deltaU << " z0 = " << z0 << " u0 = " << u0 << endmsg;
  }

  if (counterV != 0){
      graphV = new TGraphErrors(counterV, Zc_V, Xc_V, ErrZc_V, ErrXc_V);
      graphV->Fit("lineV", "R");
      fit_lineV = graphV->GetFunction("lineV");

      double deltaV = v0 - fit_lineV->Eval(z0);
      double chi2V = fit_lineV->GetChisquare();
      double ndofV = fit_lineV->GetNDF();

      for (int i = 0; i < counterV; i++){
          info() << "z = " << Zc_V[i] << ", Err_z = " << ErrZc_V[i] << ", v = " << Xc_V[i] << ", Err_v = " << ErrXc_V[i] << endmsg;
      }
      delta.push_back(TMath::Abs(deltaV));
      info() << "Blob fit line: chi2V = " << chi2V << " , ndof = " << ndofV << " , deltaV = " << deltaV << " z0 = " << z0 << " v0 = " << v0 << endmsg;
  }

  //double mag = TMath::Sqrt(delta[0]*delta[0] + delta[1]*delta[1]);
  double mag = -1.0;
  if (!delta.empty()){
      mag = delta[0];
      std::sort(delta.begin(), delta.end());
  }
  info() << "Minimum distance from two views: " << mag << endmsg;

  return mag;
}


//------------------------------------------------------------------------------
// Fit the clusters in the blob to a straight line with fixed points in the x-view
//  interaction vertex (z0, x0)
//  position of the centroid (zc, xc)
//    Return the probability to 1-P(ndof, chi2)
//------------------------------------------------------------------------------
/*
double CCProtonPi0::CalProbFitLine(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const
{
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator c;

    //-- Fill a vector of points with the (z,x,0, energy) with each cluster in the blob
    std::vector<LineFit::Point> data_X;
    std::vector<LineFit::Point> data_U;
    std::vector<LineFit::Point> data_V;
    for (c = clusters.begin(); c != clusters.end(); ++c){
      double x, z, energy;
      if ((*c)->view() == Minerva::IDCluster::X ){   // Only X-view
        x = (*c)->position();
        z = (*c)->z();
        energy = (*c)->energy();

        data_X.push_back(LineFit::Point(z,x,0.0, energy));
      }

      if ((*c)->view() == Minerva::IDCluster::U ){   // Only U-view
        x = (*c)->position();
        z = (*c)->z();
        energy = (*c)->energy();

        data_U.push_back(LineFit::Point(z,x,0.0, energy));
      }


      if ((*c)->view() == Minerva::IDCluster::V ){   // Only V-view
        x = (*c)->position();
        z = (*c)->z();
        energy = (*c)->energy();

        data_V.push_back(LineFit::Point(z,x,0.0, energy));
      }

    }


    //-- Build the Function of the straight line to fit data.
    double x0 = vert.x(), z0 = vert.z(), y0 = vert.y();
    double xc = direction.x(), zc = direction.z(), yc = direction.y();
    double slope = (xc - x0)/(zc - z0);
    TF1 *fit_line = new TF1("fit_line", "[0]+[1]*x", 5000, 15000); //--- approx z position of the vertexes and maxz point of the blobs
    double b = x0 - (slope*z0);
    fit_line->SetParameter(0, b);
    fit_line->SetParameter(1, slope);

    double xs = start_point.x(), ys = start_point.y(), zs = start_point.z();

    double distanceX = TMath::Abs(xs - fit_line->Eval(zs));


    //-- Calculate chi2
    double chi2_X = 0.0;
    for (std::vector<LineFit::Point>::const_iterator i = data_X.begin(); i != data_X.end(); ++i){
        double x_i = i->GetX();
        double y_i = i->GetY();
        //double w_i = i->GetWeight();

        double delta = y_i - fit_line->Eval(x_i);
        //chi2 += (delta*delta)/fit_line->Eval(x_i);
	chi2_X += (delta*delta)/(17*17+17*17);    // This is only considering the widht of a plane, in X-direction the error should be different

        //info() << "Debugging chi2 " << "yi = " << y_i << " , fit_line->Eval(x_i) = " << fit_line->Eval(x_i) << " , delta = " << delta << " , chi2_X = " << chi2_X << endmsg;
    }

    //-- Calculate chi2

    double chi2_U = 0.0;
    for (std::vector<LineFit::Point>::const_iterator i = data_U.begin(); i != data_U.end(); ++i){
        double x_i = i->GetX();
        double y_i = i->GetY();
        //double w_i = i->GetWeight();

        double delta = y_i - fit_line->Eval(x_i);
        //chi2 += (delta*delta)/fit_line->Eval(x_i);
	chi2_U += (delta*delta)/(17*17+17*17);    // This is only considering the widht of a plane, in X-direction the error should be different

        //info() << "Debugging chi2 " << "yi = " << y_i << " , fit_line->Eval(x_i) = " << fit_line->Eval(x_i) << " , delta = " << delta << " , chi2_U = " << chi2_U << endmsg;
    }

    //-- Calculate chi2

    double chi2_V = 0.0;
    for (std::vector<LineFit::Point>::const_iterator i = data_V.begin(); i != data_V.end(); ++i){
        double x_i = i->GetX();
        double y_i = i->GetY();
        //double w_i = i->GetWeight();

        double delta = y_i - fit_line->Eval(x_i);
        //chi2 += (delta*delta)/fit_line->Eval(x_i);
	chi2_V += (delta*delta)/(17*17+17*17);    // This is only considering the widht of a plane, in X-direction the error should be different

        //info() << "Debugging chi2 " << "yi = " << y_i << " , fit_line->Eval(x_i) = " << fit_line->Eval(x_i) << " , delta = " << delta << " , chi2_V = " << chi2_V << endmsg;
    }

    //int ndof_X = data_X.size(), ndof_U = data_U.size(), ndof_V = data_V.size();
    //double p_value = TMath::Prob(chi2, ndof);

    //info() << "Blob: chi2_X = " << chi2_X << " , ndof = " << ndof_X << " , chi2/ndof = " << chi2_X/ndof_X << endmsg;
    //info() << "Blob: chi2_U = " << chi2_U << " , ndof = " << ndof_U << " , chi2/ndof = " << chi2_U/ndof_U << endmsg;
    //info() << "Blob: chi2_V = " << chi2_V << " , ndof = " << ndof_V << " , chi2/ndof = " << chi2_V/ndof_V << endmsg;

    //return chi2_X/ndof_X;

    info() << "Distance in the X-view between start point and line Vertex-Centroid d = " << distanceX << " mm" << endmsg;     // this sould be close to 0.

    return distanceX;
}
*/
//-----------------------------------------------------------------------------
// Calculate the mean energy of the 1st half of the clusters
//-----------------------------------------------------------------------------
double CCProtonPi0::MeanEnergyClusters_1stHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const{

  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  std::vector< std::pair<double, double> > c_dis_ener;
  double dx, dz, energy;

  //--- Loop Over clusters
  for (c = clusters.begin(); c != clusters.end(); ++c){
    switch((*c)->view()){
    case Minerva::IDCluster::X:
      dx = (*c)->position() - vert.x();
      break;

    case Minerva::IDCluster::U:
      dx = (*c)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
      break;

    case Minerva::IDCluster::V:
      dx = (*c)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
      break;

    default:
      throw MinervaException("Unknown cluster view");
    }

    dz = (*c)->z() - vert.z();
    double distance = sqrt( pow(dx,2) + pow(dz,2) );
    energy = (*c)->energy();

    c_dis_ener.push_back(std::make_pair(distance, energy));

  }//-- End loop over clusters

  std::sort(c_dis_ener.begin(), c_dis_ener.end());

  //------- Obtain the Mean Energy of the 1st Half of the clusters -----
  double E_total = 0.0;
  double N_clusters = c_dis_ener.size()/2;
  for (int i = 0; i < N_clusters; i++){
    E_total += c_dis_ener[i].second;
  }
  double E_mean = E_total/N_clusters;

  return E_mean;
}

//-----------------------------------------------------------------------------
// Calculate the mean energy of the 2nd half of the clusters
//-----------------------------------------------------------------------------
double CCProtonPi0::MeanEnergyClusters_2ndHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const{

  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  std::vector< std::pair<double, double> > c_dis_ener;
  double dx, dz, energy;

  //--- Loop Over clusters
  for (c = clusters.begin(); c != clusters.end(); ++c){
    switch((*c)->view()){
    case Minerva::IDCluster::X:
      dx = (*c)->position() - vert.x();
      break;

    case Minerva::IDCluster::U:
      dx = (*c)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
      break;

    case Minerva::IDCluster::V:
      dx = (*c)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
      break;

    default:
      throw MinervaException("Unknown cluster view");
    }

    dz = (*c)->z() - vert.z();
    double distance = sqrt( pow(dx,2) + pow(dz,2) );
    energy = (*c)->energy();

    c_dis_ener.push_back(std::make_pair(distance, energy));

  }//-- End loop over clusters

  std::sort(c_dis_ener.begin(), c_dis_ener.end());

  //------- Obtain the Mean Energy of the 2nd Half of the clusters -----
  double E_total = 0.0;
  double N_clusters = c_dis_ener.size();
  int mid_clusters = N_clusters/2;
  for (int i = mid_clusters; i < N_clusters; i++){
    E_total += c_dis_ener[i].second;
  }
  double E_mean = E_total/N_clusters;

  return E_mean;
}

//-----------------------------------------------------------------------------
// Calculate the median energy of the 2nd half of the clusters
//-----------------------------------------------------------------------------
double CCProtonPi0::MedianEnergyClusters_2ndHalf(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const{

  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  std::vector< std::pair<double, double> > c_dis_ener;
  double dx, dz, energy;

  //--- Loop Over clusters
  for (c = clusters.begin(); c != clusters.end(); ++c){
    switch((*c)->view()){
    case Minerva::IDCluster::X:
      dx = (*c)->position() - vert.x();
      break;

    case Minerva::IDCluster::U:
      dx = (*c)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
      break;

    case Minerva::IDCluster::V:
      dx = (*c)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
      break;

    default:
      throw MinervaException("Unknown cluster view");
    }

    dz = (*c)->z() - vert.z();
    double distance = sqrt( pow(dx,2) + pow(dz,2) );
    energy = (*c)->energy();

    c_dis_ener.push_back(std::make_pair(distance, energy));

  }//-- End loop over clusters

  std::sort(c_dis_ener.rbegin(), c_dis_ener.rend());  //-- sort is decending order

  //------- Obtain the Median Energy of the 2nd Half of the clusters -----
  double E_median = 0.0;
  int N_clusters = c_dis_ener.size()/2;
  if (N_clusters % 2 == 0){
    E_median = (c_dis_ener[N_clusters/2 - 1 ].second + c_dis_ener[N_clusters/2].second)/2;
  }else{
    int position = N_clusters/2 + N_clusters%2;
    E_median = c_dis_ener[position].second;
  }

  //for (int i = 0; i < N_clusters; i++){
  //  E_total += c_dis_ener[i].second;
  //}
  //double frac_E = E_median/E_total;

  //return frac_E;


  return E_median;
}

//-----------------------------------------------------------------------------
// Store a pair vector with Distance to the interaction vertex and the Energy
// of the clusters c_dis_energy(distance, energy)
//-----------------------------------------------------------------------------
std::vector<std::pair<double, double> > CCProtonPi0::EnergyDistanceClusters(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const{

  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  std::vector< std::pair<double, double> > c_dis_ener;
  double dx, dz, energy;

  //--- Loop Over clusters
  for (c = clusters.begin(); c != clusters.end(); ++c){
    switch((*c)->view()){
    case Minerva::IDCluster::X:
      dx = (*c)->position() - vert.x();
      break;

    case Minerva::IDCluster::U:
      dx = (*c)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
      break;

    case Minerva::IDCluster::V:
      dx = (*c)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
      break;

    default:
      throw MinervaException("Unknown cluster view");
    }

    dz = (*c)->z() - vert.z();
    double distance = sqrt( pow(dx,2) + pow(dz,2) );
    energy = (*c)->energy();

    c_dis_ener.push_back(std::make_pair(distance, energy));

  }//-- End loop over clusters

  std::sort(c_dis_ener.begin(), c_dis_ener.end());

  return c_dis_ener;
}


//------------------------------------------------------------------------------
// Calculate the angle of the cone that contains the blob.
//      Find the cluster with max angle between the direction of the blob
//      With the max angle, then calculate the max width
//------------------------------------------------------------------------------
double CCProtonPi0::CalculateConeAngle(SmartRef<Minerva::IDBlob> blob, Gaudi::XYZPoint vert) const
{
  Gaudi::XYZVector direction = blob->direction();
  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  double dx, Dx, dz, Dz;
  double max_angle = 0.0, max_width = 0.0;
  for (c = clusters.begin(); c != clusters.end(); ++c){

    switch((*c)->view()){
    case Minerva::IDCluster::X:
      dx = (*c)->position() - vert.x();
      Dx = direction.x();
      break;

    case Minerva::IDCluster::U:
      dx = (*c)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
      Dx = m_mathTool->calcUfromXY(direction.x(), direction.y());
      break;

    case Minerva::IDCluster::V:
      dx = (*c)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
      Dx = m_mathTool->calcVfromXY(direction.x(), direction.y());
      break;

    default:
      throw MinervaException("Unknown cluster view");
    }

    dz = (*c)->z() - vert.z();
    Dz = direction.z();

    double moduled = sqrt( pow(dx,2) + pow(dz,2) ); // Distance from Vertex to the cluster
    double moduleD = sqrt( pow(Dx,2) + pow(Dz,2) ); // Distance from Vertex to the centroid

    if (moduleD == 0 || moduled == 0) continue;
    double angle = fabs(acos(Dx/moduleD) - acos(dx/moduled) );

    if (angle > max_angle){
      max_angle = angle;
      max_width = moduled*TMath::Sin(max_angle);  //-- mm?
    }

  }

  /*
  //------- Save Data
  if (Blob_ID == 1){
    event->SetDoubleData("Blob1_MaxWidth", max_witdth);
    event->SetDoubleData("Blob1_cone_angle", max_angle);
  }
  if (Blob_ID == 2){
    event->SetDoubleData("Blob2_MaxWidth", max_witdth);
    event->SetDoubleData("Blob2_cone_angle", max_angle);
  }
  */


  return max_width;
  //return max_angle*180.0/3.14159;
}





//------------------------------------------------------------------------------
// Find the maximum z cluster in the shower
//------------------------------------------------------------------------------
double CCProtonPi0::FindShowerMaxZ(SmartRef<Minerva::IDBlob> blob) const
{

  double temp_max = -5;
  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  SmartRefVector<Minerva::IDCluster>::iterator c;

  for ( c = clusters.begin(); c !=clusters.end(); ++c){

    const double z = (*c)->z();
    if (z >= temp_max) temp_max = z;
  }

  return temp_max;

}
//------------------------------------------------------------------------------
// set the track prong Geant4 truth information
//------------------------------------------------------------------------------
void CCProtonPi0::setTrackProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const
{
    std::vector<const Minerva::TG4Trajectory*> parent;
    Minerva::TG4Trajectories* mc_trajectories = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
    for(Minerva::TG4Trajectories::iterator iter = mc_trajectories->begin(); iter != mc_trajectories->end(); ++iter) {
        if( (*iter)->GetParentId() == 0 ) parent.push_back( (*iter) );
    }

    std::vector< std::pair<const Minerva::TG4Trajectory*,double> > trajectories;
    std::vector<int> ntraj(10,-1);
    std::vector<int> nprim(10,-1);
    std::vector<int> pdg(10,-1);
    std::vector<int> inside_od(10,-1);
    std::vector<double> traj_Px(10,-1);
    std::vector<double> traj_Py(10,-1);
    std::vector<double> traj_Pz(10,-1);
    std::vector<double> traj_E(10,-1);
    std::vector<double> traj_P(10,-1);
    std::vector<double> traj_PSelf(10,-1);
    std::vector<double> traj_theta(10,-1);
    std::vector<double> traj_phi(10,-1);
    std::vector<double> end_trajMomentum(10,-1);
    std::vector<double> end_trajX(10,-1);
    std::vector<double> end_trajY(10,-1);
    std::vector<double> end_trajZ(10,-1);
    for(unsigned int p = 0; p < prongs.size(); p++) {
        SmartRef<Minerva::Prong> prong = prongs[p];
        Minerva::TrackVect tracks = prong->minervaTracks();
        for(unsigned int trk = 0; trk < tracks.size(); trk++) {
            double other_energy = 0.0;
            std::map<const Minerva::TG4Trajectory*,double> trajMap = TruthMatcher->getTG4Trajectories(tracks[trk],other_energy);
            std::map<const Minerva::TG4Trajectory*,double>::iterator it;
            if( trajMap.empty() ) continue;

            std::vector< std::pair<const Minerva::TG4Trajectory*,double> > parentVect;
            for(it = trajMap.begin(); it != trajMap.end(); it++) {
                const Minerva::TG4Trajectory* tmp = (*it).first;    
                double fraction = (*it).second;

                if( tmp->GetParentId() == 0 ){
                    parentVect.push_back( std::make_pair(tmp,fraction) );
                }else {
                    int  attempts = 0;
                    bool foundPrimary = false;
                    while( !foundPrimary ) {

                        bool foundMatch = false;
                        Gaudi::XYZPoint start( tmp->GetInitialPosition().x(), tmp->GetInitialPosition().y(), tmp->GetInitialPosition().z() );

                        for(Minerva::TG4Trajectories::iterator tj = mc_trajectories->begin(); tj != mc_trajectories->end(); ++tj) {
                            Gaudi::XYZPoint end( (*tj)->GetFinalPosition().x(), (*tj)->GetFinalPosition().y(), (*tj)->GetFinalPosition().z() );
                            if( end == start ) { 
                                tmp = (*tj); foundMatch = true; 
                                break; 
                            }
                        }

                        if( !foundMatch ) break;
                        else attempts++;

                        if( foundMatch && attempts == 10 ) break;
                        if( tmp->GetParentId() == 0 ) { 
                            parentVect.push_back( std::make_pair(tmp,fraction) ); 
                            foundPrimary = true; 
                            break; 
                        }
                    } // end of while loop
                } // end of else
            } // end loop over traj 
            std::vector< std::pair<const Minerva::TG4Trajectory*,double> > updateParentVect;
            for(unsigned int i = 0; i < parent.size(); i++) {
                double energy = 0.0;
                for(unsigned int j = 0; j < parentVect.size(); j++) {
                    if( parentVect[j].first->GetTrackId() == parent[i]->GetTrackId() ){
                        energy += parentVect[j].second; 
                    }
                }
                updateParentVect.push_back( std::make_pair(parent[i],energy) );
            }

            const Minerva::TG4Trajectory* tmp_traj = NULL;
            double frac_energy = 0.0; 
            for(unsigned int i = 0; i < updateParentVect.size(); i++) {
                if( updateParentVect[i].second < frac_energy ) continue;
                frac_energy = updateParentVect[i].second; 
                tmp_traj    = updateParentVect[i].first; 
            }
            trajectories.push_back( std::make_pair(tmp_traj,frac_energy) );
        }

        if( trajectories.empty() ) return;

        const Minerva::TG4Trajectory* traj = NULL;
        double max = 0;
        for(unsigned int i = 0; i < trajectories.size(); i++) {
            if( trajectories[i].second < max ) continue;
            max  = trajectories[i].second;
            traj = trajectories[i].first;
        }

        if( !traj ) continue;

        Gaudi::LorentzVector traj_4p = traj->GetInitialMomentum();
        double PSelf = sqrt(traj_4p.Px()*traj_4p.Px() +
                traj_4p.Py()*traj_4p.Py() +
                traj_4p.Pz()*traj_4p.Pz());

        int primary = 0;
        if( traj->GetProcessName().find("Primary") != std::string::npos ) primary = 1; 

        ntraj[p]      = int(trajectories.size());
        nprim[p]      = primary;
        pdg[p]        = traj->GetPDGCode();
        traj_Px[p]    = traj_4p.Px();
        traj_Py[p]    = traj_4p.Py();
        traj_Pz[p]    = traj_4p.Pz();
        traj_E[p]     = traj_4p.E();
        traj_P[p]     = traj_4p.P();
        traj_PSelf[p] = PSelf;
        traj_theta[p] = traj_4p.Theta();
        traj_phi[p]   = traj_4p.Phi();

        const Minerva::TG4Trajectory* endTraj = trajectories[trajectories.size()-1].first;
        if( !endTraj ) continue;

        const Minerva::TG4TrajectoryPoints& points = endTraj->GetTrajectoryPoints();

        if( !points.empty() ) {
            Minerva::TG4TrajectoryPoint* endPoint = *points.rbegin();
            if( endPoint ) {
                double end_trajPx = endPoint->GetMomentum().px();
                double end_trajPy = endPoint->GetMomentum().py();
                double end_trajPz = endPoint->GetMomentum().pz();
                end_trajMomentum[p]  = sqrt(    end_trajPx*end_trajPx + 
                        end_trajPy*end_trajPy + 
                        end_trajPz*end_trajPz );

                end_trajX[p] = endPoint->GetPosition().x();
                end_trajY[p] = endPoint->GetPosition().y();
                end_trajZ[p] = endPoint->GetPosition().z();
                Gaudi::XYZPoint point(end_trajX[p],end_trajY[p],end_trajZ[p]);
                inside_od[p] = int(m_OuterDetector->isInside(point));
            } 
        }
    } 

    //-- fill
    if( prongs[0]->filtertaglist()->filterTagExists("PrimaryProton") ) {
        neutrino->setContainerIntData("ntrajProtonProng",ntraj);
        neutrino->setContainerIntData("trajProtonProngPrimary",nprim);
        neutrino->setContainerIntData("trajProtonProngPDG",pdg);
        neutrino->setContainerDoubleData("trajProtonProngPx",traj_Px);
        neutrino->setContainerDoubleData("trajProtonProngPy",traj_Py);
        neutrino->setContainerDoubleData("trajProtonProngPz",traj_Pz);
        neutrino->setContainerDoubleData("trajProtonProngEnergy",traj_E);
        neutrino->setContainerDoubleData("trajProtonProngMomentum",traj_P);
        neutrino->setContainerDoubleData("trajProtonProngPSelf",traj_PSelf);
        neutrino->setContainerDoubleData("trajProtonTheta",traj_theta);
        neutrino->setContainerDoubleData("trajProtonPhi",traj_phi);
        neutrino->setContainerIntData("isProtonInsideOD",inside_od);
        neutrino->setContainerDoubleData("endProtonTrajMomentum",end_trajMomentum);
        neutrino->setContainerDoubleData("endProtonTrajXPosition",end_trajX);
        neutrino->setContainerDoubleData("endProtonTrajYPosition",end_trajY);
        neutrino->setContainerDoubleData("endProtonTrajZPosition",end_trajZ);
    } else if( prongs[0]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
        neutrino->setIntData("ntrajMuonProng",ntraj[0]);
        neutrino->setIntData("trajMuonProngPrimary",nprim[0]);
        neutrino->setIntData("trajMuonProngPDG",pdg[0]);
        neutrino->setDoubleData("trajMuonProngPx",traj_Px[0]);
        neutrino->setDoubleData("trajMuonProngPy",traj_Py[0]);
        neutrino->setDoubleData("trajMuonProngPz",traj_Pz[0]);
        neutrino->setDoubleData("trajMuonProngEnergy",traj_E[0]);
        neutrino->setDoubleData("trajMuonProngMomentum",traj_P[0]);
        neutrino->setDoubleData("trajMuonProngPSelf",traj_PSelf[0]);
        neutrino->setDoubleData("trajMuonTheta",traj_theta[0]);
        neutrino->setDoubleData("trajMuonPhi",traj_phi[0]);
        neutrino->setIntData("isMuonInsideOD",inside_od[0]);
        neutrino->setDoubleData("endMuonTrajMomentum",end_trajMomentum[0]);
        neutrino->setDoubleData("endMuonTrajXPosition",end_trajX[0]);
        neutrino->setDoubleData("endMuonTrajYPosition",end_trajY[0]);
        neutrino->setDoubleData("endMuonTrajZPosition",end_trajZ[0]);
    }

    return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  Pi0 Reconstruction Functions 
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//==============================================================================
//  VertexBlob 
//==============================================================================
void CCProtonPi0::VertexBlob(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent) const
{
    SmartRefVector<Minerva::IDCluster> unusedClusters;
    unusedClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

    double evis = 0.0;
    double energy = 0.0;

    // Find Unused Clusters inside Vertex Blob Radius 
    SmartRefVector<Minerva::IDCluster> VertexBlobClusters = m_extraEnergyTool->getIDClustersInRadiusFromVertex( unusedClusters, vertex, m_vertex_blob_radius); 

    if (!VertexBlobClusters.empty()) {

        Minerva::IDBlob* VertexSphereBlob = new Minerva::IDBlob();
        m_blobUtils->insertIDClusters( VertexBlobClusters, VertexSphereBlob, Minerva::IDBlob::VertexBlobPatRec );

        debug()<<"Creating VertexBlob with "<<VertexSphereBlob->nclusters()<<" clusters"<<endmsg;
        addObject( event, VertexSphereBlob );
        m_hitTagger->applyColorTag( VertexSphereBlob, m_Color_VertexBlob );

        // Get VertexBlob Visible & Calibrated Energy
        evis = VertexSphereBlob->energy();
        energy = m_extraEnergyTool->getCorrectedIDEnergy( VertexBlobClusters, 1.0);

        debug()<<"Vertex Evis = "<<evis<<endmsg;
        debug()<<"Vertex Energy = "<<energy<<endmsg;
        if (truthEvent && m_study_unused_energy){
            SaveTruthUnusedClusterEnergy_NearVertex(truthEvent, VertexBlobClusters);
        }
    }

    event->setDoubleData("vertex_blob_evis", evis);
    event->setDoubleData("vertex_blob_energy", energy);
}

//==============================================================================
//  FilterInSphereClusters()
//==============================================================================
SmartRefVector<Minerva::IDCluster> CCProtonPi0::FilterInSphereClusters( Minerva::PhysicsEvent *event, const SmartRefVector<Minerva::IDCluster>& clusters, const double sphereRadius) const
{

    // Get Interaction Vertex
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();
    const Gaudi::XYZPoint& vtx_position = vertex->position();
    const double x0 = vtx_position.X();
    const double y0 = vtx_position.Y();
    const double z0 = vtx_position.Z();
    double u0 = m_mathTool->calcUfromXY(x0,y0);
    double v0 = m_mathTool->calcVfromXY(x0,y0);

    SmartRefVector<Minerva::IDCluster> sphereClusters;
    SmartRefVector<Minerva::IDCluster>::const_iterator c;

    // Loop over all INPUT clusters
    // Find the clusters stays inside the sphere with given radius
    // Fill sphereClusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {
        const double energy = (*c)->energy();
        const double z      = (*c)->z();
        const double t      = (*c)->position();
        /* Skip low-energy clusters, normal MIP cluster has 1.7*2 MeV = 3.4 MeV */
        if (fSkipLowEnergyClusterVtxEnergy && energy < fThresholdVertexEnergy) continue;

        /* Tranverse position: X,U,or V */
        double t0  = 0.0;
        switch ((*c)->view()) {
            case Minerva::IDCluster::X:
                t0 = x0;
                break;

            case Minerva::IDCluster::U:
                t0 = u0;
                break;

            case Minerva::IDCluster::V:
                t0 = v0;
                break;

            default:
                warning() << "Invalid view" << endmsg;
                exit(1);
        }

        // Distance from this cluster to the primary vertex
        double radius = sqrt(pow(z-z0,2) + pow(t-t0,2));
        if (radius < sphereRadius) sphereClusters.push_back(*c);
    }

    return sphereClusters;
}

//==============================================================================
//  ConeBlobs
//==============================================================================
bool CCProtonPi0::ConeBlobs(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent ) const
{
    //--------------------------------------------------------------------------
    // Loop over ALL UNUSED Clusters and fill the usableClusters
    //      FillUsableClusters may reject some clusters see implementation
    //--------------------------------------------------------------------------
    SmartRefVector<Minerva::IDCluster> usableClusters;
    FillUsableClusters(usableClusters, event, truthEvent);

    //--------------------------------------------------------------------------
    // Analyze usableClusters using AngleScan Class
    //--------------------------------------------------------------------------
    // Initialize Bool Variables
    bool isAngleScan        = false;
    bool isBlobsRecovered   = false;
    bool isBlobsRecovered_Direction = false;
    bool isBlobsRecovered_invMass = false;
    bool isBlobsRecovered_SmallAngle = false;
    bool isBlobsRecovered_SearchView_U = false;
    bool isBlobsRecovered_SearchView_V = false;

    // Get Vertex Position
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Create AngleScan Object
    AngleScan angleScanAlg(usableClusters,vtx_position);
    angleScanAlg.AllowUVMatchWithMoreTolerance(m_AllowUVMatchWithMoreTolerance);
    angleScanAlg.SetUVMatchTolerance(m_UVMatchTolerance);
    angleScanAlg.AllowSmallConeAngle(true);
    angleScanAlg.DoReco();

    event->setIntData("anglescan_ncandx", angleScanAlg.GetNxCandidate());   // Number of Shower Candidates in X
    event->setIntData("anglescan_ncand",  angleScanAlg.GetNCandidate());    // Number of Shower Candidates

    // Initialize foundBlobs
    std::vector<Minerva::IDBlob*> foundBlobs;
    std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
    std::vector<Minerva::IDBlob*>::const_iterator b;
    for (b = angleScanBlobs.begin(); b != angleScanBlobs.end(); ++b){
        if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
        foundBlobs.push_back(*b);
    }

    isAngleScan = (foundBlobs.size() == 2);
    debug()<<"foundBlobs.size() = "<<foundBlobs.size()<<endmsg;
    event->setIntData("anglescan_nfoundBlobs", foundBlobs.size() ); // Number of Found Blobs

    std::vector<Minerva::IDBlob*> GoodBlobs;
    Check_Quality_of_Blobs(foundBlobs, GoodBlobs, event, truthEvent);
    debug() << "The number of Good Blobs is: " << GoodBlobs.size() << endmsg;

    bool Are2BlobsGood = (GoodBlobs.size() == 2);
    bool isBestInvMass = false;

    if (!Are2BlobsGood){
      if (m_TrytoRecover_1Shower && GoodBlobs.size() == 1){
	       bool isBlobGood = Save_1ShowerInfo(GoodBlobs, event);
	       if (isBlobGood){
	          if (truthEvent){
	             Save_1ShowerTruthMatch(GoodBlobs, truthEvent);
	          }

	       if (m_recoverSingleShower_SmallAngle){
	           isBlobsRecovered = RecoverSingleShower_SmallAngle(usableClusters, GoodBlobs, event);
	           isBlobsRecovered_SmallAngle = isBlobsRecovered;
	       }

	  if (!isBlobsRecovered && m_recoverSingleShower_SearchView){
	    // Start Searching Showers in U View
	    isBlobsRecovered = RecoverSingleShower_View_U(usableClusters, GoodBlobs, event);
	    isBlobsRecovered_SearchView_U = isBlobsRecovered;

	    // If U View Fails, Start Searching Showers in V View
	    if (!isBlobsRecovered){
	      isBlobsRecovered = RecoverSingleShower_View_V(usableClusters, GoodBlobs, event);
	      isBlobsRecovered_SearchView_V = isBlobsRecovered;
	    }
	  }

	}
   }//-- End of 1 GoodBlob


      //-----------------------------------------------------------------------
      // Here is trying to find the best 2 from all of the blobs in the event
      //-----------------------------------------------------------------------
    if (GoodBlobs.size() > 2){
	       Save_AllShowersInfo(GoodBlobs, event);
	       if (truthEvent){
	          Save_AllShowersTruthMatch( GoodBlobs, truthEvent);
	       }
	       isBestInvMass = RecoverShowers_InvMass(GoodBlobs, event);
     }

    }//-- End of Are2BlobsGood 
    
    // Mark Pi0Blobs with color for Arachne Scan
    MarkFoundPi0Blobs(foundBlobs); 


    // Recover Blobs if Angle Scan can not find exactly 2 Blobs
    if ( !isAngleScan) {
        if (m_TrytoRecover_1Shower && foundBlobs.size() == 1){
            bool isShowerGood = Save_1ShowerInfo(foundBlobs,event);  
            if (isShowerGood){
                if (truthEvent && m_study_shower_recovery){
                    Save_1ShowerTruthMatch(foundBlobs, truthEvent);
                }

                if (m_recoverSingleShower_SmallAngle){
                    isBlobsRecovered = RecoverSingleShower_SmallAngle(usableClusters, foundBlobs, event);
                    isBlobsRecovered_SmallAngle = isBlobsRecovered;
                }

                if (!isBlobsRecovered && m_recoverSingleShower_SearchView){
                    // Start Searching Showers in U View
                    isBlobsRecovered = RecoverSingleShower_View_U(usableClusters, foundBlobs, event);
                    isBlobsRecovered_SearchView_U = isBlobsRecovered;

                    // If U View Fails, Start Searching Showers in V View
                    if (!isBlobsRecovered){
                        isBlobsRecovered = RecoverSingleShower_View_V(usableClusters, foundBlobs, event);
                        isBlobsRecovered_SearchView_V = isBlobsRecovered;
                    }
                }
            }
        }        

        // --------------------------------------------------------------------
        // Try to recover showers with different methods
        //  1) Check direction of 3 Showers, 
        //          if only 2 of them have good direction, keep the event
        //  2) If all 3 showers have good direction
        //      Get best invariant mass combination
        // --------------------------------------------------------------------
        if (m_TrytoRecover_3Shower && foundBlobs.size() == 3){
            bool areAllShowersGood = Save_3ShowerInfo(foundBlobs,event);  
            if (areAllShowersGood && truthEvent && m_study_shower_recovery){
                Save_3ShowerTruthMatch(foundBlobs, truthEvent);
            }

            debug()<<"Trying to Recover 3 Shower Events"<<endmsg;  
            if (!areAllShowersGood && m_recoverShower_Direction){
                debug()<<"Not All Showers direction is Good"<<endmsg;
                debug()<<"Trying to recover 2 GOOD showers"<<endmsg;
                isBlobsRecovered = RecoverShowers_Direction(foundBlobs, event);
                isBlobsRecovered_Direction = isBlobsRecovered;
            }

            if (!isBlobsRecovered && m_recoverShower_invMass){
                isBlobsRecovered = RecoverShowers_InvMass(foundBlobs, event);
                isBlobsRecovered_invMass = isBlobsRecovered;
            }
        }
    }

    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered",isBlobsRecovered);
    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered_small_angle",isBlobsRecovered_SmallAngle);
    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered_search_view_U",isBlobsRecovered_SearchView_U);
    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered_search_view_V",isBlobsRecovered_SearchView_V);
    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered_direction",isBlobsRecovered_Direction);
    event->filtertaglist()->setOrAddFilterTag( "is_blobs_recovered_invMass",isBlobsRecovered_invMass);

    // ------------------------------------------------------------------------
    //  We have 2 Reconstructed EM Showers if either of them true
    //          Process Found Blobs
    // ------------------------------------------------------------------------
    //if (isAngleScan || isBlobsRecovered) {
    if (Are2BlobsGood || isBestInvMass){
       // processBlobs(event,foundBlobs);  //Roger Commentted
       // ODActivity(event,foundBlobs);    //Rorger Commentted
        processBlobs(event,GoodBlobs);
        ODActivity(event,GoodBlobs);

        // Save FoundBlobs to Global Variables (Energetic Blob is Blob1)
        if( GoodBlobs[0]->energy() >= GoodBlobs[1]->energy()){
            m_Pi0Blob1 = GoodBlobs[0];
            m_Pi0Blob2 = GoodBlobs[1];     //Change foundBlobs by GoodBlobs
        }else{
            m_Pi0Blob1 = GoodBlobs[1];
            m_Pi0Blob2 = GoodBlobs[0];
        }

        //-- Here I should test the relation of the length of the shower vs the energy
        //-- X/Xo = (1/ln2)ln(E/Ec)
        //-- Xo = 50.31 cm   Ec = 97.56 MeV

        return true;
    }else{
        // Clear GoodBlobs Vector
        std::vector<Minerva::IDBlob*>::iterator b;
        for (b = GoodBlobs.begin(); b != GoodBlobs.end(); ++b) { 
            delete *b;
        }
        GoodBlobs.clear();

        return false;
    }
 }
    
double CCProtonPi0::GetTotalpe(SmartRefVector<Minerva::IDCluster> clusters) const
{
  double total_pe = 0.0;
  for (SmartRefVector<Minerva::IDCluster>::iterator it_clus = clusters.begin(); it_clus != clusters.end(); it_clus++ ){
    total_pe += (*it_clus)->pe();
  }

  return total_pe;
}

void CCProtonPi0::GetBlobDirection(Minerva::PhysicsEvent *event, Minerva::IDBlob *blob) const
{

  const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();
  double X0 = vtx_position.x();
  double Y0 = vtx_position.y();
  double Z0 = vtx_position.z();

  double vtx_U = (0.5)*(X0 - 1.7321*Y0);
  double vtx_V = (0.5)*(X0 + 1.7321*Y0);

  //std::vector<double> pos_centroid = {0.0, 0.0, 0.0, 0.0};    //  {X, U, V, Z}
  double pos_centroid[4] = {0.0, 0.0, 0.0, 0.0};    //  {X, U, V, Z}

  SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
  double total_blob_pe = GetTotalpe(clusters);

  for (SmartRefVector<Minerva::IDCluster>::iterator it_clus = clusters.begin(); it_clus != clusters.end(); it_clus++){

    double clus_position = (*it_clus)->position();
    Minerva::IDCluster::View clus_view = (*it_clus)->view();

    double w = (*it_clus)->pe()/total_blob_pe;
    if (clus_view == Minerva::IDCluster::X ){
      double X = clus_position - X0;
      pos_centroid[0] += X*w;
    }
    if (clus_view == Minerva::IDCluster::U ){
      double U = clus_position - vtx_U;
      pos_centroid[1] += U*w;
    }
    if (clus_view == Minerva::IDCluster::V ){
      double V = clus_position - vtx_V;
      pos_centroid[2] += V*w;
    }
    double Z = (*it_clus)->z() - Z0;
    pos_centroid[3] += Z*w;
  }

  debug() << "The position of the centroid is: "
	  << pos_centroid[0] << " X, "
	  << pos_centroid[1] << " U, "
	  << pos_centroid[2] << " V, "
	  << pos_centroid[3] << " Z, " << endmsg;

  // Due that the directional cosine are defined only in the R3 space (X,U,V,Z) ---> (X,Y,Z)
  double centroid_X = pos_centroid[0];
  double centroid_Y = (-1.7321/3)*(pos_centroid[1] - pos_centroid[2]); // Call calYformUV() in the mathTool
  double centroid_Z = pos_centroid[3];

  double norm_centroid = sqrt(pow(centroid_X, 2) + pow(centroid_Y, 2) + pow(centroid_Z, 2));

  debug() << "The position of the centroid is: "
	  << centroid_X << " X, "
	  << centroid_Y << " Y, "
	  << centroid_Z << " Z, "
	  << norm_centroid << " is the norm " << endmsg;

  // Set the direction of the blob
  Gaudi::XYZVector direction;
  direction.SetX(centroid_X/norm_centroid);
  direction.SetY(centroid_Y/norm_centroid);
  direction.SetZ(centroid_Z/norm_centroid);

  blob->setDirection(direction);

}


bool CCProtonPi0::AreBlobsDirectionGood(Minerva::PhysicsEvent *event) const
{
    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<<"Passed NULL m_Pi0Blob"<<endmsg;
        return false;
    }

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Check Gamma1 Quality
    bool goodPosition1  = m_idHoughBlob->GetStartPosition(m_Pi0Blob1, vtx_position, true );
    bool goodDirection1 = m_idHoughBlob->GetDirection(m_Pi0Blob1, vtx_position );

    if (!goodPosition1){
      debug() << " Failing fitting the start Point in Blob 1" << endmsg;
    }

    if (!goodDirection1){
      GetBlobDirection(event, m_Pi0Blob1);
      goodDirection1 = true;
    }

    bool isGoodBlob1 = false;
    if (goodPosition1 && goodDirection1) isGoodBlob1 = true;

    // Check Gamma2 Quality
    bool goodPosition2  = m_idHoughBlob->GetStartPosition(m_Pi0Blob2, vtx_position, true );
    bool goodDirection2 = m_idHoughBlob->GetDirection(m_Pi0Blob2, vtx_position );

    if (!goodPosition2){
      debug() << " Failing fitting the start Point in Blob 2" << endmsg;
    }

    if (!goodDirection2){
      GetBlobDirection(event, m_Pi0Blob2);
      goodDirection2 = true;
    }

    bool isGoodBlob2 = false;
    if (goodPosition2 && goodDirection2) isGoodBlob2 = true;

    if (isGoodBlob1 && isGoodBlob2)  return true;
    else return false;
}

void CCProtonPi0::setBlobData(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction *truthEvent) const
{
    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<< "Pi0Blobs are NULL!"<<endmsg;
    }

    // Calculate dEdX for Blobs
    Calculate_dEdx(event, m_Pi0Blob1, 1);
    Calculate_dEdx(event, m_Pi0Blob2, 2);

    // Calculate Min Separation from Vertex
    double blob1_minsep = CalcMinBlobSeparation(m_Pi0Blob1, event);
    double blob2_minsep = CalcMinBlobSeparation(m_Pi0Blob2, event);

    const Gaudi::XYZVector& blob1_direction = m_Pi0Blob1->direction();
    const Gaudi::XYZVector& blob2_direction = m_Pi0Blob2->direction();
    double gamma1_theta = blob1_direction.theta()*180/M_PI;
    double gamma2_theta = blob2_direction.theta()*180/M_PI;

    debug() << "gamma1_blob_theta = " << gamma1_theta << endmsg;
    debug() << "gamma2_blob_theta = " << gamma2_theta << endmsg;

    event->setIntData("gamma1_blob_nclusters", m_Pi0Blob1->nclusters());
    event->setIntData("gamma1_blob_ndigits", m_Pi0Blob1->getAllDigits().size());
    event->setDoubleData("gamma1_blob_minsep", blob1_minsep);
    event->setDoubleData("gamma1_blob_energy", m_Pi0Blob1->energy());
    event->setDoubleData("gamma1_blob_theta", gamma1_theta);

    event->setIntData("gamma2_blob_nclusters", m_Pi0Blob2->nclusters());
    event->setIntData("gamma2_blob_ndigits", m_Pi0Blob2->getAllDigits().size());
    event->setDoubleData("gamma2_blob_minsep", blob2_minsep);
    event->setDoubleData("gamma2_blob_energy", m_Pi0Blob2->energy());
    event->setDoubleData("gamma2_blob_theta", gamma2_theta);

    // Truth Match
    if(truthEvent){
        SaveTruthClusterEnergy_FoundBlobs(event, truthEvent);
    }

    // Shower Energy Study
    if (truthEvent && m_study_shower_energy){
        // Gamma 1
        SaveBlobDigitEnergy(event, m_Pi0Blob1, 1);
        SaveBlobTrueEvisBySubDetector(event, m_Pi0Blob1, 1);
        SaveSCALHits(event, m_Pi0Blob1, 1);
        SaveSCALHits_Improved(event, m_Pi0Blob1, 1);
        SaveSCAL_minZ_Info(event,m_Pi0Blob1,1); 

        // Gamma 2
        SaveBlobDigitEnergy(event, m_Pi0Blob2, 2);
        SaveBlobTrueEvisBySubDetector(event, m_Pi0Blob2, 2);
        SaveSCALHits(event, m_Pi0Blob2, 2);
        SaveSCALHits_Improved(event, m_Pi0Blob2, 2);
        SaveSCAL_minZ_Info(event,m_Pi0Blob2,2); 
    }
}

//==============================================================================
//  HoughBlob
//==============================================================================
StatusCode CCProtonPi0::HoughBlob(Minerva::PhysicsEvent *event, SmartRefVector<Minerva::IDCluster> idClusters, std::vector<Minerva::IDBlob*>& outBlobs) const
{
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    SmartRefVector<Minerva::IDCluster> idClusViewX, idClusViewU, idClusViewV;
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;

    if ( !m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewX, Minerva::IDCluster::X, 0 )) return StatusCode::FAILURE;
    m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewU, Minerva::IDCluster::U, 0 );
    m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewV, Minerva::IDCluster::V, 0 );

    double r, theta, spX, spZ;
    std::vector<Minerva::IDBlob*> idBlobs;

    for ( int i = 0; i < 3; i++){

        Minerva::IDBlob *idBlob = new Minerva::IDBlob;
        SmartRefVector<Minerva::IDCluster> idHTSeed;
        Gaudi::XYZPoint ref(0,0,0);// energetic cluster position

        if ( !m_idHoughTool->GetReference( idClusViewX, ref ) ) break;
        debug() << " Pass Get Reference" << endmsg;

        if ( !m_idHoughTool->Hough2D( idClusViewX, r, theta, ref, vtx_position ) ) break;
        debug() << " Pass Get Hough2d" << endmsg;

        if ( !m_idHoughBlob->Create2dHTSeed( idClusViewX, idHTSeed, r, theta, ref, spX, spZ ) ) break;
        debug() << " Pass Get Create2dHT " << endmsg;

        Gaudi::XYZPoint startpoint(spX, 0, spZ);
        Gaudi::XYZVector direction(-1/tan(theta*CLHEP::pi/180),0,1);

        if ( !m_idHoughBlob->PseudoCone(idHTSeed, idClusViewX, direction, vtx_position ) ) continue;
        debug() << " Pass PseudoCone " << endmsg;

        if ( !m_idHoughBlob->XUVMatch(idHTSeed, idClusViewU, idClusViewV, 50 ) ) continue;
        debug() << " Pass XUV Match - First " << endmsg;

        if ( !m_idHoughBlob->isPhoton(idHTSeed, vtx_position) ) continue;

        idBlob->add( idHTSeed );
        debug() << " Added idBlob with size " << idHTSeed.size() << endmsg;

        idBlobs.push_back(idBlob);
    }

    if ( idBlobs.empty() ) return StatusCode::FAILURE;

    std::vector<Minerva::IDBlob*>::iterator itBlob;
    for (itBlob = idBlobs.begin(); itBlob != idBlobs.end(); ++itBlob ){
        Minerva::IDBlob *blobTemp = new Minerva::IDBlob;
        SmartRefVector<Minerva::IDCluster> idSeed = (*itBlob)->clusters();
        if ( !m_idHoughBlob->XUVMatch( idSeed, idClusViewU, idClusViewV, 120 ) ) continue;
        blobTemp->add(idSeed);
        double Uclusters = blobTemp->nclusters(Minerva::IDCluster::U);
        double Vclusters = blobTemp->nclusters(Minerva::IDCluster::V);
        debug() << *itBlob << " " << Uclusters << " U clusters; " << Vclusters << " V clusters " << endmsg;
        (*itBlob)->clear();
        if ( Uclusters != .0 || Vclusters != .0 ) {
            debug() << *itBlob << " inserting to FinalBlobs " << endmsg;
            (*itBlob)->add(idSeed);
            m_idHoughBlob->GetStartPosition( *itBlob, vtx_position, true );
            m_idHoughBlob->GetDirection( *itBlob, vtx_position );
            outBlobs.push_back(*itBlob);
        }
    }

    for ( it_clus = idClusViewX.begin(); it_clus != idClusViewX.end(); it_clus++ )
        m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );

    for ( it_clus = idClusViewU.begin(); it_clus != idClusViewU.end(); it_clus++ )
        m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );

    for ( it_clus = idClusViewV.begin(); it_clus != idClusViewV.end(); it_clus++ )
        m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );


    debug() << " Hough Transform is done! " << endmsg;

    return StatusCode::SUCCESS;
}

//==============================================================================
//  processBlobs
//==============================================================================
void CCProtonPi0::processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const
{
  int count = 0;
  std::vector<Minerva::IDBlob*>::iterator itBlob = idBlobs.begin();
  for (; itBlob != idBlobs.end(); ++itBlob){

    (*itBlob)->setHistory( Minerva::IDBlob::Unused );
    (*itBlob)->setPatRecHistory( Minerva::IDBlob::Pi0IDBlobPatRec );

    if ( (*itBlob)->nclusters( Minerva::IDCluster::X ) > 0 )  count++;
    if ( (*itBlob)->nclusters( Minerva::IDCluster::U ) > 0 )  count++;
    if ( (*itBlob)->nclusters( Minerva::IDCluster::V ) > 0 )  count++;

    debug() << " Storing Photon blob (key: " << (*itBlob)->key() << ") with " << (*itBlob)->nclusters() << " clusters " << endmsg;
    addObject( event, *itBlob );
  }

  if ( count != 6 ){
    event->filtertaglist()->setOrAddFilterTag( "is_twoDBlob", true );
  }

  debug() << "pi0 candidate" << endmsg;
  debug() << " photon 1 is blob: " << idBlobs[0]->key() << endmsg;
  debug() << " photon 2 is blob: " << idBlobs[1]->key() << endmsg;
  idBlobs[0]->setIntData("Photon1", true);
  idBlobs[1]->setIntData("Photon2", true);

  std::pair<int,double> result1 = OneParLineFitBlob(idBlobs[0], event);
  std::pair<int,double> result2 = OneParLineFitBlob(idBlobs[1], event);

  double g1distance = TwoParLineFitBlobVtxDistance(idBlobs[0], event);
  double g2distance = TwoParLineFitBlobVtxDistance(idBlobs[1], event);

  event->setIntData("g1blob_1ParFit_ndof",result1.first);
  event->setIntData("g2blob_1ParFit_ndof",result2.first);
  event->setDoubleData("g1blob_1ParFit_fval", result1.second);
  event->setDoubleData("g2blob_1ParFit_fval", result2.second);
  event->setDoubleData("g1blob_2ParFit_vtx_distance", g1distance);
  event->setDoubleData("g2blob_2ParFit_vtx_distance", g2distance);

}

//==============================================================================
//  ApplyAttenuationCorrection
//==============================================================================
void CCProtonPi0::ApplyAttenuationCorrection(Minerva::IDBlob* blob) const
{
    const double x0 = blob->startPoint().x();
    const double y0 = blob->startPoint().y();
    const double z0 = blob->startPoint().z();

    const double sx = blob->direction().x();
    const double sy = blob->direction().y();
    const double sz = blob->direction().z();

    if (std::abs(sz) < std::cos(89.0*TMath::DegToRad())) {
        debug() << "Cannot make attenuation correction for high-angle gamma" << endmsg;
        return;
    }

    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
            c != clusters.end(); ++c) {
        const double z = (*c)->z();
        const double t = (z-z0)/sz;
        const double x = x0 + sx*t;
        const double y = y0 + sy*t;

        m_AttenuationCorrectionTool->calibrateCluster(*c,Gaudi::XYZPoint(x,y,z));
    }

    blob->updateEnergy();

    debug() << "AttenuationCorrection (after):  " << blob->energy() << endmsg;

}

//==============================================================================
//  CalcMinBlobSeparation
//==============================================================================
double CCProtonPi0::CalcMinBlobSeparation(const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const
{
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();
    const double x0 = vtx_position.X();
    const double y0 = vtx_position.Y();
    const double z0 = vtx_position.Z();
    double u0 = m_mathTool->calcUfromXY(x0,y0);
    double v0 = m_mathTool->calcVfromXY(x0,y0);

    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    double dmin = 1e6;
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin(); c != clusters.end(); ++c) {
      const double energy = (*c)->energy();
      const double z      = (*c)->z();
      const double t      = (*c)->position();
      /* Skip low-energy clusters, normal MIP cluster has 1.7*2 MeV = 3.4 MeV */
      if (energy < 1.0) continue;

        /* Tranverse position: X,U,or V */
        double t0  = 0.0; 
        switch ((*c)->view()) {

            case Minerva::IDCluster::X:
                t0 = x0;
                break;

            case Minerva::IDCluster::U:
                t0 = u0;
                break;

            case Minerva::IDCluster::V:
                t0 = v0;
                break;

            default:
                warning() << "Invalid view" << endmsg;
                exit(1);

        }
        /* Distance from this cluster to the fitted vertex */
        double distance = sqrt(pow(z-z0,2) + pow(t-t0,2));
        dmin = std::min(dmin,distance);
    }

    return dmin;
}

//==============================================================================
//  ODActivity
//==============================================================================
StatusCode CCProtonPi0::ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const
{
    // get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();

    double od_energy[6] = {0,0,0,0,0,0}; 
    int od_nclusters[6] = {0,0,0,0,0,0}; 
    double od_time[6] = {0,0,0,0,0,0};
    double od_pe[6] = {0,0,0,0,0,0};
    int frame_low = 1000;
    int frame_high = -1000;
    Minerva::BarID barId_z_low;
    Minerva::BarID barId_z_high;

    int story_low = 1000;
    int story_high = -1000;
    Minerva::BarID barId_t_low;
    Minerva::BarID barId_t_high;

    Gaudi::XYZPoint bar_pos;

    SmartRefVector<Minerva::ODCluster> odClusters = event->select<Minerva::ODCluster>( "All", "All" );
    if ( odClusters.empty() ) return StatusCode::FAILURE;

    SmartRefVector<Minerva::ODCluster>::iterator itOD = odClusters.begin();

    for ( ; itOD != odClusters.end(); itOD++ ){
        const Minerva::BarID barID = (*itOD)->barid();
        int frame = barID.frame();

        if (frame < frame_low) {
            frame_low = frame;
            barId_z_low = barID;
        }
        if (frame > frame_high) {
            frame_high = frame;
            barId_z_high = barID;
        }

        int story = barID.story();

        if (story < story_low) {
            story_low = story;
            barId_t_low = barID;
        }
        if (story > story_high) {
            story_high = story;
            barId_t_high = barID;
        }

        double energy = (*itOD)->energy();
        int tower = barID.tower();
        if (tower>=1 && tower<=6) {
            od_nclusters[tower-1]++;
            od_energy[tower-1] += energy;   
            od_time[tower-1] += (*itOD)->time()*(*itOD)->pe(); // od time charge weigthed 
            od_pe[tower-1] += (*itOD)->pe();
        } else {
            warning() << "invalid tower number: " << tower << endmsg;
        }

    }

    debug() << " Pass loop over ODClusters" << endmsg;

    double maxODE = -1;
    int od_mostEnergeticTower = 0;
    std::vector<double> energyVector, timeVector, nclustersVector, distanceVector;
    std::vector<double> timeBlobODVector, timeBlobMuonVector, timeblobVector;

    for (int i=0; i<6; i++) {
        if ( od_pe[i] != 0 ) od_time[i] = od_time[i]/od_pe[i];
        energyVector.push_back(od_energy[i]);
        timeVector.push_back(od_time[i]);
        nclustersVector.push_back(od_nclusters[i]);
        if (od_energy[i]>maxODE) {
            maxODE = od_energy[i];
            od_mostEnergeticTower = i+1;
        }
    }

    //IDBlobs and ODTower
    Minerva::IDCluster::View view;
    if ( od_mostEnergeticTower == 2 || od_mostEnergeticTower == 5 ) view =  Minerva::IDCluster::X;
    else if ( od_mostEnergeticTower == 3 || od_mostEnergeticTower == 6 ) view =  Minerva::IDCluster::U;
    else if ( od_mostEnergeticTower == 1 || od_mostEnergeticTower == 4 ) view =  Minerva::IDCluster::V;
    else {
        warning() << " invalid tower number on IDBlobTower " << od_mostEnergeticTower << endmsg;
        return StatusCode::FAILURE;
    }

    std::vector<Minerva::IDBlob*>::iterator itBlob = idBlobs.begin();

    for ( ; itBlob != idBlobs.end(); itBlob++ ){

        double min_position = 10000, max_position = -10000;
        SmartRefVector<Minerva::IDCluster> idClusters =  (*itBlob)->clusters();
        SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusters.begin();

        for ( ; itClus !=  idClusters.end(); itClus++ ){

            if ( (*itClus)->view() != view ) continue;
            if ( od_mostEnergeticTower <= 3 ) {
                if ( (*itClus)->position() < min_position )
                    min_position = (*itClus)->position();
            }
            else {
                if ( (*itClus)->position() > max_position )
                    max_position = (*itClus)->position();
            }
        }

        if ( od_mostEnergeticTower < 3 ) distanceVector.push_back( fabs(-1060-min_position) );
        else distanceVector.push_back( fabs(1060+max_position) );
        debug() << " ID Blob time " << (*itBlob)->time() << "; Tower time " << od_time[od_mostEnergeticTower-1]
            << " Muon time " << muon_position.T() << endmsg; 
        timeblobVector.push_back( (*itBlob)->time() );
        timeBlobODVector.push_back( (*itBlob)->time()-od_time[od_mostEnergeticTower-1] );
        timeBlobMuonVector.push_back( (*itBlob)->time()- muon_position.T());
    }

    event->setDoubleData( "od_upstreamFrame", frame_low );
    event->setDoubleData( "od_downstreamFrame", frame_high );
    event->setDoubleData( "od_upstreamFrame_z", m_odDet->getDeODFrame(barId_z_low)->getZCenter() );
    event->setDoubleData( "od_downstreamFrame_z", m_odDet->getDeODFrame(barId_z_high)->getZCenter() );

    event->setDoubleData( "od_highStory", story_high );
    event->setDoubleData( "od_lowStory", story_low );
    bar_pos = m_odDet->getDeODFrame(barId_t_high)->getBarPos(barId_t_high);
    event->setDoubleData( "od_highStory_t", sqrt( bar_pos.x()*bar_pos.x() + bar_pos.y()*bar_pos.y()) );
    bar_pos = m_odDet->getDeODFrame(barId_t_low)->getBarPos(barId_t_low);
    event->setDoubleData( "od_lowStory_t", sqrt( bar_pos.x()*bar_pos.x() + bar_pos.y()*bar_pos.y()) );

    event->setDoubleData( "od_maxEnergy", maxODE );
    event->setIntData( "od_energeticTower", od_mostEnergeticTower );

    event->setContainerDoubleData( "od_distanceBlobTower", distanceVector );
    event->setContainerDoubleData( "od_towerEnergy", energyVector );
    event->setContainerDoubleData( "od_towerNClusters", nclustersVector );
    event->setContainerDoubleData( "od_towerTime", timeVector );
    event->setContainerDoubleData( "od_idBlobTime", timeblobVector );
    event->setContainerDoubleData( "od_towerTimeBlobOD", timeBlobODVector );
    event->setContainerDoubleData( "od_towerTimeBlobMuon", timeBlobMuonVector );

    return StatusCode::SUCCESS;
}

//==============================================================================
//  Calculate_dEdx
//==============================================================================
void CCProtonPi0::Calculate_dEdx( Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob,  unsigned int blob_number) const
{
    //Sanity Check
    if ( blob == NULL ){
        warning()<<"Passed NULL Blob"<<endmsg;
    }

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    debug() << "Calculate dEdx for blob: " << blob_number <<  " (forward,backward) ("
        << (blob->direction().Z() > 0) << "," << (blob->direction().Z() < 0) << ")"
        << endmsg;
    std::map<int, std::vector<double> > planeClusterEnergyMap;
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    debug()<<"Vertex to increasingDistanceFromVertex = ("<<vtx_position.X()<<","<<vtx_position.Y()<<","<<vtx_position.Z()<<")"<<endmsg;

    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(vtx_position));

    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
            c != clusters.end(); ++c) {
        const Minerva::PlaneID& planeId = (*c)->planeid();
        unsigned int plane  = planeId.plane();
        unsigned int module = planeId.module();
        const unsigned int planeNo = 2*module + plane;

        planeClusterEnergyMap[planeNo].push_back((*c)->energy());
    }


    std::vector<int>    blob_cluster_occupancy;
    std::vector<double> blob_cluster_energy;

    double total     = 0.0;
    double total1    = 0.0; // exclude the first plane
    bool doublet     = false;
    bool empty_plane = false;

    for ( std::map<int, std::vector<double> >::iterator pl = planeClusterEnergyMap.begin(); pl != planeClusterEnergyMap.end(); ++pl) {
      const int plane = pl->first;
      const std::vector<double>& clusterEnergies = pl->second;

      std::map<int, std::vector<double> >::iterator next_pl = pl; ++next_pl;
      const int next_plane = next_pl->first;

      const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
      blob_cluster_occupancy.push_back(clusterEnergies.size());
      blob_cluster_energy.push_back(total_plane_energy);

      total += total_plane_energy;
      if (pl != planeClusterEnergyMap.begin()) total1 += total_plane_energy;

      if (clusterEnergies.size() > 1)     doublet = true;

      if (next_pl != planeClusterEnergyMap.end() && std::abs(next_plane-plane) > 1) empty_plane = true;

      if (std::distance(planeClusterEnergyMap.begin(),pl) > 4) break;
    }

    // Less complicated logic for cluster energies from the shower back since
    // I'm not making plot for the NIM paper!
    std::vector<double> blob_rev_cluster_energy; // cluster energies from the shower back (rev. for reverse)
    std::map<int, std::vector<double> >::reverse_iterator pl;

    for ( std::map<int, std::vector<double> >::reverse_iterator pl = planeClusterEnergyMap.rbegin(); pl != planeClusterEnergyMap.rend(); ++pl) {
      const std::vector<double>& clusterEnergies = pl->second;
      const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
      blob_rev_cluster_energy.push_back(total_plane_energy);
    }

    const Gaudi::XYZVector& blob_direction = blob->direction();
    double blob_theta = blob_direction.theta();
    double cos_blob_theta = TMath::Cos(blob_theta);

    if (blob_number == 1) {
        event->setIntData("g1dedx_nplane", planeClusterEnergyMap.size());
        event->setIntData("g1dedx_doublet", (int) doublet);
        event->setIntData("g1dedx_empty_plane", (int) empty_plane);
	event->setDoubleData("g1dedx_cos_theta", cos_blob_theta);
        event->setDoubleData("g1dedx_total", total);
        if (cos_blob_theta != 0) event->setDoubleData("g1dedx", total/(blob_cluster_energy.size()*cos_blob_theta));
        event->setDoubleData("g1dedx_total1", total1);
        event->setDoubleData("g1dedx1", total1/std::max((int)blob_cluster_energy.size()-1,1));
        event->setContainerIntData("g1dedx_cluster_occupancy", blob_cluster_occupancy);
        event->setContainerDoubleData("g1dedx_cluster_energy", blob_cluster_energy);
        event->setContainerDoubleData("g1dedx_rev_cluster_energy", blob_rev_cluster_energy);
    }

    if (blob_number == 2) {
        event->setIntData("g2dedx_nplane", planeClusterEnergyMap.size());
        event->setIntData("g2dedx_doublet", (int) doublet);
        event->setIntData("g2dedx_empty_plane", (int) empty_plane);
	event->setDoubleData("g2dedx_cos_theta", cos_blob_theta);
        event->setDoubleData("g2dedx_total", total);
	if (cos_blob_theta != 0) event->setDoubleData("g2dedx", total/(blob_cluster_energy.size()*cos_blob_theta));
        event->setDoubleData("g2dedx_total1", total1);
        event->setDoubleData("g2dedx1", total1/std::max((int)blob_cluster_energy.size()-1,1));
        event->setContainerIntData("g2dedx_cluster_occupancy", blob_cluster_occupancy);
        event->setContainerDoubleData("g2dedx_cluster_energy", blob_cluster_energy);
        event->setContainerDoubleData("g2dedx_rev_cluster_energy", blob_rev_cluster_energy);
    }
}

//==============================================================================
//  PreFilterPi0
//==============================================================================
bool CCProtonPi0::PreFilterPi0(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Visible Energy Limits for the PreFilter [MeV]
    // Other = Tracker + ECAL + HCAL
    const double max_Evis_Target = 20;
    const double min_Evis_Other = 50; // energy smaller than  50 Mev must be ignored - Can't be a Pi0
    const double max_Evis_Other = 50000; // No restrictions on max energy, therefore this high value

    const double inter_mod_z   = 45.0; // From geometry spreadsheet
    const double inter_plane_z = inter_mod_z/2.0;
    const double z0 = vtx_position.Z();

    // Use only UNUSED clusters
    SmartRefVector<Minerva::IDCluster> unusedClusters;
    unusedClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    //--------------------------------------------------------------------------
    // Find Near-Vertex Clusters
    //--------------------------------------------------------------------------
    double nearvtx_total = 0.0;

    // Save all near-vertex clusters
    SmartRefVector<Minerva::IDCluster> nearVertexClusters;
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = unusedClusters.begin();
            iter_cluster != unusedClusters.end();
            ++iter_cluster ){

        const double cluster_z = (*iter_cluster)->z();
        if (std::abs(cluster_z-z0) < 5*inter_plane_z){ 
            nearVertexClusters.push_back(*iter_cluster);
        }
    }

    debug() << "Total Number of near-vertex clusters: " << nearVertexClusters.size() << endmsg;

    // Save Visible Energy of near-vertex clusters for each view
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = nearVertexClusters.begin();
            iter_cluster != nearVertexClusters.end(); 
            ++iter_cluster ){

        const double energy = (*iter_cluster)->energy();

        nearvtx_total += energy;
    }

    //--------------------------------------------------------------------------
    // Save Visible Energy from UNUSED Clusters for each sub detector
    //      Target, Tracker, ECAL, HCAL, TOTAL
    //--------------------------------------------------------------------------
    double ntgtEvis = 0.0;
    double trkrEvis = 0.0;
    double ecalEvis = 0.0;
    double hcalEvis = 0.0;
    double total = 0.0;

    SmartRefVector<Minerva::IDCluster> ecalClusters;
    SmartRefVector<Minerva::IDCluster> hcalClusters;
    SmartRefVector<Minerva::IDCluster> otherClusters;

    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = unusedClusters.begin();
            iter_cluster != unusedClusters.end(); 
            ++iter_cluster ){

        const double energy = (*iter_cluster)->energy();
        Minerva::IDCluster::Subdet subdet = (*iter_cluster)->subdet();

        // Total Visible Energy in complete MINERvA Detector
        total += energy;

        // Visible Energy in SubDetectors
        if (subdet == Minerva::IDCluster::NuclTargs){ 
            ntgtEvis += energy;
        }else if (subdet == Minerva::IDCluster::Tracker){ 
            trkrEvis += energy;
            otherClusters.push_back(*iter_cluster);
        }else if (subdet == Minerva::IDCluster::ECAL){
            ecalEvis += energy;
            ecalClusters.push_back(*iter_cluster);
            otherClusters.push_back(*iter_cluster);
        }else if (subdet == Minerva::IDCluster::HCAL){
            hcalEvis += energy;
            hcalClusters.push_back(*iter_cluster);
            otherClusters.push_back(*iter_cluster);
        }else{ 
            debug() <<"Cluster SubDetector does not found!"<<endmsg;
        }
    } // loop over Unused Clusters

    // Visible Energy except Target Region
    const double otherevis = trkrEvis + ecalEvis + hcalEvis;

    event->setDoubleData("preFilter_evis_nearvtx", nearvtx_total);
    event->setDoubleData("preFilter_evis_total", total);
    event->setDoubleData("preFilter_evis_NuclearTarget", ntgtEvis);
    event->setDoubleData("preFilter_evis_Tracker", trkrEvis);
    event->setDoubleData("preFilter_evis_ECAL", ecalEvis);
    event->setDoubleData("preFilter_evis_HCAL", hcalEvis);
    event->setDoubleData("preFilter_evis_TotalExceptNuclearTarget", otherevis);

    // ------------------------------------------------------------------------
    // Get Cluster Truth Info
    // ------------------------------------------------------------------------
    if (truthEvent && m_study_unused_energy){ 
        SaveTruthUnusedClusterEnergyInsideDetector(truthEvent, ecalClusters, hcalClusters, otherClusters); 
    }    

    //--------------------------------------------------------------------------
    // preFilter_Result
    // 0 = Passes Filter
    // 1 = Rejected by Target Filter
    // 2 = Rejected by Max Other
    // 3 = Rejected by Min Other
    //--------------------------------------------------------------------------

    // Check Target Region
    if ( ntgtEvis > max_Evis_Target){
        debug() <<"Evis(target) bigger than max value!"<<endmsg;
        debug() <<"Evis(target) = "<<ntgtEvis<<" max =  "<<max_Evis_Target<<endmsg;
        event->setIntData("preFilter_Result",1);
        return false;
    }

    // Check Other Regions (Tracker + ECAL + HCAL)
    if ( otherevis > max_Evis_Other){
        debug() <<"Evis(other) bigger than max value!"<<endmsg;
        debug() <<"Evis(other) = "<<otherevis<<" max =  "<<max_Evis_Other<<endmsg;
        event->setIntData("preFilter_Result",2);
        event->setDoubleData("preFilter_rejectedEnergy",otherevis);
        return false;
    }else if( otherevis < min_Evis_Other ){
        debug() <<"Evis(other) lower than min value!"<<endmsg;
        debug() <<"Evis(other) = "<<otherevis<<" min =  "<<min_Evis_Other<<endmsg;
        event->setIntData("preFilter_Result",3);
        event->setDoubleData("preFilter_rejectedEnergy",otherevis);
        return false;
    }

    // Passed all Tests
    debug() <<"Evis(target) = "<<ntgtEvis<<" Evis(other) = "<<otherevis<<endmsg;
    debug() <<"max_Target = "<<max_Evis_Target<<" min_Other =  "<<min_Evis_Other<<" max_Other =  "<<max_Evis_Other<<endmsg;
    event->setIntData("preFilter_Result",0);

    return true;

}

/* Do a line fit for blob digits in the X view. Force the fitted line to
   go through the primary vertex
   */
std::pair<int,double> CCProtonPi0::OneParLineFitBlob(const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const
{
    /* Get all X blob digits via clusters since clusters have 'view' info */
    SmartRefVector<Minerva::IDDigit> digits;
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
            c != clusters.end(); ++c) {
        if ((*c)->view() != Minerva::IDCluster::X) continue;
        const SmartRefVector<Minerva::IDDigit>& cdigits = (*c)->digits();
        std::copy(cdigits.begin(), cdigits.end(), std::back_inserter(digits));
    }
    // Filter digits in the X view here
    std::vector<LineFit::Point> data;
    double total = 0.0;
    for (SmartRefVector<Minerva::IDDigit>::const_iterator d = digits.begin();
            d != digits.end(); ++d) {

        if ((*d)->pe() < 3.0) continue;

        const Minerva::DePlane* plane = m_InnerDetector->getDePlane((*d)->stripid());
        if (!plane) continue;

        total += (*d)->normEnergy();

        const double z = plane->getZCenter();
        const double x = plane->getTPos((*d)->stripid());
        const double w = (*d)->normEnergy();
        data.push_back(LineFit::Point(z,x,0.0,w));

    }

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();
    const double x0 = vtx_position.X();
    const double z0 = vtx_position.Z();

    /* Track slope near the vertex to calculate x(z) */
    SmartRef<Minerva::Track> muonTrack = m_MuonProng->minervaTracks().front();
    const double tkx = muonTrack->nearestNode(z0)->state().ax();

    const double zmin  = -20.0;      /* Scan range in Z +/- 2 cm around the vertex */
    const double zmax  = +20.0;
    const double zstep = 5.0;        /* Scan step */

    std::vector<double> fval_vector; /* fit chi2 of the 9 scanned points */
    double fval_min = 1.e12;         /* The minimum chi2 from the scanned points */
    unsigned int ndof = digits.size();
    for (double dz = zmin; dz <= zmax; dz += zstep) {
        const double z_i = z0 + dz;
        const double x_i = tkx*(z_i - z0) + x0;
        /*
           std::cout << "\tScan i: " << (int)(dz-zmin)/zstep << " "
           << z_i << " " << x_i
           << std::endl;
           */
        TFitterMinuit* minuit = new TFitterMinuit;
        OneParLineFit fittingFunction(data);
        fittingFunction.SetFixedPoint(z_i,x_i);
        minuit->SetMinuitFCN(&fittingFunction);

        /* Set paramter (index,name,init,err,range) */
        minuit->SetParameter(0,"kx",1.0,0.1,-1000,1000);

        minuit->SetPrintLevel(-1);
        minuit->CreateMinimizer();

        int err = minuit->Minimize();
        if (err > 0) {
            std::cerr << "Fitting error: " << errno << std::endl;
            fval_vector.push_back(1.e9);
            continue;
        }

        //const double kfit = minuit->GetParameter(0);
        const double fval = fittingFunction.GetFCN();

        fval_min = std::min(fval,fval_min);
    }
    std::cout << "\tLine fitting: chi2/ndof/total: " << fval_min << "/" << ndof
        << "/" << total 
        << std::endl;

    return std::make_pair<int,double>(ndof,fval_min/total);
}

double CCProtonPi0::TwoParLineFitBlobVtxDistance(const Minerva::IDBlob* blob, Minerva::PhysicsEvent *event) const
{

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();
    const double x0 = vtx_position.Z();
    const double y0 = vtx_position.X();

    /* Get all X blob digits via clusters since clusters have 'view' info */
    std::vector<LineFit::Point> data;
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
            c != clusters.end(); ++c) {

        if ((*c)->view() != Minerva::IDCluster::X) continue;

        // use cluster (x,z) position for now. Could try to use digits later
        const double x = (*c)->z();         // not a typo: z is x in 2D
        const double y = (*c)->position();  // tranverse position is y in 2D
        const double w = (*c)->energy();    // use cluster energy as weight

        data.push_back(LineFit::Point(x,y,0.0,w));

    }

    TwoParLineFit fittingFunc(data);
    TFitterMinuit* minuit = new TFitterMinuit;
    minuit->SetMinuitFCN(&fittingFunc);
    minuit->SetParameter(0,"k",  1.0,  5.0,-1000,1000);
    minuit->SetParameter(1,"b",100.0,500.0,-1.e6,1.e6);

    minuit->SetPrintLevel(-1);
    minuit->CreateMinimizer();

    int err = minuit->Minimize();

    if (err > 0) {
        std::cerr << "TwoParLineFit:: fitting error: " << errno << std::endl;
    }

    // Parameters from the fit in the form a*x + b*y + c = 0
    const double a = minuit->GetParameter(0);
    const double b = -1.0;
    const double c = minuit->GetParameter(1);

    // distance from the fitted line to the muon vertex |a*x0 + b*y0 + c|/sqrt(a*a+b*b)
    const double distance = std::abs(a*x0 + b*y0 + c)/std::sqrt(a*a + b*b);
    std::cout << "TwoParLineFit: distance to vtx: " << distance << std::endl;

    return distance;
}



double CCProtonPi0::calcDistance(  double x1, double y1, double z1,
        double x2, double y2, double z2) const
{
    double d;
    d = sqrt(   (x1-x2)*(x1-x2) + 
            (y1-y2)*(y1-y2) + 
            (z1-z2)*(z1-z2));
    return d;   
}

void CCProtonPi0::saveMichelProngToNTuple(Minerva::PhysicsEvent* event, Minerva::Prong &michelProng, bool isVertex, bool isLargeSearch) const
{


    event->setDoubleData("michelProng_energy",michelProng.getDoubleData("energy"));
    event->setDoubleData("michelProng_time_diff",michelProng.getDoubleData("time_diff"));
    event->setDoubleData("michelProng_end_Z",michelProng.getDoubleData("edz"));
    event->setDoubleData("michelProng_begin_Z",michelProng.getDoubleData("bgz"));
    

if (isVertex){
        if (isLargeSearch){
            event->setDoubleData("vtx_michelProng_Large_distance",michelProng.getDoubleData("distance"));
            event->setDoubleData("vtx_michelProng_Large_energy",michelProng.getDoubleData("energy"));
            event->setDoubleData("vtx_michelProng_Large_time_diff",michelProng.getDoubleData("time_diff"));
            event->setDoubleData("vtx_michelProng_Large_end_Z",michelProng.getDoubleData("edz"));
            event->setDoubleData("vtx_michelProng_Large_begin_Z",michelProng.getDoubleData("bgz"));
        }else{
            event->setDoubleData("vtx_michelProng_distance",michelProng.getDoubleData("distance"));
            event->setDoubleData("vtx_michelProng_energy",michelProng.getDoubleData("energy"));
            event->setDoubleData("vtx_michelProng_time_diff",michelProng.getDoubleData("time_diff"));
            event->setDoubleData("vtx_michelProng_end_Z",michelProng.getDoubleData("edz"));
            event->setDoubleData("vtx_michelProng_begin_Z",michelProng.getDoubleData("bgz"));
        }
    }else{
        event->setDoubleData("track_michelProng_distance",michelProng.getDoubleData("distance"));
        event->setDoubleData("track_michelProng_energy",michelProng.getDoubleData("energy"));
        event->setDoubleData("track_michelProng_time_diff",michelProng.getDoubleData("time_diff"));
        event->setDoubleData("track_michelProng_end_Z",michelProng.getDoubleData("edz"));
        event->setDoubleData("track_michelProng_begin_Z",michelProng.getDoubleData("bgz"));
    }
}

void CCProtonPi0::saveMichelProngTruth(Minerva::GenMinInteraction* truthEvent, Minerva::Prong &michelProng, bool isVertex, bool isLargeSearch) const
{
    SmartRefVector<Minerva::IDCluster> michelClusters = michelProng.getAllIDClusters();

    DigitVectorTruthInfo michelInfo;
    michelInfo.ParseTruth(michelClusters, fTrajectoryMap);

    if (isVertex){
        if (isLargeSearch){
            truthEvent->setIntData("vtx_michel_large_evis_most_pdg", michelInfo.GetMostEvisPdg());
            truthEvent->setDoubleData("vtx_michel_large_evis_total_truth", michelInfo.GetTotalTruthEnergy());
        }else{
            truthEvent->setIntData("vtx_michel_evis_most_pdg", michelInfo.GetMostEvisPdg());
            truthEvent->setDoubleData("vtx_michel_evis_total_truth", michelInfo.GetTotalTruthEnergy());
        }
    }else{
        truthEvent->setIntData("track_michel_evis_most_pdg", michelInfo.GetMostEvisPdg());
        truthEvent->setDoubleData("track_michel_evis_total_truth", michelInfo.GetTotalTruthEnergy());
    }
}

//---------------------------------------------------------------------------------------
// create short anchored tracks
//---------------------------------------------------------------------------------------
bool CCProtonPi0::createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex, bool make_primary_short_tracks ) const
{
    StatusCode sc;
    bool createdTracks = false;

    //-- get the time slice
    Minerva::TimeSlice* timeSlice = getTimeSlice( vertex );
    short slice = timeSlice ? timeSlice->sliceNumber() : -1;
    debug() << "   The slice associated with the event is = " << slice << endmsg;

    //-- anchored short track parameters
    double searchDist = 500.*CLHEP::mm;   //!-- maximum cluster search distance in a view  
    double startDist  = 500.*CLHEP::mm;   //!-- maximum cluster starting distance in a view  
    double gap        = 500.*CLHEP::mm;   //!-- maximum cluster allowed search gap in a view 
    double range3D    = 250.*CLHEP::mm;   //!-- maximum 3D distance between track endpoint and cluster 3D position
    double range2D    = 200.*CLHEP::mm;   //!-- maximum 2D distance between track endpoint and cluster 2D position
    bool   rawClus    = false;  //!-- used the default cluster list method for retreiving clusters
    bool   applyMS    = false;  //!-- apply the multiple scattering correction when fitting the track
    bool   force3D    = true;   //!-- enforced the x,y,z distance cuts

    //-- create a container to store tracks
    std::vector< Minerva::Track* > trackCont;
    trackCont.clear();

    //-- run the anchored short tracker
    bool can_make_anchor_tracks = true;
    while( can_make_anchor_tracks ) {
        Minerva::IDClusterVect idclusters = getClusters(event);
        debug() << "  The number of IDClusters = " << idclusters.size() << endmsg;

        if( idclusters.empty() ){ 
            can_make_anchor_tracks = false; 
            break; 
        }

        //-- creates short tracks
        std::vector<Minerva::Track*> *tracks = new std::vector<Minerva::Track*>; 
        sc = m_anchoredTracker->createAnchoredTracks(idclusters,vertex,tracks,rawClus,searchDist,startDist,gap,
                applyMS,force3D,range3D,range2D);
        debug() << "  The number of anchor tracks created = " << tracks->size() << endmsg; 

        if( !sc || tracks->empty() ){ 
            tracks->clear(); 
            delete tracks; 
            can_make_anchor_tracks = false; 
            break; 
        }

        //-- store tracks in the container
        for(std::vector<Minerva::Track*>::iterator trk = tracks->begin(); trk != tracks->end(); trk++) { 
            SmartRef<Minerva::Track> inTrack(*trk);
            addObject(event,(*trk));
            setTrackDirection(*trk,vertex);
            (*trk)->setTimeSlice(slice);
            (*trk)->setPatRecHistory( Minerva::Track::FourHitPatRec );
            (*trk)->setIntData( "createdTrack",1 );
            vertex->addOutgoingTrack(inTrack);

            Minerva::Vertex* stopPointVertex = new Minerva::Vertex();
            stopPointVertex->setIncomingTrack(inTrack);
            stopPointVertex->setPosition((*trk)->lastState().position());
            stopPointVertex->setTimeSlice(slice);
            stopPointVertex->setVertexFlags(Minerva::Vertex::StopPoint);
            addObject(event,stopPointVertex);
            createdTracks = true;
            debug() << "	The new stop vertex's position = " << stopPointVertex->position() << endmsg;

            trackCont.push_back( (*trk) );
        }
    } //-- end loop for making anchored tracks 

    //-- run the vertex energy short tracker
    bool can_make_vertex_tracks = true;
    while( can_make_vertex_tracks ) {
        Minerva::IDClusterVect idclusters = getClusters(event);
        debug() << "  The number of IDClusters = " << idclusters.size() << endmsg;

        if( idclusters.empty() ) { can_make_vertex_tracks = false; break; }

        std::vector<Minerva::Track*> *tracks = new std::vector<Minerva::Track*>; 
        Gaudi::XYZVector zDirectionVector(0, 0, 1); 
        Gaudi::XYZPoint vertexPosition = vertex->position();
        sc = m_vertexEnergyStudyTool->reconstructVertexActivity(event,vertexPosition,zDirectionVector,idclusters,tracks); 
        debug() << "  The number of vertex short tracks created = " << tracks->size() << endmsg;

        if( !sc || tracks->empty() ) { tracks->clear(); delete tracks; can_make_vertex_tracks = false; break; }

        //-- store tracks in the container
        for(std::vector<Minerva::Track*>::iterator trk = tracks->begin(); trk != tracks->end(); trk++) {
            SmartRef<Minerva::Track> inTrack(*trk);
            addObject(event,(*trk));
            setTrackDirection(*trk,vertex);
            (*trk)->setTimeSlice(slice);
            (*trk)->setPatRecHistory( Minerva::Track::NotAssigned1 );
            (*trk)->setIntData( "createdTrack",1 );
            vertex->addOutgoingTrack(inTrack);

            Minerva::Vertex* stopPointVertex = new Minerva::Vertex();
            stopPointVertex->setIncomingTrack(inTrack);
            stopPointVertex->setPosition((*trk)->lastState().position());
            stopPointVertex->setTimeSlice(slice);
            stopPointVertex->setVertexFlags(Minerva::Vertex::StopPoint);
            addObject(event,stopPointVertex);
            createdTracks = true;
            debug() << "	The new stop vertex's position = " << stopPointVertex->position() << endmsg;

            trackCont.push_back( (*trk) );
        }
    } //-- end loop for making vertex short tracks

    debug() << "	the number of tracks to add to store = " << trackCont.size() << endmsg;

    //-- add tracks to store and vertex object
    for(std::vector<Minerva::Track*>::iterator track = trackCont.begin(); track != trackCont.end(); track++) {

        //-- make primary prongs for the created primary tracks
        if( make_primary_short_tracks ) {

            //-- make new primary track prong 
            Minerva::Prong* newProng = new Minerva::Prong("CCProtonPi0PrimaryTrackProng");
            debug() << "  Made a new Prong with creation signature : " << newProng->creationSignature() << endmsg;

            //-- add track and source vertex to the prong
            (*track)->setHistory(Minerva::Track::Used);
            newProng->add((*track));
            newProng->setSourceVertex(vertex);

            //-- add prong to the store and promote to a primary prong
            addObject(event,newProng);
            event->promoteProngToPrimary(SmartRef<Minerva::Prong>(newProng));
        }
    }

    //-- refit vertex
    if( createdTracks ) m_vertexFitter->fit(vertex); 

    //-- store vertex for the primary short_tracks
    if( make_primary_short_tracks && createdTracks ) {

        std::vector<double> fit_vtx;
        fit_vtx.push_back( vertex->position().x() );
        fit_vtx.push_back( vertex->position().y() );
        fit_vtx.push_back( vertex->position().z() ); 

        event->setContainerDoubleData("fit_vtx",fit_vtx);

        SmartRef<Minerva::Vertex> vtx(vertex);
        event->setInteractionVertex( vtx );
        Minerva::ProngVect prongs = event->primaryProngs();
        for(unsigned int p = 0; p < prongs.size(); p++) prongs[p]->setSourceVertex(vertex);
    }

    debug() << "  The return vertex's position = " << vertex->position() << endmsg;   

    return createdTracks;
}

//------------------------------------------------------------------------------------
//  grab a list of clusters for creating a short track
//------------------------------------------------------------------------------------
Minerva::IDClusterVect CCProtonPi0::getClusters( Minerva::PhysicsEvent* event ) const
{
    Minerva::IDClusterVect retClusters; 
    retClusters.clear();

    StatusCode sc;

    //-- get interaction vertex
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

    //-- get the clusters
    Minerva::IDClusterVect idClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    if( idClusters.size() > 200 ){ 
        idClusters.clear(); 
        return retClusters; 
    }
    debug() << "  the number of idClusters before anchoring = " << idClusters.size() << endmsg;

    //-- create new anchored ID tracks to retrieve an axis for the cone
    std::vector<Minerva::Track*>* tracks = new std::vector<Minerva::Track*>;
    sc = m_anchoredTracker->createAnchoredTracks(idClusters,vertex,tracks,true,6000,800,1000,false,false,50,850);

    //-- sanity checks for the number of tracks
    debug() << "   the number of tracks = " << tracks->size() << ", nclusters = " << idClusters.size() << endmsg;
    if( !sc || tracks->empty() ){ 
        tracks->clear(); 
        delete tracks; 
        return retClusters; 
    }

    //-- set track direction
    std::vector<Minerva::Track*>::iterator trk;
    for(trk = tracks->begin(); trk != tracks->end(); trk++){ 
        setTrackDirection(*trk,vertex);
    }

    //-- retreive the track with the closest first state position to the vertex
    Minerva::Track* track = *tracks->begin();
    if( tracks->size() > 1 ) {
        double min_separation = 9999.9;
        for(trk = tracks->begin(); trk != tracks->end(); trk++) {
            double separation = MathTool->distanceBetween(vertex->position(),(*trk)->firstState().position());
            if( separation < min_separation ){ 
                track = *trk; 
                min_separation = separation; 
            }
        }
    }

    //-- store the track's clusters 
    Minerva::IDClusterVect clusTmp;
    Minerva::Track::NodeContainer::iterator node;
    for(node = track->nodes().begin(); node != track->nodes().end(); node++){
        clusTmp.push_back( (*node)->idcluster() );
    }

    //-- cone angle
    double coneAngle = 30*CLHEP::degree;

    //-- create cone
    Gaudi::XYZPoint firstState = track->firstState().position();
    double theta = track->theta();
    double phi   = track->phi();

    if( track->direction() == Minerva::Track::Backward ) {
        theta = Gaudi::Units::pi - theta;
        phi   = Gaudi::Units::twopi - phi;
    }

    Gaudi::Polar3DVector axis(1,theta,phi);
    Cone cone(firstState,axis,coneAngle);
    debug() << "    cone's firstState position = " << firstState << ", axis = " << axis << endmsg;  

    //-- clean up tracks
    for( trk = tracks->begin(); trk != tracks->end(); trk++){ 
        (*trk)->reset(); 
        delete *trk;  
    }
    tracks->clear();
    delete tracks;

    //-- get fill container with clusters inside of cone
    Minerva::IDClusterVect clusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    debug() << "   the number of cone candidate clusters = " << clusters.size() << endmsg;

    for(unsigned int clus = 0; clus < clusters.size(); clus++) {
        if( m_coneUtilsTool->isInsideCone(clusters[clus],cone) ) retClusters.push_back( clusters[clus] ); 
        else { 
            for(unsigned int tmp = 0; tmp < clusTmp.size(); tmp++){ 
                if( clusTmp[tmp] == clusters[clus] ) retClusters.push_back( clusters[clus] );
            }
        }
    }

    clusters.clear();
    clusTmp.clear();

    debug() << "   the number of clusters return from coning = " << retClusters.size() << endmsg;

    return retClusters;
}

//----------------------------------------------------------------------------------
// set the created track direction
//----------------------------------------------------------------------------------
void CCProtonPi0::setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const
{
    int backward = 0, forward = 0;
    Minerva::Track::NodeContainer nodes = track->nodes();
    for(unsigned int node = 0; node < nodes.size(); node++) {
        if( nodes[node]->z() < vertex->position().z() )      backward++;
        else if( nodes[node]->z() > vertex->position().z() ) forward++;
    }
    if( backward > forward )      track->setDirection(Minerva::Track::Backward);
    else if( forward > backward ) track->setDirection(Minerva::Track::Forward);
    return;
}


bool CCProtonPi0::ShouldReconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    //--------------------------------------------------------------------------
    // Do NOT Reconstruct the Event -
    //      if the ReconstructEvent is False - tagTruth() Decides This
    //--------------------------------------------------------------------------    
    if( truthEvent ){
        if ( !(truthEvent->filtertaglist()->isFilterTagTrue("ReconstructEvent")) ){
            debug() << "TagTruth() Marked Event as Do NOT Reconstruct!" <<endmsg;
            return false; 
        }
    }

    //--------------------------------------------------------------------------
    // Do NOT Analyze the Event - if Event has a Bad Object
    //--------------------------------------------------------------------------
    if( event->filtertaglist()->isFilterTagTrue( AnaFilterTags::BadObject() ) ) { 
        warning() << "Found an event flagged with a BadObject!" << endmsg;
        return false; 
    }

    return true;
}

// MAKE CUT - if NO or NULL Interaction Vertex
bool CCProtonPi0::hasEventVertex(Minerva::PhysicsEvent * event) const
{
    if( (event->hasInteractionVertex()) && event->interactionVertex() != NULL ){
        return true;        
    }else{
        debug() << "The event has NO or NULL vertex!" << endmsg;
        event->setIntData("Cut_Vertex_None",1);
        return false;
    }
}

/*
   MAKE CUT - if Interaction Vertex is NOT in Reconstructable Volume
   "Vertex is NOT in Reconstructable Volume" means we can not run vertex-anchored
   short tracker with a meaningful result outside of that volume.
   Only the events that pass this CUT are used in vertex-anchored short tracker
   */
bool CCProtonPi0::vertexInRecoVolume(Minerva::PhysicsEvent *event) const
{
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    Gaudi::XYZPoint vtx_position = vertex->position();

    if ( FiducialPointTool->isFiducial(vtx_position, m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ)){
        return true;
    }else{
        debug() <<"Interaction Vertex is NOT in reconstructable volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_Vertex_Not_Reconstructable",1);
        return false;
    }
}

bool CCProtonPi0::vertexInFiducialVolume(Minerva::PhysicsEvent *event) const
{
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    Gaudi::XYZPoint vtx_position = vertex->position();

    if( FiducialPointTool->isFiducial( vtx_position, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ){
        return true;
    }else{
        debug() <<"Interaction Vertex is NOT in fiducial volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_Vertex_Not_Fiducial",1);
        return false;
    }
}

void CCProtonPi0::RefitVertex_Using_AnchoredShortTracks(Minerva::PhysicsEvent *event) const
{
    // Get Current Vertex
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    bool make_primary_short_tracks = true;

    debug()<<"Vertex Position Before Fit = "<<vertex->position()<<endmsg;

    // Create Anchored Short Tracks
    bool createdTracks = createdAnchoredShortTracks(event, vertex, make_primary_short_tracks);
    if (createdTracks){
        debug()<<"Succesfully created Anchored Short Tracks"<<endmsg;
    }
    event->setIntData("nTracks", (int)vertex->getOutgoingTracks().size() );
    debug()<<"The vertex has " << vertex->getOutgoingTracks().size()  << " outgoing primary tracks!" << endmsg;

    // Set vertex fit information
    double fit_chi2 = -1;
    if( vertex->hasDoubleData("fit_chi2") ) {
        fit_chi2 = vertex->getDoubleData("fit_chi2"); 
        event->setDoubleData("vtx_fit_chi2",fit_chi2);
        debug()<<"\tThe vertex fit chi2 = " << fit_chi2 << endmsg;  
    }

    int fit_converged = -1;
    if( vertex->hasIntData("fit_converged") ) {
        fit_converged = vertex->getIntData("fit_converged");
        event->setIntData("vtx_fit_converged",fit_converged);
        debug()<<"\tThe vertex fit convergence = " << fit_converged << endmsg;
    }

    // re-set the vertex position for vertex fit failures for long-long tracks combination
    if( fit_converged == 0 ) {
        unsigned int longtracker = 0;
        SmartRefVector<Minerva::Track> trks = vertex->getOutgoingTracks();

        for(unsigned int t = 0; t < trks.size(); t++) {
            Minerva::Track::PatRecHistory pat = trks[t]->patRecHistory();
            if( pat == Minerva::Track::LongPatRec3View || pat == Minerva::Track::LongPatRec2View ) longtracker++;
        }

        if( longtracker == trks.size() ) {
            m_vertexFitter->fit(vertex);
            event->setInteractionVertex(vertex);
            Minerva::ProngVect prongs = event->primaryProngs();
            for(unsigned int p = 0; p < prongs.size(); p++){
                prongs[p]->setSourceVertex(vertex);
            }
        }
    }

    std::vector<double> fit_vtx;
    fit_vtx.push_back( vertex->position().x() );
    fit_vtx.push_back( vertex->position().y() );
    fit_vtx.push_back( vertex->position().z() );

    event->setContainerDoubleData("fit_vtx",fit_vtx);

    debug()<<"Vertex Position After Fit = "<<vertex->position()<<endmsg;
}

void CCProtonPi0::SetVertexCount(Minerva::PhysicsEvent * event) const
{
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();

    // Get Vertex Count
    SmartRefVector<Minerva::Vertex> allVertices  = event->select<Minerva::Vertex>( "All","StartPoint" );
    debug()<<"N(vertices) = "<<allVertices.size()<<endmsg;
    event->setIntData("vtx_total_count",allVertices.size());

    SmartRefVector<Minerva::Vertex>::iterator iter_primary_vtx = allVertices.end();
    SmartRefVector<Minerva::Vertex>::iterator iter_vtx;
    for ( iter_vtx = allVertices.begin(); iter_vtx != allVertices.end(); ++iter_vtx ){
        if (*iter_vtx == vertex) {
            iter_primary_vtx = iter_vtx;
            break;
        }
    }

    event->setIntData("vtx_primary_index", std::distance(allVertices.begin(),iter_primary_vtx));
    if (iter_primary_vtx != allVertices.end()) allVertices.erase(iter_primary_vtx);
    event->setIntData("vtx_secondary_count", allVertices.size()); /* Number of secondary vertices */
}

bool CCProtonPi0::hasEventMinosMatchedMuon(Minerva::PhysicsEvent *event) const
{
    bool foundMuon = MuonUtils->findMuonProng( event, m_MuonProng, m_MuonParticle );

    if ( !foundMuon ){
        debug() << "Did not find a muon prong!" << endmsg;
        event->setIntData("Cut_Muon_None",1);
        return false;
    }else if(m_MuonProng == NULL || m_MuonParticle == NULL){
        debug()<<"Muon Prong or Muon Particle is NULL!"<<endmsg;
        event->setIntData("Cut_Muon_None",1);
        return false;
    }else{
        debug() << "m_MuonProng and m_MuonParticle is Saved!" << endmsg;
        return true; 
    }
}

// MAKE CUT - if Muon has positive charge (AntiMuon)
bool CCProtonPi0::isMuonChargeNegative(Minerva::PhysicsEvent *event) const
{
    int charge = -99;
    MuonUtils->muonCharge(m_MuonProng,charge);
    if(charge == 1){
        debug()<<"AntiMuon Contamination"<<endmsg;
        event->setIntData("Cut_Muon_Charge",1);
        return false;
    }else{
        return true;
    }
}

void CCProtonPi0::tagPrimaryMuon(Minerva::PhysicsEvent *event) const
{
    m_MuonProng->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
    m_MuonParticle->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );

    m_hitTagger->applyColorTag(m_MuonProng, m_Color_muonProng);
    event->setTime( m_recoTimeTool->prongBestTime(m_MuonProng) );
}

bool CCProtonPi0::ImprovedMichel_EventHasMichel(Minerva::PhysicsEvent *event) const
{
    debug()<<"START: Searching Events with ImprovedMichelTool"<<endmsg;

    m_michelTool_Improved->Reset();
    m_michelTool_Improved->FindMichel(event->time() + 200.0);
    m_michelTool_Improved->ApplyQualityCuts();

    int nMichel = m_michelTool_Improved->GetMichelCount(); 
    bool foundMichel = nMichel > 0;

    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMichel", true);
    }

    debug()<<"nMichel = "<<nMichel<<endmsg;

    return foundMichel;
}

bool CCProtonPi0::VertexHasMichels(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent) const
{
    Minerva::Prong vtx_michel_prong;
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();

    bool foundMichel = m_michelTool_Default->findMichel( vertex, vtx_michel_prong );

    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        if (m_study_michel_electron){
            saveMichelProngToNTuple(event, vtx_michel_prong, true, false);
            if(truthEvent) saveMichelProngTruth(truthEvent, vtx_michel_prong, true, false);
        }
        event->setIntData("Cut_Vertex_Michel_Exist",1);
    }

    // Search For Michel Electron in a Larger Search Volume
    Minerva::Prong vtx_michel_prong_large;
    bool foundMichel_LargeBox = m_michelTool_Large->findMichel( vertex, vtx_michel_prong_large);

    if (foundMichel_LargeBox){
        debug()<<"Found a Michel Electron with Large Search Volume"<<endmsg;
        if (m_study_michel_electron){
            saveMichelProngToNTuple(event, vtx_michel_prong_large, true, true);
            if(truthEvent) saveMichelProngTruth(truthEvent, vtx_michel_prong_large, true, true);
        }
        event->setIntData("Cut_Vertex_Large_Michel_Exist",1);
    }

    bool foundMichel_Improved = ImprovedMichel_HasMichelAtPoint(vertex);
    event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_VertexHasMichel", foundMichel_Improved);
    if (foundMichel_Improved){
        event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMatchedMichel", true);
    }

    // Return value depends on Default Michel Search
    if (foundMichel) return true;
    else return false;
}

bool CCProtonPi0::TrackEndPointHasMichels(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent) const
{
    // Get All Prongs to search for end point michels
    ProngVect primaryProngs = event->primaryProngs();
    bool foundMichel_EndPoint = false;
    bool foundMichel_secEndPoint = false;
    bool foundMichel_Improved_EndPoint = false;
    bool foundMichel_Improved_secEndPoint = false;

    for (ProngVect::iterator iterProngs = primaryProngs.begin(); iterProngs != primaryProngs.end(); ++iterProngs) {

        // Skip Muon Prong
        if ( (*iterProngs) == m_MuonProng ) {
            debug()<<"Muon prong skipping!"<<endmsg;
            continue;
        }

        //Make an endpoint vertex blob for the hadron prong - currently this depends on the prong order!!!
        debug()<<"Make endpoint blob"<<endmsg;
        SmartRef<Minerva::Vertex> endpoint_vtx;
        m_objectAssociator->getVertex_fromTrackBack( endpoint_vtx, (*iterProngs)->minervaTracks().back() );
        if (!endpoint_vtx) {
            warning()<<"Could not find a back vertex for this prong!"<<endmsg;
            continue;
        }        

        IDClusterVect StopPointClusters = event->select<Minerva::IDCluster>("Unused","!XTalkCandidate");
        if (!StopPointClusters.empty()) {
            Minerva::IDBlob* StopPointBlob = new Minerva::IDBlob;
            debug()<<"CreateIDBlobs"<<endmsg;
            m_stopPointBlobTool->createIDBlobs(StopPointClusters, StopPointBlob, endpoint_vtx);
            m_hitTagger->applyColorTag(StopPointBlob, m_Color_endPointVertex);
            addObject(event, StopPointBlob);        
            endpoint_vtx->addIDBlob(StopPointBlob);
            (*iterProngs)->add(StopPointBlob);
            (*iterProngs)->setDoubleData("endpointE", m_caloUtils->applyCalConsts(StopPointBlob, "Default") );
        }

        //Search for endpoint michels
        if (!foundMichel_Improved_EndPoint){
            foundMichel_Improved_EndPoint = ImprovedMichel_HasMichelAtPoint(endpoint_vtx);
            event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EndPointHasMichel", foundMichel_Improved_EndPoint);
            if (foundMichel_Improved_EndPoint){
                event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMatchedMichel", true);
            }
        }

        // Skip Track if we already found the Michel
        if (foundMichel_EndPoint) continue;

        Minerva::Prong michelProng;
        foundMichel_EndPoint = m_michelTool_Default->findMichel( endpoint_vtx, michelProng );

        if (foundMichel_EndPoint){
            debug()<<"Found an End Point Michel Electron!"<<endmsg;
            if (m_study_michel_electron){
                saveMichelProngToNTuple(event, michelProng, false, false);
                if(truthEvent) saveMichelProngTruth(truthEvent, michelProng, false, false);
            }
            event->setIntData("Cut_EndPoint_Michel_Exist",1);
        }

        // Search for secondary Michels
        if ( (*iterProngs)->minervaTracks().size() > 1) {
            debug()<<"Searching for secondary end point michels..."<<endmsg;
            for (SmartRefVector<Minerva::Track>::iterator iterTracks = (*iterProngs)->minervaTracks().begin(); iterTracks != (*iterProngs)->minervaTracks().end() - 1; ++iterTracks) {
                SmartRef<Minerva::Vertex> current_end_vtx;
                m_objectAssociator->getVertex_fromTrackBack( current_end_vtx, *iterTracks );

                // Search with MichelToolImproved
                if (!foundMichel_Improved_secEndPoint){
                    foundMichel_Improved_secEndPoint = ImprovedMichel_HasMichelAtPoint(current_end_vtx);
                    event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_secEndPointHasMichel", foundMichel_Improved_secEndPoint);
                    if (foundMichel_Improved_secEndPoint){
                        event->filtertaglist()->setOrAddFilterTag( "ImprovedMichel_EventHasMatchedMichel", true);
                    }
                }

                // Skip Secondary Track if we already found the Michel
                if (foundMichel_secEndPoint) continue;

                foundMichel_secEndPoint = m_michelTool_Default->findMichel( current_end_vtx, michelProng);
                if (foundMichel_secEndPoint) {
                    debug()<<"Found a Secondary End Point Michel Electron!"<<endmsg;
                    if (m_study_michel_electron){
                        saveMichelProngToNTuple(event, michelProng, false, false);
                        if(truthEvent) saveMichelProngTruth(truthEvent, michelProng, false, false);
                    }
                    event->setIntData("Cut_secEndPoint_Michel_Exist",1);
                }
            }
        }
    }//end loop over primary prongs

    return foundMichel_EndPoint || foundMichel_secEndPoint;
}

bool CCProtonPi0::ImprovedMichel_HasMichelAtPoint(SmartRef<Minerva::Vertex> point) const
{
    const Gaudi::XYZPoint pos = point->position(); 
    double x = pos.x(); 
    double y = pos.y(); 
    double z = pos.z(); 
    return ImprovedMichel_HasMichelAtPoint(x, y, z);
}

bool CCProtonPi0::ImprovedMichel_HasMichelAtPoint(double x, double y, double z) const
{
    const double max_distance = 500.0; // in mm
    return m_michelTool_Improved->Match(x, y, z, max_distance, max_distance);
}

void CCProtonPi0::ColorUnusedIDClusters(Minerva::PhysicsEvent *event) const
{
    SmartRefVector<IDCluster> idClusters = event->select<IDCluster>("Used:Unused","!XTalkCandidate");

    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_c = idClusters.begin(); iter_c != idClusters.end(); ++iter_c) 
    {
        if ( (int)(*iter_c)->history() == Minerva::IDCluster::Unused ){
            m_hitTagger->applyColorTag((*iter_c), m_Color_clusterUnused);
        }
    }
}


void CCProtonPi0::SaveEventTime(Minerva::PhysicsEvent *event) const
{

    // Calculate dead time
    m_getDeadTimeTool->setDeadTimeTool( getDeadTimeTable() );

    // get dead time data
    std::map<const char*, int> intDeadData;
    std::map<const char*, double> doubleDeadData;
    DeadTimeTool->totalDead( event, intDeadData, doubleDeadData );

    debug() << "Size of dead time data: " << intDeadData.size() << endmsg;
    int dead = 0, udead = 0, ddead = 0, tdead = 0;
    for( std::map<const char*, int>::iterator it_map = intDeadData.begin(); it_map != intDeadData.end(); it_map++ ) {
        debug() << "  " << (*it_map).first << ": " << (*it_map).second << endmsg;
        if((std::string)(*it_map).first=="UDEAD")
            udead += (*it_map).second;
        if((std::string)(*it_map).first=="DDEAD")
            ddead += (*it_map).second;
        if((std::string)(*it_map).first=="TDEAD")
            tdead += (*it_map).second;
        if((std::string)(*it_map).first=="DEAD")
            dead += (*it_map).second; 
    }

    double event_time = event->time();

    event->setDoubleData("time", event_time);
    event->setIntData( "dead", dead );
    event->setIntData( "udead", udead );
    event->setIntData( "ddead", ddead );
    event->setIntData( "tdead", tdead );

}


void CCProtonPi0::GetMuonExtraEnergy(Minerva::PhysicsEvent *event) const
{
    const double extraEnergyCylinderUpstreamLength = -200.0*CLHEP::mm;
    const double extraEnergyCylinderDownstreamLength = 50.0*CLHEP::mm;
    const double extraEnergyCylinderRadius = 50*CLHEP::mm;
    const double extraEnergyLowerTimeWindow = 25.0*CLHEP::ns;
    const double extraEnergyUpperTimeWindow = 25.0*CLHEP::ns;
    const double extraEnergyPECut = 0.0;

    //  -- Muon Blob
    Minerva::IDBlob *muonBlob = new Minerva::IDBlob;	
    //SmartRef<Minerva::IDBlob> muonBlob;	
    double muon_blob_energy = 0;

    SmartRefVector<Minerva::IDCluster> extraMuonClusters
        = m_extraEnergyTool->getExtraIDClustersCylinder( m_MuonProng,
                extraEnergyCylinderUpstreamLength,
                extraEnergyCylinderDownstreamLength,
                extraEnergyCylinderRadius,
                IExtraEnergyTool::k_all,
                IExtraEnergyTool::k_inside,
                extraEnergyLowerTimeWindow,
                extraEnergyUpperTimeWindow,
                extraEnergyPECut,
                IExtraEnergyTool::k_clusters_unused,
                -100 ); /* -100 means clusters in all modules! */

    debug() << " Extra muon contain " << extraMuonClusters.size() << " clusters" << endmsg;
    std::vector<BlobSeeds::seedCandVector> VectorSeeds;
    m_blobSeedingTool->makeViewSortedBlobSeeds( extraMuonClusters, VectorSeeds );
    std::vector<BlobSeeds::seedCandVector>::iterator itSeeds = VectorSeeds.begin();
    for ( ; itSeeds != VectorSeeds.end(); itSeeds++ ){
        BlobSeeds::seedCandVector::iterator itSeed = (*itSeeds).begin();
        for ( ; itSeed != (*itSeeds).end(); itSeed++ ){
            debug() << " Seed with Energy " << (*itSeed)->seedenergy << " " << ((*itSeed)->IdClustVec).size()
                <<  endmsg;
            if ( (*itSeed)->seedenergy < 30 )	{ 
                muonBlob->add( (*itSeed)->IdClustVec );
                debug() << " Blob size " << muonBlob->nclusters() << endmsg; 
            }
        }
    }

    if ( muonBlob->nclusters() ) {
        muonBlob->setHistory( Minerva::IDBlob::Unused );
        muonBlob->setPatRecHistory( Minerva::IDBlob::IsolatedIDBlobPatRec );
        double dummy_trkrevis = 0.0;
        double dummy_ecalevis = 0.0;
        double dummy_hcalevis = 0.0;
        double dummy_scalevis = 0.0;
        m_idHoughBlob->getBlobEnergyTime_Old( muonBlob, muon_blob_energy,
                dummy_trkrevis, dummy_ecalevis, dummy_hcalevis, dummy_scalevis);

        debug()<< "Adding Muon blob with " << muonBlob->nclusters() << " clusters; energy = "  << muon_blob_energy
            << endmsg;
        addObject( event, muonBlob );
    }

    event->setDoubleData( "Extra_Energy_Muon", muon_blob_energy );

    DigitVectorTruthInfo info;
    info.ParseTruth(muonBlob->clusters(),fTrajectoryMap);
    double evis = info.GetEdepByPdg(111);
    event->setDoubleData("pi0_evis_muon_blob", evis);

}

void CCProtonPi0::FillTrajectoryMap() const
{
    /* The trajectory map must be emptied for every new physics event */
    if (fTrajectoryMap.empty()) {
        Minerva::TG4Trajectories* trajectories(NULL);
        if (exist<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default)) {
            trajectories = get<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default);
        }

        if (!trajectories) return;
        debug() << "Dump TrajectoryInfo (pdg, trackId, parentId)" << endmsg;
        for (Minerva::TG4Trajectories::iterator t = trajectories->begin();
                t != trajectories->end(); ++t) {
            Minerva::TG4Trajectory* traj = *t;
            debug() << "\tTrajectoryInfo"
                << std::setw(15) << traj->GetPDGCode()
                << std::setw(10) << traj->GetTrackId()
                << std::setw(10) << traj->GetParentId()
                << endmsg;
            fTrajectoryMap[traj->GetTrackId()] = traj;
        }
    }

    if (fTrajectoryMap.empty()) {
        warning() << "TrajectoryMap empty " << endmsg;
        return;
    }
}

void CCProtonPi0::DiscardFarTracks(Minerva::PhysicsEvent *event) const
{
    // ------------------------------------------------------------------------
    // Tracks on Primary Vertex
    // ------------------------------------------------------------------------
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   
    SmartRefVector<Minerva::Track> tracks = vertex->getTracks();

    // Get Muon Track
    SmartRef<Minerva::Track> muonTrack = m_MuonProng->minervaTracks().front();

    std::size_t trackMultiplicity = 1; // We already have Muon Track
    SmartRefVector<Minerva::Track> farTracks;
    for (SmartRefVector<Minerva::Track>::iterator t = tracks.begin(); t != tracks.end(); ++t) {
        if (muonTrack && (*t == muonTrack)) {
            debug()<<"\tSkipping the muon track"<<endmsg;
            continue;
        }

        /* Not sure which end of the track is closest to the primary vertex.
           Take the shorter of the distances to both track ends */
        Gaudi::XYZPoint firstPos = (*t)->firstNode()->position();
        Gaudi::XYZPoint lastPos  = (*t)->lastNode()->position();

        const double distance = std::min((firstPos-vertex->position()).R(), (lastPos-vertex->position()).R());
        debug()<<"\tDistance to primary vertex: "<<distance<<endmsg;

        /* Count tracks close to the primary vertex and use it
           as an updated track multiplicity */
        if (distance < 50.0) ++trackMultiplicity;
        if (distance > 50.0) farTracks.push_back(*t);
    }

    event->setIntData("nTracks_Close", trackMultiplicity);
    event->setIntData("nTracks_Far",farTracks.size());

    // ------------------------------------------------------------------------
    // Tracks on Secondary Vertices
    // ------------------------------------------------------------------------
    const Gaudi::XYZPoint& pos = vertex->position();

    // Get All Vertices in Event
    SmartRefVector<Minerva::Vertex> allVertices  = event->select<Minerva::Vertex>( "All","StartPoint" );

    // Remove Primary Vertex from All Vertices
    SmartRefVector<Minerva::Vertex>::iterator primaryVertexIterator = allVertices.end();
    for (SmartRefVector<Minerva::Vertex>::iterator vtx = allVertices.begin(); vtx != allVertices.end(); ++vtx ){
        if (*vtx == vertex) {
            primaryVertexIterator = vtx;
            break;
        }
    }
    if (primaryVertexIterator != allVertices.end()) allVertices.erase(primaryVertexIterator);

    SmartRefVector<Minerva::Track> tracksToClear;
    std::vector<int> multiplicities;
    // Loop over Secondary Vertices -- Primary is removed from allVertices 
    for (SmartRefVector<Minerva::Vertex>::iterator vtx = allVertices.begin(); vtx != allVertices.end(); ++vtx ){

        SmartRefVector<Minerva::Track> tracks = (*vtx)->getTracks();

        /* For unclear reason, Jose projects and discards later only
           the first track (outgoing from secondary vertices) */
        Gaudi::XYZPoint projectedVertex;
        SmartRef<Minerva::Track> firstTrack = tracks.front();
        m_trackPropagator->propagate(firstTrack, pos.z(), projectedVertex);
        double deviation = (pos-projectedVertex).R();
        if ( deviation < 50 ) tracksToClear.push_back(firstTrack);

        multiplicities.push_back(tracks.size());
    }

    /* Put them in the same container before deleting */
    std::copy(tracksToClear.begin(),tracksToClear.end(),
            std::back_inserter(farTracks));

    event->setIntData("nTracks_Discarded", farTracks.size());
    event->setContainerIntData("nTracks_Secondary_Vtx", multiplicities);

    if (!farTracks.empty()) {
        debug()<<"Discarding "<<farTracks.size()<<" tracks, to get photons "<<endmsg;

        // Before 
        std::size_t oldsize = farTracks.size();
        SmartRefVector<Minerva::IDCluster> old_clusters; 
        old_clusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

        Minerva::EventMgr *mgr = getEventMgr(allVertices.front());
        discardObject(mgr,farTracks);

        // After
        std::size_t newsize = farTracks.size();
        SmartRefVector<Minerva::IDCluster> new_clusters; 
        new_clusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

        debug()<<"old unusedClusters = "<<old_clusters.size()<<" new unusedClusters = "<<new_clusters.size()<<endmsg;
        debug()<<"old farTracks = "<<oldsize<<" new farTracks = "<<newsize<<endmsg; 
    }
}

void CCProtonPi0::SaveTruthUnusedClusterEnergyInsideDetector(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> ecalClusters, SmartRefVector<Minerva::IDCluster> hcalClusters, SmartRefVector<Minerva::IDCluster> otherClusters) const
{
    // Unused Visible Energy inside ECAL
    DigitVectorTruthInfo ecalInfo;
    ecalInfo.ParseTruth(ecalClusters, fTrajectoryMap, 0);

    truthEvent->setDoubleData("ecal_unused_evis_total_norm", ecalInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("ecal_unused_evis_total_truth", ecalInfo.GetTotalTruthEnergy());
    truthEvent->setDoubleData("ecal_unused_evis_pizero", ecalInfo.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("ecal_unused_evis_piplus", ecalInfo.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("ecal_unused_evis_piminus", ecalInfo.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("ecal_unused_evis_muon", ecalInfo.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("ecal_unused_evis_proton", ecalInfo.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("ecal_unused_evis_neutron", ecalInfo.GetEdepByPdg(PDG::neutron));

    // Unused Visible Energy inside HCAL
    DigitVectorTruthInfo hcalInfo;
    hcalInfo.ParseTruth(hcalClusters, fTrajectoryMap, 0);

    truthEvent->setDoubleData("hcal_unused_evis_total_norm", hcalInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("hcal_unused_evis_total_truth", hcalInfo.GetTotalTruthEnergy());
    truthEvent->setDoubleData("hcal_unused_evis_pizero", hcalInfo.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("hcal_unused_evis_piplus", hcalInfo.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("hcal_unused_evis_piminus", hcalInfo.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("hcal_unused_evis_muon", hcalInfo.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("hcal_unused_evis_proton", hcalInfo.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("hcal_unused_evis_neutron", hcalInfo.GetEdepByPdg(PDG::neutron));

    // Unused Visible Energy inside Tracker + ECAL + HCAL
    DigitVectorTruthInfo otherInfo;
    otherInfo.ParseTruth(otherClusters, fTrajectoryMap, 0);

    truthEvent->setDoubleData("other_unused_evis_total_norm", otherInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("other_unused_evis_total_truth", otherInfo.GetTotalTruthEnergy());
    truthEvent->setDoubleData("other_unused_evis_pizero", otherInfo.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("other_unused_evis_piplus", otherInfo.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("other_unused_evis_piminus", otherInfo.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("other_unused_evis_muon", otherInfo.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("other_unused_evis_proton", otherInfo.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("other_unused_evis_neutron", otherInfo.GetEdepByPdg(PDG::neutron));
}

void CCProtonPi0::SaveTruthClusterEnergy_FoundBlobs(Minerva::PhysicsEvent* event, Minerva::GenMinInteraction *truthEvent) const
{
    debug()<<" 1 ProcessName"<<endmsg;
    // Get Blob 1 Clusters
    SmartRefVector<Minerva::IDCluster> blob1Clusters = m_Pi0Blob1->clusters();
    DigitVectorTruthInfo blob1Info;
    blob1Info.ParseTruth(blob1Clusters, fTrajectoryMap, 22);

    truthEvent->setIntData("blob1_evis_most_pdg",       blob1Info.GetMostEvisPdg());
    truthEvent->setIntData("blob1_evis_most_TrackID",   blob1Info.GetMostEvisId());
    truthEvent->setIntData("blob1_evis_parent_pdg",     blob1Info.GetParentEvisPdg());
    truthEvent->setIntData("blob1_evis_parent_TrackID", blob1Info.GetParentEvisId());
    truthEvent->setDoubleData("blob1_evis_total_norm",  blob1Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("blob1_evis_total_truth", blob1Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("blob1_evis_pizero",      blob1Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("blob1_evis_piplus",      blob1Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("blob1_evis_piminus",     blob1Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("blob1_evis_muon",        blob1Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("blob1_evis_proton",      blob1Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("blob1_evis_neutron",     blob1Info.GetEdepByPdg(PDG::neutron));

    // Get Blob 2 Clusters
    debug()<<"Gamma 2 ProcessName"<<endmsg;
    SmartRefVector<Minerva::IDCluster> blob2Clusters = m_Pi0Blob2->clusters();
    DigitVectorTruthInfo blob2Info;
    blob2Info.ParseTruth(blob2Clusters, fTrajectoryMap, 22);

    truthEvent->setIntData("blob2_evis_most_pdg",       blob2Info.GetMostEvisPdg());
    truthEvent->setIntData("blob2_evis_most_TrackID",   blob2Info.GetMostEvisId());
    truthEvent->setIntData("blob2_evis_parent_pdg",     blob2Info.GetParentEvisPdg());
    truthEvent->setIntData("blob2_evis_parent_TrackID", blob2Info.GetParentEvisId());
    truthEvent->setDoubleData("blob2_evis_total_norm",  blob2Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("blob2_evis_total_truth", blob2Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("blob2_evis_pizero",      blob2Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("blob2_evis_piplus",      blob2Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("blob2_evis_piminus",     blob2Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("blob2_evis_muon",        blob2Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("blob2_evis_proton",      blob2Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("blob2_evis_neutron",     blob2Info.GetEdepByPdg(PDG::neutron));

    // Total Captured Visible Energy ==> Blob1 + Blob2
    double total_captured_evis_total_norm = blob1Info.GetTotalNormEnergy() + blob2Info.GetTotalNormEnergy();
    double total_captured_evis_total_truth = blob1Info.GetTotalTruthEnergy() + blob2Info.GetTotalTruthEnergy();
    double total_captured_evis_pizero = blob1Info.GetEdepByPdg(PDG::pi0) + blob2Info.GetEdepByPdg(PDG::pi0);

    truthEvent->setDoubleData("total_captured_evis_total_norm", total_captured_evis_total_norm);
    truthEvent->setDoubleData("total_captured_evis_total_truth", total_captured_evis_total_truth);
    truthEvent->setDoubleData("total_captured_evis_pizero", total_captured_evis_pizero);

    // Total Visible Energy for Pi0 from ALL Clusters
    SmartRefVector<Minerva::IDCluster> allClusters = event->select<Minerva::IDCluster>("All","!LowActivity&!XTalkCandidate");
    DigitVectorTruthInfo allClustersInfo;
    allClustersInfo.ParseTruth(allClusters, fTrajectoryMap, 0);

    truthEvent->setDoubleData("allClusters_evis_pizero", allClustersInfo.GetEdepByPdg(PDG::pi0));
}



void CCProtonPi0::SaveTruthUnusedClusterEnergy_NearVertex(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> vertexClusters) const
{
    DigitVectorTruthInfo vertexInfo;
    vertexInfo.ParseTruth(vertexClusters, fTrajectoryMap, 0);

    truthEvent->setIntData("vertex_unused_evis_most_pdg", vertexInfo.GetMostEvisPdg());
    truthEvent->setDoubleData("vertex_unused_evis_total_norm", vertexInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("vertex_unused_evis_total_truth", vertexInfo.GetTotalTruthEnergy());
    truthEvent->setDoubleData("vertex_unused_evis_pizero", vertexInfo.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("vertex_unused_evis_piplus", vertexInfo.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("vertex_unused_evis_piminus", vertexInfo.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("vertex_unused_evis_muon", vertexInfo.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("vertex_unused_evis_proton", vertexInfo.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("vertex_unused_evis_neutron", vertexInfo.GetEdepByPdg(PDG::neutron));
    truthEvent->setDoubleData("vertex_unused_evis_gamma", vertexInfo.GetEdepByPdg(PDG::gamma));
}

void CCProtonPi0::SaveTruthUnusedClusterEnergy_Rejected(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> clusters) const
{
    DigitVectorTruthInfo RejectedInfo;
    RejectedInfo.ParseTruth(clusters, fTrajectoryMap, 0);

    truthEvent->setIntData("Rejected_unused_evis_most_pdg", RejectedInfo.GetMostEvisPdg());
    truthEvent->setDoubleData("Rejected_unused_evis_total_norm", RejectedInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("Rejected_unused_evis_total_truth", RejectedInfo.GetTotalTruthEnergy());
}

void CCProtonPi0::SaveTruthUnusedClusterEnergy_Dispersed(Minerva::GenMinInteraction *truthEvent, SmartRefVector<Minerva::IDCluster> clusters) const
{
    DigitVectorTruthInfo dispersedInfo;
    dispersedInfo.ParseTruth(clusters,fTrajectoryMap, 0);

    truthEvent->setIntData("dispersed_unused_evis_most_pdg", dispersedInfo.GetMostEvisPdg());
    truthEvent->setDoubleData("dispersed_unused_evis_total_norm", dispersedInfo.GetTotalNormEnergy());
    truthEvent->setDoubleData("dispersed_unused_evis_total_truth", dispersedInfo.GetTotalTruthEnergy());
    truthEvent->setDoubleData("dispersed_unused_evis_pizero", dispersedInfo.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("dispersed_unused_evis_piplus", dispersedInfo.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("dispersed_unused_evis_piminus", dispersedInfo.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("dispersed_unused_evis_muon", dispersedInfo.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("dispersed_unused_evis_proton", dispersedInfo.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("dispersed_unused_evis_neutron", dispersedInfo.GetEdepByPdg(PDG::neutron));
    truthEvent->setDoubleData("dispersed_unused_evis_gamma", dispersedInfo.GetEdepByPdg(PDG::gamma));
}

void CCProtonPi0::Check_Quality_of_Blobs(std::vector<Minerva::IDBlob*> &foundBlobs, std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent) const
{

  info() << " Checking the Quality of the found Blobs... " << foundBlobs.size() << endmsg;
  std::vector<int> shower_ndigits;
  std::vector<int> shower_nclusters;
  std::vector<double> shower_energy;
  std::vector<double> shower_dist_startpoint;
  std::vector<int> shower_isGood;
  std::vector<int> shower_pdg;
  std::vector<double> shower_score;

  const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

  bool isGood;
  int i = 0;
  std::vector<Minerva::IDBlob*>::iterator blob_iterator;
  for (blob_iterator = foundBlobs.begin(); blob_iterator != foundBlobs.end(); ++blob_iterator){
    info() << " blob number: " << (*blob_iterator) << endmsg;

    shower_ndigits.push_back( (*blob_iterator)->getAllDigits().size() );
    shower_nclusters.push_back( (*blob_iterator)->nclusters());
    shower_energy.push_back( (*blob_iterator)->energy());

    isGood = isShowerGood((*blob_iterator), event);
    if (isGood){
        const Gaudi::XYZVector& shower_direction = (*blob_iterator)->direction();
        double shower_nplanes = Calculate_BlobdEdX(event, (*blob_iterator), 1);
	      double shower_length = 1.7*shower_nplanes/TMath::Abs(TMath::Cos(shower_direction.theta()));
        double shower_dist_vtx = CalcMinBlobSeparation((*blob_iterator), event);
        int n_clusters = (*blob_iterator)->nclusters();
        if (shower_dist_vtx >= 150){
            double Blob_dist_fitline = CalProbFitLine((*blob_iterator), vtx_position);
            shower_dist_startpoint.push_back(Blob_dist_fitline);
	          //if (Blob_dist_fitline >= 0.0 && Blob_dist_fitline <= 35.00){  //-- approx the width of a module, this affects too much the 1pi0 events, I'm removing it
            if (Blob_dist_fitline >= 0.0){
            	  const double SENTINEL = -9.9;
            	  std::vector<double> gamma_evis_v(5,SENTINEL);
            	  std::vector<double> gamma_energy_v(5,SENTINEL);
            	  double gamma_energy = m_idHoughBlob->getBlobEnergyTime_New((*blob_iterator), gamma_evis_v, gamma_energy_v);

            	  double shower_dEdX = ((*blob_iterator)->energy())/shower_length;
            	  double score = CalculatePhotonScore(gamma_energy, shower_dEdX);
            	  info() << "Length = " << shower_length << " , dEdX = " << shower_dEdX << " , Energy = " << gamma_energy << " , Photon Score z = " << score << endmsg;
            	  if (score <= 100){
            		    shower_isGood.push_back(1);
            		    GoodBlobs.push_back((*blob_iterator));
                    shower_score.push_back(score);
                    if (truthEvent){
                      SmartRefVector<Minerva::IDCluster> showerClusters = (*blob_iterator)->clusters();
                      DigitVectorTruthInfo showerInfo;
                      showerInfo.ParseTruth(showerClusters, fTrajectoryMap, 22);
                      int pdg = showerInfo.GetMostEvisPdg();
                      shower_pdg.push_back( pdg );
                      info() << "PDG = " << pdg << ", score = " << score << endmsg;
                    }else{
                      shower_pdg.push_back(0);
                    }
            	  } // End of if (score)
            } // End of if Blob_dist_fitline
        } // End shower dist vtx
    } // End isGood
    else shower_isGood.push_back(0);
  }

  event->setContainerIntData("shower_ndigits", shower_ndigits);
  event->setContainerIntData("shower_nclusters", shower_nclusters);
  event->setContainerDoubleData("shower_energy", shower_energy);
  event->setContainerDoubleData("shower_dist_fitline", shower_dist_startpoint);
  event->setContainerIntData("shower_isGood", shower_isGood);
  event->setContainerDoubleData("shower_score", shower_score);
  event->setContainerIntData("shower_pdg", shower_pdg);
  event->setIntData("Good_Blobs", GoodBlobs.size() );

  }

double CCProtonPi0::CalculatePhotonScore(double energy, double dEdX) const
{

  f1_photon->SetParameter(0, -1.240e-18); f1_photon->SetParameter(1, 5.177e-14); f1_photon->SetParameter(2, -8.145e-11);
  f1_photon->SetParameter(3, 6.085e-8); f1_photon->SetParameter(4, -2.393e-5); f1_photon->SetParameter(5, 0.0102);
  f1_photon->SetParameter(6, 1.868);

  f1_photon1->SetParameter(0, 0.0042); f1_photon1->SetParameter(1, 1.868);

  //--- Sigma dEdX
  std1_photon->SetParameter(0, -5.22577e-18); std1_photon->SetParameter(1, 2.09294e-14); std1_photon->SetParameter(2, -2.99908e-11);
  std1_photon->SetParameter(3, 1.73717e-8); std1_photon->SetParameter(4, -2.22225e-6); std1_photon->SetParameter(5, 0.000195544);
  std1_photon->SetParameter(6, 1.0712);

  std1_photon1->SetParameter(0, 0.00123659); std1_photon1->SetParameter(1, 0.901444);

  double z;
  if (energy <= 1475)  z = (dEdX - f1_photon->Eval(energy))/std1_photon->Eval(energy);
  else z = (dEdX - f1_photon1->Eval(energy))/std1_photon1->Eval(energy);

  return z;
}

int CCProtonPi0::Calculate_BlobdEdX(Minerva::PhysicsEvent *event, const Minerva::IDBlob* blob, unsigned int blob_number) const
{

    //Sanity Check
    if (blob == NULL){
        warning() << "Passed NULL Blob" << endmsg;
    }

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    debug() << "Calculate dEdX for Good Blob: " << blob_number << "(forward, backward) (" << (blob->direction().z() > 0) << "," << (blob->direction().z() < 0) << ")" << endmsg;

    std::map<int, std::vector<double> > planeClusterEnergyMap;
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    debug() << "Vertex to increasingDistanceFromVertex = (" << vtx_position.X() << "," << vtx_position.Y() << "," << vtx_position.Z() << ")" << endmsg;

    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(vtx_position));

    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin(); c != clusters.end(); ++c)
    {
        const Minerva::PlaneID& planeId = (*c)->planeid();
        unsigned int plane = planeId.plane();
        unsigned module = planeId.module();
        const unsigned int planeNo = 2*module + plane;

        planeClusterEnergyMap[planeNo].push_back((*c)->energy());
    }

    std::vector<int>    blob_cluster_occupancy;
    std::vector<double> blob_cluster_energy;

    double total = 0.0;
    double total1 = 0.0; // exclude the first plane
    bool doublet = false;
    bool empty_plane = false;

    for (std::map<int, std::vector<double> >::iterator pl = planeClusterEnergyMap.begin(); pl != planeClusterEnergyMap.end(); ++pl)
    {
        const int plane = pl->first;
        const std::vector<double>& clusterEnergies = pl->second;

        std::map<int, std::vector<double> >::iterator next_pl = pl; ++next_pl;
        const int next_plane = next_pl->first;

        const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
        blob_cluster_occupancy.push_back(clusterEnergies.size());
        blob_cluster_energy.push_back(total_plane_energy);

        total += total_plane_energy;
        if (pl != planeClusterEnergyMap.begin() ) total1 += total_plane_energy;

        if (clusterEnergies.size() > 1) doublet = true;

        if (next_pl != planeClusterEnergyMap.end() && std::abs(next_plane - plane) > 1) empty_plane = true;

        if (std::distance(planeClusterEnergyMap.begin(), pl) > 4) break;
    }


    std::vector<double> blob_rev_cluster_energy;
    std::map<int, std::vector<double> >::reverse_iterator pl;
    for (std::map<int, std::vector<double> >::reverse_iterator pl = planeClusterEnergyMap.rbegin(); pl != planeClusterEnergyMap.rend(); ++pl)
    {
        const std::vector<double>& clusterEnergies = pl->second;
        const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
        blob_rev_cluster_energy.push_back(total_plane_energy);
    }

    //Gaudi::XYZVector direction = blob->direction();

    //double blob_cos_theta = direction.Z();

    int nplanes = planeClusterEnergyMap.size();
    //double blob_length = TMath::Abs( (1.7*nplanes)/blob_cos_theta );

    //double blob_dEdX = total/nplanes;
    //double blob_dEdX1 = total1/blob_length;

 /*
    event->setContainerIntData("shower_dedx_nplane",    nplanes);
    event->setContainerDoubleData("shower_dedx_total",  total);
    event->setContainerDoubleData("shower_dedx",        total/blob_cluster_energy.size());
    event->setContainerDoubleData("shower_dEdX",        blob_dEdX);
    event->setContainerDoubleData("shower_dedx_total1", total1);
    event->setContainerDoubleData("shower_dedx1",       total1/std::max((int)blob_cluster_energy.size()-1, 1) );
    event->setContainerDoubleData("shower_dEdX1",       blob_dEdX1);
 */

    //return blob_dEdX;
    return nplanes;
}


void CCProtonPi0::MarkFoundPi0Blobs(std::vector<Minerva::IDBlob*> &foundBlobs) const
{
    if (foundBlobs.size() == 1){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
    }else if (foundBlobs.size() == 2){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
        m_hitTagger->applyColorTag( (foundBlobs)[1], m_Color_Gamma2Prong ); // light blue   
    }else if (foundBlobs.size() > 2){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
        m_hitTagger->applyColorTag( (foundBlobs)[1], m_Color_Gamma2Prong ); // light blue
        for(unsigned int i = 2; i < foundBlobs.size(); i++ ){
            m_hitTagger->applyColorTag( (foundBlobs)[i], m_Color_GammaOtherProng ); // orange 
        }
    }
}

//=======================================================================
//  DispersedBlob
//=======================================================================
void CCProtonPi0::DispersedBlob( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction *truthEvent ) const
{
    //  --  Dispersed blob
    double dispersed_energy = 0;

    Minerva::IDBlob* dispersedBlob = new Minerva::IDBlob();
    SmartRefVector<Minerva::IDCluster> unusedClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> dispersedClusters;
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;

    for ( it_clus = unusedClusters.begin(); it_clus != unusedClusters.end(); it_clus++){
        if ( (*it_clus)->pe()/(*it_clus)->iddigs() <= 3 ) continue;
        dispersedClusters.push_back(*it_clus);
    }

    if (!dispersedClusters.empty()) {
        m_blobUtils->insertIDClusters( dispersedClusters, dispersedBlob, Minerva::IDBlob::DispersedIDBlobPatRec );
        double dummy1 = 0.0;
        double dummy2 = 0.0;
        double dummy3 = 0.0;
        double dummy4 = 0.0;
        m_idHoughBlob->getBlobEnergyTime_Old(dispersedBlob,dispersed_energy,dummy1,dummy2,dummy3,dummy4);
        debug()<< "Adding dispersed blob with " << dispersedBlob->nclusters()
            << " clusters; Dispersed energy = "  << dispersed_energy << endmsg;
        addObject( event, dispersedBlob );

        if (truthEvent){
            SaveTruthUnusedClusterEnergy_Dispersed(truthEvent, dispersedClusters);
        }
    }

    debug()<<"Extra_Energy_Dispersed = "<<dispersed_energy<<endmsg;

    event->setDoubleData( "Extra_Energy_Dispersed", dispersed_energy );

}

int CCProtonPi0::GetMCHitPDG(const SmartRef<Minerva::MCHit> mc_hit) const
{
    TraverseHistory history(mc_hit->GetTrackId());
    history.DoTraverse(fTrajectoryMap);
    int pdg = history.GetPrimaryPdg();
    return pdg;
}


int CCProtonPi0::GetDigitPDG(const Minerva::MCIDDigit* mcdigit) const
{
    const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
    std::vector<int> all_pdg;

    if (hits.empty()){
        debug()<<"No Hits in Digit"<<endmsg;
        return -1;
    }

    // Save PDG of each Hit in Digit
    SmartRefVector<Minerva::MCHit>::const_iterator h;
    for ( h = hits.begin(); h != hits.end(); ++h){
        int hit_pdg = GetMCHitPDG(*h);             
        all_pdg.push_back(hit_pdg);
    }

    // Get Digit From Vector
    int digit_pdg;
    if (all_pdg.size() > 1) digit_pdg = GetPDGfromVector(all_pdg);
    else digit_pdg = all_pdg[0];

    return digit_pdg;
}

int CCProtonPi0::GetPDGfromVector(std::vector<int> &all_pdg) const
{
    // First Check whether all PDGs are Same
    bool isAllSame = true;

    debug()<<"PDG Vector Content"<<endmsg; 

    // Check whether all PDGs are Same
    int prev_pdg = all_pdg[0];
    debug()<<prev_pdg<<endmsg;
    for (unsigned int i = 1; i < all_pdg.size(); i++){
        debug()<<all_pdg[i]<<endmsg;
        if (prev_pdg != all_pdg[i]) isAllSame = false;
        prev_pdg = all_pdg[i];
    }

    if (isAllSame){
        debug()<<"All PDGs are same - returning unique PDG!"<<endmsg;
        return prev_pdg;
    }else{
        debug()<<"All PDGs are NOT same - returning -1"<<endmsg;
        return -1;
    }
}

void CCProtonPi0::SaveExtraEnergyLeftover(Minerva::PhysicsEvent *event) const
{
    // Get All UnusedClusters
    SmartRefVector<Minerva::IDCluster> unusedClusters_Leftover;
    unusedClusters_Leftover = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    debug()<<"N(Unused) to Get Extra Energy = "<<unusedClusters_Leftover.size()<<endmsg;

    SmartRefVector<Minerva::IDCluster> Clusters = FilterInSphereClusters(event, unusedClusters_Leftover, 300);
    double extra_energy = getClusterEnergy(Clusters);

    event->setDoubleData("Extra_Energy_Leftover", extra_energy);
}


double CCProtonPi0::getClusterEnergy( SmartRefVector<Minerva::IDCluster> clusters) const
{
    SmartRefVector<Minerva::IDCluster>::iterator c;

    double total_energy = 0.0;
    for (c = clusters.begin(); c != clusters.end(); ++c){
        total_energy += (*c)->energy();
    }

    return total_energy;
}

double CCProtonPi0::GetVertexEnergyCalConst(double evis) const
{
    double cal_const;

    if (m_ProtonParticles.size() > 0){
        cal_const = 0.91;
    }else{
        // Values determined from MC Study
        double min_cal_const = 1.82;
        double m = -0.001442;
        double c = 1.964;
        cal_const = m*evis + c;

        // We have negative slope, limit min value
        if (cal_const < min_cal_const) cal_const = min_cal_const;
    }

    return cal_const;
}

double CCProtonPi0::GetVertexEnergy(const Minerva::PhysicsEvent *event) const
{
    double vertex_evis = event->getDoubleData("vertex_blob_evis");
    double cal_const = GetVertexEnergyCalConst(vertex_evis);
    double vertex_energy = vertex_evis * cal_const;

    return vertex_energy;
}
void CCProtonPi0::FillUsableClusters(SmartRefVector<Minerva::IDCluster> &usableClusters, Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    SmartRefVector<Minerva::IDCluster> unusedClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> rejectedClusters;

    // Get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();

    double energyTracker    = 0.0;
    double energyECAL       = 0.0;
    double energyHCAL       = 0.0;
    double EnergyToAnalyze  = 0.0;

    // Loop over all UNUSED Clusters
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;
    for ( it_clus = unusedClusters.begin(); it_clus != unusedClusters.end(); ++it_clus){

        // Reject Out Time Clusters 
        bool isCluster_OutTime = std::abs( (*it_clus)->time() -  muon_position.T() ) > m_rejectedClustersTime;
        if (isCluster_OutTime){
            debug()<<"Rejecting Cluster: Out Time! -- Not saving to rejectedClusters Vector"<<endmsg;
            continue; 
        }

        // Reject Low Charge Clusters
        bool isCluster_LowCharge = (*it_clus)->pe()/(*it_clus)->iddigs() <= 3; 
        if (isCluster_LowCharge ) {
            debug()<<"Rejecting Cluster: Low Charge!"<<endmsg;
            rejectedClusters.push_back(*it_clus);
            continue;
        }

        // Reject Clusters in HCAL
        bool isCluster_HCAL = (*it_clus)->subdet() == Minerva::IDCluster::HCAL;
        if (isCluster_HCAL){
            debug()<<"Rejecting Cluster: Inside HCAL!"<<endmsg;
            rejectedClusters.push_back(*it_clus);
            continue; 
        }
       
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::Tracker ) energyTracker += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::ECAL )    energyECAL += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::HCAL )    energyHCAL += (*it_clus)->energy();

        // Fill usableClusters
        usableClusters.push_back(*it_clus);
    }


    // Save Energy Information to NTuples
    event->setDoubleData("ConeBlobs_usable_evis_Tracker", energyTracker );
    event->setDoubleData("ConeBlobs_usable_evis_ECAL", energyECAL );
    event->setDoubleData("ConeBlobs_usable_evis_HCAL", energyHCAL );

    // Check Energy to Analyze
    EnergyToAnalyze = energyTracker + energyECAL + energyHCAL;
    debug()<<"usableClusters.size() = "<<usableClusters.size()<<endmsg;
    debug()<<"Energy to analyze = "<< EnergyToAnalyze << endmsg;

    // Process Rejected Clusters
    debug()<<"rejectedClusters.size() = "<<rejectedClusters.size()<<endmsg;
    ProcessRejectedClusters(rejectedClusters, event, truthEvent);
}

void CCProtonPi0::ProcessRejectedClusters(SmartRefVector<Minerva::IDCluster> &rejectedClusters,Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    // Save Truth Information to NTuples for Rejected Clusters
    if (truthEvent && m_study_unused_energy){
        SaveTruthUnusedClusterEnergy_Rejected(truthEvent, rejectedClusters);
    }
    double rejected_energy = 0.0;

    //storing rejected id Clusters
    Minerva::IDBlob* rejectedBlob = new Minerva::IDBlob();

    if (!rejectedClusters.empty()) {
        m_blobUtils->insertIDClusters(rejectedClusters, rejectedBlob, Minerva::IDBlob::DispersedIDBlobPatRec );
        double dummy1 = 0.0;
        double dummy2 = 0.0;
        double dummy3 = 0.0;
        double dummy4 = 0.0;
        m_idHoughBlob->getBlobEnergyTime_Old(rejectedBlob, rejected_energy, dummy1, dummy2, dummy3, dummy4);
        debug()<<"rejectedBlob N(Clusters) = " << rejectedBlob->nclusters()<<endmsg;
        debug()<<"rejectedBlob Energy = "  <<rejected_energy<<endmsg;
        addObject( event, rejectedBlob );
        m_hitTagger->applyColorTag( rejectedBlob, m_Color_RejectedBlob ); // red
    
        rejected_energy = m_extraEnergyTool->getCorrectedIDEnergy( rejectedClusters, 1.0);
    }

    debug()<<"Extra_Energy_Rejected = "<<rejected_energy<<endmsg;
    event->setDoubleData( "Extra_Energy_Rejected", rejected_energy);

}


bool CCProtonPi0::RecoverShowers_InvMass( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const
{
    debug() << "Finding best shower combination using Invariant Mass" << endmsg;

    // Get Vertex Position
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    std::vector<double> masses;
    std::vector<double> good_masses;
    std::vector<Minerva::IDBlob*> bestBlobs;
    double dM = 1.e6;
    std::vector<Minerva::IDBlob*>::iterator b1;
    for (b1 = foundBlobs.begin(); b1 != foundBlobs.end(); ++b1) {
        std::vector<Minerva::IDBlob*>::iterator b2 = b1+1;

        if (b2 == foundBlobs.end()) break;

        ClusterVectorInfo clusterInfo1((*b1)->clusters());
        if (clusterInfo1.GetNx() < 2) continue;

        double e1 = 0.0;
        double dummy1 = 0.0;
        double dummy2 = 0.0;
        double dummy3 = 0.0;
        double dummy4 = 0.0;
        m_idHoughBlob->getBlobEnergyTime_Old((*b1),e1,dummy1,dummy2,dummy3,dummy4);
        // double t1 = (*b1)->time();

        bool goodPosition1  = m_idHoughBlob->GetStartPosition((*b1),vtx_position,true);
        bool goodDirection1 = m_idHoughBlob->GetDirection((*b1),vtx_position);
        Gaudi::XYZVector direction1 = (*b1)->direction();

        if (!goodPosition1 || !goodDirection1) break;

        for (; b2 != foundBlobs.end(); ++b2) {

            ClusterVectorInfo clusterInfo2((*b2)->clusters());
            if (clusterInfo2.GetNx() < 2) continue;

            double e2 = 0.0;
            m_idHoughBlob->getBlobEnergyTime_Old((*b2),e2,dummy1,dummy2,dummy3,dummy4);

            // double t2 = (*b2)->time();

            bool goodPosition2  = m_idHoughBlob->GetStartPosition((*b2),vtx_position,true);
            bool goodDirection2 = m_idHoughBlob->GetDirection((*b2),vtx_position);

            if (!goodPosition2 || !goodDirection2) continue;

            Gaudi::XYZVector direction2 = (*b2)->direction();

            const double cos_oangle = direction1.Dot(direction2);
            const double mgg = std::sqrt(2*e1*e2*(1-cos_oangle));
            debug() << "\tCombination (" << std::distance(foundBlobs.begin(),b1) << ","
                << std::distance(foundBlobs.begin(),b2) << "): " << mgg
                                                                    << endmsg;

            masses.push_back(mgg);

            if (std::abs(mgg-134.9) < dM) {
                debug() << "\t replace dM: " << dM << " " << " with " << std::abs(mgg-134.9) << endmsg;
                dM = std::abs(mgg-134.9);
                bestBlobs.clear();
                bestBlobs.push_back(*b1);
                bestBlobs.push_back(*b2);
            }
        }
    } 
    bool isBlobsRecovered = ( bestBlobs.size() == 2 );
    if (bestBlobs.size() == 2) foundBlobs.swap(bestBlobs);

    return isBlobsRecovered;
}

bool CCProtonPi0::Save_1ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const
{
    Minerva::IDBlob* shower1 = foundBlobs[0];

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();
    bool isShower1Good = isShowerGood(shower1, event);

    if (isShower1Good){
        if (m_study_shower_recovery){
            int nclusters = shower1->nclusters();
            double energy = shower1->energy();

            const Gaudi::XYZPoint& start_point = shower1->startPoint();
            const Gaudi::XYZVector& direction = shower1->direction();

            double theta = direction.theta()*180/M_PI;
            double dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(), 
                    start_point.X(), start_point.Y(), start_point.Z());

            debug()<<"1Shower nClusters = "<<nclusters<<endmsg;
            debug()<<"1Shower Energy = "<<energy<<endmsg;
            debug()<<"1Shower Theta = "<<theta<<endmsg;
            debug()<<"1Shower Distance to Vertex = "<<dist_vtx<<endmsg;


            event->setIntData("OneShower_nClusters", nclusters); 
            event->setDoubleData("OneShower_energy", energy);
            event->setDoubleData("OneShower_theta",theta);
            event->setDoubleData("OneShower_dist_vtx", dist_vtx);
        }
        return true;
    }else{
        debug()<<"Found 1Shower direction is NOT Good!"<<endmsg;
        return false;
    }
}

void CCProtonPi0::Save_1ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const
{
    Minerva::IDBlob* shower1 = foundBlobs[0];

    SmartRefVector<Minerva::IDCluster> shower1Clusters = shower1->clusters();
    DigitVectorTruthInfo shower1Info;
    shower1Info.ParseTruth(shower1Clusters, fTrajectoryMap, 0);

    truthEvent->setIntData("OneShower_evis_most_pdg", shower1Info.GetMostEvisPdg());
    truthEvent->setDoubleData("OneShower_evis_total_norm", shower1Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("OneShower_evis_total_truth", shower1Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("OneShower_evis_pizero", shower1Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("OneShower_evis_piplus", shower1Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("OneShower_evis_piminus", shower1Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("OneShower_evis_muon", shower1Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("OneShower_evis_proton", shower1Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("OneShower_evis_neutron", shower1Info.GetEdepByPdg(PDG::neutron));
}



bool CCProtonPi0::Save_AllShowersInfo( std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::PhysicsEvent *event) const
{
  const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

  std::vector<int> AllShowers_nclusters;
  std::vector<double> AllShowers_energy;
  std::vector<double> AllShowers_theta;
  std::vector<double> AllShowers_dist_vtx;
  std::vector<double> AllShowers_dEdX;


  std::vector<Minerva::IDBlob*>::iterator it_blob;
  for (it_blob = GoodBlobs.begin(); it_blob != GoodBlobs.end(); ++it_blob){
    int shower_nclusters = (*it_blob)->nclusters();
    double shower_energy = (*it_blob)->energy();
    const Gaudi::XYZPoint& shower_start_point = (*it_blob)->startPoint();
    const Gaudi::XYZVector& shower_direction = (*it_blob)->direction();
    double shower_theta = shower_direction.theta()*180/M_PI;
    double shower_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
					   shower_start_point.X(), shower_start_point.Y(), shower_start_point.Z());
    double shower_dedx = Calculate_BlobdEdX(event, (*it_blob), 1);

    double dEdX = TMath::Abs((shower_dedx/1.7)*TMath::Cos(shower_direction.theta()));

    AllShowers_nclusters.push_back(shower_nclusters);
    AllShowers_energy.push_back(shower_energy);
    AllShowers_theta.push_back(shower_theta);
    AllShowers_dist_vtx.push_back(shower_dist_vtx);
    AllShowers_dEdX.push_back(dEdX);

    debug()<<"AllShower shower_nClusters = "          << shower_nclusters << endmsg;
    debug()<<"AllShower shower_Energy = "             << shower_energy    << endmsg;
    debug()<<"AllShower shower_Theta = "              << shower_theta     << endmsg;
    debug()<<"AllShower shower_Distance to Vertex = " << shower_dist_vtx  << endmsg;
    debug()<<"AllShower shower_dEdX = "               << dEdX             << endmsg;
  }//-- End loop over blobs

  event->setContainerIntData("AllShower_nclusters",   AllShowers_nclusters);
  event->setContainerDoubleData("AllShower_energy",   AllShowers_energy);
  event->setContainerDoubleData("AllShower_theta",    AllShowers_theta);
  event->setContainerDoubleData("AllShower_dist_vtx", AllShowers_dist_vtx);
  event->setContainerDoubleData("AllShower_dEdX",     AllShowers_dEdX);

  return true;
}

bool CCProtonPi0::Save_3ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const
{
    Minerva::IDBlob* shower1 = foundBlobs[0];
    Minerva::IDBlob* shower2 = foundBlobs[1];
    Minerva::IDBlob* shower3 = foundBlobs[2];

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    bool isShower1Good = isShowerGood(shower1, event);
    bool isShower2Good = isShowerGood(shower2, event);
    bool isShower3Good = isShowerGood(shower3, event);

    if (isShower1Good && isShower2Good && isShower3Good){
        // Shower 1
        int s1_nclusters = shower1->nclusters();
        double s1_energy = shower1->energy();
        const Gaudi::XYZPoint& s1_start_point = shower1->startPoint();
        const Gaudi::XYZVector& s1_direction = shower1->direction();
        double s1_theta = s1_direction.theta()*180/M_PI;
        double s1_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
                s1_start_point.X(), s1_start_point.Y(), s1_start_point.Z());

        debug()<<"3Shower s1_nClusters = "<<s1_nclusters<<endmsg;
        debug()<<"3Shower s1_Energy = "<<s1_energy<<endmsg;
        debug()<<"3Shower s1_Theta = "<<s1_theta<<endmsg;
        debug()<<"3Shower s1_Distance to Vertex = "<<s1_dist_vtx<<endmsg;

        event->setIntData("ThreeShower_s1_nClusters", s1_nclusters);
        event->setDoubleData("ThreeShower_s1_energy", s1_energy);
        event->setDoubleData("ThreeShower_s1_theta",s1_theta);
        event->setDoubleData("ThreeShower_s1_dist_vtx", s1_dist_vtx);

        // Shower 2
        int s2_nclusters = shower2->nclusters();
        double s2_energy = shower2->energy();
        const Gaudi::XYZPoint& s2_start_point = shower2->startPoint();
        const Gaudi::XYZVector& s2_direction = shower2->direction();
        double s2_theta = s2_direction.theta()*180/M_PI;
        double s2_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
                s2_start_point.X(), s2_start_point.Y(), s2_start_point.Z());

        debug()<<"3Shower s2_gnClusters = "<<s2_nclusters<<endmsg;
        debug()<<"3Shower s2_gEnergy = "<<s2_energy<<endmsg;
        debug()<<"3Shower s2_gTheta = "<<s2_theta<<endmsg;
        debug()<<"3Shower s2_gDistance to Vertex = "<<s2_dist_vtx<<endmsg;

        event->setIntData("ThreeShower_s2_nClusters", s2_nclusters);
        event->setDoubleData("ThreeShower_s2_energy", s2_energy);
        event->setDoubleData("ThreeShower_s2_theta",s2_theta);
        event->setDoubleData("ThreeShower_s2_dist_vtx", s2_dist_vtx);

        // Shower 3
        int s3_nclusters = shower3->nclusters();
        double s3_energy = shower3->energy();
        const Gaudi::XYZPoint& s3_start_point = shower3->startPoint();
        const Gaudi::XYZVector& s3_direction = shower3->direction();
        double s3_theta = s3_direction.theta()*180/M_PI;
        double s3_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
                s3_start_point.X(), s3_start_point.Y(), s3_start_point.Z());

        debug()<<"3Shower s3_gnClusters = "<<s3_nclusters<<endmsg;
        debug()<<"3Shower s3_gEnergy = "<<s3_energy<<endmsg;
        debug()<<"3Shower s3_gTheta = "<<s3_theta<<endmsg;
        debug()<<"3Shower s3_gDistance to Vertex = "<<s3_dist_vtx<<endmsg;

        event->setIntData("ThreeShower_s3_nClusters", s3_nclusters);
        event->setDoubleData("ThreeShower_s3_energy", s3_energy);
        event->setDoubleData("ThreeShower_s3_theta",s3_theta);
        event->setDoubleData("ThreeShower_s3_dist_vtx", s3_dist_vtx);

        return true;
    }else{
        debug()<<"Not All Showers direction are good!"<<endmsg;
        return false;
    }
}

void CCProtonPi0::Save_AllShowersTruthMatch( std::vector<Minerva::IDBlob*> &GoodBlobs, Minerva::GenMinInteraction* truthEvent) const
{
  std::vector<int> showers_evis_most_pdg;
  std::vector<double> showers_evis_total_norm;
  std::vector<double> showers_evis_total_truth;
  std::vector<double> showers_evis_pizero;
  std::vector<double> showers_evis_piminus;
  std::vector<double> showers_evis_piplus;
  std::vector<double> showers_evis_muon;
  std::vector<double> showers_evis_proton;
  std::vector<double> showers_evis_neutron;

  std::vector<Minerva::IDBlob*>::iterator it_blob;
  for (it_blob = GoodBlobs.begin(); it_blob != GoodBlobs.end(); ++it_blob){

    SmartRefVector<Minerva::IDCluster> showerClusters = (*it_blob)->clusters();
    DigitVectorTruthInfo showerInfo;
    showerInfo.ParseTruth(showerClusters, fTrajectoryMap, 22);

    showers_evis_most_pdg.push_back( showerInfo.GetMostEvisPdg() );
    showers_evis_total_norm.push_back( showerInfo.GetTotalNormEnergy() );
    showers_evis_total_truth.push_back( showerInfo.GetTotalTruthEnergy() );
    showers_evis_pizero.push_back( showerInfo.GetEdepByPdg(PDG::pi0) );
    showers_evis_piminus.push_back( showerInfo.GetEdepByPdg(-(PDG::pi)) );
    showers_evis_piplus.push_back( showerInfo.GetEdepByPdg(PDG::pi) );
    showers_evis_muon.push_back( showerInfo.GetEdepByPdg(PDG::muon) );
    showers_evis_proton.push_back( showerInfo.GetEdepByPdg(PDG::proton) );
    showers_evis_neutron.push_back( showerInfo.GetEdepByPdg(PDG::neutron) );
  }

  truthEvent->setContainerIntData( "AllShowers_evis_most_pdg", showers_evis_most_pdg );
  truthEvent->setContainerDoubleData( "AllShowers_evis_total_norm", showers_evis_total_norm );
  truthEvent->setContainerDoubleData( "AllShowers_evis_total_truth", showers_evis_total_truth );
  truthEvent->setContainerDoubleData( "AllShowers_evis_pizero", showers_evis_pizero );
  truthEvent->setContainerDoubleData( "AllShowers_evis_piminus", showers_evis_piminus );
  truthEvent->setContainerDoubleData( "AllShowers_evis_piplus", showers_evis_piplus );
  truthEvent->setContainerDoubleData( "AllShowers_evis_muon", showers_evis_muon );
  truthEvent->setContainerDoubleData( "AllShowers_evis_proton", showers_evis_proton );
  truthEvent->setContainerDoubleData( "AllShowers_evis_neutron", showers_evis_neutron );
}

void CCProtonPi0::Save_3ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const
{
    // Shower 1
    Minerva::IDBlob* shower1 = foundBlobs[0];

    SmartRefVector<Minerva::IDCluster> shower1Clusters = shower1->clusters();
    DigitVectorTruthInfo shower1Info;
    shower1Info.ParseTruth(shower1Clusters, fTrajectoryMap, 0);

    truthEvent->setIntData("ThreeShower_s1_evis_most_pdg", shower1Info.GetMostEvisPdg());
    truthEvent->setDoubleData("ThreeShower_s1_evis_total_norm", shower1Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("ThreeShower_s1_evis_total_truth", shower1Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("ThreeShower_s1_evis_pizero", shower1Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("ThreeShower_s1_evis_piplus", shower1Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("ThreeShower_s1_evis_piminus", shower1Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("ThreeShower_s1_evis_muon", shower1Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("ThreeShower_s1_evis_proton", shower1Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("ThreeShower_s1_evis_neutron", shower1Info.GetEdepByPdg(PDG::neutron));

    // Shower 2
    Minerva::IDBlob* shower2 = foundBlobs[1];

    SmartRefVector<Minerva::IDCluster> shower2Clusters = shower2->clusters();
    DigitVectorTruthInfo shower2Info;
    shower2Info.ParseTruth(shower2Clusters, fTrajectoryMap, 0);

    truthEvent->setIntData("ThreeShower_s2_evis_most_pdg", shower2Info.GetMostEvisPdg());
    truthEvent->setDoubleData("ThreeShower_s2_evis_total_norm", shower2Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("ThreeShower_s2_evis_total_truth", shower2Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("ThreeShower_s2_evis_pizero", shower2Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("ThreeShower_s2_evis_piplus", shower2Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("ThreeShower_s2_evis_piminus", shower2Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("ThreeShower_s2_evis_muon", shower2Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("ThreeShower_s2_evis_proton", shower2Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("ThreeShower_s2_evis_neutron", shower2Info.GetEdepByPdg(PDG::neutron));

    // Shower 3
    Minerva::IDBlob* shower3 = foundBlobs[2];

    SmartRefVector<Minerva::IDCluster> shower3Clusters = shower3->clusters();
    DigitVectorTruthInfo shower3Info;
    shower3Info.ParseTruth(shower3Clusters, fTrajectoryMap, 0);

    truthEvent->setIntData("ThreeShower_s3_evis_most_pdg", shower3Info.GetMostEvisPdg());
    truthEvent->setDoubleData("ThreeShower_s3_evis_total_norm", shower3Info.GetTotalNormEnergy());
    truthEvent->setDoubleData("ThreeShower_s3_evis_total_truth", shower3Info.GetTotalTruthEnergy());
    truthEvent->setDoubleData("ThreeShower_s3_evis_pizero", shower3Info.GetEdepByPdg(PDG::pi0));
    truthEvent->setDoubleData("ThreeShower_s3_evis_piplus", shower3Info.GetEdepByPdg(PDG::pi));
    truthEvent->setDoubleData("ThreeShower_s3_evis_piminus", shower3Info.GetEdepByPdg(-(PDG::pi)));
    truthEvent->setDoubleData("ThreeShower_s3_evis_muon", shower3Info.GetEdepByPdg(PDG::muon));
    truthEvent->setDoubleData("ThreeShower_s3_evis_proton", shower3Info.GetEdepByPdg(PDG::proton));
    truthEvent->setDoubleData("ThreeShower_s3_evis_neutron", shower3Info.GetEdepByPdg(PDG::neutron));
}

bool CCProtonPi0::Save_4ShowerInfo( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const
{
  Minerva::IDBlob* shower1 = foundBlobs[0];
  Minerva::IDBlob* shower2 = foundBlobs[1];
  Minerva::IDBlob* shower3 = foundBlobs[2];
  Minerva::IDBlob* shower4 = foundBlobs[3];

  const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

  debug()<<"All of the 4 blobs founded"<<endmsg;
  bool isShower1Good = isShowerGood(shower1, event);
  debug()<<"Shower 1 is good"<<endmsg;
  bool isShower2Good = isShowerGood(shower2, event);
  debug()<<"Shower 2 is good"<<endmsg;
  bool isShower3Good = isShowerGood(shower3, event);
  debug()<<"Shower 3 is good"<<endmsg;
  bool isShower4Good = isShowerGood(shower4, event);
  debug()<<"Shower 4 is good"<<endmsg;

  if (isShower1Good){
	 if (m_study_shower_recovery){
		 int just_s1_nclusters = shower1->nclusters();
		 double just_s1_energy = shower1->energy();
		 const Gaudi::XYZPoint& just_s1_start_point = shower1->startPoint();
		 const Gaudi::XYZVector& just_s1_direction = shower1->direction();
		 double just_s1_theta = just_s1_direction.theta()*180/M_PI;
		 double just_s1_dist_vtx = calcDistance(vtx_position.X(), vtx_position.Y(), vtx_position.Z(), just_s1_start_point.X(), just_s1_start_point.Y(), just_s1_start_point.Z() );

	     debug()<<"4Shower s1_nClusters = " <<just_s1_nclusters<<endmsg;
	     debug()<<"4Shower s1_Energy = " <<just_s1_energy<<endmsg;
	     debug()<<"4Shower s1_Theta = " <<just_s1_theta<<endmsg;
	     debug()<<"4Shower s1_Distance to Vertex = " <<just_s1_dist_vtx<<endmsg;

	     event->setIntData("FourShower_s1_nClusters", just_s1_nclusters);
	     event->setDoubleData("FourShower_s1_energy", just_s1_energy);
	     event->setDoubleData("FourShower_s1_theta", just_s1_theta);
	     event->setDoubleData("FourShower_s1_dist_vtx", just_s1_dist_vtx);
	 }
  }

  if (isShower2Good){
	 if (m_study_shower_recovery){
		 int just_s2_nclusters = shower2->nclusters();
		 double just_s2_energy = shower2->energy();
		 const Gaudi::XYZPoint& just_s2_start_point = shower2->startPoint();
		 const Gaudi::XYZVector& just_s2_direction = shower2->direction();
		 double just_s2_theta = just_s2_direction.theta()*180/M_PI;
		 double just_s2_dist_vtx = calcDistance(vtx_position.X(), vtx_position.Y(), vtx_position.Z(), just_s2_start_point.X(), just_s2_start_point.Y(), just_s2_start_point.Z() );

	     debug()<<"4Shower s2_nClusters = " << just_s2_nclusters << endmsg;
	     debug()<<"4Shower s2_Energy = " << just_s2_energy << endmsg;
	     debug()<<"4Shower s2_Theta = " << just_s2_theta << endmsg;
	     debug()<<"4Shower s2_Distance to Vertex = " << just_s2_dist_vtx << endmsg;

	     event->setIntData("FourShower_s2_nClusters", just_s2_nclusters);
	     event->setDoubleData("FourShower_s2_energy", just_s2_energy);
	     event->setDoubleData("FourShower_s2_theta", just_s2_theta);
	     event->setDoubleData("FourShower_s2_dist_vtx", just_s2_dist_vtx);
	 }
  }

  if (isShower3Good){
	 if (m_study_shower_recovery){
		 int just_s3_nclusters = shower3->nclusters();
		 double just_s3_energy = shower3->energy();
		 const Gaudi::XYZPoint& just_s3_start_point = shower3->startPoint();
		 const Gaudi::XYZVector& just_s3_direction = shower3->direction();
		 double just_s3_theta = just_s3_direction.theta()*180/M_PI;
		 double just_s3_dist_vtx = calcDistance(vtx_position.X(), vtx_position.Y(), vtx_position.Z(), just_s3_start_point.X(), just_s3_start_point.Y(), just_s3_start_point.Z() );

	     debug()<<"4Shower s3_nClusters = " <<just_s3_nclusters<<endmsg;
	     debug()<<"4Shower s3_Energy = " <<just_s3_energy<<endmsg;
	     debug()<<"4Shower s3_Theta = " <<just_s3_theta<<endmsg;
	     debug()<<"4Shower s3_Distance to Vertex = " <<just_s3_dist_vtx<<endmsg;

	     event->setIntData("FourShower_s3_nClusters", just_s3_nclusters);
	     event->setDoubleData("FourShower_s3_energy", just_s3_energy);
	     event->setDoubleData("FourShower_s3_theta", just_s3_theta);
	     event->setDoubleData("FourShower_s3_dist_vtx", just_s3_dist_vtx);
	 }
  }

  if (isShower4Good){
	 if (m_study_shower_recovery){
		 int just_s4_nclusters = shower4->nclusters();
		 double just_s4_energy = shower4->energy();
		 const Gaudi::XYZPoint& just_s4_start_point = shower4->startPoint();
		 const Gaudi::XYZVector& just_s4_direction = shower4->direction();
		 double just_s4_theta = just_s4_direction.theta()*180/M_PI;
		 double just_s4_dist_vtx = calcDistance(vtx_position.X(), vtx_position.Y(), vtx_position.Z(), just_s4_start_point.X(), just_s4_start_point.Y(), just_s4_start_point.Z() );

	     debug()<<"4Shower s4_nClusters = " <<just_s4_nclusters<<endmsg;
	     debug()<<"4Shower s4_Energy = " <<just_s4_energy<<endmsg;
	     debug()<<"4Shower s4_Theta = " <<just_s4_theta<<endmsg;
	     debug()<<"4Shower s4_Distance to Vertex = " <<just_s4_dist_vtx<<endmsg;

	     event->setIntData("FourShower_s4_nClusters", just_s4_nclusters);
	     event->setDoubleData("FourShower_s4_energy", just_s4_energy);
	     event->setDoubleData("FourShower_s4_theta", just_s4_theta);
	     event->setDoubleData("FourShower_s4_dist_vtx", just_s4_dist_vtx);
	 }
  }

  if (isShower1Good && isShower2Good && isShower3Good && isShower4Good){
    if (m_study_shower_recovery){
      // Shower 1
      int s1_nclusters = shower1->nclusters();
      double s1_energy = shower1->energy();
      const Gaudi::XYZPoint& s1_start_point = shower1->startPoint();
      const Gaudi::XYZVector& s1_direction = shower1->direction();
      double s1_theta = s1_direction.theta()*180/M_PI;
      double s1_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
					 s1_start_point.X(), s1_start_point.Y(), s1_start_point.Z());

      debug()<<"4Shower s1_nClusters = " <<s1_nclusters<<endmsg;
      debug()<<"4Shower s1_Energy = " <<s1_energy<<endmsg;
      debug()<<"4Shower s1_Theta = " <<s1_theta<<endmsg;
      debug()<<"4Shower s1_Distance to Vertex = " <<s1_dist_vtx<<endmsg;

      event->setIntData("FourShower_s1_nClusters", s1_nclusters);
      event->setDoubleData("FourShower_s1_energy", s1_energy);
      event->setDoubleData("FourShower_s1_theta", s1_theta);
      event->setDoubleData("FourShower_s1_dist_vtx", s1_dist_vtx);

      // Shower 2
      int s2_nclusters = shower2->nclusters();
      double s2_energy = shower2->energy();
      const Gaudi::XYZPoint& s2_start_point = shower2->startPoint();
      const Gaudi::XYZVector& s2_direction = shower2->direction();
      double s2_theta = s2_direction.theta()*180/M_PI;
      double s2_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
					 s2_start_point.X(), s2_start_point.Y(), s2_start_point.Z());

      debug()<<"4Shower s2_nClusters = " <<s2_nclusters<<endmsg;
      debug()<<"4Shower s2_Energy = " <<s2_energy<<endmsg;
      debug()<<"4Shower s2_Theta = " <<s2_theta<<endmsg;
      debug()<<"4Shower s2_Distance to Vertex = " <<s2_dist_vtx<<endmsg;

      event->setIntData("FourShower_s2_nClusters", s2_nclusters);
      event->setDoubleData("FourShower_s2_energy", s2_energy);
      event->setDoubleData("FourShower_s2_theta", s2_theta);
      event->setDoubleData("FourShower_s2_dist_vtx", s2_dist_vtx);

      // Shower 3
      int s3_nclusters = shower3->nclusters();
      double s3_energy = shower3->energy();
      const Gaudi::XYZPoint& s3_start_point = shower3->startPoint();
      const Gaudi::XYZVector& s3_direction = shower3->direction();
      double s3_theta = s3_direction.theta()*180/M_PI;
      double s3_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
					 s3_start_point.X(), s3_start_point.Y(), s3_start_point.Z());

      debug()<<"4Shower s3_nClusters = " <<s3_nclusters<<endmsg;
      debug()<<"4Shower s3_Energy = " <<s3_energy<<endmsg;
      debug()<<"4Shower s3_Theta = " <<s3_theta<<endmsg;
      debug()<<"4Shower s3_Distance to Vertex = " <<s3_dist_vtx<<endmsg;

      event->setIntData("FourShower_s3_nClusters", s3_nclusters);
      event->setDoubleData("FourShower_s3_energy", s3_energy);
      event->setDoubleData("FourShower_s3_theta", s3_theta);
      event->setDoubleData("FourShower_s3_dist_vtx", s3_dist_vtx);

      // Shower 4
      int s4_nclusters = shower4->nclusters();
      double s4_energy = shower4->energy();
      const Gaudi::XYZPoint& s4_start_point = shower4->startPoint();
      const Gaudi::XYZVector& s4_direction = shower4->direction();
      double s4_theta = s4_direction.theta()*180/M_PI;
      double s4_dist_vtx = calcDistance( vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
					 s4_start_point.X(), s4_start_point.Y(), s4_start_point.Z());

      debug()<<"4Shower s4_nClusters = " <<s4_nclusters<<endmsg;
      debug()<<"4Shower s4_Energy = " <<s4_energy<<endmsg;
      debug()<<"4Shower s4_Theta = " <<s4_theta<<endmsg;
      debug()<<"4Shower s4_Distance to Vertex = " <<s4_dist_vtx<<endmsg;

      event->setIntData("FourShower_s4_nClusters", s4_nclusters);
      event->setDoubleData("FourShower_s4_energy", s4_energy);
      event->setDoubleData("FourShower_s4_theta", s4_theta);
      event->setDoubleData("FourShower_s4_dist_vtx", s4_dist_vtx);

    }
    return true;
  }else{
    debug()<<"Not All of the 4 showers direction are good!"<<endmsg;
    return false;
  }

}

void CCProtonPi0::Save_4ShowerTruthMatch( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::GenMinInteraction* truthEvent) const
{
  // Shower 1
  Minerva::IDBlob* shower1 = foundBlobs[0];

  SmartRefVector<Minerva::IDCluster> shower1Clusters = shower1->clusters();
  DigitVectorTruthInfo shower1Info;
  shower1Info.ParseTruth(shower1Clusters, fTrajectoryMap, 0);

  truthEvent->setIntData("FourShower_s1_evis_most_pdg", shower1Info.GetMostEvisPdg());
  truthEvent->setDoubleData("FourShower_s1_evis_total_norm", shower1Info.GetTotalNormEnergy());
  truthEvent->setDoubleData("FourShower_s1_evis_total_truth", shower1Info.GetTotalTruthEnergy());
  truthEvent->setDoubleData("FourShower_s1_evis_pizero", shower1Info.GetEdepByPdg(PDG::pi0));
  truthEvent->setDoubleData("FourShower_s1_evis_piplus", shower1Info.GetEdepByPdg(PDG::pi));
  truthEvent->setDoubleData("FourShower_s1_evis_piminus", shower1Info.GetEdepByPdg(-(PDG::pi)));
  truthEvent->setDoubleData("FourShower_s1_evis_muon", shower1Info.GetEdepByPdg(PDG::muon));
  truthEvent->setDoubleData("FourShower_s1_evis_proton", shower1Info.GetEdepByPdg(PDG::proton));
  truthEvent->setDoubleData("FourShower_s1_evis_neutron", shower1Info.GetEdepByPdg(PDG::neutron));

  // Shower 2
  Minerva::IDBlob* shower2 = foundBlobs[1];

  SmartRefVector<Minerva::IDCluster> shower2Clusters = shower2->clusters();
  DigitVectorTruthInfo shower2Info;
  shower2Info.ParseTruth(shower2Clusters, fTrajectoryMap, 0);

  truthEvent->setIntData("FourShower_s2_evis_most_pdg", shower2Info.GetMostEvisPdg());
  truthEvent->setDoubleData("FourShower_s2_evis_total_norm", shower2Info.GetTotalNormEnergy());
  truthEvent->setDoubleData("FourShower_s2_evis_total_truth", shower2Info.GetTotalTruthEnergy());
  truthEvent->setDoubleData("FourShower_s2_evis_pizero", shower2Info.GetEdepByPdg(PDG::pi0));
  truthEvent->setDoubleData("FourShower_s2_evis_piplus", shower2Info.GetEdepByPdg(PDG::pi));
  truthEvent->setDoubleData("FourShower_s2_evis_piminus", shower2Info.GetEdepByPdg(-(PDG::pi)));
  truthEvent->setDoubleData("FourShower_s2_evis_muon", shower2Info.GetEdepByPdg(PDG::muon));
  truthEvent->setDoubleData("FourShower_s2_evis_proton", shower2Info.GetEdepByPdg(PDG::proton));
  truthEvent->setDoubleData("FourShower_s2_evis_neutron", shower2Info.GetEdepByPdg(PDG::neutron));

  // Shower 3
  Minerva::IDBlob* shower3 = foundBlobs[2];

  SmartRefVector<Minerva::IDCluster> shower3Clusters = shower3->clusters();
  DigitVectorTruthInfo shower3Info;
  shower3Info.ParseTruth(shower3Clusters, fTrajectoryMap, 0);

  truthEvent->setIntData("FourShower_s3_evis_most_pdg", shower3Info.GetMostEvisPdg());
  truthEvent->setDoubleData("FourShower_s3_evis_total_norm", shower3Info.GetTotalNormEnergy());
  truthEvent->setDoubleData("FourShower_s3_evis_total_truth", shower3Info.GetTotalTruthEnergy());
  truthEvent->setDoubleData("FourShower_s3_evis_pizero", shower3Info.GetEdepByPdg(PDG::pi0));
  truthEvent->setDoubleData("FourShower_s3_evis_piplus", shower3Info.GetEdepByPdg(PDG::pi));
  truthEvent->setDoubleData("FourShower_s3_evis_piminus", shower3Info.GetEdepByPdg(-(PDG::pi)));
  truthEvent->setDoubleData("FourShower_s3_evis_muon", shower3Info.GetEdepByPdg(PDG::muon));
  truthEvent->setDoubleData("FourShower_s3_evis_proton", shower3Info.GetEdepByPdg(PDG::proton));
  truthEvent->setDoubleData("FourShower_s3_evis_neutron", shower3Info.GetEdepByPdg(PDG::neutron));

    // Shower 4
  Minerva::IDBlob* shower4 = foundBlobs[3];

  SmartRefVector<Minerva::IDCluster> shower4Clusters = shower4->clusters();
  DigitVectorTruthInfo shower4Info;
  shower4Info.ParseTruth(shower4Clusters, fTrajectoryMap, 0);

  truthEvent->setIntData("FourShower_s4_evis_most_pdg", shower1Info.GetMostEvisPdg());
  truthEvent->setDoubleData("FourShower_s4_evis_total_norm", shower1Info.GetTotalNormEnergy());
  truthEvent->setDoubleData("FourShower_s4_evis_total_truth", shower1Info.GetTotalTruthEnergy());
  truthEvent->setDoubleData("FourShower_s4_evis_pizero", shower1Info.GetEdepByPdg(PDG::pi0));
  truthEvent->setDoubleData("FourShower_s4_evis_piplus", shower1Info.GetEdepByPdg(PDG::pi));
  truthEvent->setDoubleData("FourShower_s4_evis_piminus", shower1Info.GetEdepByPdg(-(PDG::pi)));
  truthEvent->setDoubleData("FourShower_s4_evis_muon", shower1Info.GetEdepByPdg(PDG::muon));
  truthEvent->setDoubleData("FourShower_s4_evis_proton", shower1Info.GetEdepByPdg(PDG::proton));
  truthEvent->setDoubleData("FourShower_s4_evis_neutron", shower1Info.GetEdepByPdg(PDG::neutron));
}

bool CCProtonPi0::RecoverShowers_Direction( std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event) const
{
    std::vector<Minerva::IDBlob*> bestBlobs;
    std::vector<Minerva::IDBlob*>::iterator iter_b;

    for (iter_b = foundBlobs.begin(); iter_b != foundBlobs.end(); ++iter_b) {
        bool isCurrentShowerGood = isShowerGood((*iter_b), event);
        if (isCurrentShowerGood){
            bestBlobs.push_back((*iter_b));
        }
    }

    // Check Size of Best Blobs
    bool isBlobsRecovered = ( bestBlobs.size() == 2 );
    if (isBlobsRecovered) foundBlobs.swap(bestBlobs);

    return isBlobsRecovered;
}

bool CCProtonPi0::isShowerGood(Minerva::IDBlob* shower, Minerva::PhysicsEvent* event) const
{
    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    bool goodPosition  = m_idHoughBlob->GetStartPosition(shower, vtx_position, true );
    bool goodDirection = m_idHoughBlob->GetDirection(shower, vtx_position );

    if (goodPosition) {
      event->setIntData("shower_has_goodPosition", 1);
      bool isGoodPosition = true;
    }
    
    if (goodDirection){
      event->setIntData("shower_has_goodDirection", 1);
      bool isGoodDirection = true;
    }

    bool isShowerGood = goodPosition && goodDirection;

    return isShowerGood;
}

bool CCProtonPi0::RecoverSingleShower_SmallAngle(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const 
{
    debug()<<"Trying to Recover Single Shower Events using Small Angle"<<endmsg;  

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Create AngleScan Object
    AngleScan angleScanAlg(usableClusters,vtx_position);
    angleScanAlg.AllowUVMatchWithMoreTolerance(m_AllowUVMatchWithMoreTolerance);
    angleScanAlg.SetUVMatchTolerance(m_UVMatchTolerance);
    angleScanAlg.AllowSmallConeAngle(true);
    angleScanAlg.DoReco();

    debug()<<"N(foundShowers) with Small Angle = "<<angleScanAlg.GetNCandidate()<<endmsg;

    if (angleScanAlg.GetNCandidate() == 2){
        debug()<<"Replacing foundBlobs with new ones"<<endmsg;
        // --------------------------------------------------------------------
        //   Blobs found by the AngleScan are not managed, delete them here and
        //   empty the container 
        // --------------------------------------------------------------------
        for (std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
                b != foundBlobs.end(); ++b) {
            delete *b;
        }
        foundBlobs.clear();    

        std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
        for (   std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin(); b != angleScanBlobs.end(); ++b) 
        {
            if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
            foundBlobs.push_back(*b);
        }

        return true;
    }else{ 
        return false;
    }
}

bool CCProtonPi0::RecoverSingleShower_View_U(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const 
{
    debug()<<"Trying to Recover Single Shower Events using U View"<<endmsg;  

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Create AngleScan Object
    AngleScan_U angleScanAlg(usableClusters,vtx_position);
    angleScanAlg.AllowUVMatchWithMoreTolerance(m_AllowUVMatchWithMoreTolerance);
    angleScanAlg.SetUVMatchTolerance(m_UVMatchTolerance);
    angleScanAlg.DoReco();

    debug()<<"N(foundShowers) with View U = "<<angleScanAlg.GetNCandidate()<<endmsg;

    if (angleScanAlg.GetNCandidate() == 2){
        debug()<<"Replacing foundBlobs with new ones"<<endmsg;
        // --------------------------------------------------------------------
        //   Blobs found by the AngleScan are not managed, delete them here and
        //   empty the container 
        // --------------------------------------------------------------------
        for (std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
                b != foundBlobs.end(); ++b) {
            delete *b;
        }
        foundBlobs.clear();    

        std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
        for (   std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin(); b != angleScanBlobs.end(); ++b) 
        {
            if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
            foundBlobs.push_back(*b);
        }

        return true;
    }else{ 
        return false;
    }
}

bool CCProtonPi0::RecoverSingleShower_View_V(SmartRefVector<Minerva::IDCluster> &usableClusters, std::vector<Minerva::IDBlob*> &foundBlobs, Minerva::PhysicsEvent *event)const 
{
    debug()<<"Trying to Recover Single Shower Events using V View"<<endmsg;  

    const Gaudi::XYZPoint& vtx_position = event->interactionVertex()->position();

    // Create AngleScan Object
    AngleScan_V angleScanAlg(usableClusters,vtx_position);
    angleScanAlg.AllowUVMatchWithMoreTolerance(m_AllowUVMatchWithMoreTolerance);
    angleScanAlg.SetUVMatchTolerance(m_UVMatchTolerance);
    angleScanAlg.DoReco();

    debug()<<"N(foundShowers) with View V = "<<angleScanAlg.GetNCandidate()<<endmsg;

    if (angleScanAlg.GetNCandidate() == 2){
        debug()<<"Replacing foundBlobs with new ones"<<endmsg;
        // --------------------------------------------------------------------
        //   Blobs found by the AngleScan are not managed, delete them here and
        //   empty the container 
        // --------------------------------------------------------------------
        for (std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
                b != foundBlobs.end(); ++b) {
            delete *b;
        }
        foundBlobs.clear();    

        std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
        for (   std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin(); b != angleScanBlobs.end(); ++b) 
        {
            if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
            foundBlobs.push_back(*b);
        }

        return true;
    }else{ 
        return false;
    }
}

double CCProtonPi0::GetTotalExtraEnergy(const Minerva::PhysicsEvent* event) const
{
    double extra_energy = 0.0;
    extra_energy += event->getDoubleData("Extra_Energy_Dispersed");
    extra_energy += event->getDoubleData("Extra_Energy_Muon");
    extra_energy += event->getDoubleData("Extra_Energy_Rejected");

    return extra_energy;
}


Gaudi::LorentzVector CCProtonPi0::Get_Neutrino_4P(double Enu) const
{
    const double theta = MinervaUnits::numi_beam_angle_rad;

    double Px = 0.0;
    double Py = Enu*sin(theta);
    double Pz = Enu*cos(theta);

    Gaudi::LorentzVector beam_4P(Px,Py,Pz,Enu);

    return beam_4P;
}

void CCProtonPi0::SetExtraEnergyTotal(Minerva::PhysicsEvent* event) const
{
    double muon = event->getDoubleData("Extra_Energy_Muon"); 
    double rejected = event->getDoubleData("Extra_Energy_Rejected"); 
    double leftover = event->getDoubleData("Extra_Energy_Leftover");

    double total = muon + rejected + leftover;

    debug()<<"Extra_Energy_Muon = "<<muon<<endmsg;
    debug()<<"Extra_Energy_Rejected = "<<rejected<<endmsg;
    debug()<<"Extra_Energy_Leftover = "<<leftover<<endmsg;
    debug()<<"Extra_Energy_Total = "<<total<<endmsg;
    
    event->setDoubleData("Extra_Energy_Total", total);
}



std::vector<double> CCProtonPi0::EstimateShowerEndPoint(const Gaudi::LorentzVector& shower_4P, const std::vector<double> init_pos, double max_z) const
{
    // Get Angles
    double theta = shower_4P.theta(); 
    double phi = shower_4P.phi(); 

    // Calc Rz (We know Rz from plane locations)
    double Rz = abs(max_z - init_pos[2]);

    // Calc R using Rz and theta
    double R = Rz / cos(theta);
    double Rx = R * sin(theta) * cos(phi);
    double Ry = R * sin(theta) * sin(phi);

    // Get Final Positions 
    double final_x = init_pos[0] + Rx;
    double final_y = init_pos[1] + Ry;
    double final_z = max_z;

    std::vector<double> final_pos;
    final_pos.push_back(final_x);
    final_pos.push_back(final_y);
    final_pos.push_back(final_z);

    return final_pos;
}


void CCProtonPi0::SaveTrajectories(Minerva::PhysicsEvent* event) const
{
    event->setIntData("detmc_ntrajectory", fTrajectoryMap.size());

    std::map<int, Minerva::TG4Trajectory*> skimmed_trajectories;
    std::map<int, Minerva::TG4Trajectory*> tier1_trajectories;
    for (std::map<int, Minerva::TG4Trajectory*>::iterator t
             = fTrajectoryMap.begin();
         t != fTrajectoryMap.end(); ++t) {
        if (t->second->GetParentId() == 0) {
            tier1_trajectories.insert(std::make_pair(t->first, t->second));
            skimmed_trajectories.insert(std::make_pair(t->first, t->second));
        }
        else break;

            /*
        std::cout << "\t"
                  << setw(8) << t->second->GetParentId()
                  << setw(8) << t->second->GetTrackId()
                  << setw(8) << t->second->GetPDGCode()
                  << std::endl;
            */
    }

    std::map<int, Minerva::TG4Trajectory*> tier2_trajectories;
    for (std::map<int, Minerva::TG4Trajectory*>::iterator t
             = fTrajectoryMap.begin();
         t != fTrajectoryMap.end(); ++t) {
        if (tier1_trajectories.find(t->second->GetParentId()) != tier1_trajectories.end()) {

            std::string proc = t->second->GetProcessName();
            if (proc == "muIoni" || proc == "eBrem" || proc == "hIoni") continue;

            if (t->second->GetPDGCode() > 1e9) continue;
            
            tier2_trajectories.insert(std::make_pair(t->first, t->second));
            skimmed_trajectories.insert(std::make_pair(t->first, t->second));

                /*
            std::cout << "\t\t"
                      << setw(8)  << t->second->GetParentId()
                      << setw(8)  << t->second->GetTrackId()
                      << setw(15) << t->second->GetPDGCode()
                      << setw(30) << t->second->GetProcessName()
                      << std::endl;
                */
            
        }
    }

    std::map<int, Minerva::TG4Trajectory*> tier3_trajectories;
    for (std::map<int, Minerva::TG4Trajectory*>::iterator t
             = fTrajectoryMap.begin();
         t != fTrajectoryMap.end(); ++t) {
        if (tier2_trajectories.find(t->second->GetParentId()) != tier2_trajectories.end()) {
            
            std::string proc = t->second->GetProcessName();
            if (proc == "muIoni" || proc == "eBrem" || proc == "hIoni") continue;

            if (t->second->GetPDGCode() > 1e9) continue;
                    
            tier3_trajectories.insert(std::make_pair(t->first, t->second));
                // delay adding tiered 3 trajectories, see below
                /*
            std::cout << "\t\t\t"
                      << setw(8)  << t->second->GetParentId()
                      << setw(8)  << t->second->GetTrackId()
                      << setw(15) << t->second->GetPDGCode()
                      << setw(30) << t->second->GetProcessName()
                      << std::endl;
                */
        }
    }

    int ntrajectory2 = tier1_trajectories.size() + tier2_trajectories.size() + tier3_trajectories.size();
    event->setIntData("detmc_ntrajectory2", ntrajectory2);

    if (ntrajectory2 > 100) {
            // too many trajectories, try to remove neutrons from tiered 2 trajectories
            // and re-run tiered 3
        std::map<int, Minerva::TG4Trajectory*> reduced_tier2_trajectories;
        for (std::map<int, Minerva::TG4Trajectory*>::iterator t = tier2_trajectories.begin();
             t != tier2_trajectories.end(); ++t) {
            if (t->second->GetPDGCode() == 2112) continue;
            reduced_tier2_trajectories.insert(std::make_pair(t->first, t->second));
        }

        tier2_trajectories = reduced_tier2_trajectories;

    }

        // add tiered 3 trajectories
    for (std::map<int, Minerva::TG4Trajectory*>::iterator t
             = fTrajectoryMap.begin();
         t != fTrajectoryMap.end(); ++t) {
        if (tier2_trajectories.find(t->second->GetParentId()) != tier2_trajectories.end()) {
            
            std::string proc = t->second->GetProcessName();
            if (proc == "muIoni" || proc == "eBrem" || proc == "hIoni") continue;
            
            if (t->second->GetPDGCode() > 1e9) continue;

            if (skimmed_trajectories.size() <= 100) skimmed_trajectories.insert(std::make_pair(t->first, t->second));

        }
    }

    std::cout << "\tTotal tiered trajectories: " << skimmed_trajectories.size() << " " << ntrajectory2
              << std::endl;
    
    debug() << "\twrite trajectory info to the ntuple " << endmsg;

    std::vector<int> detmc_traj_pdg;
    std::vector<int> detmc_traj_id;
    std::vector<int> detmc_traj_mother;
    std::vector<int> detmc_traj_proc;
    std::vector<int> detmc_traj_status;
    
    std::vector<double> detmc_traj_x0;
    std::vector<double> detmc_traj_y0;
    std::vector<double> detmc_traj_z0;
    std::vector<double> detmc_traj_t0;
    
    std::vector<double> detmc_traj_px0;
    std::vector<double> detmc_traj_py0;
    std::vector<double> detmc_traj_pz0;
    std::vector<double> detmc_traj_E0;

    std::vector<double> detmc_traj_xf;
    std::vector<double> detmc_traj_yf;
    std::vector<double> detmc_traj_zf;
    std::vector<double> detmc_traj_tf;
    
    std::vector<double> detmc_traj_pxf;
    std::vector<double> detmc_traj_pyf;
    std::vector<double> detmc_traj_pzf;
    std::vector<double> detmc_traj_Ef;

    std::vector<double> detmc_traj_prepxf;
    std::vector<double> detmc_traj_prepyf;
    std::vector<double> detmc_traj_prepzf;
    std::vector<double> detmc_traj_preEf;

    
    for (std::map<int, Minerva::TG4Trajectory*>::iterator t = skimmed_trajectories.begin();
         t != skimmed_trajectories.end(); ++t) {
        std::string proc = t->second->GetProcessName();

        int code = -1;

        if      (proc == "hadElastic")                code = 0;
        else if (proc == "NeutronInelastic")          code = 1;
        else if (proc == "ProtonInelastic")           code = 2;
        else if (proc == "PionMinusInelastic")        code = 3;
        else if (proc == "PionPlusInelastic")         code = 4;
        else if (proc == "CHIPSNuclearCaptureAtRest") code = 5;
        else if (proc == "Decay")                     code = 6;
        else if (proc == "conv")                      code = 7;
        else if (proc == "compt")                     code = 8;
        else if (proc == "Primary")                   code = 10;
        else code = 20;

            /*
        std::cout << "\t\t\t"
                  << setw(8)  << t->second->GetParentId()
                  << setw(8)  << t->second->GetTrackId()
                  << setw(15) << t->second->GetPDGCode()
                  << setw(30) << t->second->GetProcessName()
                  << setw(5)  << code
                  << std::endl;
            */
                  
        const Minerva::TG4Trajectory* trajectory = t->second;
        
        const int pdg    = trajectory->GetPDGCode();                detmc_traj_pdg.push_back(pdg);
        const int id     = trajectory->GetTrackId();                detmc_traj_id.push_back(id);
        const int mother = trajectory->GetParentId();               detmc_traj_mother.push_back(mother);
        const int status = trajectory->GetLeavingCode();            detmc_traj_status.push_back(status);
        const int proc_enum = code;                                 detmc_traj_proc.push_back(proc_enum);
                                                                    
        const double x0 = trajectory->GetInitialPosition().X();     detmc_traj_x0.push_back(x0);
        const double y0 = trajectory->GetInitialPosition().Y();     detmc_traj_y0.push_back(y0);
        const double z0 = trajectory->GetInitialPosition().Z();     detmc_traj_z0.push_back(z0);
        const double t0 = trajectory->GetInitialPosition().T();     detmc_traj_t0.push_back(t0);
        
        const double px0 = trajectory->GetInitialMomentum().Px();   detmc_traj_px0.push_back(px0);
        const double py0 = trajectory->GetInitialMomentum().Py();   detmc_traj_py0.push_back(py0);
        const double pz0 = trajectory->GetInitialMomentum().Pz();   detmc_traj_pz0.push_back(pz0);
        const double E0  = trajectory->GetInitialMomentum().E();    detmc_traj_E0.push_back(E0);

        const double xf = trajectory->GetFinalPosition().X();       detmc_traj_xf.push_back(xf);
        const double yf = trajectory->GetFinalPosition().Y();       detmc_traj_yf.push_back(yf);
        const double zf = trajectory->GetFinalPosition().Z();       detmc_traj_zf.push_back(zf);
        const double tf = trajectory->GetFinalPosition().T();       detmc_traj_tf.push_back(tf);

        Minerva::TG4TrajectoryPointExtraInfo* extraInfo =  dynamic_cast<Minerva::TG4TrajectoryPointExtraInfo*>(trajectory->GetTrajectoryPoints().back());
        assert(extraInfo);
        
        const double pxf = extraInfo->GetMomentum().Px();           detmc_traj_pxf.push_back(pxf);
        const double pyf = extraInfo->GetMomentum().Py();           detmc_traj_pyf.push_back(pyf);
        const double pzf = extraInfo->GetMomentum().Pz();           detmc_traj_pzf.push_back(pzf);
        const double Ef  = extraInfo->GetMomentum().E();            detmc_traj_Ef.push_back(Ef);
        
        const double prepxf = extraInfo->GetPreMomentum().Px();     detmc_traj_prepxf.push_back(prepxf);
        const double prepyf = extraInfo->GetPreMomentum().Py();     detmc_traj_prepyf.push_back(prepyf);
        const double prepzf = extraInfo->GetPreMomentum().Pz();     detmc_traj_prepzf.push_back(prepzf);
        const double preEf  = extraInfo->GetPreMomentum().E();      detmc_traj_preEf.push_back(preEf);


    } // loop over trajectories


    event->setContainerIntData("detmc_traj_pdg",    detmc_traj_pdg);
    event->setContainerIntData("detmc_traj_id",     detmc_traj_id);
    event->setContainerIntData("detmc_traj_mother", detmc_traj_mother);
    event->setContainerIntData("detmc_traj_proc",   detmc_traj_proc);
    event->setContainerIntData("detmc_traj_status", detmc_traj_status);

    event->setContainerDoubleData("detmc_traj_x0", detmc_traj_x0);
    event->setContainerDoubleData("detmc_traj_y0", detmc_traj_y0);
    event->setContainerDoubleData("detmc_traj_z0", detmc_traj_z0);
    event->setContainerDoubleData("detmc_traj_t0", detmc_traj_t0);

    event->setContainerDoubleData("detmc_traj_px0", detmc_traj_px0);
    event->setContainerDoubleData("detmc_traj_py0", detmc_traj_py0);
    event->setContainerDoubleData("detmc_traj_pz0", detmc_traj_pz0);
    event->setContainerDoubleData("detmc_traj_E0",  detmc_traj_E0);

    event->setContainerDoubleData("detmc_traj_xf", detmc_traj_xf);
    event->setContainerDoubleData("detmc_traj_yf", detmc_traj_yf);
    event->setContainerDoubleData("detmc_traj_zf", detmc_traj_zf);
    event->setContainerDoubleData("detmc_traj_tf", detmc_traj_tf);

    event->setContainerDoubleData("detmc_traj_pxf", detmc_traj_pxf);
    event->setContainerDoubleData("detmc_traj_pyf", detmc_traj_pyf);
    event->setContainerDoubleData("detmc_traj_pzf", detmc_traj_pzf);
    event->setContainerDoubleData("detmc_traj_Ef",  detmc_traj_Ef);

    event->setContainerDoubleData("detmc_traj_prepxf", detmc_traj_prepxf);
    event->setContainerDoubleData("detmc_traj_prepyf", detmc_traj_prepyf);
    event->setContainerDoubleData("detmc_traj_prepzf", detmc_traj_prepzf);
    event->setContainerDoubleData("detmc_traj_preEf",  detmc_traj_preEf);
    
}



#endif

