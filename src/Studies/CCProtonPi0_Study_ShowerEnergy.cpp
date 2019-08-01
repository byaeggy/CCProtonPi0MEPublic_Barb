/*
   See CCProtonPi0.h header for Class Information
*/
#ifndef CCProtonPi0_Study_ShowerEnergy_cpp 
#define CCProtonPi0_Study_ShowerEnergy_cpp 1

#include "../CCProtonPi0.h"

// Gaudi
//#include "GaudiKernel/PhysicalConstants.h"

// Minerva Analysis Framework
//#include "Event/GenMinHeader.h"
//#include "Event/MCIDDigit.h"
//#include "Event/MCHit.h"
//#include "RecInterfaces/IFiducialPointTool.h"


void CCProtonPi0::SaveBlobDigitEnergy(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
   
    std::vector<double> all_digit_E;
    std::vector<double> pi0_digit_E;
    std::vector<double> pi_digit_E;
    std::vector<double> proton_digit_E;
    std::vector<double> neutron_digit_E;
    std::vector<double> muon_digit_E;

    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {
        
        // Loop Over all Digits
        const SmartRefVector<Minerva::IDDigit>& digits = (*c)->digits();
        SmartRefVector<Minerva::IDDigit>::const_iterator d;
        for(d = digits.begin(); d != digits.end(); ++d){
            // Get mcdigit
            const Minerva::IDDigit* digit = *d;
            const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);
            
            // Check mcdigit
            if (!mcdigit) continue;

            int pdg = GetDigitPDG(mcdigit);
            
            double mcdigit_E = 0;
            // Get Hits from mcdigit
            const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
            if (hits.empty()) continue;
            SmartRefVector<Minerva::MCHit>::const_iterator h;
            // Loop Over all MCHits
            for ( h = hits.begin(); h != hits.end(); ++h){
                double hit_E = (*h)->energy();                
                mcdigit_E = mcdigit_E + hit_E;
            } 
            
            double digit_E = digit->normEnergy();

            debug()<<"MCDigit Energy = "<<mcdigit_E<<endmsg;
            debug()<<"Digit Energy = "<<digit_E<<endmsg;
            
            all_digit_E.push_back(digit_E);
            if (pdg == 111) pi0_digit_E.push_back(digit_E);
            else if (std::abs(pdg) == 211) pi_digit_E.push_back(digit_E);
            else if (pdg == 2212) proton_digit_E.push_back(digit_E);
            else if (pdg == 2112) neutron_digit_E.push_back(digit_E);
            else if (pdg == 13) muon_digit_E.push_back(digit_E);
        }
    }
    
    if (blobID == 1){
       event->setContainerDoubleData("gamma1_blob_all_digit_E", all_digit_E);
       event->setContainerDoubleData("gamma1_blob_pi0_digit_E", pi0_digit_E);
       event->setContainerDoubleData("gamma1_blob_pi_digit_E", pi_digit_E);
       event->setContainerDoubleData("gamma1_blob_proton_digit_E", proton_digit_E);
       event->setContainerDoubleData("gamma1_blob_neutron_digit_E", neutron_digit_E);
       event->setContainerDoubleData("gamma1_blob_muon_digit_E", muon_digit_E);
    }else{
       event->setContainerDoubleData("gamma2_blob_all_digit_E", all_digit_E);
       event->setContainerDoubleData("gamma2_blob_pi0_digit_E", pi0_digit_E);
       event->setContainerDoubleData("gamma2_blob_pi_digit_E", pi_digit_E);
       event->setContainerDoubleData("gamma2_blob_proton_digit_E", proton_digit_E);
       event->setContainerDoubleData("gamma2_blob_neutron_digit_E", neutron_digit_E);
       event->setContainerDoubleData("gamma2_blob_muon_digit_E", muon_digit_E);
    }
}

bool CCProtonPi0::isHitInsideSCAL(double x, double y, double z) const
{
    // All distances in mm
    double trkr_zMin = 5810.0; 
    double trkr_zMax = 8590.0; 
    double trkr_ApothemMax = 920.0;
    double scal_ApothemMax = 1070.0;

    bool isInsideTracker = FiducialPointTool->isFiducial(x,y,z, trkr_ApothemMax, trkr_zMin, trkr_zMax); 
    if (isInsideTracker){ 
        return false;
    }else{
        return FiducialPointTool->isFiducial(x,y,z, scal_ApothemMax, trkr_zMin, trkr_zMax); 
    }
}

void CCProtonPi0::SaveBlobTrueEvisBySubDetector(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    // True Evis by Subdetector
    double trkr_true_evis = 0;
    double scal_true_evis = 0;
    double ecal_true_evis = 0;
    double hcal_true_evis = 0;

    // nHits on Tracker or Side ECAL
    double center_nHits_all = 0;
    double center_nHits_trkr = 0;
    double center_nHits_scal = 0;

    double side_nHits_all = 0;
    double side_nHits_trkr = 0;
    double side_nHits_scal = 0;
    
    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {

        const double cluster_energy = (*c)->energy();

        // Tracker & Side ECAL 
        if ( (*c)->subdet() == Minerva::IDCluster::Tracker ){

            const SmartRefVector<Minerva::IDDigit>& centerDigits = (*c)->centralDigits();
            const SmartRefVector<Minerva::IDDigit>& sideDigits = (*c)->sideEcalDigits();

            // Loop over ID Digits
            for(d = centerDigits.begin(); d != centerDigits.end(); ++d){
    
                // Get mcdigit
                const Minerva::IDDigit* digit = *d;
                const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                // Check mcdigit
                if (!mcdigit) continue;

                // Get Hits from mcdigit
                const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                if (hits.empty()) continue;
                SmartRefVector<Minerva::MCHit>::const_iterator h;
                // Loop Over all MCHits
                for ( h = hits.begin(); h != hits.end(); ++h){
                    // Get Hit Position (x,y,z)
                    double x = (*h)->StopX();    
                    double y = (*h)->StopY();    
                    double z = (*h)->StopZ(); 

                    // Get Hit Visible Energy
                    double hit_E = (*h)->energy();
                    
                    center_nHits_all++;
                    if ( isHitInsideSCAL(x,y,z) ) {
                        center_nHits_scal++;
                        scal_true_evis += hit_E;
                    }else{ 
                        center_nHits_trkr++;
                        trkr_true_evis += hit_E;
                    }
                }
            }

            // Loop over Side Digits
            for(d = sideDigits.begin(); d != sideDigits.end(); ++d){

                // Get mcdigit
                const Minerva::IDDigit* digit = *d;
                const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                // Check mcdigit
                if (!mcdigit) continue;

                // Get Hits from mcdigit
                const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                if (hits.empty()) continue;
                SmartRefVector<Minerva::MCHit>::const_iterator h;
                // Loop Over all MCHits
                for ( h = hits.begin(); h != hits.end(); ++h){
                    // Get Hit Position (x,y,z)
                    double x = (*h)->StopX();    
                    double y = (*h)->StopY();    
                    double z = (*h)->StopZ(); 

                    // Get Hit Visible Energy
                    double hit_E = (*h)->energy();

                    side_nHits_all++;
                    if ( isHitInsideSCAL(x,y,z) ) {
                        side_nHits_scal++;
                        scal_true_evis += hit_E;
                    }else{ 
                        side_nHits_trkr++;
                        trkr_true_evis += hit_E;
                    }
                }
            } 
        }

        // Downstream ECAL
        if ( (*c)->subdet() == Minerva::IDCluster::ECAL ){ 
            ecal_true_evis += cluster_energy;
        }
        
        // HCAL
        if ( (*c)->subdet() == Minerva::IDCluster::HCAL ){ 
            hcal_true_evis += cluster_energy;
        }
    
    } 
    
    debug()<<"Blob True Evis by SubDetector"<<endmsg;
    debug()<<"trkr_true_evis = "<<trkr_true_evis<<endmsg;
    debug()<<"scal_true_evis = "<<scal_true_evis<<endmsg;
    debug()<<"ecal_true_evis = "<<ecal_true_evis<<endmsg;
    debug()<<"hcal_true_evis = "<<hcal_true_evis<<endmsg;

    debug()<<"Tracker Hit Distribution"<<endmsg;
    debug()<<"center_nHits_all = "<<center_nHits_all<<endmsg;
    debug()<<"center_nHits_trkr = "<<center_nHits_trkr<<endmsg;
    debug()<<"center_nHits_scal = "<<center_nHits_scal<<endmsg;
    debug()<<"side_nHits_all = "<<side_nHits_all<<endmsg;
    debug()<<"side_nHits_trkr = "<<side_nHits_trkr<<endmsg;
    debug()<<"side_nHits_scal = "<<side_nHits_scal<<endmsg;

    if (blobID == 1){
        event->setDoubleData("gamma1_trkr_true_evis", trkr_true_evis);
        event->setDoubleData("gamma1_scal_true_evis", scal_true_evis);
        event->setDoubleData("gamma1_ecal_true_evis", ecal_true_evis);
        event->setDoubleData("gamma1_hcal_true_evis", hcal_true_evis);
        
        event->setDoubleData("gamma1_center_nHits_all", center_nHits_all);
        event->setDoubleData("gamma1_center_nHits_trkr", center_nHits_trkr);
        event->setDoubleData("gamma1_center_nHits_scal", center_nHits_scal);
        event->setDoubleData("gamma1_side_nHits_all", side_nHits_all);
        event->setDoubleData("gamma1_side_nHits_trkr", side_nHits_trkr);
        event->setDoubleData("gamma1_side_nHits_scal", side_nHits_scal);
    }else{
        event->setDoubleData("gamma2_trkr_true_evis", trkr_true_evis);
        event->setDoubleData("gamma2_scal_true_evis", scal_true_evis);
        event->setDoubleData("gamma2_ecal_true_evis", ecal_true_evis);
        event->setDoubleData("gamma2_hcal_true_evis", hcal_true_evis);
        
        event->setDoubleData("gamma2_center_nHits_all", center_nHits_all);
        event->setDoubleData("gamma2_center_nHits_trkr", center_nHits_trkr);
        event->setDoubleData("gamma2_center_nHits_scal", center_nHits_scal);
        event->setDoubleData("gamma2_side_nHits_all", side_nHits_all);
        event->setDoubleData("gamma2_side_nHits_trkr", side_nHits_trkr);
        event->setDoubleData("gamma2_side_nHits_scal", side_nHits_scal);
    }
}


void CCProtonPi0::SaveSCALHits(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    // nHits on Tracker or Side ECAL
    double trkr_nHits_true = 0;
    double trkr_nHits_reco = 0;

    double scal_nHits_true = 0;
    double scal_nHits_reco = 0;
    
    double SCAL_minZ = Find_SCAL_MinZ(blob);

    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {

        // Tracker & Side ECAL 
        if ( (*c)->subdet() == Minerva::IDCluster::Tracker ){

            double cluster_z = (*c)->z();
            const SmartRefVector<Minerva::IDDigit>& allDigits = (*c)->digits();

            // Loop over All Digits
            for(d = allDigits.begin(); d != allDigits.end(); ++d){
               
                // Get mcdigit
                const Minerva::IDDigit* digit = *d;
                const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                // Check mcdigit
                if (!mcdigit) continue;

                // Get Hits from mcdigit
                const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                if (hits.empty()) continue;
                SmartRefVector<Minerva::MCHit>::const_iterator h;
           
                // Loop Over all MCHits
                for ( h = hits.begin(); h != hits.end(); ++h){
                    // Get Hit Position (x,y,z)
                    double x = (*h)->StopX();    
                    double y = (*h)->StopY();    
                    double z = (*h)->StopZ(); 
                
                    // Count True Hits
                    if ( isHitInsideSCAL(x,y,z) ) {
                        scal_nHits_true++;
                    }else{ 
                        trkr_nHits_true++;
                    }

                    // Count Reco Hits
                    if ( cluster_z > SCAL_minZ ){
                        scal_nHits_reco++;
                    }else{
                        trkr_nHits_reco++;
                    }
                }
            }
        }
    } 
    
    debug()<<"Tracker Hit Distribution"<<endmsg;
    debug()<<"trkr_nHits_true = "<<trkr_nHits_true<<endmsg;
    debug()<<"trkr_nHits_reco = "<<trkr_nHits_reco<<endmsg;
    debug()<<"scal_nHits_true = "<<scal_nHits_true<<endmsg;
    debug()<<"scal_nHits_reco = "<<scal_nHits_reco<<endmsg;

    if (blobID == 1){
        event->setDoubleData("gamma1_trkr_nHits_true", trkr_nHits_true);
        event->setDoubleData("gamma1_trkr_nHits_reco", trkr_nHits_reco);
        event->setDoubleData("gamma1_scal_nHits_true", scal_nHits_true);
        event->setDoubleData("gamma1_scal_nHits_reco", scal_nHits_reco);
    }else{
        event->setDoubleData("gamma2_trkr_nHits_true", trkr_nHits_true);
        event->setDoubleData("gamma2_trkr_nHits_reco", trkr_nHits_reco);
        event->setDoubleData("gamma2_scal_nHits_true", scal_nHits_true);
        event->setDoubleData("gamma2_scal_nHits_reco", scal_nHits_reco);
    }
}

void CCProtonPi0::SaveSCALHits_Improved(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    // nHits on Tracker or Side ECAL
    double trkr_nHits_true = 0;
    double trkr_nHits_reco = 0;

    double scal_nHits_true = 0;
    double scal_nHits_reco = 0;
    
    double SCAL_minZ = Find_SCAL_MinZ(blob);
    bool use_minZ = Use_SCAL_minZ(blob,SCAL_minZ);

    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {

        // Tracker & Side ECAL 
        if ( (*c)->subdet() == Minerva::IDCluster::Tracker ){
            if (use_minZ){
            
                double cluster_z = (*c)->z();
                const SmartRefVector<Minerva::IDDigit>& allDigits = (*c)->digits();

                debug()<<"Using minZ = "<<SCAL_minZ<<endmsg;
                // Loop over All Digits
                for(d = allDigits.begin(); d != allDigits.end(); ++d){

                    // Get mcdigit
                    const Minerva::IDDigit* digit = *d;
                    const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                    // Check mcdigit
                    if (!mcdigit) continue;

                    // Get Hits from mcdigit
                    const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                    if (hits.empty()) continue;
                    SmartRefVector<Minerva::MCHit>::const_iterator h;

                    // Loop Over all MCHits
                    for ( h = hits.begin(); h != hits.end(); ++h){
                        // Get Hit Position (x,y,z)
                        double x = (*h)->StopX();    
                        double y = (*h)->StopY();    
                        double z = (*h)->StopZ(); 

                        // Count True Hits
                        if ( isHitInsideSCAL(x,y,z) ) {
                            scal_nHits_true++;
                        }else{ 
                            trkr_nHits_true++;
                        }

                        // Count Reco Hits
                        if ( cluster_z > SCAL_minZ ){
                            scal_nHits_reco++;
                        }else{
                            trkr_nHits_reco++;
                        }
                    }
                }
            }else{
                
                debug()<<"NOT using minZ = "<<SCAL_minZ<<endmsg;
                const SmartRefVector<Minerva::IDDigit>& centerDigits = (*c)->centralDigits();
                const SmartRefVector<Minerva::IDDigit>& sideDigits = (*c)->sideEcalDigits();

                // Loop over ID Digits
                for(d = centerDigits.begin(); d != centerDigits.end(); ++d){

                    // Get mcdigit
                    const Minerva::IDDigit* digit = *d;
                    const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                    // Check mcdigit
                    if (!mcdigit) continue;

                    // Get Hits from mcdigit
                    const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                    if (hits.empty()) continue;
                    SmartRefVector<Minerva::MCHit>::const_iterator h;
                    // Loop Over all MCHits
                    for ( h = hits.begin(); h != hits.end(); ++h){
                        // Get Hit Position (x,y,z)
                        double x = (*h)->StopX();    
                        double y = (*h)->StopY();    
                        double z = (*h)->StopZ(); 

                        trkr_nHits_reco++;
                        if ( isHitInsideSCAL(x,y,z) ) {
                            scal_nHits_true++;
                        }else{ 
                            trkr_nHits_true++;
                        }
                    }
                }

                // Loop over Side Digits
                for(d = sideDigits.begin(); d != sideDigits.end(); ++d){

                    // Get mcdigit
                    const Minerva::IDDigit* digit = *d;
                    const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);

                    // Check mcdigit
                    if (!mcdigit) continue;

                    // Get Hits from mcdigit
                    const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
                    if (hits.empty()) continue;
                    SmartRefVector<Minerva::MCHit>::const_iterator h;
                    // Loop Over all MCHits
                    for ( h = hits.begin(); h != hits.end(); ++h){
                        // Get Hit Position (x,y,z)
                        double x = (*h)->StopX();    
                        double y = (*h)->StopY();    
                        double z = (*h)->StopZ(); 

                        scal_nHits_reco++;
                        
                        if ( isHitInsideSCAL(x,y,z) ) {
                            scal_nHits_true++;
                        }else{ 
                            trkr_nHits_true++;
                        }
                    }
                } 
            }
        }
    } 

    debug()<<"Tracker Hit Distribution"<<endmsg;
    debug()<<"improved_trkr_nHits_true = "<<trkr_nHits_true<<endmsg;
    debug()<<"improved_trkr_nHits_reco = "<<trkr_nHits_reco<<endmsg;
    debug()<<"improved_scal_nHits_true = "<<scal_nHits_true<<endmsg;
    debug()<<"improved_scal_nHits_reco = "<<scal_nHits_reco<<endmsg;

    if (blobID == 1){
        event->setDoubleData("gamma1_improved_trkr_nHits_true", trkr_nHits_true);
        event->setDoubleData("gamma1_improved_trkr_nHits_reco", trkr_nHits_reco);
        event->setDoubleData("gamma1_improved_scal_nHits_true", scal_nHits_true);
        event->setDoubleData("gamma1_improved_scal_nHits_reco", scal_nHits_reco);
    }else{
        event->setDoubleData("gamma2_improved_trkr_nHits_true", trkr_nHits_true);
        event->setDoubleData("gamma2_improved_trkr_nHits_reco", trkr_nHits_reco);
        event->setDoubleData("gamma2_improved_scal_nHits_true", scal_nHits_true);
        event->setDoubleData("gamma2_improved_scal_nHits_reco", scal_nHits_reco);
    }
}

double CCProtonPi0::Find_SCAL_MinZ(Minerva::IDBlob* blob) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    double minZ = 999999;

    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {
        // Tracker & Side ECAL 
        if ( (*c)->subdet() == Minerva::IDCluster::Tracker ){

            // Get Side Digits inside the cluster
            const SmartRefVector<Minerva::IDDigit>& sideDigits = (*c)->sideEcalDigits();
            unsigned int nSideDigits = sideDigits.size();

            if ( nSideDigits > 0){
                double cluster_z = (*c)->z();
                if( cluster_z < minZ ) minZ = cluster_z; 
            }
        }
    } 
    return minZ;
}

bool CCProtonPi0::Use_SCAL_minZ(Minerva::IDBlob* blob, double minZ) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    bool isView_X = false;
    bool useSCAL_minZ = false;

    // if there is no minZ return false
    if (minZ == 999999) return false;

    // Loop over all Clusters and check minZ at X plane or NOT
    for (c = clusters.begin(); c != clusters.end(); ++c) {
        double cluster_z = (*c)->z();
        if (cluster_z == minZ){
            if ( (*c)->view() == Minerva::IDCluster::X ){ 
                debug()<<"minZ at X Plane"<<endmsg;
                isView_X = true;
            }else{
                debug()<<"minZ at U or V Plane"<<endmsg;
                isView_X = false;
                useSCAL_minZ = true;
            }
            break;
        }
    } 

    /* 
     * If minZ at X Plane
     *  Loop over all clusters again and save X View Z and N(SideDigits)
     *      for Z > minZ
     */
    if (isView_X){
        std::vector<double> x_cluster_z;
        std::vector<double> x_side_digits;

        for (c = clusters.begin(); c != clusters.end(); ++c) {

            if ( (*c)->view() == Minerva::IDCluster::X ){ 
                double cluster_z = (*c)->z();
                SmartRefVector< Minerva::IDDigit > sideDigits = (*c)->sideEcalDigits();
                // Save the ones after minZ
                if (cluster_z > minZ){
                    x_cluster_z.push_back(cluster_z);
                    x_side_digits.push_back(sideDigits.size());
                }
            }
        }    

        // Now Locate 2nd minZ
        debug()<<"Locating 2nd minZ"<<endmsg;
        // Return TRUE if there are no planes after minZ
        if (x_cluster_z.size() == 0){
            debug()<<"minZ plane was the last X plane - returning TRUE"<<endmsg;
            return true;
        }
        // Find 2nd minZ and check whether it has a SCAL Hit
        // Minimum Z in the array is the 2nd minZ
        double minZ_2 = 999999;
        double minZ_2_nDigits = -1;
        for (unsigned int i = 0; i < x_cluster_z.size(); i++){
            debug()<<"cluster_z = "<<x_cluster_z[i]<<endmsg;
            if (x_cluster_z[i] < minZ_2 ){
                minZ_2 = x_cluster_z[i];
                minZ_2_nDigits = x_side_digits[i];
            }
        }
        debug()<<"minZ_2 = "<<minZ_2<<" minZ_2_nDigits = "<<minZ_2_nDigits<<endmsg;
        
        if (minZ_2_nDigits > 0) useSCAL_minZ = true;
        else useSCAL_minZ = false;
    }

    return useSCAL_minZ;

}

void CCProtonPi0::SaveSCAL_minZ_Info(Minerva::PhysicsEvent *event, Minerva::IDBlob* blob, int blobID) const
{
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    SmartRefVector<Minerva::IDCluster>::const_iterator c;
    SmartRefVector<Minerva::IDDigit>::const_iterator d;

    double SCAL_minZ = Find_SCAL_MinZ(blob);
    double cluster_evis = -9.9;
    double cluster_nDigits = -9.9;

    // Loop over all Clusters
    for (c = clusters.begin(); c != clusters.end(); ++c) {
        double cluster_z = (*c)->z();
        
        if (cluster_z == SCAL_minZ){
            cluster_evis = (*c)->energy();
            cluster_nDigits = (*c)->iddigs();
        }
    } 
    
    debug()<<"SCAL minZ Info"<<endmsg;
    debug()<<"SCAL_minZ = "<<SCAL_minZ<<endmsg;
    debug()<<"cluster_evis = "<<cluster_evis<<endmsg;
    debug()<<"cluster_nDigits = "<<cluster_nDigits<<endmsg;

    if (blobID == 1){
        event->setDoubleData("gamma1_scal_minZ_evis", cluster_evis);
        event->setDoubleData("gamma1_scal_minZ_nDigits", cluster_nDigits);
    }else{
        event->setDoubleData("gamma2_scal_minZ_evis", cluster_evis);
        event->setDoubleData("gamma2_scal_minZ_nDigits", cluster_nDigits);
    }

}

#endif  

