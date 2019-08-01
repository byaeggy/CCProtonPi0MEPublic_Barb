#ifndef AngleScan_U_cpp
#define AngleScan_U_cpp

#undef NDEBUG
#include <cassert>
#include <cmath>

#include <TMath.h>
#include <TH1F.h>

#include <Event/IDCluster.h>
#include <Event/IDBlob.h>
#include <Event/VectorTypeDefs.h>
#include <Event/Vertex.h>

#include "AngleScan_U.h"
#include "ClusterVectorInfo.h"

namespace {
    struct greaterShower : public std::binary_function <
        SmartRefVector<Minerva::IDCluster>,
        SmartRefVector<Minerva::IDCluster>,
        bool > {

        bool operator()(const SmartRefVector<Minerva::IDCluster>& lhs,
                        const SmartRefVector<Minerva::IDCluster>& rhs) const {
            
            return lhs.size() > rhs.size();
        }

    };
}


AngleScan_U::AngleScan_U(const SmartRefVector<Minerva::IDCluster>& clusters,
                     const Gaudi::XYZPoint& vertex)
    : fUVMatchTolerance(10.0),
      fUVMatchMoreTolerance(100.0),
      fAllowUVMatchWithMoreTolerance(true)
{
    std::copy(clusters.begin(), clusters.end(), std::back_inserter(fAllClusters));
    
    fX = vertex.X();
    fY = vertex.Y();
    fZ = vertex.Z();
    fU = -fY*sqrt(3.)/2 + fX/2;
    fV =  fY*sqrt(3.)/2 + fX/2;
    
    Initialize();
}

void AngleScan_U::Initialize()
{

    std::cout << "    AngleScan_U::Initialize " << std::endl;

    ClusterVectorInfo clusterVectorInfo(fAllClusters,true,true,true);
    fXClusters = clusterVectorInfo.GetXClusters();
    fUClusters = clusterVectorInfo.GetUClusters();
    fVClusters = clusterVectorInfo.GetVClusters();

    std::cout << "      AngleScan_U:: x size = " << fXClusters.size() << std::endl;
    std::cout << "      AngleScan_U:: u size = " << fUClusters.size() << std::endl;
    std::cout << "      AngleScan_U:: v size = " << fVClusters.size() << std::endl;

        /* Copy to the working containers */
    fRemainingXClusters = fXClusters;
    fRemainingUClusters = fUClusters;
    fRemainingVClusters = fVClusters;

    assert(fRemainingXClusters.size() == fXClusters.size());
    assert(fRemainingUClusters.size() == fUClusters.size());
    assert(fRemainingVClusters.size() == fVClusters.size());
}

void AngleScan_U::BuildThetaHistogram()
{
    const unsigned int N = 90;
    const double lower = -180.0;
    const double upper = +180.0;
    fTheta = new TH1F("theta","Theta dist", N, lower, upper);
    
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = fRemainingUClusters.begin();
         c != fRemainingUClusters.end(); ++c) {
        const double dZ = (*c)->z() - fZ;
        const double dU = (*c)->position() - fU;
        const double theta = std::atan2(dU,dZ)*TMath::RadToDeg();
        const double w = (*c)->pe();
        
        fTheta->Fill(theta,w);
    }
}

void AngleScan_U::FindPeaks()
{
    const unsigned int N = 90;
        /* Detect and save lower and upper edges around peaks in the histogram */
    const double width = fTheta->GetBinWidth(1);      /* Use first bin since fixed bin size */
    for (unsigned int bin = 1; bin <= N; bin++){
        if (fTheta->GetBinContent(bin) > 15 ){        /* Peak is detected */
            int Limitbin = GetLimitBin(fTheta, bin ); /* Find the upper edge of the peak */
            
            const double lower_edge = fTheta->GetBinCenter(bin) - 0.5*width;
            const double upper_edge = fTheta->GetBinCenter(Limitbin) + 0.5*width;
            
            TVector2 peak(lower_edge,upper_edge);

            fPeaks.push_back(peak);
                        
            bin = Limitbin;  /* Finding the next peak starting from this peak's upper edge */
        }
    }

    delete fTheta;
}

void AngleScan_U::FormUShowerCand() 
{
    for (std::vector<TVector2>::const_iterator peak = fPeaks.begin();
         peak != fPeaks.end(); ++peak) {
        const double lower_edge = peak->X();
        const double upper_edge = peak->Y();
        const double zmin       = 4500.0;
        const double zmax       = 10000.0;

        SmartRefVector<Minerva::IDCluster> showerCand;
        coneView(fRemainingUClusters, showerCand, fU, lower_edge, upper_edge, zmin, zmax);
        
        if (showerCand.empty()) continue;

        if (showerCand.size() == 1 && showerCand.front()->pe() > 30) {

            fUShowerCandidates.push_back(showerCand);
            fGoodPeaks.push_back(*peak);
            
        } else if (showerCand.size() > 1) {

            fUShowerCandidates.push_back(showerCand);
            fGoodPeaks.push_back(*peak);
        }
        
    }

        // Calculate the distance from the shower candidates to the event vertex
        // Try two definitions of distance: closest and energy weighted
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
             = fUShowerCandidates.begin();
         s != fUShowerCandidates.end(); ++s) {

        SmartRefVector<Minerva::IDCluster>& ushowerCand = *s;

        double d_min = 1.e6;
        double d_weighted = 0.0;
        double total_energy = 0.0;
        for (SmartRefVector<Minerva::IDCluster>::iterator c = ushowerCand.begin();
             c != ushowerCand.end(); ++c) {

            double u = (*c)->position();
            double z = (*c)->z();
            double d = std::sqrt(std::pow(u-fU,2) + std::pow(z-fZ,2));

            if (d < d_min) {
                d_min = d;
            }

            d_weighted   += d * (*c)->energy();
            total_energy += (*c)->energy();
        }

        fUShowerClosestDistances.push_back(d_min);
        fUShowerWeightedDistances.push_back(d_weighted/total_energy);
    }
    
    std::sort(fUShowerCandidates.begin(),fUShowerCandidates.end(), greaterShower());
    
}

void AngleScan_U::FormXUVShowerCand() 
{

    std::vector<SmartRefVector<Minerva::IDCluster> > nogrowShowerCandidates;
    std::cout << "AngleScan_U::total U candidates: " << fUShowerCandidates.size() << std::endl;
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
             = fUShowerCandidates.begin();
         s != fUShowerCandidates.end(); ++s) {
        std::cout << "\tAngleScan_U::X cand: " << s->size() << std::endl;
        SmartRefVector<Minerva::IDCluster>& ushowerCand = *s;
        
        SmartRefVector<Minerva::IDCluster> showerCand = ushowerCand;
        double zmin = +1e6;
        double zmax = -1e6;
        double ztot = 0.0;
        for (SmartRefVector<Minerva::IDCluster>::iterator c = ushowerCand.begin();
             c != ushowerCand.end(); ++c) {
            zmin  = std::min(zmin,(*c)->z());
            zmax  = std::max(zmax,(*c)->z());
            ztot += (*c)->z();
        }

        if (zmax-zmin < 50.0) {
            const double zcenter = ztot/ushowerCand.size();
            std::cout << "\t AngleScan_U::Candidate with small z extend at: " << zcenter << std::endl;
            SmartRefVector<Minerva::IDCluster> xclusters_tmp;
            SmartRefVector<Minerva::IDCluster> vclusters_tmp;
            fRemainingXClusters.swap(xclusters_tmp);
            fRemainingVClusters.swap(vclusters_tmp);

            for (SmartRefVector<Minerva::IDCluster>::iterator c = xclusters_tmp.begin();
                 c != xclusters_tmp.end(); ++c) {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan_U::adding a U cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingXClusters.push_back(*c);
            }
                        
            for (SmartRefVector<Minerva::IDCluster>::iterator c = vclusters_tmp.begin();
                 c != vclusters_tmp.end(); ++c) {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan_U::adding a V cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingVClusters.push_back(*c);
            }
            
            fShowerCandidates.push_back(showerCand);

        } else {
            std::cout << "\t AngleScan_U::Candidate with large z extent " << std::endl;
            addClustersToBlob(ushowerCand,fRemainingXClusters,fRemainingVClusters,showerCand,
                              fUVMatchTolerance);
            
            if (showerCand.size() > ushowerCand.size()) fShowerCandidates.push_back(showerCand);
            else if (ushowerCand.size() >= 3) nogrowShowerCandidates.push_back(ushowerCand);
            else {
                std::cout << "\t AngleScan_U::Throw away U shower candidate" << std::endl;
            }
            
        }
    }

    std::cout << "AngleScan_U::Shower candidates: " << fShowerCandidates.size() << std::endl;
    std::cout << "AngleScan_U::Nogrow candidates: " << nogrowShowerCandidates.size() << std::endl;

    if (fAllowUVMatchWithMoreTolerance) {
        for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
                 = nogrowShowerCandidates.begin();
             s != nogrowShowerCandidates.end(); ++s) {
            SmartRefVector<Minerva::IDCluster>& ushowerCand = *s; 
            SmartRefVector<Minerva::IDCluster>  showerCand = ushowerCand;
            
            addClustersToBlob(ushowerCand,fRemainingXClusters,fRemainingVClusters,showerCand,
                              fUVMatchMoreTolerance);
            
            fShowerCandidates.push_back(showerCand);
        }
    }
    
    std::cout << "AngleScan_U::Shower candidates(10x): " << fShowerCandidates.size() << std::endl;
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s = fShowerCandidates.begin();
         s != fShowerCandidates.end(); ++s) {
        std::cout << "\tAngleScan_U:: cand: " << std::distance(fShowerCandidates.begin(),s) << " "
                  << fRemainingXClusters.size() << " " << fRemainingVClusters.size()
                  << std::endl;
        completeView(fRemainingXClusters, *s, fX);
        completeView(fRemainingVClusters, *s, fV);
    }

    std::cout << "AngleScan_U::Final showers: " << fShowerCandidates.size() << std::endl;
    std::sort(fShowerCandidates.begin(), fShowerCandidates.end(), greaterShower());
    
    std::copy(fRemainingUClusters.begin(),fRemainingUClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingXClusters.begin(),fRemainingXClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingVClusters.begin(),fRemainingVClusters.end(), std::back_inserter(fRemainingClusters));

}

void AngleScan_U::DoReco()
{
    std::cout << "AngleScan_U::DoReco() " << std::endl;
    BuildThetaHistogram();
    FindPeaks();
    FormUShowerCand();
    FormXUVShowerCand();
}

int AngleScan_U::GetLimitBin(const TH1F *hMax, int n_bin )const
{

    int max=n_bin, x = 1, count = 0;
    while ( count < 1 ){
        if ( hMax->GetBinContent(n_bin + x) > 0 ) max = n_bin + x;
        else count++;
        x++;
    }
    
    return max;

}
void AngleScan_U::addClustersToBlob(SmartRefVector<Minerva::IDCluster>& xshowerCand,
                                  SmartRefVector<Minerva::IDCluster>& uclusters,
                                  SmartRefVector<Minerva::IDCluster>& vclusters,
                                  SmartRefVector<Minerva::IDCluster>& showerCand,
                                  double epsilon)
{

    for (SmartRefVector<Minerva::IDCluster>::iterator c = xshowerCand.begin();
         c != xshowerCand.end(); ++c) {
        Minerva::IDCluster* cluster_x = *c;
        double min = 1e3;
        SmartRefVector<Minerva::IDCluster>::iterator ucluster = uclusters.end();
        SmartRefVector<Minerva::IDCluster>::iterator vcluster = vclusters.end();
        for (SmartRefVector<Minerva::IDCluster>::iterator itU = uclusters.begin();
             itU != uclusters.end(); ++itU) {

            if (std::abs( cluster_x->z() - (*itU)->z() ) > 50.0 ) continue;

            for (SmartRefVector<Minerva::IDCluster>::iterator itV = vclusters.begin();
                 itV != vclusters.end(); ++itV) {

                if ( std::abs( cluster_x->z() - (*itV)->z() ) > 50.0 ) continue;

                double delta = std::abs((*itU)->tpos1()+(*itU)->tpos2()+        /* |u+v-x| */
                                        (*itV)->tpos1()+(*itV)->tpos2()-
                                        cluster_x->tpos1()-cluster_x->tpos2());
                if ( delta < min ) { 
                    min = delta;
                    ucluster = itU;
                    vcluster = itV;
                }
            }
        }

        if (min <= epsilon && (ucluster != uclusters.end() && vcluster != vclusters.end())) {
            showerCand.push_back(*ucluster);
            showerCand.push_back(*vcluster);

            uclusters.erase(ucluster);
            vclusters.erase(vcluster);
        }
        
    }
    
}

void AngleScan_U::completeView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                             SmartRefVector<Minerva::IDCluster>& showerCand,
                             double vtxT)
{

    if (unusedViewClusters.empty()) return;
    
    double z_min = 10000;
    double z_max = -10000;
    double angle_min = 180;
    double angle_max = -180;

    SmartRefVector<Minerva::IDCluster> viewShowerCand;
    
    SmartRefVector<Minerva::IDCluster>::iterator it_view = unusedViewClusters.begin();
    for (SmartRefVector<Minerva::IDCluster>::iterator c = showerCand.begin();
         c != showerCand.end(); ++c){
        if ( (*c)->view() == (*it_view)->view() ) viewShowerCand.push_back(*c);
        if ( (*c)->view() == Minerva::IDCluster::U && (*c)->z() < z_min ) z_min = (*c)->z();
        if ( (*c)->view() == Minerva::IDCluster::U && (*c)->z() > z_max ) z_max = (*c)->z();
    }

    if (viewShowerCand.empty()) {
        std::cout << "\tAngleScan_U::completeView: no view cluster" << std::endl;
        return;
    }

    for ( it_view = viewShowerCand.begin(); it_view != viewShowerCand.end(); ++it_view){
        const double z   = (*it_view)->z() - fZ;
        const double u   = (*it_view)->position() - vtxT;
        const double ang = std::atan2(u,z)*TMath::RadToDeg();
        
        if ( ang >= angle_max ) angle_max = ang;
        if ( ang <= angle_min ) angle_min = ang;
        
    }

    z_min = z_min - 100;
    z_max = z_max + 100;
    angle_max = angle_max + 10.0;
    angle_min = angle_min - 10.0;
    
        /* Move clusters between (angle_min,angle_max) and (z_min,z_max) from
           'unusedClusters' to blob */
    std::cout << "\tAngleScan_U::completeView:oldsize: " << showerCand.size() << std::endl;
    coneView(unusedViewClusters, showerCand, vtxT, angle_min, angle_max, z_min, z_max );
    std::cout << "\tAngleScan_U::completeView:newsize: " << showerCand.size() << std::endl;
 
}

void AngleScan_U::coneView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                         SmartRefVector<Minerva::IDCluster>& showerCand,
                         double vtxT,
                         double min_angle, double max_angle,
                         double zmin, double zmax)
{
    SmartRefVector<Minerva::IDCluster> tmpClusters;
    tmpClusters.swap(unusedViewClusters);
    for (SmartRefVector<Minerva::IDCluster>::iterator c = tmpClusters.begin();
         c != tmpClusters.end(); ++c) {
        const double dZ = (*c)->z() - fZ;
        const double dU = (*c)->position() - vtxT;
        const double theta = std::atan2(dU,dZ)*TMath::RadToDeg();
        const double z_c = (*c)->z();
        
        if ((min_angle <= theta && theta <= max_angle) &&
            (zmin < z_c && z_c < zmax)) {
            showerCand.push_back(*c);
        } else {
            unusedViewClusters.push_back(*c);
        }
    }
}

const std::vector<TVector2>& AngleScan_U::GetPeaks() const {
    return fPeaks;
}

const std::vector<TVector2>& AngleScan_U::GetGoodPeaks() const {
    return fGoodPeaks;
}

const std::vector<double>& AngleScan_U::GetUShowerClosestDistances() const {
    return fUShowerClosestDistances;
}

const std::vector<double>& AngleScan_U::GetUShowerWeightedDistances() const {
    return fUShowerWeightedDistances;
}

unsigned int AngleScan_U::GetNuCandidate() const {
    return fUShowerCandidates.size();
}

unsigned int AngleScan_U::GetNCandidate() const {
    return fShowerCandidates.size();
}

const std::vector<SmartRefVector<Minerva::IDCluster> >&
AngleScan_U::GetUShowerCandVector() const {
    return fUShowerCandidates;
}

const std::vector<SmartRefVector<Minerva::IDCluster> >&
AngleScan_U::GetShowerCandVector() const {
    return fShowerCandidates;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetXClusters() const {
    return fXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetUClusters() const {
    return fUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetVClusters() const {
    return fVClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetUnusedClusters() const {
    return fRemainingClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetUnusedXClusters() const {
    return fRemainingXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetUnusedUClusters() const {
    return fRemainingUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_U::GetUnusedVClusters() const {
    return fRemainingVClusters;
}

std::vector<Minerva::IDBlob*> AngleScan_U::GetShowers() 
{
    std::vector<Minerva::IDBlob*> finalBlobs;
    std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s;
    
    for ( s = fShowerCandidates.begin(); s != fShowerCandidates.end(); ++s) {
        Minerva::IDBlob* newBlob = new Minerva::IDBlob;
        newBlob->add(*s);
        finalBlobs.push_back(newBlob);
    }
    
    return finalBlobs;
}

void AngleScan_U::SetUVMatchTolerance(double epsilon) {
    fUVMatchTolerance = epsilon;
}

void AngleScan_U::SetUVMatchMoreTolerance(double big_epsilon) {
    fUVMatchMoreTolerance = big_epsilon;
}

void AngleScan_U::AllowUVMatchWithMoreTolerance(bool b) {
    fAllowUVMatchWithMoreTolerance = b;
}

#endif

