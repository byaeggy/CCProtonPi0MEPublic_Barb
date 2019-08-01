#ifndef AngleScan_V_cpp
#define AngleScan_V_cpp

#undef NDEBUG
#include <cassert>
#include <cmath>

#include <TMath.h>
#include <TH1F.h>

#include <Event/IDCluster.h>
#include <Event/IDBlob.h>
#include <Event/VectorTypeDefs.h>
#include <Event/Vertex.h>

#include "AngleScan_V.h"
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


AngleScan_V::AngleScan_V(const SmartRefVector<Minerva::IDCluster>& clusters,
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

void AngleScan_V::Initialize()
{

    std::cout << "    AngleScan_V::Initialize " << std::endl;

    ClusterVectorInfo clusterVectorInfo(fAllClusters,true,true,true);
    fXClusters = clusterVectorInfo.GetXClusters();
    fUClusters = clusterVectorInfo.GetUClusters();
    fVClusters = clusterVectorInfo.GetVClusters();

    std::cout << "      AngleScan_V:: x size = " << fXClusters.size() << std::endl;
    std::cout << "      AngleScan_V:: u size = " << fUClusters.size() << std::endl;
    std::cout << "      AngleScan_V:: v size = " << fVClusters.size() << std::endl;

        /* Copy to the working containers */
    fRemainingXClusters = fXClusters;
    fRemainingUClusters = fUClusters;
    fRemainingVClusters = fVClusters;

    assert(fRemainingXClusters.size() == fXClusters.size());
    assert(fRemainingUClusters.size() == fUClusters.size());
    assert(fRemainingVClusters.size() == fVClusters.size());
}

void AngleScan_V::BuildThetaHistogram()
{
    const unsigned int N = 90;
    const double lower = -180.0;
    const double upper = +180.0;
    fTheta = new TH1F("theta","Theta dist", N, lower, upper);
    
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = fRemainingVClusters.begin();
         c != fRemainingVClusters.end(); ++c) {
        const double dZ = (*c)->z() - fZ;
        const double dV = (*c)->position() - fV;
        const double theta = std::atan2(dV,dZ)*TMath::RadToDeg();
        const double w = (*c)->pe();
        
        fTheta->Fill(theta,w);
    }
}

void AngleScan_V::FindPeaks()
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

void AngleScan_V::FormVShowerCand() 
{
    for (std::vector<TVector2>::const_iterator peak = fPeaks.begin();
         peak != fPeaks.end(); ++peak) {
        const double lower_edge = peak->X();
        const double upper_edge = peak->Y();
        const double zmin       = 4500.0;
        const double zmax       = 10000.0;

        SmartRefVector<Minerva::IDCluster> showerCand;
        coneView(fRemainingVClusters, showerCand, fV, lower_edge, upper_edge, zmin, zmax);
        
        if (showerCand.empty()) continue;

        if (showerCand.size() == 1 && showerCand.front()->pe() > 30) {

            fVShowerCandidates.push_back(showerCand);
            fGoodPeaks.push_back(*peak);
            
        } else if (showerCand.size() > 1) {

            fVShowerCandidates.push_back(showerCand);
            fGoodPeaks.push_back(*peak);
        }
        
    }

        // Calculate the distance from the shower candidates to the event vertex
        // Try two definitions of distance: closest and energy weighted
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
             = fVShowerCandidates.begin();
         s != fVShowerCandidates.end(); ++s) {

        SmartRefVector<Minerva::IDCluster>& vshowerCand = *s;

        double d_min = 1.e6;
        double d_weighted = 0.0;
        double total_energy = 0.0;
        for (SmartRefVector<Minerva::IDCluster>::iterator c = vshowerCand.begin();
             c != vshowerCand.end(); ++c) {

            double v = (*c)->position();
            double z = (*c)->z();
            double d = std::sqrt(std::pow(v-fV,2) + std::pow(z-fZ,2));

            if (d < d_min) {
                d_min = d;
            }

            d_weighted   += d * (*c)->energy();
            total_energy += (*c)->energy();
        }

        fVShowerClosestDistances.push_back(d_min);
        fVShowerWeightedDistances.push_back(d_weighted/total_energy);
    }
    
    std::sort(fVShowerCandidates.begin(),fVShowerCandidates.end(), greaterShower());
    
}

void AngleScan_V::FormXUVShowerCand() 
{
    std::vector<SmartRefVector<Minerva::IDCluster> > nogrowShowerCandidates;
    std::cout << "AngleScan_V::total V candidates: " << fVShowerCandidates.size() << std::endl;
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
             = fVShowerCandidates.begin();
         s != fVShowerCandidates.end(); ++s) {
        std::cout << "\tAngleScan_V::V cand: " << s->size() << std::endl;
        SmartRefVector<Minerva::IDCluster>& vshowerCand = *s;
        
        SmartRefVector<Minerva::IDCluster> showerCand = vshowerCand;
        double zmin = +1e6;
        double zmax = -1e6;
        double ztot = 0.0;
        for (SmartRefVector<Minerva::IDCluster>::iterator c = vshowerCand.begin();
             c != vshowerCand.end(); ++c) {
            zmin  = std::min(zmin,(*c)->z());
            zmax  = std::max(zmax,(*c)->z());
            ztot += (*c)->z();
        }

        if (zmax-zmin < 50.0) {
            const double zcenter = ztot/vshowerCand.size();
            std::cout << "\t AngleScan_V::Candidate with small z extend at: " << zcenter << std::endl;
            SmartRefVector<Minerva::IDCluster> xclusters_tmp;
            SmartRefVector<Minerva::IDCluster> uclusters_tmp;
            fRemainingXClusters.swap(xclusters_tmp);
            fRemainingUClusters.swap(uclusters_tmp);

            for (SmartRefVector<Minerva::IDCluster>::iterator c = xclusters_tmp.begin();
                 c != xclusters_tmp.end(); ++c) {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan_V::adding a X cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingXClusters.push_back(*c);
            }
                        
            for (SmartRefVector<Minerva::IDCluster>::iterator c = uclusters_tmp.begin();
                 c != uclusters_tmp.end(); ++c) {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan_V::adding a U cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingUClusters.push_back(*c);
            }
            
            fShowerCandidates.push_back(showerCand);

        } else {
            std::cout << "\t AngleScan_V::Candidate with large z extent " << std::endl;
            addClustersToBlob(vshowerCand,fRemainingXClusters,fRemainingUClusters,showerCand,
                              fUVMatchTolerance);
            
            if (showerCand.size() > vshowerCand.size()) fShowerCandidates.push_back(showerCand);
            else if (vshowerCand.size() >= 3) nogrowShowerCandidates.push_back(vshowerCand);
            else {
                std::cout << "\t AngleScan_V::Throw away V shower candidate" << std::endl;
            }
            
        }
    }

    std::cout << "AngleScan_V::Shower candidates: " << fShowerCandidates.size() << std::endl;
    std::cout << "AngleScan_V::Nogrow candidates: " << nogrowShowerCandidates.size() << std::endl;

    if (fAllowUVMatchWithMoreTolerance) {
        for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s
                 = nogrowShowerCandidates.begin();
             s != nogrowShowerCandidates.end(); ++s) {
            SmartRefVector<Minerva::IDCluster>& vshowerCand = *s; 
            SmartRefVector<Minerva::IDCluster>  showerCand = vshowerCand;
            
            addClustersToBlob(vshowerCand,fRemainingXClusters,fRemainingUClusters,showerCand,
                              fUVMatchMoreTolerance);
            
            fShowerCandidates.push_back(showerCand);
        }
    }
    
    std::cout << "AngleScan_V::Shower candidates(10x): " << fShowerCandidates.size() << std::endl;
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s = fShowerCandidates.begin();
         s != fShowerCandidates.end(); ++s) {
        std::cout << "\tAngleScan_V:: cand: " << std::distance(fShowerCandidates.begin(),s) << " "
                  << fRemainingXClusters.size() << " " << fRemainingUClusters.size()
                  << std::endl;
        completeView(fRemainingXClusters, *s, fX);
        completeView(fRemainingUClusters, *s, fU);
    }

    std::cout << "AngleScan_V::Final showers: " << fShowerCandidates.size() << std::endl;
    std::sort(fShowerCandidates.begin(), fShowerCandidates.end(), greaterShower());
    
    std::copy(fRemainingVClusters.begin(),fRemainingVClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingXClusters.begin(),fRemainingXClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingUClusters.begin(),fRemainingUClusters.end(), std::back_inserter(fRemainingClusters));

}

void AngleScan_V::DoReco()
{
    std::cout << "AngleScan_V::DoReco() " << std::endl;
    BuildThetaHistogram();
    FindPeaks();
    FormVShowerCand();
    FormXUVShowerCand();
}

int AngleScan_V::GetLimitBin(const TH1F *hMax, int n_bin )const
{

    int max=n_bin, x = 1, count = 0;
    while ( count < 1 ){
        if ( hMax->GetBinContent(n_bin + x) > 0 ) max = n_bin + x;
        else count++;
        x++;
    }
    
    return max;

}
void AngleScan_V::addClustersToBlob(SmartRefVector<Minerva::IDCluster>& xshowerCand,
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

void AngleScan_V::completeView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
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
        std::cout << "\tAngleScan_V::completeView: no view cluster" << std::endl;
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
    std::cout << "\tAngleScan_V::completeView:oldsize: " << showerCand.size() << std::endl;
    coneView(unusedViewClusters, showerCand, vtxT, angle_min, angle_max, z_min, z_max );
    std::cout << "\tAngleScan_V::completeView:newsize: " << showerCand.size() << std::endl;
 
}

void AngleScan_V::coneView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
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
        const double dV = (*c)->position() - vtxT;
        const double theta = std::atan2(dV,dZ)*TMath::RadToDeg();
        const double z_c = (*c)->z();
        
        if ((min_angle <= theta && theta <= max_angle) &&
            (zmin < z_c && z_c < zmax)) {
            showerCand.push_back(*c);
        } else {
            unusedViewClusters.push_back(*c);
        }
    }
}

const std::vector<TVector2>& AngleScan_V::GetPeaks() const {
    return fPeaks;
}

const std::vector<TVector2>& AngleScan_V::GetGoodPeaks() const {
    return fGoodPeaks;
}

const std::vector<double>& AngleScan_V::GetVShowerClosestDistances() const {
    return fVShowerClosestDistances;
}

const std::vector<double>& AngleScan_V::GetVShowerWeightedDistances() const {
    return fVShowerWeightedDistances;
}

unsigned int AngleScan_V::GetNuCandidate() const {
    return fVShowerCandidates.size();
}

unsigned int AngleScan_V::GetNCandidate() const {
    return fShowerCandidates.size();
}

const std::vector<SmartRefVector<Minerva::IDCluster> >&
AngleScan_V::GetVShowerCandVector() const {
    return fVShowerCandidates;
}

const std::vector<SmartRefVector<Minerva::IDCluster> >&
AngleScan_V::GetShowerCandVector() const {
    return fShowerCandidates;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetXClusters() const {
    return fXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetUClusters() const {
    return fUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetVClusters() const {
    return fVClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetUnusedClusters() const {
    return fRemainingClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetUnusedXClusters() const {
    return fRemainingXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetUnusedUClusters() const {
    return fRemainingUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan_V::GetUnusedVClusters() const {
    return fRemainingVClusters;
}

std::vector<Minerva::IDBlob*> AngleScan_V::GetShowers() 
{
    std::vector<Minerva::IDBlob*> finalBlobs;
    std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s;

    for ( s= fShowerCandidates.begin(); s != fShowerCandidates.end(); ++s) {
        Minerva::IDBlob* newBlob = new Minerva::IDBlob;
        newBlob->add(*s);
        finalBlobs.push_back(newBlob);
    }

    return finalBlobs;
}

void AngleScan_V::SetUVMatchTolerance(double epsilon) {
    fUVMatchTolerance = epsilon;
}

void AngleScan_V::SetUVMatchMoreTolerance(double big_epsilon) {
    fUVMatchMoreTolerance = big_epsilon;
}

void AngleScan_V::AllowUVMatchWithMoreTolerance(bool b) {
    fAllowUVMatchWithMoreTolerance = b;
}


#endif

