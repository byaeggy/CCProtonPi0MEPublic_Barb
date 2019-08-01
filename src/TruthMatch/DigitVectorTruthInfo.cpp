#undef NDEBUG
#include <cassert>
#include <set>
#include <iterator>
#include <iostream>


#include <Event/IDCluster.h>
#include <Event/IDDigit.h>
#include <Event/MCIDDigit.h>
#include <Event/MCHit.h>
#include <Event/TG4Trajectory.h>


#include "DigitVectorTruthInfo.h"
#include "TraverseHistory.h"

DigitVectorTruthInfo::DigitVectorTruthInfo()
    : fTotalNormEnergy(0.0),
      fDataNormEnergy(0.0),
      fMCNormEnergy(0.0),
      fMCXtalkNormEnergy(0.0),
      fTotalTruthEnergy(0.0),
      fSharedTruthEnergy(0.0)
{}

double DigitVectorTruthInfo::GetTotalNormEnergy() const {
    return fTotalNormEnergy;
}
double DigitVectorTruthInfo::GetDataNormEnergy() const {
    return fDataNormEnergy;
}
double DigitVectorTruthInfo::GetMCNormEnergy() const {
    return fMCNormEnergy;
}
double DigitVectorTruthInfo::GetMCXtalkNormEnergy() const {
    return fMCXtalkNormEnergy;
}
double DigitVectorTruthInfo::GetTotalTruthEnergy() const {
    return fTotalTruthEnergy;
}
double DigitVectorTruthInfo::GetSharedTruthEnergy() const {
    return fSharedTruthEnergy;
}

double DigitVectorTruthInfo::GetDataFraction() const {
    return fDataNormEnergy/fTotalNormEnergy;
}
double DigitVectorTruthInfo::GetXtalkFraction() const {
    return fMCXtalkNormEnergy/fMCNormEnergy;
}
double DigitVectorTruthInfo::GetSharedFraction() const {
    return fSharedTruthEnergy/fTotalTruthEnergy;
}

double DigitVectorTruthInfo::GetEdepByPdg(int pdg) const {
    std::map<int,double>::const_iterator target = fPdgEdepMap.find(pdg);
    if (target != fPdgEdepMap.end()) return target->second;

    return 0.0; /* If the particle does not exist, return 0.0 */
}

double DigitVectorTruthInfo::GetEdepByTrackId(int trackId) const {
    std::map<int,double>::const_iterator target = fNextToPrimaryIdEdepMap.find(trackId);
    if (target != fNextToPrimaryIdEdepMap.end()) return target->second;

    return 0.0; /* If trackId does not exist */
}

/* Return the pdg of the primary particle contributing most to this digit vector */
int DigitVectorTruthInfo::GetMostEvisPdg() const {

    if (fPdgEdepMap.empty()) return -1;
    
    double max = -1.e6;
    int pdg = -1;
    for (std::map<int,double>::const_iterator p = fPdgEdepMap.begin(); p != fPdgEdepMap.end(); ++p) {
        if (p->second > max) {
            max = p->second;
            pdg = p->first;
        }
    }

    assert(pdg != -1);

    return pdg;
}

int DigitVectorTruthInfo::GetMostEvisId() const{

  if (fIdEdepMap.empty()) return -1;

  double max = -1.e6;
  int id = -1;
  for (std::map<int,double>::const_iterator p = fIdEdepMap.begin(); p != fIdEdepMap.end(); ++p) {
    if (p->second > max) {
      max = p->second;
      id = p->first;
    }
  }

    assert(id != -1);

    return id;
  
}


int DigitVectorTruthInfo::GetParentEvisPdg() const {

    if (fParentPdgEdepMap.empty()) return -1;
    
    double max = -1.e6;
    int pdg = -1;
    for (std::map<int,double>::const_iterator p = fParentPdgEdepMap.begin(); p != fParentPdgEdepMap.end(); ++p) {
        if (p->second > max) {
            max = p->second;
            pdg = p->first;
        }
    }

    assert(pdg != -1);

    return pdg;
}

int DigitVectorTruthInfo::GetParentEvisId() const{

  if (fParentIdEdepMap.empty()) return -1;

  double max = -1.e6;
  int id = -1;
  for (std::map<int,double>::const_iterator p = fParentIdEdepMap.begin(); p != fParentIdEdepMap.end(); ++p) {
    if (p->second > max) {
      max = p->second;
      id = p->first;
    }
  }

    assert(id != -1);

    return id;
  
}


const std::map<int,double>& DigitVectorTruthInfo::GetPdgEdepMap() const {
    return fPdgEdepMap;
}

const std::map<int,double>& DigitVectorTruthInfo::GetNextToPrimaryIdEdepMap() const {
    return fNextToPrimaryIdEdepMap;
}

// Modified version to work with clusters instead of digits
void DigitVectorTruthInfo::ParseTruth(const SmartRefVector<Minerva::IDCluster>& clusters, 
				      const std::map<int, Minerva::TG4Trajectory*>& trajectories, int goalcode)
{

    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin(); c != clusters.end(); ++c) {
        SmartRefVector<Minerva::IDDigit> digits = (*c)->digits();
        
    
        for (SmartRefVector<Minerva::IDDigit>::const_iterator d = digits.begin(); d != digits.end(); ++d) {
            const double norm_energy = (*d)->normEnergy();                    // norm energy of this digit
            const double fraction    = (*c)->getDigitOwnershipFraction(*d);   // fraction of the digit in the cluster
                                                                              // use correct fraction to avoid double counting
            const double digit_energy = norm_energy * fraction;
            fTotalNormEnergy += digit_energy;
            
            const Minerva::IDDigit* digit = *d;
            const Minerva::MCIDDigit* mcdigit = dynamic_cast<const Minerva::MCIDDigit*>(digit);
            if (!mcdigit) {
                fDataNormEnergy += digit_energy;
                continue;
            }
            
                /* This is an MC digit, including MC x-talk digit */
            fMCNormEnergy += digit_energy;
            
            const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
            
            if (hits.empty()) {
                fMCXtalkNormEnergy += digit_energy;
                continue;
            }
            
                /* Handle non x-talk, MC digits */
                /* Implementation note: This part could be replaced by HitVectorTruthInfo */
            double total = 0.0;
            std::set<int> primaryContributors;
            for (SmartRefVector<Minerva::MCHit>::const_iterator h = hits.begin(); h != hits.end(); ++h) {
                double hit_energy = (*h)->energy();
                hit_energy *= fraction;                // all MC hits inside the digit must share the same
                                                       // fraction
                total += hit_energy;
                
                TraverseHistory history((*h)->GetTrackId());
		if (goalcode == 0){
		  history.DoTraverse(trajectories);		  
		}else{
		  history.DoTraverse(trajectories);
		}
                int primaryPdg = history.GetPrimaryPdg();
		int primaryId  = history.GetPrimaryId();
		int parentPdg  = history.GetParentPdg();
		int parentId   = history.GetParentId();
                int nextToPrimaryId = history.GetNextToPrimaryId();
                
		if (goalcode != 0){
		  if (history.IsPresent(goalcode)){
		    primaryPdg = goalcode;
		  }

		  fPdgEdepMap[primaryPdg] += hit_energy;
		  fIdEdepMap[primaryId] += hit_energy;
		  fParentPdgEdepMap[parentPdg] += hit_energy;
		  fParentIdEdepMap[parentId] += hit_energy;
		}
		/*
		if (goalcode != 0  && false){
		  if (history.IsPresent(goalcode)){
		    
		    //for (std::map<int, double>::const_iterator p = fPdgEdepMap.begin(); p != fPdgEdepMap.end(); ++p){
		      
		    //}
		  }
		  }*/
 
		//                fPdgEdepMap[primaryPdg] += hit_energy;
                fNextToPrimaryIdEdepMap[nextToPrimaryId] += hit_energy;
                
                primaryContributors.insert(primaryPdg);
            
            }
            
            fTotalTruthEnergy += total;
        
            if (primaryContributors.size() > 1) fSharedTruthEnergy += total;
            
        } // loop over cluster's digits
 
    }
   
}

void DigitVectorTruthInfo::PrintMap() const
{
    std::cout << "Primary particle PDG | energy deposit " << std::endl;
    for (std::map<int,double>::const_iterator p = fPdgEdepMap.begin();
         p != fPdgEdepMap.end(); ++p) {
        std::cout << "\t--> " << p->first << " " << p->second << std::endl;
    }
}
