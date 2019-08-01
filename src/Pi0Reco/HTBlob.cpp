#include "HTBlob.h"
#include <cmath>
#include <limits>

#include "Event/IDCluster.h"
#include "Event/IDBlob.h"

#include "MinervaUtils/IMinervaMathTool.h"
#include "MinervaUtils/MinervaObjectSort.h"

#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"

using namespace MinervaSort;

DECLARE_TOOL_FACTORY( HTBlob );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
HTBlob::HTBlob( const std::string& type, const std::string& name, const IInterface* parent ) : 
    MinervaHistoTool( type, name, parent ) 
{
    debug() << "Instantiating HTBlob..." << endmsg;
    declareInterface<IHoughBlob>(this);

    //dEdx
    declareProperty( "NumberPlanesdEdX",              m_planesdEdx = 4 );

    // Blob vertex and Event vertex
    declareProperty( "MinDistanceStripPhoton",        m_minDistanceStripPhoton = 75* CLHEP::mm );
    declareProperty( "MaxDisntaceModulePhoton",       m_minDistanceModulePhoton = 50* CLHEP::mm );

    // Calibration Constants - Latest Update 11/05/2015
    m_kT = 1.326; // This kT is only for kS calculation - Actual kT is evis dependent
    m_kE = 3.2367; 
    m_kS_X = (4*m_kE - 2*m_kT)/2;
    m_kS_UV = (4*m_kE - 3*m_kT);
    m_kH = 12.650;
}

//=============================================================================
/// Destructor
//=============================================================================
HTBlob::~HTBlob() {}


//=============================================================================
// Initialize
//=============================================================================
StatusCode HTBlob::initialize() 
{

    debug() << "Initializing HTBlob..." << endmsg;
    StatusCode sc = this->MinervaHistoTool::initialize();
    if( sc.isFailure() ) { return Error( "Failed to initialize!", sc ); }

    try {
        m_mathTool = tool<IMinervaMathTool>("MinervaMathTool");
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: MinervaMathTool!" << endmsg;
        return StatusCode::FAILURE;
    }

    return sc;
}

//=============================================================================
// GetClusters
//=============================================================================
StatusCode HTBlob::GetViewClustersAboveThreshold( SmartRefVector<Minerva::IDCluster> &idClusterVec,
        SmartRefVector<Minerva::IDCluster> &idClusterView,
        Minerva::IDCluster::View view, double pecut ) const
{

    if( !idClusterVec.size() ) { 
        debug() << " ALERT: Input Vector of Clusters to the Get_Clusters tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }

    SmartRefVector<Minerva::IDCluster> ClusTemp = idClusterVec;
    SmartRefVector<Minerva::IDCluster>::iterator itClus = ClusTemp.begin();
    idClusterVec.clear();

    for ( ; itClus != ClusTemp.end(); itClus++ ){

        if ( (*itClus)->view()==view && (*itClus)->pe()/(*itClus)->iddigs()>3 && (*itClus)->pe()>pecut  ) {
            idClusterView.push_back(*itClus);
        }
        else idClusterVec.push_back(*itClus);

    }

    return StatusCode::SUCCESS;

}



//======================================================================
//  XUVMatch
//=======================================================================
StatusCode HTBlob::XUVMatch(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorU,
        SmartRefVector<Minerva::IDCluster> &ClusVectorV, double match) const
{
    debug() << " == HTBlob::XUVMatch " << endmsg;

    SmartRefVector<Minerva::IDCluster> ClusTemp = Seed;
    SmartRefVector<Minerva::IDCluster>::iterator itClusX;
    double zmin = 10000, zmax = 0;

    debug() << " Seed size = " << Seed.size()  << endmsg;

    for ( itClusX = ClusTemp.begin(); itClusX != ClusTemp.end(); itClusX++ ){
        double z  = (*itClusX)->z();
        if ( z < zmin ) zmin = z;
        if ( z > zmax ) zmax = z;
    }
    debug() << " Finding clusters with macth and zmin " << zmin << "; zmax " << zmax << endmsg;

    for ( itClusX = ClusTemp.begin(); itClusX != ClusTemp.end(); itClusX++ ) {
        if ( (*itClusX)->view() != Minerva::IDCluster::X ) continue;
        XUVMatch( *itClusX, Seed, ClusVectorU, ClusVectorV, zmin, zmax, match );
    }

    debug() << " == HTBlob::XUVMatch - NO MORE SEED CLUSTERS to Match, leaving seed with size: " << Seed.size()
        << endmsg << endmsg;

    return StatusCode::SUCCESS;

}

//======================================================================
//  XUVMatch overload
//=======================================================================
StatusCode HTBlob::XUVMatch( SmartRef<Minerva::IDCluster> Cluster, SmartRefVector<Minerva::IDCluster> &Seed,
        SmartRefVector<Minerva::IDCluster> &ClusVectorU,
        SmartRefVector<Minerva::IDCluster> &ClusVectorV,
        double zmin, double zmax, double match) const
{
    debug() << " HTBlob::XUVMatch  - Overload with match = " << match <<endmsg;

    SmartRefVector<Minerva::IDCluster>::iterator itClusU, itClusV;
    SmartRef<Minerva::IDCluster> U, V;
    double dmin = 1000, distance;

    debug() << " MATCH, X cluster, pe " << Cluster->pe() << "; z " << Cluster->z() << "; position " << Cluster->position()
        << "; sum pos " << Cluster->position()+Cluster->tpos1()+Cluster->tpos2() << endmsg;

    for ( itClusU = ClusVectorU.begin(); itClusU != ClusVectorU.end(); itClusU++ ){

        if ( fabs( Cluster->z() - (*itClusU)->z() ) > 50 ) continue;

        for ( itClusV = ClusVectorV.begin(); itClusV != ClusVectorV.end(); itClusV++ ) {

            if ( fabs( Cluster->z() - (*itClusV)->z() ) > 50 ) continue;

            distance =  Cluster->position()+Cluster->tpos1()+Cluster->tpos2();
            distance -= ((*itClusU)->position()+(*itClusU)->tpos1()+(*itClusU)->tpos2());
            distance -= ((*itClusV)->position()+(*itClusV)->tpos1()+(*itClusV)->tpos2());
            distance = fabs(distance);

            if ( distance < dmin) {
                debug() << " MATCH Cand. " << " U, pe " << (*itClusU)->pe() << "; z " << (*itClusU)->z()
                    << "; position " << (*itClusU)->position()
                    << "; sum pos " << (*itClusU)->position()+(*itClusU)->tpos1()+(*itClusU)->tpos2() << endmsg;

                debug() << " V, pe " << (*itClusV)->pe() << "; z " << (*itClusV)->z()
                    << "; position " << (*itClusV)->position()
                    << "; sum pos " << (*itClusV)->position()+(*itClusV)->tpos1()+(*itClusV)->tpos2()
                    << endmsg;

                dmin = distance;
                U = *itClusU;
                V = *itClusV;
            }

        }
    }

    double efmatch = fabs(zmax-zmin) < 190 ? match*2 : match; // High angles  or shorts?

    if ( dmin <= efmatch && Cluster->z() >= zmin && Cluster->z() <= zmax ) {
        debug() << " FOUND MATCH " << U << " " << V << " Match " << efmatch << endmsg;
        SmartRefVector<Minerva::IDCluster>::iterator itU, itV;

        itU = remove(ClusVectorU.begin(),ClusVectorU.end(),U); // move elements to the end to can erase
        itV = remove(ClusVectorV.begin(),ClusVectorV.end(),V);
        ClusVectorU.erase(itU,ClusVectorU.end()); // found it in doxygen
        ClusVectorV.erase(itV,ClusVectorV.end()); // found it in doxygen
        Seed.push_back(U);
        Seed.push_back(V);
    }

    return StatusCode::SUCCESS;

}

//======================================================================
//  AddClustersInsideCone
//=======================================================================
StatusCode HTBlob::AddClusterInsideCone(SmartRef<Minerva::IDCluster> UnuCluster, std::vector<Minerva::IDBlob*> &idBlobs, 
        Gaudi::XYZPoint vert ) const
{
    debug() << " HTBlob::AddClusterInsideCone " << endmsg;
    debug() << " Cluster view = " << UnuCluster->view() << "; z = " << UnuCluster->z()
        << "; position = " << UnuCluster->position() << "; pe = " << UnuCluster->pe() << endmsg;

    std::vector<Minerva::IDBlob*> BlobsTemp = idBlobs; idBlobs.clear();
    std::vector<Minerva::IDBlob*>::iterator itBlob;
    double amin = 1000, angle;
    int count = -1, marker = 0;

    for (itBlob = BlobsTemp.begin(); itBlob != BlobsTemp.end(); itBlob++ ){

        count++; 
        if ( (*itBlob)->direction().z() > 0 && ((*itBlob)->startPoint().z() - 25) > UnuCluster->z() ) continue;
        if ( (*itBlob)->direction().z() < 0 && ((*itBlob)->startPoint().z() + 25) < UnuCluster->z() ) continue;
        if ( !Angle( UnuCluster, (*itBlob)->direction(), vert, angle ) ) continue;
        debug() << " Blob " << count << "; angle = " << angle << "; amin = " << amin << endmsg;
        if ( angle < amin ) {
            amin = angle;
            marker = count;
        }

    }

    debug() << " amin = " << amin << " marker = " << marker  << endmsg;

    if ( amin < 0.174) { //10 degrees

        itBlob = BlobsTemp.begin();
        itBlob = itBlob + marker;
        (*itBlob)->add(UnuCluster);

    }

    idBlobs = BlobsTemp;

    return StatusCode::SUCCESS;
}

//======================================================================
//  PseudoCone
//=======================================================================
StatusCode HTBlob::PseudoCone(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorX,
        Gaudi::XYZVector direction, Gaudi::XYZPoint vert ) const
{

    debug() << " HTBlob::PseudoCone, clusters with Angles < 0.06 will be include in the seed " << endmsg;

    SmartRefVector<Minerva::IDCluster> ClusTemp = ClusVectorX; ClusVectorX.clear();
    SmartRefVector<Minerva::IDCluster>::iterator itClusX;

    double angle;

    for ( itClusX = ClusTemp.begin(); itClusX != ClusTemp.end(); itClusX++ ){
        if ( Angle( *itClusX, direction, vert, angle ) ) { 
            if ( angle < 0.06 && (*itClusX)->z() > vert.z() ) Seed.push_back(*itClusX); // must be carefull with backward showers
            else ClusVectorX.push_back(*itClusX);
        } else ClusVectorX.push_back(*itClusX);
    }

    debug() << endmsg;

    return StatusCode::SUCCESS;

}

//======================================================================
//  Angle
//=======================================================================
StatusCode HTBlob::Angle( SmartRef<Minerva::IDCluster> Cluster, Gaudi::XYZVector direction, Gaudi::XYZPoint vert, double &angle ) const
{

    if ( direction.x() == -9999 || vert.x() == -9999 ) return StatusCode::FAILURE;
    double dx, dz, Dx, Dz;

    switch(Cluster->view())
    {
        case Minerva::IDCluster::X:
            dx = Cluster->position() - vert.x();
            Dx = direction.x();
            break;

        case Minerva::IDCluster::U:
            dx = Cluster->position() - m_mathTool->calcUfromXY(vert.x(),vert.y());
            Dx = m_mathTool->calcUfromXY(direction.x(),direction.y());
            break;

        case Minerva::IDCluster::V:
            dx = Cluster->position() - m_mathTool->calcVfromXY(vert.x(),vert.y());
            Dx = m_mathTool->calcVfromXY(direction.x(),direction.y());
            break;

        default:
            throw MinervaException("Unknown cluster view");

    }

    dz = Cluster->z() - vert.z();
    Dz = direction.z();

    double moduled = sqrt( pow(dx,2) + pow(dz,2) );
    double moduleD = sqrt( pow(Dx,2) + pow(Dz,2) );

    dx = dx/moduled; dz = dz/moduled;
    Dx = Dx/moduleD; Dz = Dz/moduleD;

    angle = acos(fabs(dx*Dx+dz*Dz));

    debug() << " pe = " << Cluster->pe() << "; z = " << Cluster->z() << "; pos = "  << Cluster->position()
        << " Angle " << angle << endmsg;

    return StatusCode::SUCCESS;

}

//=============================================================================
// Create2dHTBlob
//=============================================================================
StatusCode HTBlob::Create2dHTSeed( SmartRefVector<Minerva::IDCluster> &idClusterView,
        SmartRefVector<Minerva::IDCluster> &HT2dClusters, 
        double r, double theta, Gaudi::XYZPoint ref, double &spX, double &spZ ) const
{

    debug() << " HTtool::Create2dHTSeed " << endmsg;

    double rmin, rmax, x, z, zmin = 10000, Total_e = 0;
    SmartRefVector<Minerva::IDCluster> ClusTemp = idClusterView;
    SmartRefVector<Minerva::IDCluster>::iterator itClus = ClusTemp.begin();
    idClusterView.clear();

    debug() << " Will study " << ClusTemp.size() << " clusters " << endmsg;

    debug() << " Seed with, r: " << r << ", theta = " << theta << ";contains these clusters: " << endmsg;

    for ( ; itClus != ClusTemp.end(); itClus++ ){
        z = (*itClus)->z() - ref.z();
        x = (*itClus)->tpos1() - ref.x();
        rmin = x*sin(theta*CLHEP::pi/180) + z*cos(theta*CLHEP::pi/180);
        x = (*itClus)->tpos2() - ref.x();
        rmax = x*sin(theta*CLHEP::pi/180) + z*cos(theta*CLHEP::pi/180);

        if ( fabs ( 2*r - rmin - rmax ) <= 90 ){
            if ( (*itClus)->z()< zmin ) {
                zmin = (*itClus)->z();
                spZ  = (*itClus)->z();
                spX  = (*itClus)->position();
            }
            debug() << " pe = " << (*itClus)->pe() << "; z = " << (*itClus)->z() << "; pos = "  << (*itClus)->position()
                << endmsg;
            HT2dClusters.push_back(*itClus);
            Total_e += (*itClus)->energy();

            continue;
        }

        idClusterView.push_back(*itClus);

    }

    debug() << " Total energy comming from seed = " << Total_e << "; this energy must be bigger than 19"
        << endmsg << endmsg;

    if ( Total_e < 19 ) {
        idClusterView.insert(idClusterView.end(),HT2dClusters.begin(),HT2dClusters.end());
        return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;

}
//=============================================================================
// GetDirection - need the startPoint
//=============================================================================
StatusCode HTBlob::GetDirection( Minerva::IDBlob *idBlob ) const
{
    debug() << " HTtool::GetDirection " << endmsg;

    Gaudi::XYZPoint vertex = idBlob->startPoint();
    GetDirection( idBlob, vertex );
    return StatusCode::SUCCESS;

}

//=============================================================================
//  GetDirection - Cesar Sotelo's idea 
//  Calculate direction using every cluster weighted by Energy 
//  and distance inverse weighted
//=============================================================================
bool HTBlob::GetDirection( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vert ) const
{

    debug() << " HTBlob::GetDirection " << endmsg;

    SmartRefVector<Minerva::IDCluster> idClusters = idBlob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusters.begin();
    SmartRefVector<Minerva::IDCluster>::iterator endClus = idClusters.end();

    std::sort(itClus, endClus, ClusterSorterTime());

    double Xx = 0.0;
    double Zx = 0.0;
    double Xu = 0.0;
    double Xv = 0.0;
    double totalX = 0.0;
    double totalU = 0.0;
    double totalV = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double distance = 0.0;

    for ( ; itClus != idClusters.end(); itClus++ ){
      dz = (*itClus)->z() - vert.z();

      switch((*itClus)->view())
        {
	case Minerva::IDCluster::X:
	  dx = (*itClus)->position() - vert.x();
	  distance = sqrt(pow(dx,2)+pow(dz,2));
	  if ( distance == .0 ) distance = 25;
	  Xx += dx*(*itClus)->energy()/distance;
	  Zx += dz*(*itClus)->energy()/distance;
	  totalX += (*itClus)->energy()/distance;
	  break;
	  
	case Minerva::IDCluster::U:
	  dx = (*itClus)->position() - m_mathTool->calcUfromXY(vert.x(), vert.y());
	  distance = sqrt(pow(dx,2)+pow(dz,2));
	  if ( distance == .0 ) distance = 25;
	  Xu += dx*(*itClus)->energy()/distance;
	  totalU += (*itClus)->energy()/distance;
	  break;
	  
	case Minerva::IDCluster::V:
	  dx = (*itClus)->position() - m_mathTool->calcVfromXY(vert.x(), vert.y());
	  distance = sqrt(pow(dx,2)+pow(dz,2));
	  if ( distance == .0 ) distance = 25;
	  Xv += dx*(*itClus)->energy()/distance;
	  totalV += (*itClus)->energy()/distance;
	  break;
		
	default:
	  throw MinervaException("Unknown cluster view");
        }
      
    }

    Gaudi::XYZVector direction;

    dx = Xx/totalX;
    dz = Zx/totalX;
    dy = -9999;
    bool valid_dY = false;
    const double epsilon = std::numeric_limits<double>::epsilon();

    if ( std::abs(Xu) > epsilon && std::abs(Xv) > epsilon ) {
      dy = m_mathTool->calcYfromUV(Xu/totalU,Xv/totalV);
      valid_dY = true;
    }
    else if ( std::abs(Xu) > epsilon ){
      dy = (dx*.5 - Xu/totalU)*2/sqrt(3); 	//calcYfromXU?
      valid_dY = true;
    }
    else if ( std::abs(Xv) > epsilon ){
      dy = (Xv/totalV - dx*.5)*2/sqrt(3); 	//calcYfromXV?
      valid_dY = true;
    }

    if  ( !valid_dY ){
      debug() << " Bad direction" << endmsg;
      idBlob->setDirection(Gaudi::XYZVector(-9999,-9999,-9999));
      return false;
    }
    else {

      double mod = sqrt( pow(dx,2)+pow(dy,2)+pow(dz,2) );

      direction.SetX(dx/mod);
      direction.SetY(dy/mod);
      direction.SetZ(dz/mod);
      
    }

    debug() << " Setting direction " << direction << " Blob" << idBlob << endmsg;
    idBlob->setDirection(direction);

    return true;
}

//=============================================================================
// GetStartPosition
// is_vertex must be "true" when the vert is muon vertex
// is_vertex must be "false" when the vert is for reference
//=============================================================================
bool HTBlob::GetStartPosition( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vert, bool is_vertex ) const
{

    debug() << " HTBlob::GetStartPosition " << endmsg;

    TH2D *hU = new TH2D ( "hU", "hU", 480,4510,9990,127,-1075,1075);
    TH2D *hV = new TH2D ( "hV", "hV", 480,4510,9990,127,-1075,1075);

    if ( is_vertex ){	
        hU->Fill( vert.z(), m_mathTool->calcUfromXY(vert.x(), vert.y()), 20 );
        hV->Fill( vert.z(), m_mathTool->calcVfromXY(vert.x(), vert.y()), 20);
    }

    SmartRefVector<Minerva::IDCluster> idClusters = idBlob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusters.begin();
    Gaudi::XYZPoint pos;

    double Dx = 0.0;
    double Dz = 0.0;
    double distance = 0.0;

    double vt_x = -9999;
    double vt_u;
    double vt_v;
    double vt_z = -9999;

    double vtX = -9999;
    double vtY = -9999;
    double vtZ = -9999;


    double dis_min_x = 10000;
    double xu[2] = {0.0};
    double xv[2] = {0.0};
    double zu[2] = {0.0};
    double zv[2] = {0.0}; // to NC Pi0 events, there is no muon vertex

    int countu = 0;
    int countv = 0; // to NC Pi0 events, there is no muon vertex

    for ( ; itClus != idClusters.end(); itClus++ ){

        Dx = (*itClus)->position() - vert.x(); 
        Dz = (*itClus)->z() - vert.z();// to avoid 0

        if( (*itClus)->view()== Minerva::IDCluster::X ) {
            distance = sqrt( pow(Dx,2) + pow(Dz,2) );
            if (distance <= dis_min_x  ) {
                dis_min_x = distance;
                vt_x = (*itClus)->position();
                vt_z = (*itClus)->z();
            }
        }

        if( (*itClus)->view()== Minerva::IDCluster::U ){
            debug() <<  " StartPoint U view, pe " << (*itClus)->pe() << "; z = " << (*itClus)->z() << "; coord " << (*itClus)->position() << endmsg;
            Dx = (*itClus)->position() - m_mathTool->calcUfromXY(vert.x(),vert.y());
            distance = sqrt( pow(Dx,2) + pow(Dz,2) );

            if ( is_vertex ) {
                hU->Fill( (*itClus)->z()-12,(*itClus)->tpos1(), (*itClus)->pe()/distance );
                hU->Fill( (*itClus)->z()+12,(*itClus)->tpos2(), (*itClus)->pe()/distance );
            }
            hU->Fill( (*itClus)->z(),(*itClus)->position(), (*itClus)->pe()/distance );
            if ( countu < 2 ){
                zu[countu] = (*itClus)->z();
                xu[countu] = (*itClus)->position();
                countu++;
            }
        }

        if( (*itClus)->view()== Minerva::IDCluster::V ){
            debug() <<  " StartPoint V view, pe " << (*itClus)->pe() << "; z = " << (*itClus)->z() << "; coord " << (*itClus)->position() << endmsg;
            Dx = (*itClus)->position() -  m_mathTool->calcVfromXY(vert.x(),vert.y());
            distance = sqrt( pow(Dx,2) + pow(Dz,2) );

            if ( is_vertex ){
                hV->Fill( (*itClus)->z()-12,(*itClus)->tpos1(), (*itClus)->pe()/distance );
                hV->Fill( (*itClus)->z()+12,(*itClus)->tpos2(), (*itClus)->pe()/distance );
            }
            hV->Fill( (*itClus)->z(),(*itClus)->position(), (*itClus)->pe()/distance );
            if ( countv < 2 ){
                zv[countu] = (*itClus)->z();
                xv[countu] = (*itClus)->position();
                countv++;
            }
        }
    }

    TF1 *fU, *fV;

    double slopeu = -9999;
    double slopev = -9999;
    double bu = -9999;
    double bv = -9999;

    bool goodFit_U = false;
    bool goodFit_V = false;
    if ( hU->GetEntries() > 2 ){
        hU->Fit("pol1","Q0");
        fU = hU->GetFunction("pol1");
        bu = fU->GetParameter(0);
        slopeu = fU->GetParameter(1);
        goodFit_U = true;

        delete fU;
    }
    else if ( hU->GetEntries() == 2 ){ // to deal with 2 clusters on NCPi0
        if ( zu[0] > zu[1] ){
            slopeu = (xu[0] - xu[1]) / (zu[0] - zu[1]);
            bu = xu[1] - zu[1]*slopeu;
            goodFit_U = true;
        }
        else if (zu[0] < zu[1] ) {
            slopeu = (xu[1] - xu[0]) / (zu[1] - zu[0]);
            bu = xu[0] - zu[0]*slopeu;
            goodFit_U = true;
        }
    }

    if ( hV->GetEntries() > 2 ){
        hV->Fit("pol1","Q0");
        fV = hV->GetFunction("pol1");
        bv = fV->GetParameter(0);
        slopev = fV->GetParameter(1);
        goodFit_V = true;

        delete fV;
    }
    else if ( hV->GetEntries() == 2 ){ // to deal with 2 clusters on NCPi0
        if ( zv[0] > zv[1] ){
            slopev = (xv[0] - xv[1]) / (zv[0] - zv[1]);
            bv = xv[1] - zv[1]*slopeu;
            goodFit_V = true;
        }
        else if (zu[1] < zu[0] ) { /* Trung: why zu instead of zv? */
            slopev = (xv[1] - xv[0]) / (zv[1] - zv[0]);
            bv = xv[0] - zv[0]*slopev;
            goodFit_V = true;
        }
    }

    vtX = vt_x;
    vtZ = vt_z;
    debug() << " Startpoint, slope u " << slopeu << " slope v" << slopev << endmsg;
    if ( goodFit_U && goodFit_V ){ 	     //3D blobs
        vt_u = slopeu*vt_z + bu;
        vt_v = slopev*vt_z + bv;
        vtY = m_mathTool->calcYfromUV(vt_u,vt_v);
    }
    else if ( goodFit_U ) {                  //2D blobs 
        vt_u = slopeu*vt_z + bu;
        vtY = (vt_x*.5 - vt_u)*2/sqrt(3);    //calcYfromXU?
    }
    else if ( goodFit_V ) { 	             //2D blobs 
        vt_v = slopev*vt_z + bv;
        vtY = (vt_v - vt_x*.5)*2/sqrt(3);    //calcYfromXV?
    }

    pos.SetX(vtX); pos.SetY(vtY); pos.SetZ(vtZ);

    idBlob->setStartPoint(pos);

    debug() << " Setting StarPoint " << pos << " Blob" << idBlob << endmsg;

    delete hU;
    delete hV;

    if (goodFit_U || goodFit_V) return true;

    return false;

}

//=======================================================================
//  idBlobdEdx  4 planes
//=======================================================================
StatusCode HTBlob::idBlobdEdx( Minerva::IDBlob *idblob, double &dEdx ) const
{
    debug() << "CCPi0HoughTool::idBlobdEdx"  << endmsg;

    SmartRefVector< Minerva::IDCluster > idClusters = idblob->clusters();
    SmartRefVector< Minerva::IDCluster >::iterator it_clus = idClusters.begin();
    TH1D *h = new TH1D ("h","h", 100,0,4500); // bin size 45mm <> 1 module
    dEdx = 0;

    Gaudi::XYZPoint vert = idblob->startPoint();
    double x, z, distance, vt_x, vt_u, vt_v, vt_z;

    // event vertex per view
    vt_x =  vert.x();
    vt_u = -vert.y()*sqrt(3)/2 + vert.x()*.5; //  -ycos30 + xsin30
    vt_v =  vert.y()*sqrt(3)/2 + vert.x()*.5; //   ycos30 + xsin30
    vt_z =  vert.z();


    for ( ; it_clus != idClusters.end(); it_clus++ ){

        x = vt_x - (*it_clus)->position();
        z = vt_z - (*it_clus)->z();

        if ( (*it_clus)->view() ==  Minerva::IDCluster::U ){
            x = vt_u - (*it_clus)->position();
            z = vt_z - (*it_clus)->z();
        }

        if ( (*it_clus)->view() ==  Minerva::IDCluster::V ){
            x = vt_v - (*it_clus)->position();
            z = vt_z - (*it_clus)->z();
        }

        distance = sqrt( pow(x,2) + pow(z,2) );
        h->Fill(distance, (*it_clus)->energy());
    }

    // dEdx calculation, maybe not clear
    int binmin = h->FindFirstBinAbove();
    int binmax = (h->FindLastBinAbove()-m_planesdEdx/2);

    for ( int i = binmin; i <= binmax; i++ ){
        if ( FinddEdxPlanes(h, i, dEdx) ) break;
        dEdx = 0;
    }
    if ( dEdx == 0 ) dEdx = -999;
    else dEdx = dEdx/m_planesdEdx;

    info () << " dEdx = " << dEdx << " number planes " << m_planesdEdx << endmsg;

    delete h;

    return StatusCode::SUCCESS;

}

//=======================================================================
//  FinddEdxPlanes
//=======================================================================
StatusCode HTBlob::FinddEdxPlanes(TH1D *h, int &index, double &dEdx) const
{

    dEdx = 0;
    int count = 0;

    for ( int i = index; i <= m_planesdEdx/2; i++){
        if ( h->GetBinContent(i) > 0 ) {
            count++;
            dEdx += h->GetBinContent(i);
        }
        else { index = i; break; }
    }
    if ( count == m_planesdEdx/2 ) 	return StatusCode::SUCCESS;
    else 	return StatusCode::FAILURE;

}

//=======================================================================
//  isPhoton
//  Calculating distance from event vertex to photon, 
//  this distance must be bigger than
//  m_minDistanceStripPhoton and  m_minDistanceModulePhoton
//=======================================================================
StatusCode HTBlob::isPhoton( SmartRefVector<Minerva::IDCluster> Seed, Gaudi::XYZPoint vtX ) const
{

    debug() << " HTBlob::isPhoton, asking vtx_z = " << vtX.z()  << "; vtx_x " << vtX.x() << endmsg;

    SmartRefVector<Minerva::IDCluster>::iterator itClus = Seed.begin();
    double min_radius = 10000.0;
    Gaudi::XYZPoint upstream;

    for ( ; itClus != Seed.end(); ++itClus ) {

        if ( (*itClus)->view() != Minerva::IDCluster::X ) continue;

        double radius =sqrt( pow(vtX.x()-(*itClus)->position(),2) + pow(vtX.z() - (*itClus)->z(),2) );
        if ( radius < min_radius ){
            min_radius = radius;
            upstream.SetX( (*itClus)->position() );
            upstream.SetY( 0.0 );
            upstream.SetZ( (*itClus)->z() );
        }
    }

    if ( fabs(upstream.x()-vtX.x()) > m_minDistanceStripPhoton ||
            fabs(upstream.z()-vtX.z()) > m_minDistanceModulePhoton ) 
        return StatusCode::SUCCESS;

    else 
        return StatusCode::FAILURE;

}

//=======================================================================
//  getBlobEnergy - Old Version - Jaewon's Study 
//=======================================================================
void HTBlob::getBlobEnergyTime_Old( Minerva::IDBlob *idblob, double &energy,
        double& tracker_evis,
        double& ecal_evis,
        double& hcal_evis,
        double& scal_evis) const
{

    debug() << "HTBlob::getBlobEnergy"  << endmsg;
    SmartRefVector< Minerva::IDCluster > idClusters = idblob->clusters();
    SmartRefVector< Minerva::IDCluster >::iterator it_clus = idClusters.begin();
    SmartRefVector< Minerva::IDDigit >::iterator it_dig;

    double m_scalefactor = 1.326;
    double m_calibrationTracker = 1.0;
    double m_calibrationECal = 2.341;
    double m_calibrationHCal = 9.54;

    energy = 0;
    double time = 0, total_pe = 0, factor = 2;

    tracker_evis = 0.0;
    ecal_evis    = 0.0;
    hcal_evis    = 0.0;
    scal_evis    = 0.0;

    for ( ; it_clus != idClusters.end(); it_clus++ ){
        time += (*it_clus)->time()*(*it_clus)->pe();
        total_pe += (*it_clus)->pe();
        const double cluster_energy = (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::ECAL ) {
            energy    += cluster_energy*m_scalefactor*m_calibrationECal;
            ecal_evis += cluster_energy;
        }

        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::HCAL ) {
            energy += cluster_energy*m_scalefactor*m_calibrationHCal;
            hcal_evis   += cluster_energy;
        }

        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::Tracker ) {

            if ( (*it_clus)->view() == Minerva::IDCluster::X ) factor = 2;
            else factor = 4;

            SmartRefVector< Minerva::IDDigit > idDigits = (*it_clus)->centralDigits();
            SmartRefVector< Minerva::IDDigit > sideDigits = (*it_clus)->sideEcalDigits();

            for ( it_dig = idDigits.begin(); it_dig != idDigits.end(); it_dig++ ){
                energy += (*it_dig)->normEnergy()*m_scalefactor*m_calibrationTracker;
                tracker_evis += (*it_dig)->normEnergy();
            }

            for ( it_dig = sideDigits.begin(); it_dig != sideDigits.end(); it_dig++ ){ //Jaewon formula DocDB 7950
                energy += (*it_dig)->normEnergy()*m_scalefactor*(factor*m_calibrationECal-1);
                scal_evis += (*it_dig)->normEnergy();
            }

        }
    }

    if ( total_pe > 0 )  time = time/total_pe;
    else time = 0;

    idblob->setTime(time);
    debug() << " Setting time " << time << " Id blob " << *idblob << endmsg;
}

//=======================================================================
//  getBlobEnergy -- New Method 
//=======================================================================
double HTBlob::getBlobEnergyTime_New( Minerva::IDBlob *idblob, std::vector<double>& evis_v, std::vector<double>& energy_v) const
{
    debug()<<"Enter HTBlob::getBlobEnergy() -- New Method"<<endmsg;
    SmartRefVector< Minerva::IDCluster > idClusters = idblob->clusters();
    SmartRefVector< Minerva::IDCluster >::iterator it_clus = idClusters.begin();
    SmartRefVector< Minerva::IDDigit >::iterator it_dig;

    double time = 0;
    double total_pe = 0;

    bool isView_X = false;
    double tracker_evis = 0.0;
    double ecal_evis = 0.0;
    double hcal_evis = 0.0;
    double scal_evis_X = 0.0;
    double scal_evis_UV = 0.0;

    // Loop Over all Clusters and Collect all Visible Energy
    for ( ; it_clus != idClusters.end(); it_clus++ ){
        time += (*it_clus)->time()*(*it_clus)->pe();
        total_pe += (*it_clus)->pe();

        const double cluster_energy = (*it_clus)->energy();

        // Tracker: Central and SideECAL
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::Tracker ) {

            if ( (*it_clus)->view() == Minerva::IDCluster::X ) isView_X = true;
            else isView_X = false;

            SmartRefVector< Minerva::IDDigit > idDigits = (*it_clus)->centralDigits();
            SmartRefVector< Minerva::IDDigit > sideDigits = (*it_clus)->sideEcalDigits();

            for ( it_dig = idDigits.begin(); it_dig != idDigits.end(); it_dig++ ){
                tracker_evis += (*it_dig)->normEnergy();
            }

            for ( it_dig = sideDigits.begin(); it_dig != sideDigits.end(); it_dig++ ){
                if(isView_X) scal_evis_X += (*it_dig)->normEnergy();
                else scal_evis_UV += (*it_dig)->normEnergy();
            }
        }
       
        // Downstream ECAL
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::ECAL ) {
            ecal_evis += cluster_energy;
        }

        // Downstream HCAL
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::HCAL ) {
            hcal_evis   += cluster_energy;
        }
    }

    // Fill evis_v
    evis_v[0] = tracker_evis;
    evis_v[1] = ecal_evis;
    evis_v[2] = scal_evis_X;
    evis_v[3] = scal_evis_UV;
    evis_v[4] = hcal_evis;

    // Get Shower Energy from Visible Energy
    double total_energy = getShowerEnergy(evis_v,energy_v);
  
    // Debugging 
    debug()<<"tracker evis = "<<evis_v[0]<<" energy = "<<energy_v[0]<<endmsg;
    debug()<<"ecal evis = "<<evis_v[1]<<" energy = "<<energy_v[1]<<endmsg;
    debug()<<"scal_X evis = "<<evis_v[2]<<" energy = "<<energy_v[2]<<endmsg;
    debug()<<"scal_UV evis = "<<evis_v[3]<<" energy = "<<energy_v[3]<<endmsg;
    debug()<<"hcal evis = "<<evis_v[4]<<" energy = "<<energy_v[4]<<endmsg;

    // Get Time
    if ( total_pe > 0 )  time = time/total_pe;
    else time = 0;

    idblob->setTime(time);
    debug() << " Setting time " << time << " Id blob " << *idblob << endmsg;
    
    debug()<<"Exit HTBlob::getBlobEnergy() -- New Method"<<endmsg;

    return total_energy;
}

/*
 *  Get Calorimetric Energy from Visible Energy
 *      Uses new kT, new kE, SAME kS and kH          
 */
double HTBlob::getShowerEnergy(std::vector<double>& evis_v, std::vector<double>& energy_v) const
{
    // Get Visible Energy in each subdetector
    double trkr_evis = evis_v[0];
    double ecal_evis = evis_v[1];
    double scal_evis_X = evis_v[2];
    double scal_evis_UV = evis_v[3];
    double hcal_evis = evis_v[4];

    // Get kT from evis inside Tracker
    double kT = get_kT(trkr_evis);
    
    // Get Actual Energy in each subdetector
    double trkr_energy = trkr_evis * kT;
    double ecal_energy = ecal_evis * m_kE;
    double scal_energy_X = scal_evis_X * m_kS_X;
    double scal_energy_UV = scal_evis_UV * m_kS_UV;
    double hcal_energy = hcal_evis * m_kH;

    // Fill energy_v
    energy_v[0] = trkr_energy;
    energy_v[1] = ecal_energy;
    energy_v[2] = scal_energy_X;
    energy_v[3] = scal_energy_UV;
    energy_v[4] = hcal_energy;

    double total_energy = trkr_energy + ecal_energy + scal_energy_X + scal_energy_UV + hcal_energy;

    return total_energy;
}

double HTBlob::get_kT(double evis) const
{
    const double evis_limit = 697.77; // MeV
    const double kT_default = 1.326;
    double kT;

    if ( evis < evis_limit){
        const double m = -0.0002336;
        const double c = 1.489;
        kT = m*evis + c;
    }else{ 
        kT = kT_default; 
    }

    return kT;
}

//=======================================================================
//  Finalize
//=======================================================================
StatusCode HTBlob::finalize()
{
    debug() << "Finalizing HTBlob..." << endmsg;
    StatusCode sc = this->MinervaHistoTool::finalize();
    if( sc.isFailure() ) { return Error( "Failed to finalize!", sc ); }
    return sc;
}

