#include "HTtool.h"

#include "Event/IDCluster.h"

#include <cmath>
#include "TMath.h"
#include "TH2D.h"

DECLARE_TOOL_FACTORY( HTtool );

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
HTtool::HTtool( const std::string& type, const std::string& name, const IInterface* parent ) : 
  MinervaHistoTool( type, name, parent ) 
{
    debug() << "Instantiating HTtool..." << endmsg;
    declareInterface<IHoughTool>(this);
	
}

//=============================================================================
/// Destructor
//=============================================================================
HTtool::~HTtool() {}


//=============================================================================
// Initialize
//=============================================================================
StatusCode HTtool::initialize() 
{
 
    debug() << "Initializing HTtool..." << endmsg;
    StatusCode sc = this->MinervaHistoTool::initialize();
    if( sc.isFailure() ) { return Error( "Failed to initialize!", sc ); }
    
    return sc;
}



//=============================================================================
// GetReference - Clusters
//=============================================================================
StatusCode HTtool::GetReference( SmartRefVector<Minerva::IDCluster> idClusterView, Gaudi::XYZPoint &ref ) const
{
    if( !idClusterView.size() ) {
        warning() << " Input Vector of Clusters to the GetReference tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }
    
    SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusterView.begin();
    double maxpe = 0, total_pe = 0;
    
    for ( ; itClus != idClusterView.end(); itClus++ ){
        
        total_pe += (*itClus)->pe();
        if ( (*itClus)->pe() > maxpe  ) {
            maxpe = (*itClus)->pe();
            ref.SetX((*itClus)->position());
            ref.SetY(0);
            ref.SetZ((*itClus)->z());
        }
        
    }

    if ( total_pe < 20 ) {
        debug() << " ALERT: No enough energy to use Hough Transformation " << endmsg;
        return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
}

//=============================================================================
// GetReference - Seeds
//=============================================================================
StatusCode HTtool::GetReference( BlobSeeds::seedCandVector Seeds, Gaudi::XYZPoint &ref ) const
{
    if( !Seeds.size() ) {
        warning() << " Input Vector of Seeds to the GetReference tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }

    BlobSeeds::seedCandVector::iterator itSeed = Seeds.begin();
    double maxpe = 0, total_pe = 0;
    
    for ( ; itSeed != Seeds.end(); itSeed++ ){
        
        total_pe += (*itSeed)->seedpe;
        if ( (*itSeed)->seedpe > maxpe  ) {
            maxpe = (*itSeed)->seedpe;
            ref.SetX( (*itSeed)->seedcoordcentroid );
            ref.SetY( 0 );
            ref.SetZ( (*itSeed)->seedzcentroid );
        }
        
    }

    debug() << " Setting Ref point " << ref << endmsg;
  
    if ( total_pe < 20 ) {
        debug() << " ALERT: No enough energy to use Hough Transformation " << endmsg;
        return StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
}

//=============================================================================
// Hough2D - Clusters
//=============================================================================
StatusCode HTtool::Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta,
                            const Gaudi::XYZPoint &ref ) const
{
	
    debug() << "Working in ConeScanCreator::Hough2D " << endmsg;
    
    if( !idClusterVec.size() ) { 
        debug() << " ALERT: Input Vector of Clusters to the Hough2D tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }
	
    TH2D *h = new TH2D("h","Hough Space",46,-92,92,275,-5500,5500); // less granulated because lot of hits
    SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusterVec.begin();
    
    for ( ; itClus != idClusterVec.end(); itClus++ ) FillHough1Cluster(*itClus, h, ref);
    
    int x,y,z, maxbin;
    if ( h->GetEntries() ){
        maxbin = h->GetMaximumBin(x, y, z);
        r      = h->GetYaxis()->GetBinCenter(y);
        theta  = h->GetXaxis()->GetBinCenter(x);
    }
    else {
        debug() << " ALERT: Hough space is empty " << endmsg;
        delete h;
        return StatusCode::FAILURE;
    }
    
    delete h;
    return StatusCode::SUCCESS;
    
}

//=============================================================================
// Hough2D - Clusters Overload with vertex
//=============================================================================
StatusCode HTtool::Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta,
                            const Gaudi::XYZPoint &ref, const Gaudi::XYZPoint& vert ) const
{
	
    debug() << "Working in ConeScanCreator::Hough2D " << endmsg;
    
    if( !idClusterVec.size() ) { 
        debug() << " ALERT: Input Vector of Clusters to the Hough2D tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }
	
    TH2D *h = new TH2D("h","Hough Space",46,-92,92,275,-5500,5500); // less granulated because lot of hits
    SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusterVec.begin();
    
    double x = vert.x()-ref.x();
    double z = vert.z()-ref.z();
    double radius, theta_var, maxpe = 0;
    debug() << " Ref " << ref << " vert " << vert << endmsg; 
  
    for ( ; itClus != idClusterVec.end(); itClus++ ) {
        if ( (*itClus)->pe() > maxpe ) maxpe = (*itClus)->pe();
        FillHough1Cluster(*itClus, h, ref);
    }

    debug() << " Hough2D we found maxpe for vertex " << maxpe << endmsg;

        //must be filled with the vertex direction and weighted by max pe
    for ( int i = h->GetXaxis()->GetFirst(); i <= h->GetXaxis()->GetLast(); i++){
        theta_var = h->GetXaxis()->GetBinCenter(i)*CLHEP::pi/180; // angle radians
        radius = x*sin(theta_var) + z*cos(theta_var);
        h->Fill( theta_var*180/CLHEP::pi, radius, maxpe / 2 ); // must be filled in degrees
    }

    int xbin,ybin,zbin, maxbin;
    if ( h->GetEntries() ) {
        maxbin = h->GetMaximumBin(xbin, ybin, zbin);
        r      = h->GetYaxis()->GetBinCenter(ybin);
        theta  = h->GetXaxis()->GetBinCenter(xbin);
    }
    else {
        debug() << " ALERT: Hough space is empty " << endmsg;
        delete h;
        return StatusCode::FAILURE;
    }

    delete h;
    return StatusCode::SUCCESS;
    
}

//=============================================================================
// Hough2D - Seeds
//=============================================================================
StatusCode HTtool::Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta,
                            const Gaudi::XYZPoint& ref ) const
{
	
    debug() << "Working in ConeScanCreator::Hough2D - Seeds" << endmsg;
    
    if( !Seeds.size() ) { 
        debug() << " ALERT: Input Vector of Seeds to the Hough2D tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }
    
    TH2D *h = new TH2D("h","Hough Space",90,-90,90,275,-5500,5500);
    BlobSeeds::seedCandVector::iterator itSeed = Seeds.begin();
    
    for ( ; itSeed != Seeds.end(); itSeed++ ) FillHough1Cluster(*itSeed, h, ref);
    
    int x,y,z, maxbin;
    if ( h->GetEntries() ) {
        maxbin = h->GetMaximumBin(x, y, z);
        r      = h->GetYaxis()->GetBinCenter(y);
        theta  = h->GetXaxis()->GetBinCenter(x);
    }
    else {
        debug() << " ALERT: Hough space is empty " << endmsg;
        delete h;
        return StatusCode::FAILURE;
    }

    delete h;
    return StatusCode::SUCCESS;
    
}

//=============================================================================
// Hough2D - Seeds Overload
//=============================================================================
StatusCode HTtool::Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta, const Gaudi::XYZPoint& ref,
                            const Gaudi::XYZPoint& vert ) const
{
	
    debug() << "Working in ConeScanCreator::Hough2D " << endmsg;
    
    if( !Seeds.size() ) { 
        debug() << " ALERT: Input Vector of Seeds to the Hough2D tool is empty " << endmsg;
        return StatusCode::FAILURE;
    }
    
	//TH2D *h = new TH2D("h","Hough Space",45,-90,90,275,-5500,5500);
    TH2D *h = new TH2D("h","Hough Space",92,-92,92,260,-5500,5500); // more granulated because few points?

    double x = vert.x()-ref.x();
    double z = vert.z()-ref.z();
    double radius, theta_var, maxpe = 0;
    
    BlobSeeds::seedCandVector::iterator itSeed = Seeds.begin();
    for ( ; itSeed != Seeds.end(); itSeed++ ) {
        if ( (*itSeed)->seedpe > maxpe ) maxpe = (*itSeed)->seedpe;
        FillHough1Cluster( *itSeed, h, ref );
    }
    
    debug() << " Hough2D we found maxpe for vertex " << maxpe << endmsg;

        //must be filled with the vertex direction and weighted by max pe
    for ( int i = h->GetXaxis()->GetFirst(); i <= h->GetXaxis()->GetLast(); i++){
        theta_var = h->GetXaxis()->GetBinCenter(i)*CLHEP::pi/180; // angle radians
        radius = x*sin(theta_var) + z*cos(theta_var);
        h->Fill( theta_var*180/CLHEP::pi, radius, maxpe ); // must be filled in degrees
    }
    
    int xbin, ybin, zbin, maxbin;
    if ( h->GetEntries() ) {
        maxbin = h->GetMaximumBin(xbin, ybin, zbin);
        r      = h->GetYaxis()->GetBinCenter(ybin);
        theta  = h->GetXaxis()->GetBinCenter(xbin);
        debug() << " Hough2D we found theta " << theta << " radius " << r << endmsg;
    }
    else {
        debug() << " ALERT: Hough space is empty " << endmsg;
        delete h;

        return StatusCode::FAILURE;
    }
    
    delete h;

    return StatusCode::SUCCESS;

}

//=============================================================================
// FillHough1Cluster - Cluster
//=============================================================================
StatusCode HTtool::FillHough1Cluster( SmartRef<Minerva::IDCluster> idCluster, TH2D *h, const Gaudi::XYZPoint& ref) const
{
    double x = idCluster->position()-ref.x();
    double z = idCluster->z()-ref.z();
    double r, theta;
    debug() << " Fill Hough; pos " << idCluster->position() << "; z clus " << idCluster->z()
            << "; x " << x << "; z " << z << "; pe " << idCluster->pe() << endmsg;
    
    for ( int i = h->GetXaxis()->GetFirst(); i <= h->GetXaxis()->GetLast(); i++){
        
        theta = h->GetXaxis()->GetBinCenter(i)*CLHEP::pi/180; // angle radians
        r = x*sin(theta) + z*cos(theta);
        h->Fill( theta*180/CLHEP::pi, r, idCluster->pe() ); // must be filled in degree
        
    }

    return StatusCode::SUCCESS;
}

//=============================================================================
// FillHough1Cluster - Seed
//=============================================================================
StatusCode HTtool::FillHough1Cluster( BlobSeeds::seedCandidate *Seed, TH2D *h, const Gaudi::XYZPoint& ref) const
{
    double x = Seed->seedcoordcentroid-ref.x();
    double z = Seed->seedzcentroid-ref.z();
    double r, theta;
    debug() << " Fill Hough; pos " << Seed->seedcoordcentroid << "; x " << x << "; z clus "
            << Seed->seedzcentroid << "; z " << z << "; pe " <<  Seed->seedpe << endmsg;
    
    for ( int i = h->GetXaxis()->GetFirst(); i <= h->GetXaxis()->GetLast(); i++){
        
        theta = h->GetXaxis()->GetBinCenter(i)*CLHEP::pi/180; // angle radians
        r = x*sin(theta) + z*cos(theta);
        h->Fill( theta*180/CLHEP::pi, r, Seed->seedpe ); // must be filled in degree
        
    }
    
    return StatusCode::SUCCESS;
}


//=============================================================================
// Finalize
//=============================================================================
StatusCode HTtool::finalize() 
{
    
    debug() << "Finalizing HTtool..." << endmsg;
    StatusCode sc = this->MinervaHistoTool::finalize();
    if( sc.isFailure() ) { return Error( "Failed to finalize!", sc ); }
    return sc;
    
}
