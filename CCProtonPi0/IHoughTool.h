#ifndef IHOUGHTOOL_H 
#define IHOUGHTOOL_H 1

// Include files
// -------------
//inheritance
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/SmartRefVector.h"

// includes
#include "BlobFormation/BlobSeedDef.h" //for seedCandVector (could be forwarded?)
#include "GaudiKernel/Point3DTypes.h"

// forwards
namespace Minerva {
  class IDCluster;
}
static const InterfaceID IID_IHoughTool( "IHoughTool", 1, 0 );

class IHoughTool : virtual public IAlgTool {

  public:

    /// Return the interface ID
    static const InterfaceID& interfaceID() { return IID_IHoughTool; }

    virtual StatusCode initialize() = 0;
    virtual StatusCode finalize()   = 0;

    //todo: these clustervecs should be passed by (const?) reference
		virtual StatusCode GetReference( SmartRefVector<Minerva::IDCluster> idClusterView, Gaudi::XYZPoint &ref ) const = 0;
    virtual StatusCode GetReference( BlobSeeds::seedCandVector Seeds, Gaudi::XYZPoint &ref ) const = 0;

		virtual StatusCode Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta, const Gaudi::XYZPoint& ref ) const = 0;
		virtual StatusCode Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta, const Gaudi::XYZPoint& ref, const Gaudi::XYZPoint& vert ) const = 0;
    virtual StatusCode Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta, const Gaudi::XYZPoint& ref ) const = 0;
    virtual StatusCode Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta, const Gaudi::XYZPoint& ref, const Gaudi::XYZPoint& vert ) const = 0;


};
#endif 
