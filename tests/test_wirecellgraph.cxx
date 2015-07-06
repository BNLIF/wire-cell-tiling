#include "WireCellTiling/WireCellGraph.h"

#include "WireCellNav/ExampleWires.h"
#include "WireCellNav/GeomDataSource.h"

#include "WireCellUtil/Testing.h"

#include "WireCellIface/IWireGeometry.h"

using namespace WireCell;
using namespace std;

int main()
{
    IWireGeometry* wires = make_example_wires(10*units::mm);
    GeomDataSource gds;
    gds.use_wires(*wires);
    
    
    WireCell::Graph::Filler wcg(gds);

    cerr << "#cells: " << wcg.cellVertex.size() << " =?= " <<  wcg.cell_index.collection.size() << endl;
    cerr << "#wires: " << wcg.wireVertex.size() << " =?= " <<  wcg.wire_index.collection.size() << endl;

    Assert(wcg.cellVertex.size() == wcg.cell_index.collection.size());
    Assert(wcg.wireVertex.size() == wcg.wire_index.collection.size());
									 

}
