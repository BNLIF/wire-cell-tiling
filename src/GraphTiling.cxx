#include "WireCellTiling/GraphTiling.h"

#include <iterator>		// std::back_inserter
#include <algorithm>		// std::set_difference

using namespace WireCell;
using namespace WireCell::Graph;

GraphTiling::GraphTiling(const IWireDatabase& wdb)
    : wcgraph(wdb)
{
}
GraphTiling::~GraphTiling()
{
}

const WireCell::GeomCellSet& GraphTiling::get_cells() const
{
    return wcgraph.cellset;
}


GeomWireSelection GraphTiling::wires(const GeomCell& cell) const
{
    Vertex v = wcgraph.cell_vertex(cell);
    Type::adjacency_iterator it, end;
    boost::tie(it,end) = boost::adjacent_vertices(v, wcgraph.graph);
    GeomWireSelection ret;
    for (; it != end; ++it) {
	const Property& p = wcgraph.graph[*it];
	if (p.vtype == Property::wire) {
	    ret.push_back(wcgraph.wire_index.collection[p.index]);
	}
    }
    return ret;
}

GeomCellSelection GraphTiling::cells(const GeomWire& wire) const
{
    Vertex v = wcgraph.wire_vertex(wire);
    Type::adjacency_iterator it, end;
    boost::tie(it,end) = boost::adjacent_vertices(v, wcgraph.graph);
    GeomCellSelection ret;
    for (; it != end; ++it) {
	const Property& p = wcgraph.graph[*it];
	if (p.vtype == Property::cell) {
	    ret.push_back(wcgraph.cell_index.collection[p.index]);
	}
    }
    return ret;

}

const GeomCell* GraphTiling::cell(const GeomWireSelection& wires) const
{
    // this method returns the first cell associated with one wire in
    // the selection which is also associated with the other wires in
    // the selection.

    GeomCellSelection cellsel = this->cells(*wires[0]);
    if (!cellsel.size()) {
	return 0;
    }

    GeomWireSelection others(wires.begin()+1, wires.end());

    for (auto cit = cellsel.begin(); cit != cellsel.end(); ++cit) {
	const GeomCell* thecell = *cit;
	GeomWireSelection wiresel = this->wires(*thecell);

	GeomWireSelection result;
	std::set_difference(others.begin(), others.end(),
			    wiresel.begin(), wiresel.end(),
			    std::back_inserter(result));
			    
	if (0 == result.size()) {
	    return thecell;
	}
    }

    return 0;
}
