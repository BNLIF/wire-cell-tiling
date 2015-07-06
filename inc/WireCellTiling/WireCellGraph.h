#ifndef WIRECELLTILING_WIRECELLGRAPH
#define WIRECELLTILING_WIRECELLGRAPH

#include "WireCellData/GeomCell.h"

#include "WireCellIface/IWireDatabase.h"
#include "WireCellUtil/IndexedSet.h"

#include <boost/graph/adjacency_list.hpp>
#include <vector>


namespace WireCell {


    namespace Graph {


	/// Each vertex of the graph has a "type" and an index into a
	/// vector of objects of corresponding type.:
	/// - cell indexes into an array of GeomCell
	/// - wire indexes into an array of GeomWire
	/// - point index into an array of Point 
	struct Property {
	    enum VertexType_t { unknown=0, cell, wire, point};
	    VertexType_t vtype;
	    int index;
	    Property(VertexType_t vt=unknown, int ind=-1) : vtype(vt), index(ind) {}

	    bool operator<(const Property& rhs) const {
		if (vtype < rhs.vtype) return true;
		return index < rhs.index;
	    }

	};
	/// The edges are implicitly of these four types:
	/// - wirecell :: cell <--> wire
	/// - neighbor :: cell <--> cell
	/// - boundary :: point <--> point
	/// - internal :: cell <--> point

	typedef boost::adjacency_list<boost::setS, boost::vecS, 
				      boost::undirectedS,
				      Property > Type;	
	typedef boost::graph_traits<Type>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Type>::edge_descriptor Edge;


	class Filler {
	public:
	    Filler(const IWireDatabase& wdb);

	    /// BGL does not apparel support lookups by property, so DIY. 
	    mutable std::map<Property, Vertex> cellVertex, wireVertex;

	    Vertex cell_vertex(const GeomCell& cell) const;
	    Vertex wire_vertex(const GeomWire& wire) const;

	    /// The wire db used.  The "point" type vertices of the
	    /// graph index into the vector returned by
	    /// wdb.get_wires()
	    const IWireDatabase& wdb;

	    /// Owning collection of geom cells
	    WireCell::GeomCellSet cellset;

	    /// Pertinent points.
	    mutable WireCell::IndexedSet<WireCell::Point> point_index;
	    /// The cells
	    mutable WireCell::IndexedSet<const WireCell::GeomCell*> cell_index;
	    /// The wires
	    mutable WireCell::IndexedSet<const WireCell::GeomWire*> wire_index;

	    Type graph;
	    

	    

	private:
	    void generate();
	    void record(const WireCell::PointVector& pcell,
			const WireCell::GeomWireSelection& uvw_wires);
	    const GeomWireSelection& index_wires(WireCell::WirePlaneType_t wiretype);
	};



    } // namespace WireCell::Graph

} // namespace WireCell
#endif
