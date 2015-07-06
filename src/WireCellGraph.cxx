#include "WireCellTiling/WireCellGraph.h"

using namespace WireCell::Graph;

#include <iostream>
using namespace std;		// debugging

Filler::Filler(const IWireDatabase& wdb)
    : wdb(wdb)
{
    this->generate();
}



Vertex Filler::cell_vertex(const GeomCell& cell) const
{
    Property prop(Property::cell, cell_index(&cell));
    return cellVertex[prop];
}
Vertex Filler::wire_vertex(const GeomWire& wire) const
{
    Property prop(Property::wire, wire_index(&wire));
    return wireVertex[prop];
}



// static
// WireCell::Vector find_center(const WireCell::PointVector& pcell)
// {
//     WireCell::Point center;
//     size_t npos = pcell.size();
//     for (int ind=0; ind<npos; ++ind) {
// 	center = center + pcell[ind];
//     }
//     float norm = 1.0/npos;
//     return  norm * center;
// }


void WireCell::Graph::Filler::record(const WireCell::PointVector& pcell,
				     const WireCell::GeomWireSelection& uvw_wires)
{
    /// fixme: beware of duplicate entries!
    auto itok = cellset.insert(GeomCell(cellset.size(), pcell));
    const GeomCell* thecellptr = &(*itok.first);
    Property cell_prop(Property::cell, cell_index(thecellptr));
    Vertex cv = boost::add_vertex(cell_prop, graph);
    cellVertex[cell_prop] = cv;

    // corners
    size_t npos = pcell.size();
    std::vector<Vertex> corner_vertices;

    // wire vertices and wire-cell edges
    for (int ind=0; ind<uvw_wires.size(); ++ind) {
	Property prop(Property::wire, wire_index(uvw_wires[ind]));
	Vertex w = boost::add_vertex(prop, graph);
	wireVertex[prop] = w;
	boost::add_edge(cv, w, graph);
    }
    
    // corner vertices and internal edges
    for (int ind=0; ind < npos; ++ind) {
	Property prop(Property::point, point_index(pcell[ind]));
        Vertex v = boost::add_vertex(prop, graph);
	corner_vertices.push_back(v);

	// center to corner internal edges
	boost::add_edge(cv, v, graph);
    }

    // cell boundary edges
    for (int ind=0; ind<npos; ++ind) {
	boost::add_edge(corner_vertices[ind], corner_vertices[(ind+1)%npos], graph);
    }
    
    // done except nearest neighbors
}


const WireCell::GeomWireSelection& WireCell::Graph::Filler::index_wires(WireCell::WirePlaneType_t wiretype)
{
    const GeomWireSelection& pwires = wdb.wires_in_plane(wiretype);
    for (auto it = pwires.begin(); it != pwires.end(); ++it) {
	wire_index(*it);
    }
    return pwires;
}


void WireCell::Graph::Filler::generate()
{
    // This is Xin's cell definition algorithm sung to the tune of
    // graph minor.

    const WireCell::GeomWireSelection& u_wires = index_wires(kUwire);
    const WireCell::GeomWireSelection& v_wires = index_wires(kVwire);
    const WireCell::GeomWireSelection& w_wires = index_wires(kYwire);

    float pitch_u = wdb.pitch(kUwire);
    float pitch_v = wdb.pitch(kVwire);
    float pitch_w = wdb.pitch(kYwire);
    
    std::pair<int,int> box_ind[4] = { // allows loops over a box of indices
        std::pair<int,int>(0,0),
        std::pair<int,int>(0,1),
        std::pair<int,int>(1,1),
        std::pair<int,int>(1,0)
    };

    // pack it up and send it out
    GeomWireSelection uvw_wires(3);

    for (int u_ind=0; u_ind < u_wires.size(); ++u_ind) {
        const WireCell::GeomWire& u_wire = *u_wires[u_ind];
	uvw_wires[0] = &u_wire;
        float dis_u_wire = wdb.wire_dist(u_wire);
        float dis_u[2] = { dis_u_wire - pitch_u, dis_u_wire + pitch_u }; // half-line minmax

        for (int v_ind=0; v_ind < v_wires.size(); ++v_ind) {
            const WireCell::GeomWire& v_wire = *v_wires[v_ind];
	    uvw_wires[1] = &v_wire;
            float dis_v_wire = wdb.wire_dist(v_wire);
            float dis_v[2] = { dis_v_wire - pitch_v, dis_v_wire + pitch_v };


            std::vector<Vector> puv(4);
            float dis_puv[4];
            for (int ind=0; ind<4; ++ind) {
                // fixme: we are not handling the case where one
                // of these crossing points are outside the wire
                // plane boundary.
                wdb.crossing_point(dis_u[box_ind[ind].first], dis_v[box_ind[ind].second], kUwire, kVwire, puv[ind]);
                dis_puv[ind] = wdb.wire_dist(puv[ind], kYwire);
            }

            for (int w_ind=0; w_ind < w_wires.size(); ++w_ind) {
                const WireCell::GeomWire& w_wire = *w_wires[w_ind];
		uvw_wires[2] = &w_wire;
                float dis_w_wire = wdb.wire_dist(w_wire);
                float dis_w[2] = { dis_w_wire - pitch_w, dis_w_wire + pitch_w };

                WireCell::PointVector pcell;

                for (int ind=0; ind<4; ++ind) {
                    if (dis_w[0] <= dis_puv[ind] && dis_puv[ind] < dis_w[1]) {
                        pcell.push_back(puv[ind]);
                    }
                }

                if (!pcell.size()) { 
                    continue;
                }

                for (int ind=0; ind<4; ++ind) {
                    {           // fresh scope
                        WireCell::Vector pointvec;
                        if (wdb.crossing_point(dis_u[box_ind[ind].first], dis_w[box_ind[ind].second], 
                                               kUwire, kYwire, pointvec)) 
                            {
                                float disv = wdb.wire_dist(pointvec, kVwire);
                                if (dis_v[0] <= disv && disv < dis_v[1]) {
                                    pcell.push_back(pointvec);
                                }
                            }
                    }

                    {           // fresh scope
                        WireCell::Vector pointvec;
                        if (wdb.crossing_point(dis_v[box_ind[ind].first], dis_w[box_ind[ind].second], 
                                               kVwire, kYwire, pointvec)) 
                            {
                                float disu = wdb.wire_dist(pointvec, kUwire);
                                if (dis_u[0] <= disu && disu < dis_u[1]) {
                                    pcell.push_back(pointvec);
                                }
                            }
                    }

                }

                this->record(pcell, uvw_wires);

            } // W wires
        } // v wires
    } // u wires
} // generate()




