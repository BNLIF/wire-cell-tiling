#ifndef WIRECELL_GRAPHTILING_H
#define WIRECELL_GRAPHTILING_H

#include "WireCellIface/ICellTiling.h"
#include "WireCellIface/ICellGeometry.h"
#include "WireCellTiling/WireCellGraph.h"

namespace WireCell {
    /**
     *  The GraphTiling tiles the 2D wire plane and provides fast
     *  queries such as:
     *
     *  - all wires associated with a cell
     *  - all cells associated with a wire
     *  - the cell associated with a set of wires
     *  - the cells nearest a given cell
     */
    class GraphTiling : public ICellTiling,
			public ICellGeometry
    { 
	WireCell::Graph::Filler wcgraph;

    public:
	GraphTiling(const WireCell::IWireDatabase& wdb);
	virtual ~GraphTiling();

	/// ICellTiling
	GeomWireSelection wires(const GeomCell& cell) const;
	GeomCellSelection cells(const GeomWire& wire) const;
	const GeomCell* cell(const GeomWireSelection& wires) const;
	//GeomCellSelection neighbors(const GeomCell& cell) const;

	/// ICellGeometry
	const WireCell::GeomCellSet& get_cells() const;

    };

}
#endif
