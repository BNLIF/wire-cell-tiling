#ifndef WIRECELL_TILINGBASE_H
#define WIRECELL_TILINGBASE_H


#include "WireCellNav/GeomDataSource.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"

namespace WireCell {

    /** WireCell::TilingBase - base class for providing a tiling of 2D
     * space with cells based on planes of parallel wires.
     *
     * A Tiling defines a number of "cells" which completely fill a
     * region of 2D space.  These cells are expected to be related to
     * a number of sets of parallel wires.  Each set of parallel wires
     * fall into one of three wire planes called U, V and Y.  Y wire
     * planes run parallel to the Y-axis.  U and V wires are at some
     * angle w.r.t. the Y-axis.
     *
     * The methods of the Tiling expose relationships between a cell
     * and the associated wires.
     */

    class TilingBase {
    public:

	virtual ~TilingBase();

	/// Must return all wires associated with the given cell
	virtual GeomWireSelection wires(const GeomCell& cell) const = 0;
	
	/// Must return all cells associated with the given wire
	virtual GeomCellSelection cells(const GeomWire& wire) const = 0;

	/// Must the one cell associated with the collection of wires or 0.
	virtual const GeomCell* cell(const GeomWireSelection& wires) const = 0;

    };

}

#endif
