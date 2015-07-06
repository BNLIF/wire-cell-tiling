#ifndef WIRECELL_BOGUSTILING_H
#define WIRECELL_BOGUSTILING_H

#include "WireCellIface/ICellTiling.h"

namespace WireCell {
    /**
     *  A bogus tiling class that doesn't do anything.
     */
    class BogusTiling : public ICellTiling { 
    public:
	BogusTiling();
	virtual ~BogusTiling();

	GeomWireSelection wires(const GeomCell& cell) const;
	GeomCellSelection cells(const GeomWire& wire) const;
	virtual GeomCell* cell(const GeomWireSelection& wires) const;

    };
}
#endif
