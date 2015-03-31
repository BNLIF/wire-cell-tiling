#ifndef WIRECELL_BOGUSTILING_H
#define WIRECELL_BOGUSTILING_H

#include "WireCellTiling/TilingBase.h"

namespace WireCell {
    /**
     *  A bogus tiling class that doesn't do anything.
     */
    class BogusTiling : public TilingBase { 
    public:
	BogusTiling();
	virtual ~BogusTiling();

	WireCell::WireSelection wires(const WireCell::Cell& cell) const;
	WireCell::CellSelection cells(const WireCell::Wire& wire) const;
	virtual WireCell::Cell* cell(const WireCell::WireSelection& wires) const;

    };
}
#endif
