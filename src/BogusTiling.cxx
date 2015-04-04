#include "WireCellTiling/BogusTiling.h"
#include "WireCellData/GeomCell.h"
#include "WireCellData/GeomWire.h"

using namespace WireCell;

BogusTiling::BogusTiling()
{
}

BogusTiling::~BogusTiling()
{
}

GeomWireSelection BogusTiling::wires(const GeomCell& cell) const
{
    return GeomWireSelection();
}
	
GeomCellSelection BogusTiling::cells(const GeomWire& wire) const
{
    return GeomCellSelection();
}

GeomCell* BogusTiling::cell(const GeomWireSelection& wires) const
{
    return 0;
}
