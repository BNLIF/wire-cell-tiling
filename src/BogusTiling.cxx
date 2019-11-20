#include "WCPTiling/BogusTiling.h"
#include "WCPData/GeomCell.h"
#include "WCPData/GeomWire.h"

using namespace WCP;

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
