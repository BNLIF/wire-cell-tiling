#include "WireCellTiling/BogusTiling.h"
#include "WireCellData/Cell.h"
#include "WireCellData/Wire.h"

using namespace WireCell;

BogusTiling::BogusTiling()
{
}

BogusTiling::~BogusTiling()
{
}

WireSelection BogusTiling::wires(const Cell& cell) const
{
    return WireSelection();
}
	
CellSelection BogusTiling::cells(const Wire& wire) const
{
    return CellSelection();
}

Cell* BogusTiling::cell(const WireSelection& wires) const
{
    return 0;
}
