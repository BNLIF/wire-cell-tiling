#ifndef WIRECELL_TILEMAKER_H
#define WIRECELL_TILEMAKER_H

#include "WireCellNav/GeomDataSource.h"

#include "WireCellData/GeomWireCellMap.h"

#include "WireCellIface/ICellTiling.h"

namespace WireCell {

    /** WireCellTiling::TileMaker - tiling using Michael Mooney's algorithm.

	This class is a transliterated copy of the tile generation
	code from the original monolithic CellMaker (+Plotter) app.
     */
    class TileMaker : public ICellTiling { 
    public:
	TileMaker(const WireCell::GeomDataSource& geom);
	virtual ~TileMaker();

	// base API

	/// Must return all wires associated with the given cell
	GeomWireSelection wires(const GeomCell& cell) const;
	
	/// Must return all cells associated with the given wire
	GeomCellSelection cells(const GeomWire& wire) const;

	/// Returns the one cell associated with the collection of wires or 0.
	virtual GeomCell* cell(const GeomWireSelection& wires) const;

    private:

	// Our connection to the wire geometry
	const GeomDataSource& geo;
	// What we make
	GeomCellSet cellset;
	GeomWireMap wiremap;
	GeomCellMap cellmap;
	
	// Cache some values between methods
	GeomWireSelection Uwires;
	GeomWireSelection Vwires;
	GeomWireSelection Ywires;
	double maxHeight;
	double UdeltaY, VdeltaY;
	double wirePitchU, wirePitchV, wirePitchY;
	double firstYwireZval, angleUrad, angleVrad;
	double firstYwireUoffsetYval, firstYwireVoffsetYval;
	double leftEdgeOffsetZval, rightEdgeOffsetZval;
	double UspacingOnWire, VspacingOnWire;

	void constructCells();
	void constructCellChain(double wireZval, double YvalOffsetU, double YvalOffsetV);
	void constructCell(double YwireZval, double UwireYval, double VwireYval);
	bool formsCell(double UwireYval, double VwireYval);
	std::vector<std::pair<double,double> > getCellVertices(double YwireZval, double UwireYval, double VwireYval);

	int getUwireID(double Yval, double Zval);
	int getVwireID(double Yval, double Zval);
	int getYwireID(double Zval);


    };

}
#endif
