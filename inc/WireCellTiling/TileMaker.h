#ifndef WIRECELL_TILEMAKER_H
#define WIRECELL_TILEMAKER_H

#include "WireCellTiling/TilingBase.h"

#include "WireCellNav/GeomDataSource.h"

#include "WireCellData/WireCellMap.h"

namespace WireCell {

    /** WireCellTiling::TileMaker - tiling using Michael Mooney's algorithm.

	This class is a transliterated copy of the tile generation
	code from the original monolithic CellMaker (+Plotter) app.
     */
    class TileMaker : public TilingBase { 
	// Our connection to the wire geometry
	const WireCell::GeomDataSource& geo;
	// What we make
	WireCell::CellSet cellset;
	WireCell::WireCellMap wcmap;
	
	// Cache some values between methods
	WireCell::WireSelection Uwires;
	WireCell::WireSelection Vwires;
	WireCell::WireSelection Ywires;
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


    public:
	TileMaker(const WireCell::GeomDataSource& geom);
	virtual ~TileMaker();

	// base API

	/// Load of tiling from a .wct file.  If null filename is given, will generate.
	bool load(const char* filename);

	/// Save current tiling the given filename.
	bool save(const char* filename);

	/// Must return all wires associated with the given cell
	WireCell::WireSelection wires(const WireCell::Cell& cell) const;
	
	/// Must return all cells associated with the given wire
	WireCell::CellSelection cells(const WireCell::Wire& wire) const;

    };

}
#endif
