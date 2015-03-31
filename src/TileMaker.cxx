#include "WireCellTiling/TileMaker.h"

#include <cmath>
#include <iostream>
#include <algorithm> 
using namespace WireCell;

const double epsilon = 0.0000000001;

TileMaker::TileMaker(const WireCell::GeomDataSource& geom)
    : TilingBase(), geo(geom)
{
    Uwires = geo.wires_in_plane(WireCell::kUtype);
    Vwires = geo.wires_in_plane(WireCell::kVtype);
    Ywires = geo.wires_in_plane(WireCell::kYtype);

    std::vector<float> ext = geo.extent();
    maxHeight = ext[1];
    
    wirePitchU = geo.pitch(WireCell::kUtype);
    wirePitchV = geo.pitch(WireCell::kVtype);
    wirePitchY = geo.pitch(WireCell::kYtype);

    angleUrad = geo.angle(WireCell::kUtype) / units::radian; // explicitly carry 
    angleVrad = geo.angle(WireCell::kVtype) / units::radian; // value in radians

    UdeltaY = wirePitchY/tan(angleUrad);
    VdeltaY = wirePitchY/tan(angleVrad);

    std::pair<float,float> YwireZminmax = geo.minmax(2, WireCell::kYtype);
    firstYwireZval = YwireZminmax.first + 0.5*wirePitchY;
    leftEdgeOffsetZval = rightEdgeOffsetZval = 0.0*units::cm;
    firstYwireUoffsetYval = firstYwireVoffsetYval = 0.0 * units::cm;

    std::cerr << "maxHeight=" << maxHeight/units::m << " meter " 
	      << "angleUrad=" << angleUrad << " radians, " << angleUrad * units::radian/units::degree << " degree "
	      << "angleVrad=" << angleVrad << " radians, " << angleVrad * units::radian/units::degree << " degree "
	      << "wirePitchY=" << wirePitchY << " "
	      <<std::endl;

    angleUrad = geo.pitch(WireCell::kUtype)/units::radian;
    angleVrad = geo.pitch(WireCell::kVtype)/units::radian;

    UspacingOnWire = std::abs(wirePitchU/sin(angleUrad));
    VspacingOnWire = std::abs(wirePitchV/sin(angleVrad));

    std::cerr << "Tiling..." << std::endl;
    this->constructCells();
}

TileMaker::~TileMaker()
{
}


WireSelection TileMaker::wires(const WireCell::Cell& cell) const
{
    return WireSelection();
}

CellSelection TileMaker::cells(const WireCell::Wire& wire) const
{
    return CellSelection();
}

WireCell::Cell* TileMaker::cell(const WireCell::WireSelection& wires) const
{
}


bool TileMaker::formsCell(double UwireYval, double VwireYval)
{
    bool isCell = false;

    const double tanU = tan(angleUrad);
    const double tanV = tan(angleVrad);

    double deltaY = 0;
    if (UwireYval > VwireYval) {
	deltaY = (UwireYval-UspacingOnWire/2.0)-(VwireYval+VspacingOnWire/2.0);
    }
    else {
	deltaY = (VwireYval-VspacingOnWire/2.0)-(UwireYval+UspacingOnWire/2.0);
    }

    if ((UwireYval+UspacingOnWire/2.0 < -1.0*epsilon) && (VwireYval+VspacingOnWire/2.0 < -1.0*epsilon)) {
	isCell = false;
    }
    else if ((UwireYval-UspacingOnWire/2.0 > maxHeight) && (VwireYval-VspacingOnWire/2.0 > maxHeight+epsilon)) {
	isCell = false;
    }
    else if (deltaY < epsilon) {
	isCell = true;
    }
    else if (UwireYval == VwireYval) {
	isCell = true;
    }
    else if (fabs((tanU*tanV*deltaY)/(tanU+tanV))-wirePitchY/2.0 < epsilon) {
	isCell = true;
    }
    else {
	isCell = false;
    }
    return isCell;
}


static std::pair<double,double> getIntersectionUV(double ZvalOffset, double slope1, double intercept1, double slope2, double intercept2)
{
    double Yval = (slope2*intercept1-slope1*intercept2)/(slope2-slope1);
    double Zval = ZvalOffset + (Yval-intercept1)/slope1;
    
    std::pair<double,double> intersectionPoint;
    intersectionPoint.first = Zval;
    intersectionPoint.second = Yval;
    
    return intersectionPoint;
}
static std::pair<double, double> getIntersectionY(double ZvalOffset, double YwireZval, double slope, double intercept)
{
    double Yval = slope*(YwireZval-ZvalOffset) + intercept;
    double Zval = YwireZval;

    std::pair<double,double> intersectionPoint;
    intersectionPoint.first = Zval;
    intersectionPoint.second = Yval;

    return intersectionPoint;
}

static std::pair<double,double> calcCellCenter(std::vector<std::pair<double,double> > vertices)
{
    double Zval = 0.0;
    double Yval = 0.0;

    const int numVertices = vertices.size();

    for (int ind = 0; ind < numVertices; ++ind) {
	Zval += vertices.at(ind).first;
	Yval += vertices.at(ind).second;
    }  

    Zval /= ((double) numVertices);
    Yval /= ((double) numVertices);

    std::pair<double,double> centerPoint;
    centerPoint.first = Zval;
    centerPoint.second = Yval;

    return centerPoint;
}


static bool compareOrientedVertices(std::pair<double,std::pair<double,double> > orientedVertex1,
				    std::pair<double,std::pair<double,double> > orientedVertex2)
{
    return (orientedVertex1.first > orientedVertex2.first);
}


static std::vector<std::pair<double,double> > sortVertices(std::vector<std::pair<double,double> > vertices)
{
    std::vector<std::pair<double,double> > verticesSorted;

    if (vertices.size() < 2) {
	return vertices;
    }

    std::pair<double,double> vertex, center = calcCellCenter(vertices);
    std::pair<double,std::pair<double,double> > orientedVertex;
    std::vector<std::pair<double,std::pair<double,double> > > orientedVertices;
    for (int ind = 0; ind < vertices.size(); ++ind) {
	vertex = vertices.at(ind);
	orientedVertex.first = atan2(vertex.second-center.second,vertex.first-center.first);
	orientedVertex.second = vertex;
	orientedVertices.push_back(orientedVertex);
    }

    sort(orientedVertices.begin(),orientedVertices.end(),compareOrientedVertices);

    for (int ind = 0; ind < orientedVertices.size(); ++ind) {
	verticesSorted.push_back(orientedVertices.at(ind).second);
    }

    return verticesSorted;
}

static bool 
vertexOutsideBoundary(std::pair<double,double> vertex, int edgeType, double edgeVal)
{
    switch(edgeType) {
    case 1:
	return (vertex.first < edgeVal-epsilon);
    case 2:
	return (vertex.first > edgeVal+epsilon);
    case 3:
	return (vertex.second < edgeVal-epsilon);
    case 4:
	return (vertex.second > edgeVal+epsilon);
    default:
	return false;
    }
}

static std::vector<std::pair<double,double> > 
trimEdgeCellVertices(std::vector<std::pair<double,double> > vertices, 
		     int edgeType, double edgeVal)
{
    std::vector<std::pair<double,double> > edgeCellVertices;

    if((edgeType < 1) && (edgeType > 4)) {
	return vertices;
    }
  
    int numVertOutside = 0;
    for (int ind = 0; ind < vertices.size(); ++ind) {
	if (vertexOutsideBoundary(vertices.at(ind),edgeType,edgeVal) == true) {
	    numVertOutside++;
	}
    }

    if (numVertOutside == 0) {
	return vertices;
    }


    std::pair<double,double> vertex;
    for(int ind = 0; ind < vertices.size(); ++ind) {
	int otherind = (ind+1) % vertices.size();
    
	double Zval1 = vertices.at(ind).first;
	double Zval2 = vertices.at(otherind).first;
	double Yval1 = vertices.at(ind).second;
	double Yval2 = vertices.at(otherind).second;

	if(Zval1 != Zval2) {
	    double slope = (Yval2-Yval1)/(Zval2-Zval1);
	    double intercept = Yval1-slope*Zval1;
	    
	    if ((edgeType == 1) || (edgeType == 2)) {
		if (((edgeVal > Zval1+epsilon) && (edgeVal < Zval2-epsilon)) ||
		    ((edgeVal > Zval2+epsilon) && (edgeVal < Zval1-epsilon))) {
		    vertex.first = edgeVal;
		    vertex.second = slope*edgeVal+intercept;
		    edgeCellVertices.push_back(vertex);
		}
	    }
	    else if((edgeType == 3) || (edgeType == 4)) {
		if(((edgeVal > Yval1+epsilon) && (edgeVal < Yval2-epsilon)) || 
		   ((edgeVal > Yval2+epsilon) && (edgeVal < Yval1-epsilon))) {
		    vertex.first = (edgeVal-intercept)/slope;
		    vertex.second = edgeVal;
		    edgeCellVertices.push_back(vertex);
		}
	    }
	}
	else if((edgeType == 3) || (edgeType == 4)) {
	    if(((edgeVal > Yval1+epsilon) && (edgeVal < Yval2-epsilon)) || 
	       ((edgeVal > Yval2+epsilon) && (edgeVal < Yval1-epsilon))) {
		vertex.first = Zval1;
		vertex.second = edgeVal;
		edgeCellVertices.push_back(vertex);
	    }
	}
    } // loop over vertices
    
    for (int ind = 0; ind < vertices.size(); ++ind) {
	if(vertexOutsideBoundary(vertices.at(ind),edgeType,edgeVal) == false) {
	    edgeCellVertices.push_back(vertices.at(ind));
	}
    }

    return sortVertices(edgeCellVertices);
}

std::vector<std::pair<double,double> > TileMaker::getCellVertices(double YwireZval, double UwireYval, double VwireYval)
{
    std::vector<std::pair<double,double> > cellVertices;

    const double U1slope = 1.0/tan(angleUrad);
    const double U2slope = U1slope;
    const double U1intercept = UwireYval - UspacingOnWire/2.0;
    const double U2intercept = UwireYval + UspacingOnWire/2.0;

    const double V1slope = 1.0/tan(angleVrad);
    const double V2slope = V1slope;
    const double V1intercept = VwireYval - VspacingOnWire/2.0;
    const double V2intercept = VwireYval + VspacingOnWire/2.0;

    const double Y1Zval = YwireZval - wirePitchY/2.0;
    const double Y2Zval = YwireZval + wirePitchY/2.0;

    std::pair<double,double> U1V1intersectionPoint = getIntersectionUV(YwireZval,U1slope,U1intercept,V1slope,V1intercept);
    std::pair<double,double> U1V2intersectionPoint = getIntersectionUV(YwireZval,U1slope,U1intercept,V2slope,V2intercept);
    std::pair<double,double> U2V1intersectionPoint = getIntersectionUV(YwireZval,U2slope,U2intercept,V1slope,V1intercept);
    std::pair<double,double> U2V2intersectionPoint = getIntersectionUV(YwireZval,U2slope,U2intercept,V2slope,V2intercept);

    double UVminZval = U1V1intersectionPoint.first;
    double UVmaxZval = U1V1intersectionPoint.first;
    double UVminYval = U1V1intersectionPoint.second;
    double UVmaxYval = U1V1intersectionPoint.second;

    if (U1V2intersectionPoint.first < UVminZval+epsilon) {
	UVminZval = U1V2intersectionPoint.first;
    }
    if (U1V2intersectionPoint.first > UVmaxZval-epsilon) {
	UVmaxZval = U1V2intersectionPoint.first;
    }
    if (U1V2intersectionPoint.second < UVminYval+epsilon) {
	UVminYval = U1V2intersectionPoint.second;
    }
    if (U1V2intersectionPoint.second > UVmaxYval-epsilon) {
	UVmaxYval = U1V2intersectionPoint.second;
    }

    if (U2V1intersectionPoint.first < UVminZval+epsilon) {
	UVminZval = U2V1intersectionPoint.first;
    }
    if (U2V1intersectionPoint.first > UVmaxZval-epsilon) {
	UVmaxZval = U2V1intersectionPoint.first;
    }
    if (U2V1intersectionPoint.second < UVminYval+epsilon) {
	UVminYval = U2V1intersectionPoint.second;
    }
    if (U2V1intersectionPoint.second > UVmaxYval-epsilon) {
	UVmaxYval = U2V1intersectionPoint.second;
    }

    if (U2V2intersectionPoint.first < UVminZval+epsilon) {
	UVminZval = U2V2intersectionPoint.first;
    }
    if (U2V2intersectionPoint.first > UVmaxZval-epsilon) {
	UVmaxZval = U2V2intersectionPoint.first;
    }
    if (U2V2intersectionPoint.second < UVminYval+epsilon) {
	UVminYval = U2V2intersectionPoint.second;
    }
    if (U2V2intersectionPoint.second > UVmaxYval-epsilon) {
	UVmaxYval = U2V2intersectionPoint.second;
    }

    std::pair<double,double> Y1U1intersectionPoint = getIntersectionY(YwireZval,Y1Zval,U1slope,U1intercept);
    std::pair<double,double> Y1U2intersectionPoint = getIntersectionY(YwireZval,Y1Zval,U2slope,U2intercept);
    std::pair<double,double> Y1V1intersectionPoint = getIntersectionY(YwireZval,Y1Zval,V1slope,V1intercept);
    std::pair<double,double> Y1V2intersectionPoint = getIntersectionY(YwireZval,Y1Zval,V2slope,V2intercept);
    std::pair<double,double> Y2U1intersectionPoint = getIntersectionY(YwireZval,Y2Zval,U1slope,U1intercept);
    std::pair<double,double> Y2U2intersectionPoint = getIntersectionY(YwireZval,Y2Zval,U2slope,U2intercept);
    std::pair<double,double> Y2V1intersectionPoint = getIntersectionY(YwireZval,Y2Zval,V1slope,V1intercept);
    std::pair<double,double> Y2V2intersectionPoint = getIntersectionY(YwireZval,Y2Zval,V2slope,V2intercept);

    if((Y1U1intersectionPoint.second < UVmaxYval+epsilon) && (Y1U1intersectionPoint.second > UVminYval+epsilon) && (Y1U1intersectionPoint.first < UVmaxZval-epsilon) && (Y1U1intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y1U1intersectionPoint);
    }
    if((Y1U2intersectionPoint.second < UVmaxYval-epsilon) && (Y1U2intersectionPoint.second > UVminYval+epsilon) && (Y1U2intersectionPoint.first < UVmaxZval-epsilon) && (Y1U2intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y1U2intersectionPoint);
    }
    if((Y1V1intersectionPoint.second < UVmaxYval-epsilon) && (Y1V1intersectionPoint.second > UVminYval+epsilon) && (Y1V1intersectionPoint.first < UVmaxZval-epsilon) && (Y1V1intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y1V1intersectionPoint);
    }
    if((Y1V2intersectionPoint.second < UVmaxYval-epsilon) && (Y1V2intersectionPoint.second > UVminYval+epsilon) && (Y1V2intersectionPoint.first < UVmaxZval-epsilon) && (Y1V2intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y1V2intersectionPoint);
    }

    if((Y2U1intersectionPoint.second < UVmaxYval-epsilon) && (Y2U1intersectionPoint.second > UVminYval+epsilon) && (Y2U1intersectionPoint.first < UVmaxZval-epsilon) && (Y2U1intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y2U1intersectionPoint);
    }
    if((Y2U2intersectionPoint.second < UVmaxYval-epsilon) && (Y2U2intersectionPoint.second > UVminYval+epsilon) && (Y2U2intersectionPoint.first < UVmaxZval-epsilon) && (Y2U2intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y2U2intersectionPoint);
    }
    if((Y2V1intersectionPoint.second < UVmaxYval-epsilon) && (Y2V1intersectionPoint.second > UVminYval+epsilon) && (Y2V1intersectionPoint.first < UVmaxZval-epsilon) && (Y2V1intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y2V1intersectionPoint);
    }
    if((Y2V2intersectionPoint.second < UVmaxYval-epsilon) && (Y2V2intersectionPoint.second > UVminYval+epsilon) && (Y2V2intersectionPoint.first < UVmaxZval-epsilon) && (Y2V2intersectionPoint.first > UVminZval+epsilon)) {
	cellVertices.push_back(Y2V2intersectionPoint);
    }

    if((U1V1intersectionPoint.first >= Y1Zval-epsilon) && (U1V1intersectionPoint.first <= Y2Zval+epsilon)) {
	cellVertices.push_back(U1V1intersectionPoint);
    }
    if((U1V2intersectionPoint.first >= Y1Zval-epsilon) && (U1V2intersectionPoint.first <= Y2Zval+epsilon)) {
	cellVertices.push_back(U1V2intersectionPoint);
    }
    if((U2V1intersectionPoint.first >= Y1Zval-epsilon) && (U2V1intersectionPoint.first <= Y2Zval+epsilon)) {
	cellVertices.push_back(U2V1intersectionPoint);
    }
    if((U2V2intersectionPoint.first >= Y1Zval-epsilon) && (U2V2intersectionPoint.first <= Y2Zval+epsilon)) {
	cellVertices.push_back(U2V2intersectionPoint);
    }

    std::vector<std::pair<double,double> > cellVerticesSorted = sortVertices(cellVertices);

    const int numYwires = Ywires.size();

    if(UVminZval < (firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval)+epsilon) {
	cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,1,firstYwireZval-0.5*wirePitchY+leftEdgeOffsetZval);
    }
    else if(UVmaxZval > (firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval)-epsilon) {
	cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,2,(firstYwireZval+(numYwires-0.5)*wirePitchY-rightEdgeOffsetZval));
    }

    if(UVminYval < epsilon) {
	cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,3,0.0);
    }
    else if(UVmaxYval > maxHeight-epsilon) {
	cellVerticesSorted = trimEdgeCellVertices(cellVerticesSorted,4,maxHeight);
    }

    return cellVerticesSorted;
}


static double calcCellArea(std::vector<std::pair<double,double> > vertices)
{ 
    double cellArea = 0.0;
    const int numVertices = vertices.size();

    int prev = numVertices-1;
    for (int ind = 0; ind < numVertices; ++ind) {
	cellArea += (vertices.at(prev).first+vertices.at(ind).first)*(vertices.at(prev).second-vertices.at(ind).second); 
	prev = ind;
    }
    cellArea /= 2.0;

    return cellArea;
}


int TileMaker::getUwireID(double Yval, double Zval)
{
    return round((Zval/tan(angleUrad) + maxHeight - (firstYwireZval/tan(angleUrad)) - firstYwireUoffsetYval - Yval)/UspacingOnWire);
}
int TileMaker::getVwireID(double Yval, double Zval)
{
    return round((Zval/tan(angleVrad) - firstYwireZval/tan(angleVrad) - firstYwireVoffsetYval + Yval)/VspacingOnWire);
}

int TileMaker::getYwireID(double Zval)
{
    return round((Zval-firstYwireZval)/wirePitchY);
}


void TileMaker::constructCell(double YwireZval, double UwireYval, double VwireYval)
{
    std::vector<std::pair<double,double> > vertices = getCellVertices(YwireZval,UwireYval,VwireYval);
    if(vertices.size() > 2) {
	return;
    }

    WireCell::Cell cell;
    cell.ident = cellset.size();
    for (size_t ind=0; ind<vertices.size(); ++ind) {
	std::pair<double,double> v = vertices[ind];
	cell.boundary.push_back(Point(0, v.second, v.first));
    }
    std::pair<double,double> center = calcCellCenter(vertices);
    cell.center = Point(0, center.second, center.first);

    cell.area = calcCellArea(vertices);
    cellset.push_back(cell);
    
    WireSelection ws;
    ws.push_back(Uwires[getUwireID(UwireYval,YwireZval)]);
    ws.push_back(Vwires[getVwireID(VwireYval,YwireZval)]);
    ws.push_back(Ywires[getYwireID(YwireZval)]);
    cellmap[&cellset.back()] = ws;

    // fixme: need to fill up the wcmap
    //  cell.UwireID = getUwireID(UwireYval,YwireZval);
    //  cell.VwireID = getVwireID(VwireYval,YwireZval);
    //  cell.YwireID = getYwireID(YwireZval);

    //cell.hitType = kNoHit;
}



//vector<Cell> 
void TileMaker::constructCellChain(double wireZval, double YvalOffsetU, double YvalOffsetV)
{ 
    int numUcrosses = std::ceil(((UdeltaY-UspacingOnWire)/2.0+YvalOffsetU)/UspacingOnWire)+1;
    int numVcrosses = std::ceil((maxHeight-(VdeltaY+VspacingOnWire)/2.0-YvalOffsetV)/VspacingOnWire)+1;


    for (int indU = 0; indU < numUcrosses; indU++) {
	bool flag1 = false, flag2 = false;
	
	for (int indV=0; indV < numVcrosses && !flag2; ++indV) {

	    if (formsCell(YvalOffsetU-indU*UspacingOnWire,YvalOffsetV+indV*VspacingOnWire)) {
		flag1 = true;
		constructCell(wireZval,YvalOffsetU-indU*UspacingOnWire,YvalOffsetV+indV*VspacingOnWire);
	    }
	    else if (flag1 == true) {
		flag2 = true;
	    }
	}
    }
 
    //return cellChain;
}


//vector<Cell> 
void TileMaker::constructCells()
{
    double Zval = firstYwireZval;
    double Uoffset = maxHeight-firstYwireUoffsetYval;
    double Voffset = firstYwireVoffsetYval;

    while(Uoffset < maxHeight-((UspacingOnWire-UdeltaY)/2.0)-epsilon) {
	std::cerr << Uoffset << " " << maxHeight-((UspacingOnWire-UdeltaY)/2.0)-epsilon << " " << UspacingOnWire << std::endl;
	Uoffset += UspacingOnWire;
    }
    while(Voffset > ((VspacingOnWire+VdeltaY)/2.0)+epsilon) {
	Voffset -= VspacingOnWire;
    }

    const int numYwires = Ywires.size();
    for (int ind = 0; ind < numYwires; ++ind) {
	std::cerr << "Constructing cell chain " << ind << " " << Zval << " " << Uoffset << " " << Voffset << std::endl; 
	constructCellChain(Zval,Uoffset,Voffset); 
	Zval += wirePitchY;
	Uoffset += UdeltaY;
	Voffset += VdeltaY;

	while(Uoffset < maxHeight-((UspacingOnWire-UdeltaY)/2.0)-epsilon) {
	    Uoffset += UspacingOnWire;
	}
	while(Voffset < ((VdeltaY-VspacingOnWire)/2.0)-epsilon) {
	    Voffset += VspacingOnWire;
	}
    }

    std::cerr << "Filling wire-cell mesh" << std::endl;

    CellMap::iterator it, done = cellmap.end();
    for (it=cellmap.begin(); it != done; ++it) {
	const Cell* cell = it->first;
	const WireSelection& wires = it->second;
	for (size_t ind=0; ind < wires.size(); ++ind) {
	    const Wire* wire = wires[ind];
	    wiremap[wire].push_back(cell);
	}
    }
}


