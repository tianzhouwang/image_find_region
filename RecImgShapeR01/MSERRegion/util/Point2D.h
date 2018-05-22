#ifndef _POINT2D_H_
#define _POINT2D_H_

//#include <iostream.h>

/**
 * Primitive 2D point class used as input for the LineParamEstimator.
 *
 * Author: Ziv Yaniv (zivy@cs.huji.ac.il)
 */
class Point2D {
public:
	Point2D(double px, double py) : x(px), y(py) {}
	double x;
	double y;
	Point2D(){x=0;y=0;};
	int Selected;
	double R,G,B;

	double NewX,NewY;
	double Weight;
};

/*inline ostream &operator<<(ostream &output,const Point2D &pnt)
{
//	output<<pnt.x<<' '<<pnt.y;
	return(output);
}*/

#endif //_POINT2D_H_