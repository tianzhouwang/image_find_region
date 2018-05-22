#ifndef Split_And_Get_Rect_Rotation_Angle
#define Split_And_Get_Rect_Rotation_Angle
#include <vector>
using namespace std;
//#include "Point2D.h"
//#include "mregion.h"

//将一个点旋转到新的角度
void RotPointToNewPt(double OldX, double OldY,
					 double Dx  , double Dy  ,
					 double&NewX, double&NewY )
{
	NewX = OldX *   Dx  + OldY * Dy;
	NewY = OldX * (-Dy) + OldY * Dx;
}


//点是否在旋转矩形内
int PointInRange(int PtX,int PtY, int RectCentX, int RectCentY,double Dx,double Dy,double &RWidth,double&RHeight)
{
	double X,Y, NewX,NewY;
	X = PtX - RectCentX;
	Y = PtY - RectCentY;
	RotPointToNewPt(X,  Y,   Dx, Dy,  NewX,  NewY ); 

	if( 
		(fabs(NewX*2)<RWidth) &&
		(fabs(NewY*2)<RHeight)
  	 )
	{
		return 1;
	}

	return 0;
}

//得到凸包的投影宽高
void GetRegionProjHeight(Region& R,double Dx,double Dy,double &RWidth,double&RHeight)
{
	int i;
    
	vector<double> XVec,YVec;

	Loopi(R.ConvexHullPtVec.size())
	{ 
		double X,Y, NewX,NewY;
		X = R.ConvexHullPtVec[i].x - R.GeoX;
		Y = R.ConvexHullPtVec[i].y - R.GeoY;

		RotPointToNewPt(X,  Y,   Dx, Dy,  NewX,  NewY );
		
		XVec.push_back( NewX);
		YVec.push_back( NewY);
	}  

	RHeight = ( *max_element(YVec.begin(),YVec.end())  ) - 
		      ( *min_element(YVec.begin(),YVec.end())  ); 

	RWidth =  ( *max_element(XVec.begin(),XVec.end())  ) - 
		      ( *min_element(XVec.begin(),XVec.end())  ); 
}


//得到区域凸包向矩形中心的投影 
void GetRegionConvexHullProj(Region& R,int RectCentX, int RectCentY,double Dx,double Dy,					 
						     double &XMin,double&YMin,
						     double &XMax,double&YMax)
{
	int i;
    
	vector<double> XVec,YVec;
	
	Loopi(R.ConvexHullPtVec.size())
	{ 
		double X,Y, NewX,NewY;
		X = R.ConvexHullPtVec[i].x - RectCentX;
		Y = R.ConvexHullPtVec[i].y - RectCentY;
		
		RotPointToNewPt(X,  Y,   Dx, Dy,  NewX,  NewY );
		
		XVec.push_back( NewX);
		YVec.push_back( NewY);
	}  
	
	YMax = ( *max_element(YVec.begin(),YVec.end())  ) ; 
	YMin = ( *min_element(YVec.begin(),YVec.end())  ) ; 
	
	XMax = ( *max_element(XVec.begin(),XVec.end())  ) ; 
	XMin = ( *min_element(XVec.begin(),XVec.end())  ) ; 
}
//得到一组对象的投影能量
double GetProjEnergy(vector<Point2D>&RectVec,double Angle)
{
	int i;
	Angle =   Angle / 180 * M_PI ;
	double Dx,Dy;
	//Angle =0;
	Dx =  cos(Angle); Dy = - sin(Angle);
	
	vector <double> XaxisVec;
	vector <double> YaxisVec;
	
	
	Loopi(RectVec.size())
	{ 
		RotPointToNewPt(RectVec[i].x, RectVec[i].y,
			    Dx, Dy, RectVec[i].NewX,  RectVec[i].NewY );
		
		XaxisVec.push_back(RectVec[i].NewX);
		YaxisVec.push_back(RectVec[i].NewY);
	}  
	
	double MinX,MaxX;
	double MinY,MaxY;
	
	MinX = *min_element(XaxisVec.begin(),XaxisVec.end());
    MinY = *min_element(YaxisVec.begin(),YaxisVec.end());
	
	MaxX = *max_element(XaxisVec.begin(),XaxisVec.end());
    MaxY = *max_element(YaxisVec.begin(),YaxisVec.end());
	
	vector <int> VoteBin;
    int BinSize = 3;
	VoteBin.resize(BinSize);
	
	double ProjHeight  = MaxY - MinY + 10; //return ProjHeight;
	double BinWidth    = ProjHeight /double(BinSize);
	
	Loopi(RectVec.size())
	{ 
		//XaxisVec.push_back(RectVec[i].NewX);
		//YaxisVec.push_back(RectVec[i].NewY);
		double BinIdx = (RectVec[i].NewY - (MinY-5))/ BinWidth;
        int MaxId = (BinSize-1);
		int Idx = BOUND(BinIdx, 0, MaxId );
		VoteBin[Idx] +=  sqrt(RectVec[i].Weight);
	}  
    
	double BinEntropyScore =0;
	
	Loopi(BinSize)
	{
		if(VoteBin[i]==0) continue;
		
		double Cnt = VoteBin[i];
		//BinEntropyScore += Cnt * log(Cnt);
		BinEntropyScore += Cnt *  (Cnt);
	}
	return BinEntropyScore/(ProjHeight);
}

//得到一直对象外接矩形的最佳投影角度
double GetRotAngle(vector<Point2D> & RectVec, double AngleStep = 5.0)
{
	double MaxEntropy =     0;
	double MaxAngle   = -60.0;
	
	double RotAngle;
	for(RotAngle = -50.0; RotAngle < 50.0; RotAngle+= AngleStep)
	{
		
		double AngleEntropy = GetProjEnergy(RectVec,RotAngle); 

		//printf("angle:%lf, AngleEntropy:%lf\n",RotAngle,AngleEntropy);
		if(MaxEntropy < AngleEntropy )
		{
			MaxEntropy = AngleEntropy;
			MaxAngle   = RotAngle;
		}
	}   
    return MaxAngle;	
}

#endif