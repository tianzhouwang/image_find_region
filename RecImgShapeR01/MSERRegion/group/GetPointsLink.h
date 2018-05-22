 
#ifndef GetPoint_Link_H
#define GetPoint_Link_H
#include "c24bitmap.h"
#include "c256bitmap.h"
//#include "region.h"
//#include "Point2D.h"
#include <math.h>

//#include "hzmerge.h"
#include "mregion.h"


#define VAL_IN_RANGE(a,b,c)  ((a) >= (b) && (a) <= (c)) 

double getRegionCompare(Region&A,Region&B)
{
	if((ColorDist(A.RR,A.GG,A.BB,B.RR,B.GG,B.BB) > 7500) )
		return -1;
	
	//float H1 = max(A.HH, B.HH);
	//float H2 = min(A.HH, B.HH);
	//if( (H1 - H2 > 60 &&  (H2 + 360 - H1) > 60 && (A.SS * A.VV >= 5) && (B.SS * B.VV >= 5) && max(A.VV, B.VV) > 50) /*|| (fabs(A.SS - B.SS) > 0.5 && max(A.VV, B.VV) > 50)*/ || fabs(A.VV - B.VV) > 50)
	//		return -1;		
	/*用于相似区域判断的Filter
	1、KMeans参数
	2、宽高比
	3、面积比
	4、距离与宽高的比
	*/
	//if(A.Mid.disTo(B.Mid)>ConstValue.EPS)
	//    return -1;
	double MinW = min(A.rwidth,B.rwidth);
	double MinH = min(A.rheight,B.rheight);
	double DisT = min(MinW, MinH);
	//double Value = DisT/sqrt((A.x-B.x)*(A.x-B.x) + (A.y-B.y)*(A.y-B.y));
	double Value = 0.25;		
	//if(Value<0.25)  return -1;
	
	double RadioW = min(double(A.rwidth)  / double(B.rwidth) , double(B.rwidth ) / double(A.rwidth) );
	double RadioH = min(double(A.rheight) / double(B.rheight), double(B.rheight) / double(A.rheight));
	double Radio  = max(RadioW, RadioH);
	if(Radio <0.35)
		return -1;
	
	//double AreaR = min(float(A.PtVec.size())/float(B.PtVec.size()),float(B.PtVec.size())/float(A.PtVec.size()));
	//if(AreaR < 0.25)
	//  return -1;
	
	if((A.Feature.StrokeMean==0)||
		(B.Feature.StrokeMean==0))
		return -1;
	
	double StkR = min(float(A.Feature.StrokeMean)/float(B.Feature.StrokeMean),
		float(B.Feature.StrokeMean)/float(A.Feature.StrokeMean));
	if(StkR < 0.6) return -1;
	// return max(Value,max(Radio,AreaR));
	return 1;
}



void GetRotBoundBox(Region&InRegion , double RotCenterX, double RotCenterY, 
					double Dx, double Dy, ImgRECT&BoundRect,int ContourStep=3)
{
	int i;
	
	BoundRect.left  = BoundRect.top    =  9999;
	BoundRect.right = BoundRect.bottom = -9999;
	
	for(i=0;i<InRegion.ContourPtVec.size(); i+=ContourStep)
	{  
		double Dv1x =  InRegion.ContourPtVec[i].x - RotCenterX;
		double Dv1y =  InRegion.ContourPtVec[i].y - RotCenterY;
		double newX = Dv1x *   Dx  + Dv1y * Dy;
		double newY = Dv1x * (-Dy) + Dv1y * Dx;
		double xx_ = RotCenterX + newX;
		double yy_ = RotCenterY + newY;	
		if(xx_ < BoundRect.left  ) BoundRect.left  = xx_;
		if(xx_ > BoundRect.right ) BoundRect.right = xx_;
		if(yy_ < BoundRect.top   ) BoundRect.top   = yy_;
		if(yy_ > BoundRect.bottom) BoundRect.bottom= yy_;
	}
	
}


void GetRotateImageRECTVec(double RotCenterX, double RotCenterY, double Dx, double Dy,
						   vector<Region> &InputRegion, vector<ImgRECT>& OutPutRectVec)
{
	int i;
	Loopi(InputRegion.size())
	{ ImgRECT BoundRect;
	  GetRotBoundBox(InputRegion[i], RotCenterX,  RotCenterY, Dx,  Dy, BoundRect);
	  OutPutRectVec.push_back(BoundRect);
	}
}  


void RotateRegion(Region&OriginPt,double Dx, double Dy, Region&InRegion, Region&OutRegion)
{
	int i,j;
	OutRegion = InRegion;
	double Dv1x,Dv1y, newX,newY;
	Loopi(InRegion.PtVec.size())
	{
		
		
		Dv1x = InRegion.PtVec[i].x - OriginPt.GeoX;
		Dv1y = InRegion.PtVec[i].y - OriginPt.GeoY;
		newX = Dv1x *   Dx  + Dv1y * Dy;
		newY = Dv1x * (-Dy) + Dv1y * Dx;
		OutRegion.PtVec[i].x = OriginPt.GeoX + newX;
		OutRegion.PtVec[i].y = OriginPt.GeoY + newY;
	}
	
	Loopi(InRegion.ContourPtVec.size())
	{
		
		Dv1x =  InRegion.ContourPtVec[i].x - OriginPt.GeoX;
		Dv1y =  InRegion.ContourPtVec[i].y - OriginPt.GeoY;
		newX = Dv1x *   Dx  + Dv1y * Dy;
		newY = Dv1x * (-Dy) + Dv1y * Dx;
		OutRegion.ContourPtVec[i].x = OriginPt.GeoX + newX;
		OutRegion.ContourPtVec[i].y = OriginPt.GeoY + newY;	   
	}
	
	Loopi(InRegion.MContours.size())
	{
		Loopj(InRegion.MContours[i].size())
		{
			Dv1x =  InRegion.MContours[i][j].x - OriginPt.GeoX;
			Dv1y =  InRegion.MContours[i][j].y - OriginPt.GeoY;
			newX = Dv1x *   Dx  + Dv1y * Dy;
			newY = Dv1x * (-Dy) + Dv1y * Dx;
			OutRegion.MContours[i][j].x = OriginPt.GeoX + newX;
			OutRegion.MContours[i][j].y = OriginPt.GeoY + newY;	
		}
		
	}
	
}

void GetOutterRect(vector<Region>& RgVec, ImgRECT & Rect_)
{
    int i;
	
	if(RgVec.size()>0)
	{
		Rect_.left   = RgVec[0].left  ;
		Rect_.right  = RgVec[0].right ;
		Rect_.top    = RgVec[0].top   ;
		Rect_.bottom = RgVec[0].bottom;
	}
    
	Loopi(RgVec.size())
	{
		if(RgVec[i].left  < Rect_.left  ) Rect_.left   = RgVec[i].left  ;
		if(RgVec[i].right > Rect_.right ) Rect_.right  = RgVec[i].right ;
		if(RgVec[i].top   < Rect_.top   ) Rect_.top    = RgVec[i].top   ;
		if(RgVec[i].bottom> Rect_.bottom) Rect_.bottom = RgVec[i].bottom;
	}
	
}

struct LkPointPair 
{
	Point2D Pt1, Pt2;
	Region*R1,*R2;
	int Idx1, Idx2;
	int PointNum;
	double MeanPtDis;
	vector<int> IdxVec;
};

bool operator<(const LkPointPair&P1,const LkPointPair&P2)
{
	//return Pt>y.fitness;
	if(P1.PointNum ==P2.PointNum)
		return P1.MeanPtDis < P2.MeanPtDis;
	
	
	return P1.PointNum > P2.PointNum;
}


double Magnitude( Point2D *Point1, Point2D *Point2 )
{
	Point2D Vector;
    Vector.x = Point2->x - Point1->x;
    Vector.y = Point2->y - Point1->y;
	return (float) sqrt( Vector.x * Vector.x + Vector.y * Vector.y );
}

double Magnitude( Region *Point1, Region *Point2 )
{
	Region Vector;
    Vector.x = Point2->x - Point1->x;
    Vector.y = Point2->y - Point1->y;
	return (float) sqrt( Vector.x * Vector.x + Vector.y * Vector.y );
}

//点到直线距离
int DistancePointLine( Point2D *Point, Point2D *LineStart, Point2D *LineEnd, double&Distance,
					  Point2D &Intersection)
{
	double LineMag;
	double U;
	LineMag = Magnitude( LineEnd, LineStart );
		  
	U = (  
		( Point->x - LineStart->x ) * ( LineEnd->x - LineStart->x )  +
		( Point->y - LineStart->y ) * ( LineEnd->y - LineStart->y ) 
		) /   ( LineMag * LineMag );
		  
		  
		  if( U < 0.0f || U > 1.0f )
			  return 0;   // closest point does not fall within the line segment
		  
		  Intersection.x = LineStart->x + U * ( LineEnd->x - LineStart->x );
		  Intersection.y = LineStart->y + U * ( LineEnd->y - LineStart->y );
		  
		  
		  Distance = Magnitude( Point, &Intersection );
		  
		  return 1;

}

//区域到直线距离
int RegionToLineDistance( Region *Point, Region *LineStart, Region *LineEnd, double&Distance,
						 Region &Intersection)
{
	double LineMag;
	double U;
	LineMag = Magnitude( LineEnd, LineStart );
		  
	U = (  
		( Point->x - LineStart->x ) * ( LineEnd->x - LineStart->x )  +
		( Point->y - LineStart->y ) * ( LineEnd->y - LineStart->y ) 
		) /   ( LineMag * LineMag );
		  
		  
		  if( U < 0.00f || U > 1.01f )
			  return 0;   // closest point does not fall within the line segment
		  
		  Intersection.x = LineStart->x + U * ( LineEnd->x - LineStart->x );
		  Intersection.y = LineStart->y + U * ( LineEnd->y - LineStart->y );
		  
		  
		  Distance = Magnitude( Point, &Intersection );
		  
		  return 1;
}

int GetPointNum(Point2D&Pt1, Point2D&Pt2,vector<Point2D> &PtVec,double&DisAverage, double distrange=40)
{
	int i;
	double dissum;
	dissum = 0;
	int PtCount;
	PtCount = 0;
	for (i=0;i<PtVec.size();i++)
	{
		double  Distance;
        Point2D Intersection;
		
		
		
		int in_range =  DistancePointLine( &PtVec[i], &Pt1, &Pt2, Distance, Intersection);
        if(!in_range)    continue;
		if(Distance>distrange)  continue;
        PtCount++;
		dissum += Distance;
	}
	
	if(PtCount==0) return 0;
	
	DisAverage = dissum/double(PtCount);
	
	return PtCount;
} 

//----------------------------------------------------------------------------------
class FastSelCoLinearRegion
{
public:
	vector< vector<double> > AngleVec;
    vector< vector<int>    > RegionSimilar;
	vector< vector<int>    > RegionLink;
	void GenRegionSimilarInfo(vector<Region>&RegionVec);
	void GenRegionNeighbourInfo(vector<Region>&RegionVec, vector< vector<int> >&RegionLink);
	void MakeNearlink(vector<Region>&RegionVec);
	int  PicWidth, PicHeight;
	void SetWidthheight(int PicWidth,int PicHeight);
};

void FastSelCoLinearRegion::SetWidthheight(int PicWidth_,int PicHeight_)
{
  PicWidth  = PicWidth_;
  PicHeight = PicHeight_;
}

int SkipInvalidRegion(Region&R)
{
	int max_size = max(R.rwidth, R.rheight);
	if( max_size/ R.Feature.StrokeMean > 10.5 && R.Feature.roundratio >0.99)
		return 1;
	return 0;
}

int RegionXYOverLap(Region & R1,Region & R2)
{
  int size1 = max(R1.rwidth,R1.rheight);
  int size2 = max(R2.rwidth,R2.rheight);
  
  if(size1/size2>3)return 0;
  if(size2/size1>3)return 0;

  int Val1 = VAL_IN_RANGE(R1.top,  R2.top - R2.Feature.StrokeMean, R2.bottom + R2.Feature.StrokeMean);  
  int Val2 = VAL_IN_RANGE(R2.top,  R1.top - R1.Feature.StrokeMean, R1.bottom + R1.Feature.StrokeMean);  
  int Val3 = VAL_IN_RANGE(R1.left, R2.left - R2.Feature.StrokeMean, R2.right + R2.Feature.StrokeMean);   
  int Val4 = VAL_IN_RANGE(R2.left, R1.left - R1.Feature.StrokeMean, R1.right + R1.Feature.StrokeMean);
  
  if(Val1||Val2||Val3||Val4)
	  return 1;

  return 0;
}


void FastSelCoLinearRegion::GenRegionNeighbourInfo(vector<Region>&RegionVec, vector< vector<int> >&RegionLink)
{
	int i,j,t;

//	int xDir[8],yDir[8];
//	xDir[0] =1; xDir[1] =1; xDir[2] =0; xDir[3] = -1; xDir[4] = -1; xDir[5] =-1; xDir[6] = 0; xDir[7] = 1;
//  yDir[0] =0; yDir[1] =1; yDir[2] =1; yDir[3] =  1; yDir[4] =  0; yDir[5] =-1; yDir[6] =-1; yDir[7] =-1;

	int xDir[12],yDir[12];
	xDir[0] =1; xDir[1] =1; xDir[2] =0; xDir[3] = -1; xDir[4] = -1; xDir[5] =-1; xDir[6] = 0; xDir[7] = 1;
    yDir[0] =0; yDir[1] =1; yDir[2] =1; yDir[3] =  1; yDir[4] =  0; yDir[5] =-1; yDir[6] =-1; yDir[7] =-1;

	xDir[8] = 2; xDir[9] =-2; xDir[10] = -2; xDir[11] =  2;
	yDir[8] = 1; yDir[9] = 1; yDir[10] = -1; yDir[11] = -1;

	double searchLenNorm[12];
    
	for(i=0 ; i < 12 ; i++)
	{
		double dx = xDir[i];
		double dy = yDir[i];
		searchLenNorm[i] = sqrt(dx * dx + dy * dy);
	}

	if(RegionLink.size()!= RegionVec.size())
	{
		RegionLink.resize(RegionVec.size());
		for(i=0;i< RegionVec.size();i++)
		{
			RegionLink[i].resize(RegionVec.size());
		}
	}

	vector<int> LabelVec;
	LabelVec.resize(PicWidth*PicHeight);

	Loopi(PicWidth*PicHeight)
		 LabelVec[i]= -1;
	
	for (i=0;i<RegionVec.size();i++)
	{
		for(j=0;j<RegionVec[i].PtVec.size();j++)
		{
			int x,y;
			x = RegionVec[i].PtVec[j].x;
			y = RegionVec[i].PtVec[j].y;
			LabelVec[y*PicWidth+x] = i;
		}
	}
	
	for(i=0; i<RegionVec.size(); i++)
	{
        //if(SkipInvalidRegion(RegionVec[i]))continue;
		//if(i==86)
		//{
		//	int debug_here;debug_here =1;
		//}
		int  max_size      = max(RegionVec[i].rwidth, RegionVec[i].rheight);
		int  stroke_length = RegionVec[i].Feature.StrokeMean * 8;
		int  size_length   = max_size  * 3/2; 
		int  search_length = size_length;//min(stroke_length, size_length); 
        

		for (j=0;j<12;j++)
		{
          int region_seach_len;
		  region_seach_len = search_length / searchLenNorm[j] ;//(xDir[j]*yDir[j]==0) ? search_length : (search_length* 3/2); 
		  int startX,startY;		  
		   
          startX = RegionVec[i].GeoX + RegionVec[i].rwidth /3 * xDir[j];
		  startY = RegionVec[i].GeoY + RegionVec[i].rheight/3 * yDir[j];
		 
          
		  for(t=0;t<region_seach_len;t++)
		  {

            int newX,newY;
			newX = startX + t * xDir[j];
			newY = startY + t * yDir[j];

			if(newY > (PicHeight-1)|| newX >(PicWidth-1) || newX <0 || newY <0 )
				break;

			int RgLabel;
			RgLabel = LabelVec[newY * PicWidth +newX];
			
			if(RgLabel==-1) continue;


			if(RgLabel!= i && RegionSimilar[i][RgLabel])
			{
                if(!SkipInvalidRegion(RegionVec[RgLabel])) 
				{
					if( RegionXYOverLap(RegionVec[i],RegionVec[RgLabel]))
					    RegionLink[i][RgLabel] = 1;
				}
				break;
			}
		  }
		}
	    
	}
	
}


void FastSelCoLinearRegion::GenRegionSimilarInfo(vector<Region>&RegionVec)
{
	int i,j;
	int RegionVecLength = RegionVec.size();
	
	RegionSimilar.resize(RegionVecLength);
	AngleVec.resize(RegionVecLength);
	
	Loopi(RegionVecLength)
	{
		RegionSimilar[i].resize(RegionVecLength);
		AngleVec[i].resize(RegionVecLength);
	}
	
	
	for(i=0;i<RegionVecLength;i++)
		for(j=i+1;j<RegionVecLength;j++)
		{
		   // if(i==86 && j==87)
		//	{int stop_here;stop_here =1;}//TEMP_DEBUG
			if( SkipNotSimilarRegion(RegionVec[i], RegionVec[j]))
			{
				RegionSimilar[i][j]=0;
				RegionSimilar[j][i]=0;
			}
			else
			{
				RegionSimilar[i][j]=1;
				RegionSimilar[j][i]=1;
			}
		}
		
		for(i=0;i<RegionVecLength;i++)
		{	 for(j=i+1;j<RegionVecLength;j++)
		{
			double Dx = RegionVec[i].GeoX - RegionVec[j].GeoX;
			double Dy = RegionVec[i].GeoY - RegionVec[j].GeoY;
			double angle;
			angle = atan2(Dx,Dy)  *180.0/3.1415926;  AngleVec[i][j] = angle;
			angle = atan2(-Dx,-Dy)*180.0/3.1415926;  AngleVec[j][i] = angle;
		}
		}
		
}
//----------------------------------------------------------------------------------------

int GetInSideRegion(FastSelCoLinearRegion & RegionCoAngleInfo,int Pt1Idx,int Pt2Idx,vector<Region> &RgVec,vector<int> &IdxVec,
					double&DisAverage, double distrange=40)
{
	int i;
	double dissum;
	dissum = 0;
	int PtCount;
	PtCount = 0;
	IdxVec.clear();
	for (i=0;i<RgVec.size();i++)
	{ 
		double  Distance;
        Region Intersection;
		Region DebugR = RgVec[i];
		
		if(getRegionCompare(RgVec[i], RgVec[Pt1Idx])==-1) continue;
        if(getRegionCompare(RgVec[i], RgVec[Pt2Idx])==-1) continue;
		
		if(i==Pt1Idx||i==Pt2Idx)
		{
			IdxVec.push_back(i);
		}
		else
		{
			double angle12 = RegionCoAngleInfo.AngleVec[Pt1Idx][Pt2Idx];
			double anglePt = RegionCoAngleInfo.AngleVec[Pt1Idx][i];
			
			double diff1 = fabs(anglePt - angle12);
			double diff2 = fabs(diff1 -180.0) ;
			
			int SkipFlag ;
			
			
			if( (diff1 < 15.0) || (diff2<15.0)) 
				SkipFlag = 0;
			else
				SkipFlag = 1;
			
			if(SkipFlag) continue;
			
			int in_range =  RegionToLineDistance( &RgVec[i], &RgVec[Pt1Idx], &RgVec[Pt2Idx], Distance, Intersection);
			if(!in_range)    continue;
			if(Distance>distrange)  continue;
			PtCount++;
			dissum += Distance;
			IdxVec.push_back(i);
		}
		
	}
	
	if(PtCount==0) return 0;
	
	DisAverage = dissum/double(PtCount);
	
	return PtCount;
} 


void SelectRegionBetweenALine(FastSelCoLinearRegion & RegionCoAngleInfo ,
							  int Pt1Idx,int Pt2Idx, vector<Region>&InRgVec,
							  vector<int> & IdxList,
							  //vector<Region>&OutRgVec,
							  int distance =40)
{
	int i,j;
	vector<int> IdxVec;
	double DisAverage;
	IdxList.clear();
	
	GetInSideRegion(RegionCoAngleInfo, Pt1Idx ,  Pt2Idx , InRgVec, IdxVec, DisAverage, distance);
	vector<Region> ObjVec;
	
	Loopi(IdxVec.size()) //可以优化,只查边界便可以
		ObjVec.push_back(InRgVec[IdxVec[i]]);
	
	ImgRECT  RectBound;
    GetOutterRect(ObjVec,  RectBound);
	
	//OutRgVec.clear();
	
	Loopi(InRgVec.size())
	{
		if(RegionInRect(InRgVec[i], RectBound))
		{
			if(ColorDist(InRgVec[i].RR,InRgVec[i].GG,InRgVec[i].BB,
				InRgVec[Pt1Idx].RR,     InRgVec[Pt1Idx].GG,     InRgVec[Pt1Idx].BB)< (40*40*3))
			{	//OutRgVec.push_back(InRgVec[i]);
				IdxList.push_back(i);
			}
		}
	} 
	
}
void GetMinPointLink(vector<Point2D>&PtVec)
{ 
	int i,j;
	
	
	LkPointPair P;
	P.Pt1  = PtVec[0];
	P.Pt2  = PtVec[1];
	P.Idx1 = i;
	P.Idx2 = j;
	double DisAverage; DisAverage =999;
	P.PointNum = GetPointNum(PtVec[0], PtVec[11], PtVec, DisAverage);
    P.MeanPtDis = DisAverage;
	
	
	vector <LkPointPair> PtPairVec;
	for(i=0;i<PtVec.size();i++)
		for (j=i+1;j<PtVec.size();j++)
		{
			
			LkPointPair P;
			P.Pt1  = PtVec[i];
			P.Pt2  = PtVec[j];
			P.Idx1 = i;
			P.Idx2 = j;
			double DisAverage; DisAverage =999;
			P.PointNum = GetPointNum(PtVec[i], PtVec[j], PtVec, DisAverage);
			P.MeanPtDis = DisAverage;
			PtPairVec.push_back(P);
		}
		
		sort(PtPairVec.begin(),PtPairVec.end());
		
		FILE* file = fopen("ptlink1.txt","wt+");
		
		for(i=0; i< PtPairVec.size(); i++)
		{
			LkPointPair P = PtPairVec[i];
			fprintf(file,"%i, %i, %i, %.3lf\n",  PtPairVec[i].Idx1,     PtPairVec[i].Idx2,
				PtPairVec[i].PointNum, PtPairVec[i].MeanPtDis);
		}
		
		fclose(file);
}


void GetRegionPointsLink(vector<Region>&RgVec,C24BitMap&CPic)
{ 
	int i,j;
	
	FastSelCoLinearRegion  RegionCoAngleInfo;
	RegionCoAngleInfo.GenRegionSimilarInfo(RgVec);
	
	FILE* file = fopen("regioninfo.txt","wt+");
    for(i=0;i<RgVec.size();i++)
	{
		int xx,yy;
		xx = RgVec[i].GeoX;
		yy = RgVec[i].GeoY;
		CPic.PenColor.R = RgVec[i].RR; 
		CPic.PenColor.G = RgVec[i].GG;
		CPic.PenColor.B = RgVec[i].BB;
		char Buffer[40];
		sprintf(Buffer,"%i",i);
		MergeTxtStr(CPic, xx, yy, 15, Buffer, RgVec[i].RR, RgVec[i].GG, RgVec[i].BB);
		CPic.SigDot(xx - 1, yy   );
		CPic.SigDot(xx  +1, yy   );
		CPic.SigDot(xx    , yy -1);
		CPic.SigDot(xx    , yy +1);
		fprintf(file,"%i,%f,%f --- %f,%f,%f\n", i,RgVec[i].GeoX,RgVec[i].GeoY, RgVec[i].RR,RgVec[i].GG,RgVec[i].BB);
	}
	
	fclose(file);
	CPic.Save("dest.bmp");
	
	LkPointPair P;
	vector<LkPointPair> PtPairVec;
	for(i=0;i<RgVec.size();i++)
		for(j=i+1;j<RgVec.size();j++)
		{
			//Point2D Pt1,Pt2;
			//Pt1.x = RgVec[i].x ; Pt1.y = RgVec[i].y ;
            //Pt2.x = RgVec[j].x ; Pt2.y = RgVec[j].y ;
			int a=0;
			if(i==17&&j==24)
			{
				a = 2000;
			}
			if(getRegionCompare(RgVec[i], RgVec[j])==-1)
				continue;
			double DisAverage;
			P.R1 = &RgVec[i];
			P.R2 = &RgVec[j];
			P.Idx1 = i;
   	        P.Idx2 = j;
			P.PointNum = GetInSideRegion( RegionCoAngleInfo, i , j , RgVec, P.IdxVec, DisAverage);//, double distrange=40);
			P.MeanPtDis = DisAverage;
			PtPairVec.push_back(P);
			//P.Pt1  = PtVec[i];
			//P.Pt2  = PtVec[j];
			//P.Idx1 = i;
			//P.Idx2 = j;
		}
		
		
		
		sort(PtPairVec.begin(),PtPairVec.end());
		
		file = fopen("ptlink1.txt","wt+");
		
		for(i=0; i< PtPairVec.size(); i++)
		{
			LkPointPair P = PtPairVec[i];
			fprintf(file,"%i, %i, %i, %.3lf\n",  PtPairVec[i].Idx1,     PtPairVec[i].Idx2,
				PtPairVec[i].PointNum, PtPairVec[i].MeanPtDis);
			fprintf(file,"########################\n");
			for(j=0;j<P.IdxVec.size();j++)
				fprintf(file,"%i,",PtPairVec[i].IdxVec[j]);
			fprintf(file,"\n--------------------------\n");
		}
		
		fclose(file);
		
		for(i=0;i<100;i++)
		{
			CPic.DrawLine(PtPairVec[i].R1->x,PtPairVec[i].R1->y,
				PtPairVec[i].R2->x,PtPairVec[i].R2->y);
		}
}

/*
typedef struct tagXYZ
{
double X, Y, Z;
}
XYZ;

  double Magnitude( XYZ *Point1, XYZ *Point2 )
  {
  XYZ Vector;
  
	Vector.X = Point2->X - Point1->X;
	Vector.Y = Point2->y - Point1->y;
	Vector.Z = Point2->Z - Point1->Z;
	
	  return (float)sqrt( Vector.X * Vector.X + Vector.Y * Vector.Y + Vector.Z * Vector.Z );
	  }
*/

//最小二乘法直线拟合
int CalculateLineKB(vector<Point2D> &PointVec, double &k, double &b, double &meanx, double&meany)
{
	//最小二乘法直线拟合
	//m_FoldList为关键点(x,y)的列表
	//拟合直线方程(Y=kX+b)
	int i;
	if(PointVec.size()==0) return 0;
	long lCount = PointVec.size();
	if(lCount<2)return FALSE;
	
	double mX,mY,mXX,mXY,n;
	mX=mY=mXX=mXY=0;
	n=lCount;
	
	Loopi(PointVec.size())
	{
		//pFold=m_FoldList->GetNext(pos);
		mX  += PointVec[i].x;
		mY  += PointVec[i].y;
		mXX += PointVec[i].x * PointVec[i].x;
		mXY += PointVec[i].x * PointVec[i].y;
	}
	
	if(mX*mX-mXX*n==0)
		return FALSE;
	
	k = (mY * mX - mXY * n) / (mX * mX - mXX * n);
	b = (mY - mX*k ) / n;
	
	meanx = mX / double( PointVec.size() );
    meany = mY / double( PointVec.size() );
	
	return 1;
}




#endif 
