#ifndef Region_Rect_Operation_Head
#define Region_Rect_Operation_Head
#include <math.h>
//#include "mregion.h"
#include "mktextchain.h"
#include "chisquare.h"

#define IMG_RECT_GEOX(Rect)  ((Rect.left + Rect.right )/2)
#define IMG_RECT_GEOY(Rect)  ((Rect.top  + Rect.bottom)/2)
#define VECT_NROM(X,Y)       sqrt((X)*(X) + (Y)*(Y))
#define VECT_DOT_PRODUCT(X1,Y1,X2,Y2)  ( (X1)*(X2) + (Y1)*(Y2) )
 

/// \brief  得到向量的统计数值
template<typename ADType>
double GetVectorMeanStdev(vector<ADType>&fvec,
						  double&minv,double&maxv,
						  double&mean,double&median,double&stdev)
{
	int i;
	double sumVal,powVal;
	sumVal = powVal = 0;
	if(fvec.size()==0) return 0;
	minv = maxv = fvec[0];
	Loopi(fvec.size())
	{
		ADType val = fvec[i];
		sumVal += val;
		powVal += val*val;
		if(maxv<val)maxv = val;
		if(minv>val)minv = val;
	}
    // nth_element(lengthVec.begin(),lengthVec.begin()+middle,lengthVec.end()); 
	mean   = sumVal/double(fvec.size());
	stdev  = sqrt(powVal /double(fvec.size()) - mean * mean);
	int middle = fvec.size()/2;
	nth_element(fvec.begin(),fvec.begin()+middle,fvec.end());
    median = fvec[middle];

//	VAL_IN_RANGE
	float Count(0);
	Loopi(fvec.size())
	{
		double ratio = double(fvec[i])/ median;
		if(VAL_IN_RANGE(ratio, 0.8 , 1.2))
          Count++;
	}

	return Count/double(fvec.size());
}


/// \brief  外接矩形连接线的角度分数
///         小于10度的角在所有角度中的占比
/// \param   
/// \return  
double GetAngleScore(vector<double>&AngleVec, double &MeanDiff)
{
  int i;
  double Cnt    =0;
  double AngSum =0;
  Loopi(AngleVec.size())
  {
	 if(AngleVec[i]<10.0)
	 {
	   Cnt++;
       AngSum += AngleVec[i];
	 }
  }

  if(Cnt!=0)
    MeanDiff = AngSum /double(Cnt);
  else
    MeanDiff = 60.0;
  return Cnt / double(AngleVec.size());
}


/// \brief  外接矩形连接线的综合分数
/// \param   
/// \return  
double ScoreTxtLine(vector<ImgRECT> & CombineRegion,vector<ImgRECT> & ChrRegion)
{
  sort(CombineRegion.begin(),CombineRegion.end(), LessByX);

  if(ChrRegion.size()*9/10>CombineRegion.size())
	  return 0;
  int i;

  vector<double> DistVec;
  vector<double> AngVec;
  vector<double> HeightVec;

   Loopi(CombineRegion.size())
   {
     HeightVec.push_back(CombineRegion[i].bottom - CombineRegion[i].top + 1) ;
   }


   Loopi(CombineRegion.size()-2)
  {
    double Dx1,Dy1;
	double Dx2,Dy2;

    Dx1 = IMG_RECT_GEOX(CombineRegion[i+1]) - IMG_RECT_GEOX(CombineRegion[i ]);
    Dy1 = IMG_RECT_GEOY(CombineRegion[i+1]) - IMG_RECT_GEOY(CombineRegion[i ]);

    Dx2 = IMG_RECT_GEOX(CombineRegion[i+2]) - IMG_RECT_GEOX(CombineRegion[i+1]);
    Dy2 = IMG_RECT_GEOY(CombineRegion[i+2]) - IMG_RECT_GEOY(CombineRegion[i+1]);
    
	double Len1,Len2,Angle;
	Len1  = VECT_NROM(Dx1,Dy1); Len2 = VECT_NROM(Dx2,Dy2);
    Angle = acos( abs( VECT_DOT_PRODUCT(Dx1,Dy1,Dx2,Dy2) )/(Len1 * Len2) );
    Angle = Angle *180.0 / 3.1415926;
	DistVec.push_back(Len1);
    DistVec.push_back(Len2); 
    AngVec.push_back(Angle);
  }

   double minv, maxv, mean, median, stdev;

   double AngMeanDiff;
   double score1 =  GetAngleScore(AngVec, AngMeanDiff);//GetVectorMeanStdev(AngVec , minv, maxv, mean, median, stdev);
   double score3 = (180.0 - AngMeanDiff) / 180.0;
   double score2 = GetVectorMeanStdev(DistVec, minv, maxv, mean, median, stdev);
   double score4 = GetVectorMeanStdev(HeightVec, minv, maxv, mean, median, stdev); 

   return score1 + score2 + score3 + score4;
  //GetVectorMeanStdev(vector<ADType>&fvec,ADType&minv,ADType&maxv,ADType&mean,ADType&median,ADType&stdev)
}


double GetObjLinkUniformity(vector<ImgRECT>&RgVec,double& MeanDistance)
{
   int i;
   vector<double> ObjDistanceVec;
   Loopi(RgVec.size()-1)
   {
	 double x1,x2,y1,y2;
	 x1 = (RgVec[i].left + RgVec[i].right )/2;
	 y1 = (RgVec[i].top  + RgVec[i].bottom)/2;

	 x2 = (RgVec[i+1].left + RgVec[i+1].right )/2;
	 y2 = (RgVec[i+1].top  + RgVec[i+1].bottom)/2;
     double dist = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)); 
	 ObjDistanceVec.push_back(dist);
	 printf("dis:%.2lf\n,",dist);
   }
   
   int NMid = ObjDistanceVec.size()/2;
   
   nth_element(ObjDistanceVec.begin(), 
	           ObjDistanceVec.begin()+NMid,
			   ObjDistanceVec.end());

   MeanDistance = ObjDistanceVec[NMid];

   double sumdif      = 0;
   double average_dis = 0;
   Loopi(ObjDistanceVec.size())
   {
	   sumdif += (ObjDistanceVec[i] - MeanDistance)*
		         (ObjDistanceVec[i] - MeanDistance)/MeanDistance;

	   average_dis+= ObjDistanceVec[i];
   }
   
   average_dis/=double(ObjDistanceVec.size());
  //sumdif
   double V1 = pochisq(0, ObjDistanceVec.size()-1);
   

     sumdif = 0;

   Loopi(ObjDistanceVec.size())
   {
	   sumdif += (ObjDistanceVec[i] - average_dis)*
		         (ObjDistanceVec[i] - average_dis)/average_dis;
   }

   double V2 = pochisq(sumdif, ObjDistanceVec.size()-1);
   
   return max(V1,V2);
}


double GetObjHeightUniformity(vector<ImgRECT>&RgVec,double& MeanDistance,int refence_Top)
{
	int i;
	vector<double> ObjDistanceVec;
	Loopi(RgVec.size() )
	{   
		double dist = (RgVec[i].top + RgVec[i].bottom)/2 - refence_Top;
		ObjDistanceVec.push_back(dist);
	}
	
	int NMid = ObjDistanceVec.size()/2;
	
	nth_element(ObjDistanceVec.begin(), 
		ObjDistanceVec.begin()+NMid,
		ObjDistanceVec.end());
	
	MeanDistance = ObjDistanceVec[NMid];
	
	double sumdif      = 0;
	double average_dis = 0;
	Loopi(ObjDistanceVec.size())
	{
		sumdif += (ObjDistanceVec[i] - MeanDistance)*
			      (ObjDistanceVec[i] - MeanDistance)/MeanDistance;
		
		average_dis+= ObjDistanceVec[i];
	}
	
	average_dis/=double(ObjDistanceVec.size());
	
	double V1 = pochisq(sumdif, RgVec.size());
	
	
	sumdif = 0;
	
	Loopi(ObjDistanceVec.size())
	{
		sumdif += (ObjDistanceVec[i] - average_dis)*
			(ObjDistanceVec[i] - average_dis)/average_dis;
	}
	double V2 = pochisq(sumdif, RgVec.size());
	
	return max(V1,V2);
}

#endif