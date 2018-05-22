 
#ifndef Region_Select_Head
#define Region_Select_Head
#include "mregion.h"


//���һ��ֵ�Ƿ��ڸ�����ֵ��Χ��
#define valueInRange(value_, min_, max_)  ((value_ >= (min_)) && (value_ <= (max_)))


/// \brief  ������������Ƿ��ཻ
///  
/// \param   
/// \return   
bool checkRectIntersect(int ax1,int ay1, int ax2,int ay2,
						int bx1,int by1, int bx2,int by2)
{   //rect A, rect B
	
    bool xOverlap = valueInRange(ax1, bx1, bx2) || valueInRange(bx1, ax1, ax2);
    bool yOverlap = valueInRange(ay1, by1, by2) || valueInRange(by1, ay1, ay2);	
    return xOverlap && yOverlap;
}

/// \brief  �Ƿ�����ĳ����
/// ����������:  �ʻ���ȱ������,�ʻ��ڲ��Ŀն�������,
///              �ն���ռ�����������
/// \param   
/// \return   
int SkipRegion(Region&R)
{
  if(R.Feature.StrokeStd/R.Feature.StrokeMean>0.7)return 1;
  if(R.Feature.InnerHoles>9)return 1;
  if(R.Feature.roundratio>0.99)return 1;
  return 0;
}


/// \brief  ������������Ƿ�����
///  
/// \param   
/// \return   
int SkipNotSimilarRegion(Region&R1, Region&R2)
{
	double Cdis = ColorDist(R1.RR,R1.GG,R1.BB, R2.RR,R2.GG,R2.BB);
		  if(Cdis>40*40*3)
			  return 1;
		  
    double StrokeR1 = R1.Feature.StrokeMean;
    double StrokeR2 = R2.Feature.StrokeMean;

    if(StrokeR1==0||StrokeR2==0)return 1;
    StrokeR1/=StrokeR2;
    if(StrokeR1 >2 ||StrokeR1<0.5)return 1;

    return 0;
}


/// \brief  ͨ���������������������Աȣ�������������
///  
/// \param   
/// \return   
void MarkInvalidRegionByCrossValidate(vector<Region> & RegionVec)
{ 
  
  int i,j;
  double Weight;
   
  Loopi(RegionVec.size())
	  RegionVec[i].IsValid = 0;

  Weight = 1.0 + 1.5 +1.0 + 2.0 + 2.5;
  for(i=0;i<RegionVec.size();i++)
  {
	  if(SkipRegion(RegionVec[i]) ) continue;
      
	  
	  for (j=i+1;j<RegionVec.size();j++)
	  {
		  if(SkipRegion(RegionVec[j]) ) continue;
		  if( SkipNotSimilarRegion(RegionVec[i], RegionVec[j]))continue;

	
		  Region R1,R2;
		  R1 = RegionVec[i];
		  R2 = RegionVec[j];

		  int region_interset = checkRectIntersect(R1.left,R1.top, R1.right, R1.bottom,
						                          R2.left,R2.top, R2.right, R2.bottom);

	      if(!region_interset)

		  {double maxWidth   = max(RegionVec[i].rwidth  , RegionVec[j].rwidth );
		  double maxHeight  = max(RegionVec[i].rheight , RegionVec[j].rheight);

		  double maxWH      = max(maxWidth ,  maxHeight);

		  double sz1 = max(RegionVec[i].rwidth,RegionVec[i].rheight);
		  double sz2 = max(RegionVec[j].rwidth,RegionVec[j].rheight);
          double minWH = min(sz1,sz2); 


		  double Dist = sqrt( (R1.GeoX - R2.GeoX) *  (R1.GeoX - R2.GeoX) + 
			                  (R1.GeoY - R2.GeoY) *  (R1.GeoY - R2.GeoY) );
          
		  if(Dist > 2*minWH)continue;
		  }

		  RegionVec[i].IsValid =  RegionVec[j].IsValid = 1;
          
		  //printf("pair:%i,%i\n",i,j);

		  /*double Space      = min(abs(RegionVec[i].left - RegionVec[j].right), abs(RegionVec[j].left - RegionVec[i].right));
		  double Distance = 
			  fabs(RegionVec[i].rwidth  - RegionVec[j].rwidth) / maxWidth             *1.0    +
			  fabs(RegionVec[i].rheight - RegionVec[j].rheight)/(maxHeight*0.5)       *1.5    +
			  Space/maxWidth                                                          *1.0    +
			  min(fabs(RegionVec[i].top- RegionVec[j].top), fabs(RegionVec[i].bottom- RegionVec[j].bottom))/(maxHeight * 0.15) * 2.5;
			  Distance = Distance / Weight;
           if(Distance<1.1)
		  {
             RegionVec[i].IsValid =  RegionVec[j].IsValid = 1;
		  }*/
		   
	  }
  }


}

#endif