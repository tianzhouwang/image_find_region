 
#ifndef  Object_Stroke_Head
#define  Object_Stroke_Head
 
#include "mregion.h"
#include "objfeature.h"
#include <algorithm>

#define ABS(x) ((x) <0 ? -(x): x)

/// \brief  ������������
///  
/// \param val  ����ĸ�����    
/// \return     ��������������
int roundf(double val)
{   
	val += 0.5;
	return int(val);
}

/// \brief    ͨ�������صı�����������Ƿ�һ����
///           �õ������������
///  
/// \param RegionVec  һ��MSER����
/// \param LabelImg   ���жϵľ��� 
/// \param PicWidth   ͼ��Ŀ��
/// \param PicHeight  ͼ��ĸ߶�
/// \return     
void  GetRegionCountorPts(vector<Region>& RegionVec,vector<int>&LabelImg,int PicWidth,int PicHeight)
{
	int i,j,t;
	
	int offsetY[] = {-1, 1,  0, 0 };// 0,  1, 1, 1};
	int offsetX[] = { 0, 0, -1, 1 };// 1, -1, 0, 1};
	
	for(t=0;t<RegionVec.size();t++)
	{     
		Loopi(RegionVec[t].PtVec.size())
		{
			Loopj(4)
			{
				int x,y, pixpos;
				x = RegionVec[t].PtVec[i].x + offsetX[j];
				y = RegionVec[t].PtVec[i].y + offsetY[j];
				
				pixpos = y * PicWidth + x;
				
				if((x<0)||x>(PicWidth-1))  continue;
				if((y<0)||y>(PicHeight-1)) continue;
				
				if(LabelImg[pixpos] != RegionVec[t].Feature.label)
				{
					RegionVec[t].ContourPtVec.push_back(RegionVec[t].PtVec[i]);
					break;
				}
			}
		}
		
	}
	
}

//****************************************************************
//---------����Ϊ����8������ҿ��ٵõ������Ƕȡ�Dx��Dy��----------
//****************************************************************

double contour_angle[256]={
	180.00000, 45.00000, 90.00000, 71.56505,135.00000, 90.00000,108.43495, 90.00000,360.00000, 18.43495, 45.00000, 45.00000, 45.00000, 45.00000, 71.56505, 63.43495,
		180.00000,135.00000,135.00000,108.43495,161.56505,135.00000,135.00000,116.56505,180.00000, 45.00000, 90.00000, 71.56505,135.00000, 90.00000,108.43495, 90.00000,
		315.00000,360.00000, 45.00000, 45.00000,180.00000, 45.00000, 90.00000, 71.56505,341.56505,360.00000, 18.43495, 26.56505,360.00000, 18.43495, 45.00000, 45.00000,
		225.00000,180.00000,135.00000, 90.00000,180.00000,135.00000,135.00000,108.43495,315.00000,360.00000, 45.00000, 45.00000,180.00000, 45.00000, 90.00000, 71.56505,
		270.00000,315.00000,180.00000, 45.00000,225.00000,180.00000,135.00000, 90.00000,315.00000,341.56505,360.00000, 18.43495,315.00000,360.00000, 45.00000, 45.00000,
		225.00000,225.00000,180.00000,135.00000,198.43495,180.00000,161.56505,135.00000,270.00000,315.00000,180.00000, 45.00000,225.00000,180.00000,135.00000, 90.00000,
		288.43495,315.00000,315.00000,360.00000,270.00000,315.00000,180.00000, 45.00000,315.00000,333.43495,341.56505,360.00000,315.00000,341.56505,360.00000, 18.43495,
		251.56505,270.00000,225.00000,180.00000,225.00000,225.00000,180.00000,135.00000,288.43495,315.00000,315.00000,360.00000,270.00000,315.00000,180.00000, 45.00000,
		225.00000,180.00000,135.00000, 90.00000,180.00000,135.00000,135.00000,108.43495,315.00000,360.00000, 45.00000, 45.00000,180.00000, 45.00000, 90.00000, 71.56505,
		198.43495,180.00000,161.56505,135.00000,180.00000,161.56505,153.43495,135.00000,225.00000,180.00000,135.00000, 90.00000,180.00000,135.00000,135.00000,108.43495,
		270.00000,315.00000,180.00000, 45.00000,225.00000,180.00000,135.00000, 90.00000,315.00000,341.56505,360.00000, 18.43495,315.00000,360.00000, 45.00000, 45.00000,
		225.00000,225.00000,180.00000,135.00000,198.43495,180.00000,161.56505,135.00000,270.00000,315.00000,180.00000, 45.00000,225.00000,180.00000,135.00000, 90.00000,
		251.56505,270.00000,225.00000,180.00000,225.00000,225.00000,180.00000,135.00000,288.43495,315.00000,315.00000,360.00000,270.00000,315.00000,180.00000, 45.00000,
		225.00000,225.00000,198.43495,180.00000,206.56505,198.43495,180.00000,161.56505,251.56505,270.00000,225.00000,180.00000,225.00000,225.00000,180.00000,135.00000,
		270.00000,288.43495,270.00000,315.00000,251.56505,270.00000,225.00000,180.00000,296.56505,315.00000,315.00000,341.56505,288.43495,315.00000,315.00000,360.00000,
		243.43495,251.56505,225.00000,225.00000,225.00000,225.00000,198.43495,180.00000,270.00000,288.43495,270.00000,315.00000,251.56505,270.00000,225.00000,180.00000};
	
	double angle_dx[256] ={
		-1.00000, 0.70711, 0.00000, 0.31623,-0.70711, 0.00000,-0.31623, 0.00000, 1.00000, 0.94868, 0.70711, 0.70711, 0.70711, 0.70711, 0.31623, 0.44721,
			-1.00000,-0.70711,-0.70711,-0.31623,-0.94868,-0.70711,-0.70711,-0.44721,-1.00000, 0.70711, 0.00000, 0.31623,-0.70711, 0.00000,-0.31623, 0.00000,
			0.70711, 1.00000, 0.70711, 0.70711,-1.00000, 0.70711, 0.00000, 0.31623, 0.94868, 1.00000, 0.94868, 0.89443, 1.00000, 0.94868, 0.70711, 0.70711,
			-0.70711,-1.00000,-0.70711, 0.00000,-1.00000,-0.70711,-0.70711,-0.31623, 0.70711, 1.00000, 0.70711, 0.70711,-1.00000, 0.70711, 0.00000, 0.31623,
			-0.00000, 0.70711,-1.00000, 0.70711,-0.70711,-1.00000,-0.70711, 0.00000, 0.70711, 0.94868, 1.00000, 0.94868, 0.70711, 1.00000, 0.70711, 0.70711,
			-0.70711,-0.70711,-1.00000,-0.70711,-0.94868,-1.00000,-0.94868,-0.70711,-0.00000, 0.70711,-1.00000, 0.70711,-0.70711,-1.00000,-0.70711, 0.00000,
			0.31623, 0.70711, 0.70711, 1.00000,-0.00000, 0.70711,-1.00000, 0.70711, 0.70711, 0.89443, 0.94868, 1.00000, 0.70711, 0.94868, 1.00000, 0.94868,
			-0.31623,-0.00000,-0.70711,-1.00000,-0.70711,-0.70711,-1.00000,-0.70711, 0.31623, 0.70711, 0.70711, 1.00000,-0.00000, 0.70711,-1.00000, 0.70711,
			-0.70711,-1.00000,-0.70711, 0.00000,-1.00000,-0.70711,-0.70711,-0.31623, 0.70711, 1.00000, 0.70711, 0.70711,-1.00000, 0.70711, 0.00000, 0.31623,
			-0.94868,-1.00000,-0.94868,-0.70711,-1.00000,-0.94868,-0.89443,-0.70711,-0.70711,-1.00000,-0.70711, 0.00000,-1.00000,-0.70711,-0.70711,-0.31623,
			-0.00000, 0.70711,-1.00000, 0.70711,-0.70711,-1.00000,-0.70711, 0.00000, 0.70711, 0.94868, 1.00000, 0.94868, 0.70711, 1.00000, 0.70711, 0.70711,
			-0.70711,-0.70711,-1.00000,-0.70711,-0.94868,-1.00000,-0.94868,-0.70711,-0.00000, 0.70711,-1.00000, 0.70711,-0.70711,-1.00000,-0.70711, 0.00000,
			-0.31623,-0.00000,-0.70711,-1.00000,-0.70711,-0.70711,-1.00000,-0.70711, 0.31623, 0.70711, 0.70711, 1.00000,-0.00000, 0.70711,-1.00000, 0.70711,
			-0.70711,-0.70711,-0.94868,-1.00000,-0.89443,-0.94868,-1.00000,-0.94868,-0.31623,-0.00000,-0.70711,-1.00000,-0.70711,-0.70711,-1.00000,-0.70711,
			-0.00000, 0.31623,-0.00000, 0.70711,-0.31623,-0.00000,-0.70711,-1.00000, 0.44721, 0.70711, 0.70711, 0.94868, 0.31623, 0.70711, 0.70711, 1.00000,
			-0.44721,-0.31623,-0.70711,-0.70711,-0.70711,-0.70711,-0.94868,-1.00000,-0.00000, 0.31623,-0.00000, 0.70711,-0.31623,-0.00000,-0.70711,-1.00000};
		
		double angle_dy[256] ={
			0.00000, 0.70711, 1.00000, 0.94868, 0.70711, 1.00000, 0.94868, 1.00000,-0.00000, 0.31623, 0.70711, 0.70711, 0.70711, 0.70711, 0.94868, 0.89443,
				0.00000, 0.70711, 0.70711, 0.94868, 0.31623, 0.70711, 0.70711, 0.89443, 0.00000, 0.70711, 1.00000, 0.94868, 0.70711, 1.00000, 0.94868, 1.00000,
				-0.70711,-0.00000, 0.70711, 0.70711, 0.00000, 0.70711, 1.00000, 0.94868,-0.31623,-0.00000, 0.31623, 0.44721,-0.00000, 0.31623, 0.70711, 0.70711,
				-0.70711, 0.00000, 0.70711, 1.00000, 0.00000, 0.70711, 0.70711, 0.94868,-0.70711,-0.00000, 0.70711, 0.70711, 0.00000, 0.70711, 1.00000, 0.94868,
				-1.00000,-0.70711, 0.00000, 0.70711,-0.70711, 0.00000, 0.70711, 1.00000,-0.70711,-0.31623,-0.00000, 0.31623,-0.70711,-0.00000, 0.70711, 0.70711,
				-0.70711,-0.70711, 0.00000, 0.70711,-0.31623, 0.00000, 0.31623, 0.70711,-1.00000,-0.70711, 0.00000, 0.70711,-0.70711, 0.00000, 0.70711, 1.00000,
				-0.94868,-0.70711,-0.70711,-0.00000,-1.00000,-0.70711, 0.00000, 0.70711,-0.70711,-0.44721,-0.31623,-0.00000,-0.70711,-0.31623,-0.00000, 0.31623,
				-0.94868,-1.00000,-0.70711, 0.00000,-0.70711,-0.70711, 0.00000, 0.70711,-0.94868,-0.70711,-0.70711,-0.00000,-1.00000,-0.70711, 0.00000, 0.70711,
				-0.70711, 0.00000, 0.70711, 1.00000, 0.00000, 0.70711, 0.70711, 0.94868,-0.70711,-0.00000, 0.70711, 0.70711, 0.00000, 0.70711, 1.00000, 0.94868,
				-0.31623, 0.00000, 0.31623, 0.70711, 0.00000, 0.31623, 0.44721, 0.70711,-0.70711, 0.00000, 0.70711, 1.00000, 0.00000, 0.70711, 0.70711, 0.94868,
				-1.00000,-0.70711, 0.00000, 0.70711,-0.70711, 0.00000, 0.70711, 1.00000,-0.70711,-0.31623,-0.00000, 0.31623,-0.70711,-0.00000, 0.70711, 0.70711,
				-0.70711,-0.70711, 0.00000, 0.70711,-0.31623, 0.00000, 0.31623, 0.70711,-1.00000,-0.70711, 0.00000, 0.70711,-0.70711, 0.00000, 0.70711, 1.00000,
				-0.94868,-1.00000,-0.70711, 0.00000,-0.70711,-0.70711, 0.00000, 0.70711,-0.94868,-0.70711,-0.70711,-0.00000,-1.00000,-0.70711, 0.00000, 0.70711,
				-0.70711,-0.70711,-0.31623, 0.00000,-0.44721,-0.31623, 0.00000, 0.31623,-0.94868,-1.00000,-0.70711, 0.00000,-0.70711,-0.70711, 0.00000, 0.70711,
				-1.00000,-0.94868,-1.00000,-0.70711,-0.94868,-1.00000,-0.70711, 0.00000,-0.89443,-0.70711,-0.70711,-0.31623,-0.94868,-0.70711,-0.70711,-0.00000,
				-0.89443,-0.94868,-0.70711,-0.70711,-0.70711,-0.70711,-0.31623, 0.00000,-1.00000,-0.94868,-1.00000,-0.70711,-0.94868,-1.00000,-0.70711, 0.00000};

using namespace std;
//ע������û��һ����Point��Region�ṹ��ԭ����
//�����˽ṹ�е�������������ǰ�����ٶ���һ�����ϣ��������ڴ�Ĺ�ϵ

struct PointDxy       //��ǵ�
{
	int x,y;           //�������
	double dx,dy;      //��ķ���
	double angle;      //��ĽǶ�
	int ifBoundPoint;  //�Ƿ�߽��
	int chain       ;  //������ֵ
	int pixIntensity;  //����ֵ
	int ptOffset    ;  //���ص㵽buffer��ƫ��
	int Label       ;  //��ͨ������
	int groupIdx    ;  //������ͨ���ڵ����
	int R,G,B;
	float strokeWidth;
};

struct PointPair
{
	PointDxy Pt1,Pt2;
	double  Distance;
	double  AngleDiff;
	double  Angle;
	int     idx1,idx2;
};
//-------------------------------------------------------------
struct PairInfoStat
{
	double meanVal, stdevVal, medianVal,minVal,maxVal, Fraction;//ƽ���Ȼ�
	double AnglemeanVal, AnglestdevVal, AnglemedianVal, AngleFraction;
	double meanAngle;      //ƽ���Ƕ�
};
//-------------------------------------------------------------

//ר�����Լ���ʻ���ȵĽṹ��
struct CharRegion   
{
	vector<PointDxy> PtVec;
	vector<PointDxy> contourPts;
	vector<PointDxy> CenterLine;
	vector<float>    strokelen;
	int indx,peak;
	int left,right,top,bottom;
	int width,height;
	int Label;               //<-----����ı�����
	double RoundRatio;       //Բ�Ķ�
	//double Average;
	double   CharEnergy;
	vector<float> feature;
	PairInfoStat ObjStrokeInfo;
};

/// \brief  ����ʻ�ƽ���������ܱʻ��е�ռ��
/// \param  PairVec    �����ڶԳƵ������
/// \param  pairInfo   ����ıʻ�������Ϣ
/// \return ��
void CptPointPairParameterLength(vector<PointPair>&PairVec, PairInfoStat&pairInfo)
{
   int i;
   vector<double> lengthVec;
   Loopi(PairVec.size())
	   lengthVec.push_back(PairVec[i].Distance);

   double sumVal,powVal;
   sumVal = powVal = 0;
   Loopi(lengthVec.size())
   {
	   sumVal += lengthVec[i];
	   powVal += lengthVec[i]*lengthVec[i];
   }
 
   pairInfo.meanVal  = sumVal/double(lengthVec.size());
   pairInfo.stdevVal = sqrt(powVal /lengthVec.size() - pairInfo.meanVal * pairInfo.meanVal);
   int middle = int(lengthVec.size()/2);
   nth_element(lengthVec.begin(),lengthVec.begin()+middle,lengthVec.end()); 
   pairInfo.medianVal =  lengthVec[middle];
   
   double cnt=0;
   
   Loopi(lengthVec.size())
   {
	   double diff = ABS((lengthVec[i] - pairInfo.medianVal)/pairInfo.medianVal);
	   if(diff<0.2)
		   cnt++;
   }
   pairInfo.Fraction = cnt / double(lengthVec.size());
}

/// \brief  ����ƽ���Ƕ����ܽǶ��е�ռ��
/// \param  PairVec    �����ڶԳƵ������
/// \param  pairInfo   ����ıʻ��Ƕ���Ϣ
/// \return ��
void CptPointPairParameterAngle(vector<PointPair>&PairVec,
								 PairInfoStat&pairInfo)
{
	int i;
	vector<double> lengthVec;
	Loopi(PairVec.size())
		lengthVec.push_back(PairVec[i].Angle);
	
	double sumVal,powVal;
	sumVal = powVal = 0;

	Loopi(lengthVec.size())
	{
		sumVal += lengthVec[i];
		powVal += lengthVec[i]*lengthVec[i];
	}
    
	pairInfo.AnglemeanVal  = sumVal/double(lengthVec.size());
	pairInfo.AnglestdevVal = sqrt(powVal /lengthVec.size() - pairInfo.meanVal * pairInfo.meanVal);
	int middle = int(lengthVec.size()/2);
	nth_element(lengthVec.begin(),lengthVec.begin()+middle,lengthVec.end()); 
	pairInfo.AnglemedianVal =  lengthVec[middle];
	
	double cnt=0;
	
	Loopi(lengthVec.size())
	{
		double diff = ABS((lengthVec[i] - pairInfo.medianVal)/pairInfo.medianVal);
		if(diff<0.2)
			cnt++;
	}
	pairInfo.AngleFraction = cnt/double(lengthVec.size());
}

//����ʻ��������
class ObjStroke
{
public:
	vector<int>       LableMap;      //������ı����ֵ
	vector<int>       contourPos;    //������ľ���λ��
	vector<PointDxy>   PtsMap;       //��������Ϣ
	int picWidth,picHeight,picSize;  //ͼ��Ŀ���
 
	void  CalculateObjStroke(CharRegion & Region,vector<PointPair>&pairVec);
	float GetAngle(float x,float y);
	void CalcuObjEnery(CharRegion & Region);
	void DispObjvecSpectrum(vector<CharRegion>&ObjVec,char*filename);
	void computePointsAngle(vector<CharRegion>&regionVec,int*Ptlbs);
	void InitRegions(vector<Region>&RegionVec,int Width,int Height);
	void InitRegions(int Width,int Height);
	vector<CharRegion> StkObjRegions;
	void  GetObjVecStrokeInfo();
    float GetCharStrokeLengthFromCenter(CharRegion & Region,int x,int y);

};

/// \brief  ����һ�������ľ�ֵ�����Сֵ����ֵ������
/// \param  fvec    ��������
/// \param  minv, maxv, mean, median, stdev ����������ֵ
/// \return ��
void GetVectorInfo(vector<float>&fvec,float&minv,float&maxv,float&mean,float&median,float&stdev)
{
	int i;
	double sumVal,powVal;
	sumVal = powVal = 0;
	if(fvec.size()==0)return;
	minv = maxv = fvec[0];
	Loopi(fvec.size())
	{
		float val = fvec[i];
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
}


/// \brief  �õ�����ʻ���Ϣ
/// \return ��
void ObjStroke::GetObjVecStrokeInfo()
{
	int i,j;
    for(i=0;i<StkObjRegions.size();i++)
	{
		StkObjRegions[i].CharEnergy =0;
		vector<PointPair> PairVec; PairVec.clear();
		CalculateObjStroke(StkObjRegions[i],PairVec);
		
		if(PairVec.size()==0) 
			continue;	

		if( StkObjRegions[i].strokelen.size()>0)
		{   
            float minv, maxv, mean, median, stdev;
			GetVectorInfo(StkObjRegions[i].strokelen, minv, maxv, mean, median, stdev);
            StkObjRegions[i].ObjStrokeInfo.meanVal   = mean; 
            StkObjRegions[i].ObjStrokeInfo.stdevVal  = stdev;
            StkObjRegions[i].ObjStrokeInfo.medianVal = median;
			StkObjRegions[i].ObjStrokeInfo.minVal    = minv;
			StkObjRegions[i].ObjStrokeInfo.maxVal    = maxv;
		}
		else
			StkObjRegions[i].ObjStrokeInfo.meanVal = 0;
	   

	}
}

void ShowObjvecSpectrum(C24BitMap&KBPic,vector<CharRegion>&ObjVec);

/// \brief  ��ʾ�ʻ��Ĳʺ�ɫ��ͼ
/// \param  ObjVec     һ������
/// \param  filename   ���ɵ�ͼ���ļ���
/// \return ��
void ObjStroke::DispObjvecSpectrum(vector<CharRegion>&ObjVec,char*filename)
{
   C24BitMap KBPic;
   KBPic.FormatF(picWidth,picHeight);
   KBPic.ClearPic(0);
  // ShowObjvecSpectrum( KBPic, ObjVec);
   KBPic.Save(filename);
}


/// \brief  ����õ������ǶȵĽǶȷ���  
/// \param  regionVec  ͼ���е�һ������
/// \param  Ptlbs      �������ͨ����Ϣ  
/// \return ��
void ObjStroke::computePointsAngle(vector<CharRegion>&regionVec,int*Ptlbs)
{
	int i,j,t;
	int offsetY[] = {-1, -1, -1,  0, 0,  1, 1, 1};
	int offsetX[] = {-1,  0,  1, -1, 1, -1, 0, 1};
	int pixoffset[8];
	int nNeighbors = 8;
	for (i=0;i<nNeighbors;i++)
	{
		pixoffset[i] = offsetX[i] + picWidth * offsetY[i];
	}
	
	Loopi(regionVec.size())
	{
		Loopj(regionVec[i].contourPts.size())
		{
			int x,y;
			x = regionVec[i].contourPts[j].x;
			y = regionVec[i].contourPts[j].y;
			int pixel_pos = y*picWidth+x;
			int pixlabel = Ptlbs[pixel_pos];
			int nbpattern = 0;
			for(t=0;t<nNeighbors;t++)
			{  
				int offset = pixoffset[t];
				int nblabel = Ptlbs[pixel_pos+offset];
				if(nblabel != pixlabel )
					nbpattern+= 1<<t;
			}
			regionVec[i].contourPts[j].angle = contour_angle[nbpattern];
			regionVec[i].contourPts[j].dx    = angle_dx[nbpattern];
			regionVec[i].contourPts[j].dy    = angle_dy[nbpattern];
		}
	}
		
}

/// \brief  ��ʼ���ʻ�����ṹ��
/// \param  Width      ͼ��Ŀ��
/// \param  Height     ͼ��ĸ߶� 
/// \return ��
void ObjStroke::InitRegions(int Width,int Height)
{
	int i,j;
	picWidth  = Width; picHeight = Height;
	picSize   = Width * Height;
	
	PtsMap.clear(); PtsMap.resize(picSize);//LableMap.clear(); contourPos.clear();
	if(LableMap.size()!=picSize)    LableMap.resize(picSize);
	if(contourPos.size()!=picSize)  contourPos.resize(picSize);	 
	
	fill(LableMap.begin()   , LableMap.end()   , -1);
    fill(contourPos.begin() , contourPos.end() ,  0);
}

/// \brief  ��ʼ���ʻ�����ṹ��
/// \param  RegionVec  ��������������
/// \param  Width      ͼ��Ŀ��
/// \param  Height     ͼ��ĸ߶� 
/// \return ��
void ObjStroke::InitRegions(vector<Region>&RegionVec,int Width,int Height)
{
	int i,j;
	picWidth  = Width; picHeight = Height;
	picSize   = Width * Height;
	
	PtsMap.clear(); PtsMap.resize(picSize);//LableMap.clear(); contourPos.clear();
	if(LableMap.size()!=picSize)    LableMap.resize(picSize);
	if(contourPos.size()!=picSize)  contourPos.resize(picSize);	 

	fill(LableMap.begin()   , LableMap.end()   , -1);
    fill(contourPos.begin() , contourPos.end() ,  0);

    //----------------��������stroke��Obj�ṹ������-------------------
	StkObjRegions.resize(RegionVec.size());
	//--------------------����������ͨ���----------------------------
    for (i=0; i<RegionVec.size(); i++)
	{
	    RegionVec[i].Feature.label =i;
        for(j=0;j<RegionVec[i].PtVec.size();j++)
		{
		   int Pos = 
			   RegionVec[i].PtVec[j].y * picWidth + RegionVec[i].PtVec[j].x; 
		   LableMap[Pos] = i;
		}
	}
    //----------------ɨ�������ǵõ��߽��---------------------------
	GetRegionCountorPts(RegionVec,LableMap,picWidth,picHeight);

	for (i=0; i<RegionVec.size(); i++)
	{
		StkObjRegions[i].Label = i;
		PointDxy Pt;
		int left, right, top, bottom;
        left = right = top = bottom =0 ;
        if(RegionVec[i].PtVec.size()>0)
		{
           left =  right  = RegionVec[i].PtVec[0].x;
		   top  =  bottom = RegionVec[i].PtVec[0].y;
		}

        for(j=0;j<RegionVec[i].PtVec.size();j++)
		{	
			Pt.x = RegionVec[i].PtVec[j].x;
			Pt.y = RegionVec[i].PtVec[j].y;
			Pt.Label = i;
			StkObjRegions[i].PtVec.push_back(Pt);
			if(Pt.x < left)left = Pt.x; if(Pt.x > right)  right  = Pt.x;
			if(Pt.y < top)  top = Pt.y; if(Pt.y > bottom) bottom = Pt.y;
		}
		
		StkObjRegions[i].left = left;  StkObjRegions[i].right  = right;
		StkObjRegions[i].top  = top;   StkObjRegions[i].bottom = bottom;
        StkObjRegions[i].width  = right  - left + 1;
		StkObjRegions[i].height = bottom - top  + 1;
		
		Loopj(RegionVec[i].ContourPtVec.size())
		{
			Pt.x = RegionVec[i].ContourPtVec[j].x;
			Pt.y = RegionVec[i].ContourPtVec[j].y;
			
			if((Pt.x == 0)||(Pt.x == (picWidth -1))) 
				continue;
			if((Pt.y == 0)||(Pt.y == (picHeight-1))) 
				continue;
			Pt.Label =i;
			StkObjRegions[i].contourPts.push_back(Pt);
	   }
	}

	computePointsAngle(StkObjRegions,&LableMap[0]);

	Loopi(StkObjRegions.size())
	{
		StkObjRegions[i].Label = i;

		Loopj(StkObjRegions[i].contourPts.size())
		{
			StkObjRegions[i].contourPts[j].Label    = i;
			StkObjRegions[i].contourPts[j].groupIdx = j;                //������ͨ�������ڵ����
			StkObjRegions[i].contourPts[j].ptOffset = 
				StkObjRegions[i].contourPts[j].y * picWidth + StkObjRegions[i].contourPts[j].x; 
			
			int Pos          = StkObjRegions[i].contourPts[j].ptOffset; //�������ʼ��ƫ��
			PtsMap[Pos]      = StkObjRegions[i].contourPts[j];          //���õ�
			contourPos[Pos]  = 1;                                    //���Ϊ����
		}
	}

   

}

/// \brief  ����ʻ����
///  �����ĵ㿪ʼ�����ĸ�����Ѱ�ҳ�����̵�������Ϊ�ʻ����
/// \param  Region     ������һ������ 
/// \param  x     ���ĵ��x
/// \param  y     ���ĵ��y
/// \return ��
float ObjStroke::GetCharStrokeLengthFromCenter(CharRegion & Region,int x,int y)
{
	int   offsetY[] = {0,  1,  1,  1 };// 0,  1, 1, 1};
    int   offsetX[] = {1,  0,  1, -1 };// 1, -1, 0, 1};
	float distlen[] = {1.0, 1.0, 1.414, 1.414};
    int i,j,t;
	int MinLen = 9999;
    float Len1,Len2;

	Len1 = Len2 = MinLen;
	 for(t=0;t<4;t++)
    {
       for(i=0;i<MinLen;i++)
	  {
	    int endX,endY;	 
		endX = x + i*offsetX[t];
		endY = y + i*offsetY[t];
		if((endX > (picWidth-1 )) ||
		   (endY > (picHeight-1)) ||
		   (endX <0) ||
		   (endY <0))
		{   Len1 = MinLen;
			break;
		}

		int LabelPos = endY * picWidth + endX;
		int EndLabel = LableMap[LabelPos];

		if(EndLabel!=Region.Label)
		{
           Len1 = distlen[t] *float(i); 
		   break;
		}
	  }
      
	   for(i=0;i<MinLen;i++)
	   {
		   int endX,endY;	 
		   endX = x - i*offsetX[t];
		   endY = y - i*offsetY[t];
		   if((endX > (picWidth-1 )) ||
			   (endY > (picHeight-1)) ||
			   (endX <0) ||
			   (endY <0))
		   {   Len2 = MinLen;
		     break;
		   }
		   
		   int LabelPos = endY * picWidth + endX;
		   int EndLabel = LableMap[LabelPos];
		   
		   if(EndLabel!=Region.Label)
		   {
			   Len2 = distlen[t] *float(i); 
			   break;
		   }
	  }

	   if( (Len1+Len2) <MinLen) 
		   MinLen = Len1 + Len2;
	}
	 return MinLen;
}

/// \brief  ����ʻ����
///  ����һ�������ڵ�����������������������ıʻ���Ϣ
/// \param  Region     ������һ������ 
/// \param  pairVec    �������ںͳ��ֶԳ�״̬�ıʻ��߶�
/// \return ��
void ObjStroke::CalculateObjStroke(CharRegion & Region,vector<PointPair>&pairVec)
{
	int i,j;

	Loopi(Region.contourPts.size())
	{
      PointDxy Pt = Region.contourPts[i];
      double step =1;
	  int maxstepSize = Region.width + Region.height +2;
	  double startX,startY;
	  int startLb;
	  startLb = Region.Label;//Region.contourPts[i].Label;
	  startX  = Region.contourPts[i].x;
	  startY  = Region.contourPts[i].y;

	  if(maxstepSize> picWidth/5) maxstepSize = picWidth/5;
     
	  for(step =1; step < maxstepSize; step++)
	  {
		 
         float nextX = roundf(startX + Pt.dx * step);
         float nextY = roundf(startY + Pt.dy * step);
         if(nextX > (picWidth -1) || nextX <1)break;
		 if(nextY > (picHeight-1) || nextY <1)break;
		 int  Posoffset;
	     Posoffset = nextY * picWidth + nextX;
         if(LableMap[Posoffset]!=startLb)
		 {
			int LbNdebug = LableMap[Posoffset]; 
			if(step == 1)break;
            step -=1;
			nextX = roundf(startX + Pt.dx * step); nextY = roundf(startY + Pt.dy * step);
			Posoffset = nextY * picWidth + nextX;

			if(LableMap[Posoffset]!=startLb)break;
			if(contourPos[Posoffset]!=1) break;
            
			int CenterX,CenterY;
			CenterX = (startX + nextX)/2;
			CenterY = (startY + nextY)/2;
			PointDxy PtT = PtsMap[Posoffset];

			if(PtT.Label!=startLb)
                break;

			if(ABS(ABS(PtT.angle - Pt.angle) -180.0) >25)
				break;

			PointDxy CenterPt;
            CenterPt.x = CenterX;  CenterPt.y = CenterY;
			Region.CenterLine.push_back(CenterPt);
            Region.strokelen.push_back(GetCharStrokeLengthFromCenter(Region,CenterPt.x,CenterPt.y));
          
			double AngleDiff = ABS(ABS(Pt.angle - PtT.angle) - 180.0);
	        if( AngleDiff>35.0)
				 break;

            PointPair P;
			P.AngleDiff =  AngleDiff;
			P.Angle = (Pt.angle+PtT.angle-180.0)/2;
			P.idx1  = i;
			P.idx2  = PtT.groupIdx;
			P.Pt1   = Region.contourPts[i];
			P.Pt2   = Region.contourPts[PtT.groupIdx]; //P.Pt2   = Region.contourPts[PtT.groupIdx];
			P.Distance = step;
			pairVec.push_back(P);
			break;
		 }
	  }
	}

}


#endif