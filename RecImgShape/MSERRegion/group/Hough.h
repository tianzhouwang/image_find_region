#ifndef Hough_Line_Dec_Head
#define Hough_Line_Dec_Head
#include <vector>
#include <math.h>
#define M_PI 3.1415926535897932384626433832795
#define DBG_HOUGHMAT_FILE 1
using namespace std;

struct Pt2D
{
  int x,y;
};

void HoughTransAcumlt(vector<Pt2D> &PtVec, int PicWidth,int PicHeight,
					  vector< vector<int> >&HoughMat,
                      float LowAngle,float HighAngle,float Accuracy,float DistStep,
					  int & AngleDomain,int& DistMatLen)//12
{
 int i,j;
 
 
 AngleDomain = (HighAngle-LowAngle)/Accuracy + 1;//角度范围
 
 
 DistMatLen = double(PicHeight*(3/2))/DistStep;

if(DistMatLen>500)
DistMatLen/=2;
/*矩阵的最大列计算(p,S)的S,S为直线到原点的距离
所以S Max＝sqrt(X^2+Y^2)取正+1*/

HoughMat.resize(AngleDomain); 
for(i=0;i<AngleDomain;i++)HoughMat[i].resize(DistMatLen);  

int CalcBMPHeight(PicHeight),CalcBMPWidth(PicWidth);
//if(CalcBMPHeight>500)CalcBMPHeight/=2;
if(CalcBMPWidth>800)CalcBMPWidth/=2;

//Loopi(PtVec.size())
for(i=0;i<PtVec.size();i++)
{
  for (int AglNote=0;AglNote<AngleDomain;AglNote++)
  {
	double Dist;
	float Angle;
	 Angle=AglNote*Accuracy+LowAngle;
     Dist= (double(PtVec[i].x) * cos(Angle/180.0*M_PI) +
		    double(PtVec[i].y) * sin(Angle/180.0*M_PI)) /DistStep;
      if((Dist>=1)&&((Dist+1)<DistMatLen))
	  {
		  HoughMat[AglNote][int(Dist-1)]+=1;
		  HoughMat[AglNote][int(Dist  )]+=3;
		  HoughMat[AglNote][int(Dist+1)]+=1;
	  }
  }

}

if(DBG_HOUGHMAT_FILE)
{
C24BitMap HoughRawPic;
HoughRawPic.FormatF(PicWidth,PicHeight);
HoughRawPic.CleanPic();

Loopi(PtVec.size())
{
	HoughRawPic.DrawCircle(PtVec[i].x,PtVec[i].y,1);
}	

HoughRawPic.Save("result\\rawdot.bmp");

FILE*file = fopen("result\\vote.txt","wt+");
for(j=0;j<HoughMat.size();j++)
{
	for(i=0;i<DistMatLen;i++)
	{
	  fprintf(file,"%2i,",HoughMat[j][i]);
	}

	fprintf(file,"\n");
}
fclose(file);
 }
 
}

#endif