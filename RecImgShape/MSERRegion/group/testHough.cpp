#include <stdio.h>
#include <stdlib.h>
#include "Hough.h"
#include "c24bitmap.h"
#include "C256BITMAP.H"
//#include "SubBmp.h"
#include "region.h"
#include "SelLinPoints.h"
#include <string>
using namespace std;

void main(int argc,char*argv[])
{
	int i,j;
	C256BitMap Pic;
	C24BitMap CPic;
	Pic.Load(argv[1]);
	CPic.FormatF(Pic.Width,Pic.Height);
	vector <ObjInfo> m_udResult;
	GetObjFromC256(Pic,m_udResult);
	vector<Pt2D>  PtVec;

	for(i=0;i<m_udResult.size();i++)
	{
		Pt2D Pt;
		Pt.x = m_udResult[i].x;
		Pt.y = m_udResult[i].y;
		PtVec.push_back(Pt);
		CPic.Rectangle(m_udResult[i].x-2,m_udResult[i].y-2,m_udResult[i].x+2,m_udResult[i].y+2);
	}
	
	/*CPic.PenColor.R =0;CPic.PenColor.G =255;CPic.PenColor.B =0;
	vector<int> Points;
	GetTriPoints(m_udResult, Points);
	
	for(i=0;i<(Points.size()/3);i++)
	{
		int Idx1,Idx2,Idx3;
		Idx1 = Points[3*i];
		Idx2 = Points[3*i+1];
		Idx3 = Points[3*i+2];
		CPic.DrawLine(m_udResult[Idx1].x,m_udResult[Idx1].y,m_udResult[Idx2].x,m_udResult[Idx2].y);
		CPic.DrawLine(m_udResult[Idx2].x,m_udResult[Idx2].y,m_udResult[Idx3].x,m_udResult[Idx3].y);
	}
	*/
	
	int PicWidth, PicHeight;
	vector< vector<int> > HoughMat;
	int  AngleDomain,  DistMatLen;
	HoughTransAcumlt( PtVec, CPic.Width, CPic.Height,HoughMat,0,180,2,15, AngleDomain,  DistMatLen);
    CPic.Save("color.bmp");

	C256BitMap GPic;
	GPic.FormatF(AngleDomain, DistMatLen);

	Loopi(AngleDomain)
		Loopj(DistMatLen)
	{
		*get_pix_color(GPic,i,j) = BOUND(HoughMat[i][j]*10,0,255);
	}


	GPic.Save("hough.bmp");
}