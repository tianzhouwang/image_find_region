#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "c24bitmap.h"
#include "c256bitmap.h"
#include "GenMser.h"
#include "region.h"
//#include "LoadFont.h"


/*
int i,j;
C24BitMap  CPic, DestPic;
C256BitMap GPic;
//CPic.Load(argv[1]);
GPic.Load("extchar.bmp");

  ConnectedComponents  ConnT;
  
	ConnT.left          = 0;
	ConnT.right         = GPic.LineWidth  - 1;
	ConnT.top           = 0;
	ConnT.bottom        = GPic.Height     - 1;
	ConnT.scan_width    = GPic.LineWidth;
	ConnT.scan_height   = GPic.Height;
	ConnT.pix_threshold = 80;
	ConnT.img_buffer_pt = GPic.Buffer;
	ConnT.Pic = &GPic;
	ConnT.alloc_space();
	ConnT.label_image(1);
	int obnum=ConnT.relabel_image(); 
  ConnT.GetResultGroup(); 
*/

//char ChrImgBuff[400000];

void ExpandPic(C256BitMap &GPic,C256BitMap&ExpPic)
{
  int i,j;
  int ExpWidth =50;
  ExpPic.FormatF(GPic.Width+ExpWidth,GPic.Height+ExpWidth);

  Loopi(ExpPic.LineWidth * ExpPic.Height)
	    ExpPic.Buffer[i] = 255;

  Loopi(GPic.Width)
	  Loopj(GPic.Height)
  {
	  *get_pix_color(ExpPic,i+ExpWidth/2,j+ExpWidth/2)=
		  *get_pix_color(GPic,i,j);
  }
}

void GetBlackPicRegion(C256BitMap &Pic,C24BitMap&CPic,vector<Region>&RegionVec)
{
	int i,j;
	int ww,hh;
	ww = Pic.Width;
	hh = Pic.Height;
	ConnectedComponents Cp;
		
	Cp.scan_width  = Pic.LineWidth;
	Cp.scan_height = Pic.Height;
	
	Cp.left  = 0;
	Cp.right = ww-1;
	Cp.top   = 0;
	Cp.bottom= hh-1; 
	Cp.pix_threshold =128;		
	Cp.alloc_space();
	
	//Cp.label_image(Pic.Buffer,1);
	Cp.label_image(Pic.Buffer,1);
	int obnum=Cp.relabel_image(); 
	Cp.GetResultGroup(); 	 
	CPic.FormatF(Pic.Width,Pic.Height);
	CPic.CleanPic(0);
	
	RegionVec.clear();
	RegionVec.resize(Cp.ResultVec.size());
	
	for( i=0;i<Cp.ResultVec.size();i++)
	{
		if(Cp.ResultVec[i].PtVec.size()< 10 )
			continue;
		
		CPic.RandPenColor();
		Loopj(Cp.ResultVec[i].PtVec.size())
		{
			CPic.SigDot(Cp.ResultVec[i].PtVec[j].x,
				CPic.Height-1-Cp.ResultVec[i].PtVec[j].y);
			
			RPoint Pt;
			Pt.x = Cp.ResultVec[i].PtVec[j].x;
			Pt.y = CPic.Height-1-Cp.ResultVec[i].PtVec[j].y;
			RegionVec[i].PtVec.push_back(Pt);
		}
	}

}

void AddRegionVec(vector<Region>&AllRegionVec,vector<Region>&RegionVec)
{
   int i;
   Loopi(RegionVec.size())
   {
	   AllRegionVec.push_back(RegionVec[i]);
   }
} 

int main(int argc, char*argv[])
{
 int i,j;
 C256BitMap Pic;
 vector<Region> RegionVec;
 int  PicWidth, PicHeight;
 vector<Region> AllRegionVec;
 vector<int> LabelVec;

	   C24BitMap CPic;

	   string fname ="";// argv[1];
	   fname = fname.substr(0,fname.length()-4);
 
	   C256BitMap GPic;
	  // GPic.FormatF(PicWidth, PicHeight);
	   GPic.Load("shapeRec.bmp");
	   GPic.Save("tmpA.bmp");
	  // memcpy(GPic.Buffer,ChrImgBuff,GPic.LineWidth*GPic.Height);
	   Pic.Clear();    
	   ExpandPic( GPic, Pic); 
	   Pic.Save("tmp2.bmp");
	   CPic.CleanPic(0);
//================================================================================
	   RegionVec.clear();
	   GetBlackPicRegion( Pic, CPic, RegionVec);
	   LabelVec.clear();
       GetObjContourColor(RegionVec, CPic, LabelVec);
	   AddRegionVec( AllRegionVec, RegionVec);
	   CPic.Save("region.bmp");
 
 

 DispReionObjvec(RegionVec,CPic.Width, CPic.Height,"char_region.bmp");
 CPic.Load("char_region.bmp");

 Loopi(RegionVec.size())
 {
    Loopj(RegionVec[i].CenterPts.size())
	 CPic.SigDot(RegionVec[i].CenterPts[j].x, RegionVec[i].CenterPts[j].y);
    
	vector<RPoint> ContourExp;
	GetRegionExpContour(RegionVec[i], ContourExp,4);
	Loopj(ContourExp.size())
	{
        CPic.SigDot(ContourExp[j].x, ContourExp[j].y);
	}
	ContourExp.clear();
 }/**/

 //PrintfRegionVec(AllRegionVec,"regionInfo.csv");

 CPic.Save("dest.bmp");
 return 1;
}
