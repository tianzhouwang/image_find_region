#ifndef GetObjContour_Head_H
#define GetObjContour_Head_H
//#include "msquare.h"
#include "mregion.h"
#include "objstrok.h"

void GetObjContour(int Width,int Height,vector<Region>&RgVec,vector<int>&LabelVec)
{
	int i,j;
	vector<int> ContourMp;
	int TotalCnt =  Width * Height;
	LabelVec.resize(TotalCnt);
    ContourMp.resize(TotalCnt);

	int ObjNum = RgVec.size();
    
	Loopi(TotalCnt)  LabelVec[i]  = -1;
	Loopi(TotalCnt)  ContourMp[i] =  0;
		 
	Loopi(RgVec.size())
	{   RgVec[i].ContourPtVec.clear();
	    RgVec[i].InflectionPts.clear();

		Loopj(RgVec[i].PtVec.size())
		{
			int xx,yy;
			xx = RgVec[i].PtVec[j].x ;
			yy = RgVec[i].PtVec[j].y ;
			int Pos;
			Pos = yy * Width + xx;
			LabelVec[Pos] = i;
		}
  }  
	
  
  MatchSquare EdgeLineTracer;
  EdgeLineTracer.PicWidth = Width;
  EdgeLineTracer.PicHeight= Height;
  EdgeLineTracer.GenerateSquareLinesMultiObj(LabelVec, Width, Height, ObjNum);
  EdgeLineTracer.RenamePointId();
  EdgeLineTracer.GetPath();
  EdgeLineTracer.GetInflectionPoint();
  int s1,s2;

  s1 = EdgeLineTracer.LinkPath.size();
  s2 = EdgeLineTracer.PathLbVec.size();

  Loopi(EdgeLineTracer.LinkPath.size())
  {
	  EdgeLineTracer.SmoothTangentDirections(EdgeLineTracer.LinkPath[i], 3);
	  Loopj(EdgeLineTracer.LinkPath[i].size())
	  {
		  float x,y;
		  RPoint Pt;
		  int Idx = EdgeLineTracer.LinkPath[i][j]-1;
		  Pt.x = EdgeLineTracer.PtVec[Idx].x;  
		  Pt.y = EdgeLineTracer.PtVec[Idx].y;
		  int Label = EdgeLineTracer.PtVec[Idx].Label;
		  RgVec[Label].ContourPtVec.push_back(Pt);
		  RgVec[Label].Feature.AngleCnt[0] = EdgeLineTracer.PathLbVec[i].VoteVec[0];
		  RgVec[Label].Feature.AngleCnt[1] = EdgeLineTracer.PathLbVec[i].VoteVec[1];
		  RgVec[Label].Feature.AngleCnt[2] = EdgeLineTracer.PathLbVec[i].VoteVec[2];
		  RgVec[Label].Feature.AngleCnt[3] = EdgeLineTracer.PathLbVec[i].VoteVec[3];
		  if(EdgeLineTracer.PtVec[Idx].InflectionPt)
             RgVec[Label].InflectionPts.push_back(Pt);
	  }
  }
   
  ObjStroke ObjStk;
  ObjStk.InitRegions(Width,Height);
  ObjStk.LableMap  = LabelVec;
  ObjStk.picWidth  = Width;
  ObjStk.picHeight = Height;

  vector<CharRegion> regionVec;
  regionVec.resize(RgVec.size());
   
  Loopi(RgVec.size())
  {
	int left, right, top, bottom;
	left = right = top = bottom =0 ;
	if(RgVec[i].PtVec.size()>0)
	 {
	  left =  right  = RgVec[i].PtVec[0].x;
	  top  =  bottom = RgVec[i].PtVec[0].y;
	}

    Loopj(RgVec[i].ContourPtVec.size())
	{
	 PointDxy Pt;
	 Pt.x = RgVec[i].ContourPtVec[j].x;
	 Pt.y = RgVec[i].ContourPtVec[j].y;
	 Pt.Label = i;
     if(Pt.x < left)left = Pt.x; if(Pt.x > right)  right  = Pt.x;
	 if(Pt.y < top)  top = Pt.y; if(Pt.y > bottom) bottom = Pt.y;

	 for(int ii = (Pt.x-1);ii<=(Pt.x+1);ii++)
		 for(int jj=(Pt.y-1);jj<=(Pt.y+1);jj++)
		 {
			  ii = BOUND(ii,0,Width  - 1);
			  jj = BOUND(jj,0,Height - 1);

			 int Pos = jj*Width +ii;
			 if(ContourMp[Pos]==0&&LabelVec[Pos]==i)
			 { ContourMp[Pos]=1;
               Pt.x = ii;
			   Pt.y = jj;
			   regionVec[i].contourPts.push_back(Pt);
			   break;
			 }
		 }
	}

	regionVec[i].left = left;  regionVec[i].right  = right;
	regionVec[i].top  = top;   regionVec[i].bottom = bottom;
	regionVec[i].width  = right  - left + 1;
	regionVec[i].height = bottom - top  + 1;
  }

  ObjStk.computePointsAngle(regionVec, &LabelVec[0]);
 
  ObjStk.contourPos = ContourMp;
  ObjStk.StkObjRegions = regionVec;
  
  Loopi(ObjStk.StkObjRegions.size())
  {
	  ObjStk.StkObjRegions[i].Label = i;
	  
	  Loopj(ObjStk.StkObjRegions[i].contourPts.size())
	  {
		  ObjStk.StkObjRegions[i].contourPts[j].Label    = i;
		  ObjStk.StkObjRegions[i].contourPts[j].groupIdx = j;                //设置联通域在组内的序号
		  ObjStk.StkObjRegions[i].contourPts[j].ptOffset = 
			  ObjStk.StkObjRegions[i].contourPts[j].y * Width + ObjStk.StkObjRegions[i].contourPts[j].x; 
		  
		  int Pos          = ObjStk.StkObjRegions[i].contourPts[j].ptOffset; //相对于起始的偏移
		  ObjStk.PtsMap[Pos]      = ObjStk.StkObjRegions[i].contourPts[j];          //设置点
		  ObjStk.contourPos[Pos]  = 1;                                    //标记为轮廓
	  }
	}

  ObjStk.GetObjVecStrokeInfo();

  Loopi(ObjStk.StkObjRegions.size())
  {
	  RgVec[i].Feature.StrokeMean = ObjStk.StkObjRegions[i].ObjStrokeInfo.meanVal;
	  RgVec[i].Feature.StrokeStd  = ObjStk.StkObjRegions[i].ObjStrokeInfo.stdevVal; 
	  RgVec[i].Feature.StkMin     = ObjStk.StkObjRegions[i].ObjStrokeInfo.minVal;
	  RgVec[i].Feature.StkMax     = ObjStk.StkObjRegions[i].ObjStrokeInfo.maxVal;
	  RgVec[i].Feature.StkMedian  = ObjStk.StkObjRegions[i].ObjStrokeInfo.medianVal;
	  //printf("%.2f,%.2f\n",RgVec[i].Feature.StrokeMean, RgVec[i].Feature.StrokeStd);
  }

}


#endif