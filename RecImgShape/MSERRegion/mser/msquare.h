#ifndef MQS_Head_H
#define MQS_Head_H

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include <vector>
#include "c24bitmap.h"
#include "c256bitmap.h"


using namespace std;

struct St_IdPoint {
	int oldId,newID;
	float x, y, z;
	int xi,yi,zi;
	int PatIdx;
	int InflectionPt;
	int Label;
	//double Val;
};

struct St_Line
{
	int Pt1Idx,Pt2Idx;
	St_IdPoint Pt1,Pt2;
	//int Label;
};

struct PathLb
{
  int Label;
  float VoteVec[4]; 
}; 

typedef  map< int, St_IdPoint> Id2PointMap;

class MatchSquare {
public:
	// Constructor and destructor.
	MatchSquare(){};
	~MatchSquare(){};
	int GetPointIdx(int x,int y,int edgnum,St_IdPoint&Pt);
	int PicWidth,PicHeight;
	static const  int m_LineTable[16][4];
	int GetLinePattern(int x1,int y1,int x2,int y2);
	C256BitMap GPic;
	vector<St_IdPoint> PtVec;
	vector<St_IdPoint> ProblemPts;
	vector<St_Line> LineVec;
	vector<int> LinkIdx;
	vector< vector<int> >LinkPath;
	vector< PathLb > PathLbVec;
	void GenerateSquareLines();
	void GetPath();
	void AddLink(int Id1,int Id2);
	void RenamePointId();
	int	 SmoothTangentDirections( vector<int>&PtIdxVec, int nIter );
	void GetInflectionPoint();
	void MakeWhiteEdge();
	int  InflectPtCnt;
	int  ObjArea;
	void GenerateSquareLinesMultiObj(vector<int>&ImgVec,int ww,int hh,int ObjNum);
};

void MatchSquare::MakeWhiteEdge()
{
	int i,j;
	Loopi( GPic.Width )
	{
		*get_pix_color(GPic,i,0)             = 255; 
		*get_pix_color(GPic,i,GPic.Height-1) = 255; 
	}
	
	Loopi(GPic.Height)
	{
		*get_pix_color(GPic,0,i)            = 255;
		*get_pix_color(GPic,GPic.Width-1,i) = 255;
	}

}

void MatchSquare::GetInflectionPoint()
{
	//double atan2(double y, double x)
	int i,j;
	InflectPtCnt = 0;
	Loopi(PtVec.size())
		PtVec[i].InflectionPt = 0;
	
	float ValArr[4];
   
    float Cnt = 0;
	Loopi(LinkPath.size())
	{   
		 Cnt=0;
		 ValArr[0] = ValArr[1] = ValArr[2] = ValArr[3] = 0;
		//printf("%i,%i\n",i,LinkPath[i].size());
		for(j=4;j<(LinkPath[i].size()-4);j++)
		{
			int Idx1, Idx3, Idx5;
			Idx3 = LinkPath[i][j  ] - 1;
			Idx1 = LinkPath[i][j-4] - 1;
			Idx5 = LinkPath[i][j+4] - 1;
			
			int xxi,yyi;
			xxi = PtVec[Idx3].x;
			yyi = PtVec[Idx3].y;
			 
			double dx1,dy1,dx2,dy2;		
			dx1 = PtVec[Idx1].x - PtVec[Idx3].x;
			dy1 = PtVec[Idx1].y - PtVec[Idx3].y;
			
			dx2 = PtVec[Idx3].x - PtVec[Idx5].x;
			dy2 = PtVec[Idx3].y - PtVec[Idx5].y;
			
			int Vote =  (atan2(dy1, dx1) + M_PI) / (M_PI/2.0);
			Vote = BOUND(Vote, 0, 3);
            ValArr[Vote] += 1.0;
            Cnt++;

			double Val = (dx1*dx2 + dy1*dy2) / (sqrt(dx1*dx1+dy1*dy1)*sqrt(dx2*dx2+ dy2*dy2));
			if(fabs(Val)<0.8)
			{
				PtVec[Idx3].InflectionPt = 1;
				InflectPtCnt ++;
			}
		
		}

			if( Cnt!=0 )
			{
			 ValArr[0] /= float(Cnt); ValArr[1] /= float(Cnt);
			 ValArr[2] /= float(Cnt); ValArr[3] /= float(Cnt);
			}
          
			PathLb PInfo;
			PInfo.Label = i;
		    PInfo.VoteVec[0] = ValArr[0]; PInfo.VoteVec[1] = ValArr[1];
			PInfo.VoteVec[2] = ValArr[2]; PInfo.VoteVec[3] = ValArr[3];

            PathLbVec.push_back(PInfo);
	}

	
}

int	MatchSquare::SmoothTangentDirections( vector<int>&PtIdxVec, int nIter=2)
{
	//return 1;
	int nPnt = PtIdxVec.size();
	
	if( nPnt <=2 ) 	 return -1;
	
	vector<St_IdPoint> tempPnt;
	tempPnt.resize(nPnt);
	for( int n=0; n<nIter; n++ )
	{
		int i;
		for( i=1; i<nPnt-1; i++ )
		{ int idx1,idx2,idx3;
		
		idx1 = PtIdxVec[i-1]  -1;
		idx2 = PtIdxVec[i]    -1;
		idx3 = PtIdxVec[i+1]  -1;
		
		tempPnt[i] = PtVec[idx2];
		
		tempPnt[i].x = 0.25 * PtVec[idx1].x + 0.5*PtVec[idx2].x + 0.25*PtVec[idx3].x;
		tempPnt[i].y = 0.25 * PtVec[idx1].y + 0.5*PtVec[idx2].y + 0.25*PtVec[idx3].y;
		}
		
		for( i=1; i<nPnt-1; i++ )
		{
			int Idx = PtIdxVec[i]-1;
			PtVec[Idx] = tempPnt[i];
			
		}
	}
	
	return 0;
}

void MatchSquare::AddLink(int Id1,int Id2)
{
	if(LinkIdx[2*Id1]==-1) 
		LinkIdx[2*Id1] = Id2;
	else
	{
		if(LinkIdx[2*Id1+1]==-1)
			LinkIdx[2*Id1+1] = Id2;
	}
}

void MatchSquare::GetPath()
{
	
	LinkIdx.resize((PtVec.size()+1)*2);
	fill(LinkIdx.begin(),LinkIdx.end(),-1);
	int i;
	Loopi(LineVec.size())
	{
		int id1,id2;
		id1 = LineVec[i].Pt1.newID;
		id2 = LineVec[i].Pt2.newID;
		AddLink(id1,id2); 
		AddLink(id2,id1);
	}
	//----------------------start to link--------------------------------
	map<int,int> PtMark;
	
	LinkPath.clear();
	
	Loopi(PtVec.size())
	{
		int PtIdx = i+1;

		if(PtMark[PtIdx])continue;//这个点已经经过了
		{
			vector<int> linkvec;
			linkvec.clear();
			//linkvec.push_back(PtIdx);
			
			while(1)
			{
				int LPt1,LPt2;
				
				LPt1 = LinkIdx[2*PtIdx];
				LPt2 = LinkIdx[2*PtIdx+1];
				
				if( PtMark[LPt1] && PtMark[LPt2])
					break;
                 
				if( LPt1==-1 || LPt2==-1)
				{   PtMark[PtIdx] =1;
				    St_IdPoint Pt;
				    Pt = PtVec[PtIdx-1];
					ProblemPts.push_back(Pt);
					linkvec.push_back(PtIdx);
					break;
				}

				if( (LPt1!=-1) && (PtMark[LPt1] ==0))
				{
					PtMark[LPt1] = 1;
					linkvec.push_back(LPt1);
					PtIdx = LPt1;
				}
				else
				{   
				   if(LPt2!=-1)
				   {
					 PtMark[LPt2] = 1;
					 linkvec.push_back(LPt2);
					 PtIdx = LPt2;
					}
				}

				 
			}
			if(linkvec.size()>5)
			LinkPath.push_back(linkvec);
		}
	}
}
const int MatchSquare::m_LineTable[16][4] =
{
	{-1,-1,-1,-1}, { 2, 3,-1,-1}, { 3, 4,-1,-1}, { 2, 4,-1,-1},
	{ 1, 4,-1,-1}, { 1, 2, 3, 4}, { 1, 3,-1,-1}, { 1, 2,-1,-1},
	{ 1, 2,-1,-1}, { 1, 3,-1,-1}, { 1, 4, 2, 3}, { 1, 4,-1,-1},
	{ 2, 4,-1,-1}, { 3, 4,-1,-1}, { 2, 3,-1,-1}, {-1,-1,-1,-1}	
};

void MatchSquare::RenamePointId()
{
	int NewId = 1;
	map<int,int> NewIdMap;
	
	int i;
	Loopi(LineVec.size())
	{
		int Id ;
		Id = LineVec[i].Pt1.oldId;
		if(NewIdMap[Id]==0)
		{
			NewIdMap[Id] = NewId;
			LineVec[i].Pt1.newID = NewId;
			PtVec.push_back(LineVec[i].Pt1);
			NewId++;
		}
		else
			LineVec[i].Pt1.newID = NewIdMap[Id];
		Id = LineVec[i].Pt2.oldId;
		if(NewIdMap[Id]==0)
		{
			NewIdMap[Id] = NewId;
			LineVec[i].Pt2.newID = NewId;
			PtVec.push_back(LineVec[i].Pt2);
			NewId++;
		}
		else
			LineVec[i].Pt2.newID = NewIdMap[Id];
	}
	//printf("%i,%i\n", PtVec.size(),LineVec.size());
}

int MatchSquare::GetLinePattern(int v1,int v2,int v3,int v4)
{
 	return (v1<<3) + (v4<<2) + (v3<<1) + (v2);
}

int MatchSquare::GetPointIdx(int x,int y,int edgnum,St_IdPoint&Pt)
{
	switch (edgnum)
	{
	case 1:
		Pt.x = float(x) +0.5;
		Pt.y = float(y);
		return y*(2*PicWidth) + 2*x;
	case 2:
		Pt.x = float(x);
		Pt.y = float(y) +0.5;
        return y*(2*PicWidth) + 2*x +1;
	case 3:
		Pt.x = float(x) +0.5;
		Pt.y = float(y) +1.0;
		return (y+1)*(2*PicWidth) + 2*(x);
	case 4:
		Pt.x = float(x) +1.0;
		Pt.y = float(y) +0.5;
        return (y)*(2*PicWidth)  + 2*(x+1) +1;
	}
}

void MatchSquare::GenerateSquareLinesMultiObj(vector<int>&ImgVec,int ww,int hh,int ObjNum)
{int i,j,t;
PicWidth  = ww;
PicHeight = hh;
vector<int> ObjAreas;
ObjAreas.resize(ObjNum);

Loopi( PicWidth )
{
	ImgVec[i                           ] = -1; 
	ImgVec[(PicHeight-1) * PicWidth + i] = -1;
}

Loopi(PicHeight)
{
	ImgVec[i*PicWidth                ] = -1;
	ImgVec[i*PicWidth +  (PicWidth-1)] = -1;
}

Loopj(PicHeight-1)
Loopi(PicWidth-1)
{
	int v1, v2, v3, v4;
	int Base = j * PicWidth;
	
	v1 = ImgVec[Base + i ];
	v2 = ImgVec[Base + i + PicWidth ];
	v3 = ImgVec[Base + i + PicWidth + 1 ];
	v4 = ImgVec[Base + i + 1 ];
	
	
	if((v1==v2)&&(v2==v3)&&(v3==v4))continue;
	
	int KindNum; //how much color
	KindNum = 0;
	int Label = -1;   
	if(v1!=-1)
	{	Label = v1;
		if(v2!= v1 &&  (v2!=-1)) ImgVec[Base + i + PicWidth ]     = -1;
		if(v3!= v1 &&  (v3!=-1)) ImgVec[Base + i + PicWidth + 1 ] = -1;
		if(v4!= v1 &&  (v4!=-1)) ImgVec[Base + i + 1 ]            = -1;    
	}
	 if(v4!=-1)
	  {
		 Label = v4;
		 if(v2!= v4 &&  (v2!=-1)) ImgVec[Base + i + PicWidth ]     = -1;
		 if(v3!= v4 &&  (v3!=-1)) ImgVec[Base + i + PicWidth + 1 ] = -1;		  
	  }

	 if(v2!=-1)
	 {	 Label = v2;
	     if(v3!= v2 &&  (v3!=-1)) ImgVec[Base + i + PicWidth + 1 ] = -1;   
	 }
	 if(v3!=-1)
		 Label = v3;

	  v1 = (v1 == Label) ? 1:0;
	  v2 = (v2 == Label) ? 1:0;
	  v3 = (v3 == Label) ? 1:0;
 	  v4 = (v4 == Label) ? 1:0;
     
	  int PatIdx = GetLinePattern(v1,v2,v3,v4);
	  
	  for(t=0;t<4;t+=2)
	  {
		  if(m_LineTable[PatIdx][t]==-1)
			  break;
		  St_IdPoint Pt;
          
		  St_Line  LineSeg;
		  Pt.Label = Label;
		  //LineSeg.Label = 
		  Pt.oldId    = GetPointIdx(i,j,m_LineTable[PatIdx][t ], Pt);
		  Pt.xi =i;Pt.yi =j; Pt.PatIdx = m_LineTable[PatIdx][t ];
		  LineSeg.Pt1 = Pt;  
		  
		  Pt.oldId    = GetPointIdx(i,j,m_LineTable[PatIdx][t+1], Pt);
		  Pt.xi =i;Pt.yi =j; Pt.PatIdx = m_LineTable[PatIdx][t+1];
		  LineSeg.Pt2 = Pt;
		  
		  LineVec.push_back(LineSeg);
	}
  }

}

void MatchSquare::GenerateSquareLines()
{int i,j,t;

ObjArea =0;
Loopj(GPic.Height)
Loopi(GPic.Width)
{
	int v1,v2,v3,v4;

	v1 = *get_pix_color(GPic, i  ,  j)>128? 0:1;
	v2 = *get_pix_color(GPic, i  ,j+1)>128? 0:1;
	v3 = *get_pix_color(GPic, i+1,j+1)>128? 0:1;
	v4 = *get_pix_color(GPic, i+1,  j)>128? 0:1;

	if(v1)
		ObjArea++;
	
	int vv[4];
	vv[0] =  *get_pix_color(GPic, i, j);
	vv[1] =  *get_pix_color(GPic, i, j+1);
	vv[2] =  *get_pix_color(GPic, i+1, j+1);
	vv[3] =  *get_pix_color(GPic, i+1, j);
	
	int PatIdx = GetLinePattern(v1,v2,v3,v4);
	
	for(t=0;t<4;t+=2)
	{
		if(m_LineTable[PatIdx][t]==-1)
			break;
		St_IdPoint Pt;
		St_Line  LineSeg;
		Pt.oldId    = GetPointIdx(i,j,m_LineTable[PatIdx][t ], Pt);
		Pt.xi =i;Pt.yi =j; Pt.PatIdx = m_LineTable[PatIdx][t ];
		LineSeg.Pt1 = Pt;  
		
		Pt.oldId    = GetPointIdx(i,j,m_LineTable[PatIdx][t+1], Pt);
		Pt.xi =i;Pt.yi =j; Pt.PatIdx = m_LineTable[PatIdx][t+1];
		LineSeg.Pt2 = Pt;
	
		LineVec.push_back(LineSeg);
	}
}

}



#endif