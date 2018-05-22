
#ifndef _FillRg_CDIB_H
#define _FillRg_CDIB_H
#include <vector>
#include <stdlib.h>
#include <string>
#include <stack>
//#include "mser.h"
#include "c24bitmap.h"
#include "c256bitmap.h"
#include "mregion.h"

using namespace std;



class FillFindRegion
{
public:
	C256BitMap Pic;
	
	//bool IfRot;
	FillFindRegion(C256BitMap&Picture);
	~FillFindRegion();
	//PointRng FindRegion(int x,int y);
	bool IfSimilar(int pixVal);
	bool FillSeed(int x,int y);
	
	//bool IsValid(int x,int y);
	bool checkFill(int x,int y);
	
	stack<RPoint> CalCuStack;
	
	vector<RPoint> PtVec;
	int fill_threshold;  //<----填充的阈值
	
	int Left,Right,Top,Bottom;
	
	vector<char*>Matrix;
	vector<char>  Buffer;
	void DrawMatrix(char*filename);
};

void FillFindRegion::DrawMatrix(char*filename)
{
	
}

/*bool FillFindRegion::IsValid(int x,int y)
{
if(x<0||x>=Pic.Width)return false;
if(y<0||y>=Pic.Height)return false;
return true;
}*/

//----------------------------------------------------对所找区域进行种子填充
bool FillFindRegion::FillSeed(const int X,const int Y)
{
	PtVec.clear();
	while(!CalCuStack.empty())CalCuStack.pop();
	RPoint TmpPt;
	TmpPt.x=X;
	TmpPt.y=Y;
	int x,y;
	Matrix[X][Y]=true;
	CalCuStack.push(TmpPt);
	PtVec.push_back(TmpPt);
	while(!CalCuStack.empty())
	{
		TmpPt=CalCuStack.top();
		CalCuStack.pop();
		x=TmpPt.x;
		y=TmpPt.y;
		if(checkFill(x-1,y))
			if(!Matrix[x-1][y])
			{
				if(Left>(x-1))Left=x-1;
				Matrix[x-1][y]=true;
				TmpPt.x=x-1;
				TmpPt.y=y;
				CalCuStack.push(TmpPt);
				PtVec.push_back(TmpPt);
			}
			if(checkFill(x+1,y))
				if(!Matrix[x+1][y])
				{
					if(Right<(x+1))Right=x+1;
					Matrix[x+1][y]=true;
					TmpPt.x=x+1;
					TmpPt.y=y;
					CalCuStack.push(TmpPt);
					PtVec.push_back(TmpPt);
				}
				if(checkFill(x,y-1))
					if(!Matrix[x][y-1])
					{
						if(Top>(y-1))Top=y-1;
						Matrix[x][y-1]=true;
						TmpPt.x=x;
						TmpPt.y=y-1;
						CalCuStack.push(TmpPt);
						PtVec.push_back(TmpPt);
					}
					if(checkFill(x,y+1))
						if(!Matrix[x][y+1])
						{
							if(Bottom<(y+1))Bottom=y+1;
							Matrix[x][y+1]=true;
							TmpPt.x=x;
							TmpPt.y=y+1;
							CalCuStack.push(TmpPt);
							PtVec.push_back(TmpPt);
						}
	}
	return true;
}


FillFindRegion::FillFindRegion(C256BitMap&Picture) 
{
	int i; 
	Pic = Picture;
	Buffer.resize(Pic.Width*Pic.Height);
	Matrix.resize(Pic.Width);
	
	Loopi(Pic.Width*Pic.Height) Buffer[i]=false;
	Loopi(Pic.Width) Matrix[i] =&Buffer[Pic.Height*i];
}

FillFindRegion::~FillFindRegion()
{
	
}

bool FillFindRegion::checkFill(int x,int y)
{
	if(x<0||x>=Pic.Width)return false;
	if(y<0||y>=Pic.Height)return false;
	if(Matrix[x][y]==true)return false;
	
	int PixV = *get_pix_color(Pic,x,y);
	if(!IfSimilar(PixV))return false;
	return true;
}



bool FillFindRegion::IfSimilar(int pixVal)
{
	if(pixVal <= fill_threshold)  
		return true;
	else
		return false;
}



#endif
