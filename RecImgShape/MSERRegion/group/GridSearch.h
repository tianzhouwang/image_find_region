 

#ifndef Obj_Grid_Search_H
#define Obj_Grid_Search_H
#include "mregion.h"
#define GRID_SIZE 9
#define BOUND(x,a,b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))


// \brief ������������ṹ
class GridSearch
{
public:
  GridSearch(int W,int H);
  GridSearch(int W,int H,vector<Region>&RegionVec);
  vector<Region*>  ObjGrid[GRID_SIZE+1][GRID_SIZE+1];  
  int ImgWidth,ImgHeight,GridSzX,GridSzY;
  int MaxGridSize;
   
  void AddObj(Region&Obj);
  int  FindClosePts(Region&Obj);
  Region* AdjacentObjs[2000];

};

/// \brief  �����ٽ�����s
/// \param  Obj �����ҵ�Ԫ��
/// \return �ٽ�����ĸ���, ���ٽ������ָ�뽫������AdjacentObjs������ 
int GridSearch::FindClosePts(Region&Obj)
{
  int ObjCnt = 0;
  int i,j,t;

  int IdxX, IdxY; 
  IdxX = Obj.x / GridSzX;
  IdxY = Obj.y / GridSzY;

  for(i=IdxX-1; i<=IdxX+1; i++)
	  for(j=IdxY-1; j<=IdxY+1; j++)
	  {
       if(i<0)continue; if(j<0)continue;
       if(i>= GRID_SIZE)continue;
	   if(j>= GRID_SIZE)continue;
       
	   for(t=0;t<ObjGrid[i][j].size();t++)
	   {if(ObjCnt>=2000)break;
       AdjacentObjs[ObjCnt] =  ObjGrid[i][j][t];
	   ObjCnt++;}
	  }

  return  ObjCnt;
}

/// \brief  ��ʼ�����������
/// \param  W ͼ��Ŀ�
/// \param  H ͼ��ĸ�
/// \return ��
GridSearch::GridSearch(int W,int H)
{
 ImgWidth  = W; ImgHeight = H;
  GridSzX  = ImgWidth  / GRID_SIZE;
  GridSzY  = ImgHeight / GRID_SIZE;
  MaxGridSize = GRID_SIZE;
}

/// \brief  ���һ�����嵽������Ҷ���
/// \param  Obj ��ӵ�����
/// \return ��
void GridSearch::AddObj(Region&Obj)
{
  int IdxX, IdxY;

  GridSzX = ImgWidth  / GRID_SIZE;
  GridSzY = ImgHeight / GRID_SIZE;
  
  IdxX = Obj.x / GridSzX;
  IdxY = Obj.y / GridSzY;

  ObjGrid[IdxX][IdxY].push_back(&Obj);
}

/// \brief  ��ʼ��������Ҷ���,��������������Ԫ����ӵ�������
/// \param  W ͼ��Ŀ�
/// \param  H ͼ��ĸ�
/// \param  RegionVec ������������
/// \return ��
GridSearch::GridSearch(int W,int H,vector<Region>&RegionVec)
{   int i;
	ImgWidth  = W; ImgHeight = H;
	GridSzX  = ImgWidth  / GRID_SIZE;
	GridSzY  = ImgHeight / GRID_SIZE;
    MaxGridSize = GRID_SIZE;
	for(i=0;i<RegionVec.size();i++)
	{
      AddObj(RegionVec[i]);
	}
};

#endif