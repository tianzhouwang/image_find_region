 
#ifndef Mser_txt_Region_H
#define Mser_txt_Region_H
//#include "ObjContour.h"
#include "mser.h"
#include "regionfill.h"
#include "regionslt.h"
#include "dumpMesr.h"
#include "rgreduction.h"
#include "region.h"
//#include "kutil.h"
#define DUMP_MSER_DEBUG      0
#define DUMP_MSER_TREE_DEBUG 0
 
/// \brief  得到区域正反两面的MSER区域
void GetPosNegMesrTxtRegion(C24BitMap&CPic,vector<Region>&RegionB,vector<Region>&RegionW);
 
void TransGray(C24BitMap&Sr,C256BitMap&Ds)
{
	IJT;
	Ds.FormatF(Sr.Width,Sr.Height);
	LoopPic_ij(Ds){
		C24PixVal Tmp=get_pix_color(Sr,i,j);
		*get_pix_color(Ds,i,j)=
			float(*Tmp.r)*.30+float(*Tmp.g)*.59+float(*Tmp.b)*.11;}
}

//MSER 文字区域挑选结构体
class MserRegion
{
public:
  void InitPicture(C24BitMap&CPic);

  vector<Region> RegionW; //白色的字
  vector<Region> RegionB; //黑色的字
  int PicWidth,PicHeight;
  
  vector<uint8_t> grayscale_w;
  vector<uint8_t> grayscale_b;
  C256BitMap GPic_w, GPic_b;
  //void CompleteRegion(FillFindRegion&Region,C256BitMap&Pic,int thrd,vector<Point> &PtVec);
  void CompleteRegion(ConnectedComponents& ConnT,FillFindRegion&RgFill,C256BitMap&Pic,Region&Rg);
  void CompPictures(C256BitMap&GPic, Region*rootRegion,vector<Region>&ReginVec);
  void GenerateMser(vector<uint8_t>&grayscale, C256BitMap&GPic, vector<Region> &RgVec);
  double MinSizeRatio, MaxSizeRatio;
  int    NoiseRegionSz;
  MserRegion();
  
  void GenerateMserGrayPic(C256BitMap&GPic, vector<Region> &RgVec);
};


/// \brief  得到区域正反两面的MSER区域
/// \param  CPic         待分析的彩色图像
/// \param  RegionB      白底黑字的稳定区域
/// \param  RegionW      黑底白字的稳定区域
/// \return 无
void GetPosNegMesrTxtRegion(C24BitMap&CPic,vector<Region>&RegionB,vector<Region>&RegionW)
{
	MserRegion MserAlg;
	MserAlg.InitPicture(CPic);


	MserAlg.GenerateMser(MserAlg.grayscale_b,MserAlg.GPic_b, RegionB);//MserAlg.RegionB);
    //RegionB = MserAlg.RegionB ;
	 

    MserAlg.GenerateMser(MserAlg.grayscale_w,MserAlg.GPic_w, RegionW);//MserAlg.RegionW);
	//RegionW = MserAlg.RegionW ;
}


/// \brief  MserRegion初始化函数
/// \return 无
MserRegion::MserRegion()
{
  MinSizeRatio  = 0.0001;
  MaxSizeRatio  = 0.3;
  NoiseRegionSz = 25;
}


/// \brief  得到一个通道灰度图像的MSER稳定区域
/// \param  grayscale     灰度的裸数据矩阵
/// \param  GPic          原始灰度图像
/// \param  RgVec         稳定区域元素矩阵
/// \return 无
void MserRegion::GenerateMser(vector<uint8_t>&grayscale, C256BitMap&GPic, vector<Region> &RgVec )
{
	DumpMserRg Dmp;

	if(DUMP_MSER_DEBUG)
	{
     Dmp.OPic.FormatF(PicWidth,PicHeight);
	 Dmp.GPic =   GPic_b;//;GPic; //<------注意调试时候选择合适的+ or -	 
	}
//---------------------------------------------
    int i;
    MSER mser8(false);
     //MSER mser8(true);//八邻域

	double minsize,maxsize;
	minsize = MinSizeRatio  * double(PicWidth * PicHeight);
	maxsize = MaxSizeRatio  * double(PicWidth * PicHeight);  
	
	/*Region::GRegionMinSize = minsize;
	Region::GRegionMaxSize = maxsize;
	Region::PicWidth  = PicWidth;
	Region::PicHeight = PicHeight;*/

	mser8.GenerateMSER(&grayscale[0], PicWidth, PicHeight);    


    
	for(i=0;i<mser8.pool_.size();i++)
	{
		mser8.pool_[i].r_id = i;
		//SetRegionAreaConstrain(int GRegionMinSize_,int GRegionMaxSize_,int PicWidth_,int PicHeight_);
		mser8.pool_[i].SetRegionAreaConstrain(minsize,maxsize,PicWidth,PicHeight);
	}
	
    Region*RgPt ;
	RgPt = mser8.TopRegion;
	
    RgPt->getvariation(mser8.delta_);
    Region*NewRoot ; 

    if(DUMP_MSER_DEBUG)Dmp.PrintfRegion("res\\lintree0.txt",RgPt);

	NewRoot = cut_noiseNode(RgPt, NoiseRegionSz);
    
	if(DUMP_MSER_DEBUG)Dmp.PrintfRegion("res\\lintree0_noisecut.txt",NewRoot);

	NewRoot->update_childnum();
    NewRoot = linear_reduction(NewRoot);
	NewRoot->update_childnum();
	NewRoot =  linear_reduction(NewRoot);
	
	if(DUMP_MSER_TREE_DEBUG)
	{
     Dmp.PrintfRegion("res\\lintree.txt",NewRoot);
	 Dmp.DumpPictures("res\\lin\\",NewRoot);
     
	}

    vector<Region*> SubNodes = tree_accumulation(NewRoot);
	InsertChildNodes(NewRoot, SubNodes);
    
	if(DUMP_MSER_TREE_DEBUG)Dmp.PrintfRegion("res\\tree_acc.txt",NewRoot);

	CompPictures( GPic, NewRoot, RgVec);
	
	if(DUMP_MSER_DEBUG)Dmp.DumpPictures("res\\dmp\\",RgVec);
}

/// \brief     对MSER操作进行初始化
/// 分配正反两面MSER区域正反两面所需要的数据 
/// \param  CPic       待分析的彩图
/// \return 无
void MserRegion::InitPicture(C24BitMap&CPic)
{
  int i,j;
  PicWidth  = CPic.Width ;
  PicHeight = CPic.Height;

  TransGray(CPic, GPic_b);
  
  if(DUMP_MSER_DEBUG) GPic_b.Save("gray.bmp");
         
  GPic_w.FormatF(PicWidth, PicHeight);
  
  Loopi(GPic_w.LineWidth * GPic_w.Height)
	  GPic_w.Buffer[i] = 255 - GPic_b.Buffer[i];

  grayscale_w.resize(CPic.Width * CPic.Height);
  grayscale_b.resize(CPic.Width * CPic.Height);

  for (j=0;j<PicHeight;j++)
  {
	  //int sum = 0; 
	  
	  for (i = 0; i <PicWidth ; ++i) 
	  {	
		  grayscale_w[j*PicWidth+i]   =   *get_pix_color(GPic_w,i, j) ;//sum / depth;
		  grayscale_b[j*PicWidth+i]   =   *get_pix_color(GPic_b,i, j) ;
	  }
	}
}


/// \brief     对灰度图像生成连通域
///  
/// \param  GPic      待分析的灰度图
/// \param  RgVec     生成的MSER区域
/// \return 无
void MserRegion::GenerateMserGrayPic(C256BitMap&GPic, vector<Region> &RgVec)
{ 
	int i,j;
	PicWidth  = GPic.Width ;
	PicHeight = GPic.Height;
	
	vector<uint8_t> grayscale_pic;
	grayscale_pic.resize(GPic.Width * GPic.Height);
	
	for (j=0;j<PicHeight;j++)
	{ 	
		for (i = 0; i <PicWidth ; ++i) 
		{	
		    grayscale_pic[j*PicWidth+i] = *get_pix_color(GPic , i, j) ;//sum / depth;	 
		}
	}

	GenerateMser(grayscale_pic, GPic,  RgVec);
}
	 
/// \brief   遍历MSER结构树 
///          将一定阈值范围下的稳定区域保存到区域向量内
/// \param   GPic         原始图像
/// \param   rootRegion   MSER结构树根节点   
/// \param   ReginVec     保存完的MSER结果区域
/// \return  无
void MserRegion::CompPictures(C256BitMap&GPic, Region*rootRegion,vector<Region>&ReginVec)
{   
	int i; 
    FillFindRegion RegionFill(GPic);
	rootRegion->SaveMserRegion(ReginVec, 0.3);
    Loopi(ReginVec.size())ReginVec[i].IsValid = 1;
	//MarkInvalidRegionByCrossValidate(ReginVec);

	ConnectedComponents Cp;
// 	Cp.scan_width  = GPic.LineWidth;
//  Cp.scan_height = GPic.Height;
// 	Cp.alloc_space();

    Loopi(ReginVec.size())
	{
		if(ReginVec[i].IsValid)
		  CompleteRegion(Cp, RegionFill,GPic,ReginVec[i]);
		else
		  ReginVec[i].PtVec.clear();
	}
}


/// \brief   MSER合并完后可能有区域元素不全 
///          通过区域种子点补全mser区域
/// \param   RgFill      种子填充算法结构
/// \param   Pic         图像原图
/// \param   Rg          找到的区域
/// \return 无
void MserRegion::CompleteRegion(ConnectedComponents& ConnT,FillFindRegion& RgFill,C256BitMap&Pic,Region&Rg)
//void MserRegion::CompleteRegion(FillFindRegion&Region,C256BitMap&Pic,int thrd,vector<Point> &PtVec)
{
	int i,j;
	RgFill.fill_threshold = Rg.level_;
   
	for (i=Rg.left;i<=Rg.right;i++)
	{
       for(j=Rg.top;j<=Rg.bottom;j++)
		   RgFill.Matrix[i][j] = 0;
	}

	
	if(Rg.PtVec.size()>0)
	{
		int x, y;
		x = Rg.PtVec[0].x;
		y = Rg.PtVec[0].y;
		
		RgFill.FillSeed(x,y);
		
		if(RgFill.PtVec.size()>Rg.PtVec.size())
			Rg.PtVec = RgFill.PtVec;

	}
	return;
//---------------------------------------------
 /*	ConnT.left          = Rg.left;
	ConnT.right         = Rg.right;
	ConnT.top           = Rg.top;
	ConnT.bottom        = Rg.bottom;
	ConnT.pix_threshold = Rg.level_;

	//ConnT.label_image(1);
	int obnum=ConnT.relabel_image(); 
    ConnT.GetResultGroup(); 
 
Cp.left  = 0;
 Cp.right = ww-1;
 Cp.top   = 0;
 Cp.bottom= hh-1; 
 Cp.pix_threshold =128;

  

 Cp.alloc_space();
 Cp.label_image(Pic.Buffer,1);
 int obnum=Cp.relabel_image(); 
  Cp.GetResultGroup(); */
} 

#endif

