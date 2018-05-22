 
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
 
/// \brief  �õ��������������MSER����
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

//MSER ����������ѡ�ṹ��
class MserRegion
{
public:
  void InitPicture(C24BitMap&CPic);

  vector<Region> RegionW; //��ɫ����
  vector<Region> RegionB; //��ɫ����
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


/// \brief  �õ��������������MSER����
/// \param  CPic         �������Ĳ�ɫͼ��
/// \param  RegionB      �׵׺��ֵ��ȶ�����
/// \param  RegionW      �ڵװ��ֵ��ȶ�����
/// \return ��
void GetPosNegMesrTxtRegion(C24BitMap&CPic,vector<Region>&RegionB,vector<Region>&RegionW)
{
	MserRegion MserAlg;
	MserAlg.InitPicture(CPic);


	MserAlg.GenerateMser(MserAlg.grayscale_b,MserAlg.GPic_b, RegionB);//MserAlg.RegionB);
    //RegionB = MserAlg.RegionB ;
	 

    MserAlg.GenerateMser(MserAlg.grayscale_w,MserAlg.GPic_w, RegionW);//MserAlg.RegionW);
	//RegionW = MserAlg.RegionW ;
}


/// \brief  MserRegion��ʼ������
/// \return ��
MserRegion::MserRegion()
{
  MinSizeRatio  = 0.0001;
  MaxSizeRatio  = 0.3;
  NoiseRegionSz = 25;
}


/// \brief  �õ�һ��ͨ���Ҷ�ͼ���MSER�ȶ�����
/// \param  grayscale     �Ҷȵ������ݾ���
/// \param  GPic          ԭʼ�Ҷ�ͼ��
/// \param  RgVec         �ȶ�����Ԫ�ؾ���
/// \return ��
void MserRegion::GenerateMser(vector<uint8_t>&grayscale, C256BitMap&GPic, vector<Region> &RgVec )
{
	DumpMserRg Dmp;

	if(DUMP_MSER_DEBUG)
	{
     Dmp.OPic.FormatF(PicWidth,PicHeight);
	 Dmp.GPic =   GPic_b;//;GPic; //<------ע�����ʱ��ѡ����ʵ�+ or -	 
	}
//---------------------------------------------
    int i;
    MSER mser8(false);
     //MSER mser8(true);//������

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

/// \brief     ��MSER�������г�ʼ��
/// ������������MSER����������������Ҫ������ 
/// \param  CPic       �������Ĳ�ͼ
/// \return ��
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


/// \brief     �ԻҶ�ͼ��������ͨ��
///  
/// \param  GPic      �������ĻҶ�ͼ
/// \param  RgVec     ���ɵ�MSER����
/// \return ��
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
	 
/// \brief   ����MSER�ṹ�� 
///          ��һ����ֵ��Χ�µ��ȶ����򱣴浽����������
/// \param   GPic         ԭʼͼ��
/// \param   rootRegion   MSER�ṹ�����ڵ�   
/// \param   ReginVec     �������MSER�������
/// \return  ��
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


/// \brief   MSER�ϲ�������������Ԫ�ز�ȫ 
///          ͨ���������ӵ㲹ȫmser����
/// \param   RgFill      ��������㷨�ṹ
/// \param   Pic         ͼ��ԭͼ
/// \param   Rg          �ҵ�������
/// \return ��
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

