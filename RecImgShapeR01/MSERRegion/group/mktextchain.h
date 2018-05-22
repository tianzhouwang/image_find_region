 

#ifndef make_text_line_h
#define make_text_line_h
#include <vector>
#include <algorithm>
//#include "mregion.h"
using namespace std;
typedef unsigned char byte;
//==============================================================================
#define MIN(a,b)        (((a)<(b))   ? (a)  : (b)) 
#define MAX(a,b)        (((a)>(b))   ? (a)  : (b))
#define VAL_IN_RANGE(a,b,c)  ((a) >= (b) && (a) <= (c)) 
#define Loopi(k_) for(i=0;i<(k_);i++)
#define DEBUG_OUT 0

//�м亯��
//ImgRECT MergRectBox(ImgRECT&Rect1,ImgRECT&Rect2)//�ϲ���Ӿ���



/// \brief  ��ȡRegion�����е�����ͳ�����Ϣ,������Ӿ����������Է���   
/// \param  ObjVec    �������Ԫ������
/// \param  RectVec   �������Ӿ�������
/// \return ��
//===============================================================================
void transObjsToRects(vector<Region>&ObjVec,vector<ImgRECT>&RectVec)
{
	int i; RectVec.clear();
	Loopi(ObjVec.size())
	{
		ImgRECT rect;
		rect.right  = ObjVec[i].right  ;                                         
		rect.top    = ObjVec[i].top    ;                                         
		rect.bottom = ObjVec[i].bottom ;                                         
		rect.left   = ObjVec[i].left   ;                                         
		rect.width  = ObjVec[i].right  - ObjVec[i].left + 1;                                         
		rect.height = ObjVec[i].bottom - ObjVec[i].top  + 1;                                         
		RectVec.push_back(rect);
	}
}
//===============================================================================
#define CHR_WH_RATIO(char_box) ( float(char_box.height) / float(char_box.width) )
#define CHR_BOX_GAP(rect1  ,rect2)    ((rect1).left   - (rect2).right )
#define CHR_BOX_STEP(rect1 ,rect2)    ((rect1).right  - (rect2).left  )

#define CHR_BOX_GAP_V(rect1 ,rect2)   ((rect1).top    - (rect2).bottom)
#define CHR_BOX_STEP_V(rect1,rect2)   ((rect1).bottom - (rect2).top   )

/// \brief  ��������������Ϣ���ƾ��εĳ�����Ϣ  
/// \param  rect  ��������Ϣ�ľ��ο�
/// \return ��
void getRectWH(ImgRECT&rect) 
{ 
	rect.width  = rect.right  - rect.left + 1;
	rect.height = rect.bottom - rect.top  + 1;
}

/// \brief �ϲ���ѡ��Ӿ���
/// \param Rect1   ���ϲ�����1
/// \param Rect2   ���ϲ�����2
/// \return  �ϲ��Ժ�ľ���
ImgRECT MergRectBox(ImgRECT&Rect1,ImgRECT&Rect2)
{
	ImgRECT rRect; 
	rRect.right  = MAX( Rect1.right , Rect2.right );
	rRect.top    = MIN( Rect1.top   , Rect2.top   );
	rRect.bottom = MAX( Rect1.bottom, Rect2.bottom);
	rRect.left   = MIN( Rect1.left  , Rect2.left  );
	getRectWH(rRect) ;
	return rRect;
}

// \brief ���εĺ���������������
bool LessByX(const ImgRECT&ob1,const ImgRECT&ob2) { return ob1.left < ob2.left;}
bool LessByY(const ImgRECT&ob1,const ImgRECT&ob2) { return ob1.top  < ob2.top ;}

// \brief ��ѡ�����������������������м����з���Ϣ
class GetTextLineComptChain
{
public:   
    GetTextLineComptChain();//���ʼ��
	
	
	//------------------------------------------------------------------
	//�ϲ������ص�����
	int  CombOverlapRectX(vector<ImgRECT>&InCompts,vector<ImgRECT>&OutRects,int nOverlap);
	//�ϲ������ص�����
	int  CombOverlapRectY(vector<ImgRECT>&InCompts,vector<ImgRECT>&OutRects,int nOverlap);
	//�õ��߶ȵ���ֵ
	int  GetMedianHeight(vector<ImgRECT>&RectVec);
	//�õ���ȵ���ֵ
	int  GetMedianWidth(vector<ImgRECT>&RectVec);
	//�ı���ˮƽ�з�
    int  GetHorizontalTextChain(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects,int MHeight);
	int  GetTextLineHorizontal(int chrHeight,vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects);
	
	//�õ��������з�
	int  GetVerticalTextChain(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects,int Mwidth);
	int  GetTextLineVertical(int chrWidth,vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects);
	//GetValidChar(ObjComb, ObjH, Rects);
	
	bool  StrokeNumberY(ImgRECT Rect, int StrokeNumTh, int StrokeSumTh, int gap);
	bool  StrokeNumberX(ImgRECT Rect, int StrokeNumTh, int StrokeSumTh, int gap);
	//-------------------------------------------------------------------------
	void MarkValidRectH1(float hwRatio, int chrHeight,vector<ImgRECT>&RectVec,
		vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum);
	void MarkValidRectH2(float hwRatio, int chrHeight,vector<ImgRECT>&RectVec,
		vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum);
	void MarkValidRectEN( float hwRatio, int chrHeight,vector<ImgRECT>&RectVec,
		vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum);
	//-------------------------------------------------------------------------	
	void MarkValidRectH_overlap(float hwRatio, int MeanWidth,vector<ImgRECT>&RectVec,
		vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum);
	void MarkValidRectV_overlap(float hwRatio, int MeanHeight,vector<ImgRECT>&RectVec,
		vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum);
	
	void GetHwRatio(vector<ImgRECT>&RectVec, float &hwRatio);
	int  SimpleHorizontalTextLineCut( vector<ImgRECT>&RectVec,
								vector<ImgRECT>&CombRectVec,vector<ImgRECT>&OutRects);
	//<-------һ���򵥵��з�ģ��
	//------------------------------------------------------------------
	//�õ��߶ȵ��м�ֵ
	
	int  GetAignTop(vector<ImgRECT> & BoundBoxVec, int chHeight);
	vector<ImgRECT> AlignObjVec;  //<------��õ�һ���߶ȵ��м�ֵ
	
	//���ư�Χ���ε�����ֿ�
    void DispBoxPositionPic(vector<ImgRECT> &TempVec, int Width, int Height,char*filenmae);
	
	int OverLayPix;
	bool ChineseOnly;
	C256BitMap Pic;
	bool FirstLineFlag;
	double AlignTop_LowR ,  AlignTop_HighR; //�߶ȶ������ֵ����Χ
	double CNchr_WH_LowR ,  CNchr_WH_HighR; //�����ַ��������ֵ����Χ 
	double MinLowChrHeightR;  //�߶ȵ���СҪ��
	double MinCharWidthR;     //�����СҪ��   
	double ChrGapRatio;      //�ַ�����Ҫ��
	double RatioTolerance0,RatioTolerance1;
};


/// \brief �õ�����������ӿ�ĺ����з�
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param RectVec       ��ʼ������ͨ�����Ӿ�������
/// \param OutRects      �������ַ������
/// \param chrHeight     �ο����ָ߶�     
/// \return  �з��Ժ�õ������ָ���
int GetTextLineComptChain::SimpleHorizontalTextLineCut(vector<ImgRECT>&RectVec,
										 vector<ImgRECT>&CombRectVec,vector<ImgRECT>&OutRects)
{
	int i;
    vector<ImgRECT> CombinedChrVec;
	sort(RectVec.begin(), RectVec.end(), LessByX); 
    //---------------------------------------------------------
	ImgRECT Rect_;
	
	if(RectVec.size()>0)
	{
		Rect_.left   = RectVec[0].left  ;
		Rect_.right  = RectVec[0].right ;
		Rect_.top    = RectVec[0].top   ;
		Rect_.bottom = RectVec[0].bottom;
	}
    
	Loopi(RectVec.size())
	{
		if(RectVec[i].left  < Rect_.left  ) Rect_.left   = RectVec[i].left  ;
		if(RectVec[i].right > Rect_.right ) Rect_.right  = RectVec[i].right ;
		if(RectVec[i].top   < Rect_.top   ) Rect_.top    = RectVec[i].top   ;
		if(RectVec[i].bottom> Rect_.bottom) Rect_.bottom = RectVec[i].bottom;
	}
	//---------------------------------------------------------
	int MHeight = Rect_.bottom - Rect_.top +1;

	CombOverlapRectX(RectVec, CombinedChrVec, OverLayPix);
    CombRectVec = CombinedChrVec;   	//DispBoxPositionPic(CombinedChr,  Pic.Width,  Pic.Height,"step0.bmp");//<---debug

	int cHeight = MHeight; 
	cHeight = GetAignTop(CombinedChrVec,cHeight);

	vector<int> chObjDone;
	chObjDone.resize(CombinedChrVec.size());
    fill(chObjDone.begin(),chObjDone.end(),0);
	
	int Num = RectVec.size();
	int nChractNum = 0;
	int WidthSum   = 0;
	OutRects.clear();
	
	float nW,nH;
	ImgRECT rRect;  
	float hwRatio = 1.0f;
	GetHwRatio( CombinedChrVec, hwRatio);
	
	MarkValidRectH1( hwRatio, cHeight, CombinedChrVec, chObjDone, OutRects, WidthSum, nChractNum);
    MarkValidRectH2( hwRatio, cHeight, CombinedChrVec, chObjDone, OutRects, WidthSum, nChractNum);

	 
	sort(OutRects.begin(),OutRects.end(), LessByX);
	
	if(nChractNum!=0)
		return WidthSum/nChractNum;
	return -1;
};

/// \brief ���ʼ����������ֵ
///     
/// \param  
/// \return  ��
GetTextLineComptChain::GetTextLineComptChain()
{
	OverLayPix        = -4;
	FirstLineFlag     = false;
	RatioTolerance0   = 0.1;
	RatioTolerance1   = 0.2;
	AlignTop_LowR     = 0.70;  AlignTop_HighR  = 1.1;
	CNchr_WH_LowR     = 0.8 ;  CNchr_WH_HighR  = 1.2;
	MinLowChrHeightR  = 0.75;
	ChrGapRatio       = 0.3;
	MinCharWidthR = 0.75;
}

/// \brief ���ư�Χ���ε��������
///   
/// �½�һ����ʱͼƬ���������Ӿ���ͼ��  
/// \param TempVec  ��Ҫ���Ƶ���Ӿ�������
/// \param Width    ������ͼ�񱳾����
/// \param Height   ������ͼ�񱳾��߶� 
/// \param filename ����ļ���
/// \return  ��
void GetTextLineComptChain::DispBoxPositionPic(vector<ImgRECT> &TempVec, int Width, int Height,char*filenmae)
{
	C24BitMap pictemp;
	pictemp.FormatF( Width,  Height);
	pictemp.ClearPic(255);
	pictemp.SetPenColor(0,0,0);
	pictemp.FillRect = 1;

	for(int i = 0; i < TempVec.size(); i++)
	{
		
		char buff[30]; sprintf(buff,"%i",i);
        MergeTxtStr(pictemp,(TempVec[i].left+TempVec[i].right)/2,
			                (TempVec[i].top +TempVec[i].bottom)/2, 15, buff, 255,0,0); 

		pictemp.Rectangle(TempVec[i].left,TempVec[i].top, TempVec[i].right,TempVec[i].bottom );
	}
	
	pictemp.Save(filenmae);
}


/// \brief ��������ֵ��Χ�ڵĸ߶Ƚ�����������߶ȵ���ֵ
///  �����ο������ָ߶�(chHeight),���з��϶���߶ȷ�Χ(AlignTop_LowR,AlignTop_HighR)��
///  �ĺ�ѡ���򽫱��浽AlignObjVec
///  ���û���ҵ����������ĺ�ѡ����,AlignObjVec����Ϊ���ҷ���chHeight
/// \param BoundBoxVec  ����ѡ�ľ�������
/// \param chHeight     �ο����ָ߶�     
/// \return  �����Ĳο��߶�
int GetTextLineComptChain::GetAignTop(vector<ImgRECT> & BoundBoxVec,int chHeight)
{
	int i;
	float nW, nH;
	
	AlignObjVec.clear();
	vector<int> HeightVec;
	HeightVec.clear();
    Loopi(BoundBoxVec.size()) 
	{
		nW = BoundBoxVec[i].width  ; nH = BoundBoxVec[i].height ;
		
		if(VAL_IN_RANGE( nH, AlignTop_LowR * chHeight, AlignTop_HighR * chHeight ))
		{
			AlignObjVec.push_back(BoundBoxVec[i]);
			HeightVec.push_back(nH);
		}

		if(VAL_IN_RANGE(CHR_WH_RATIO(BoundBoxVec[i]),
			1.0 - RatioTolerance1, 1.0 + RatioTolerance1))
		{
           HeightVec.push_back(nH);
		}
	}
	
	if(HeightVec.size()>0)
	{   
		int MdSz = HeightVec.size()/2;
		nth_element(HeightVec.begin(),HeightVec.begin()+MdSz,HeightVec.end());
		return HeightVec[MdSz];
	}
	else
		return chHeight;
}

 
/// \brief �õ��������������з�
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param RectVec     ��ʼ������ͨ�����Ӿ�������
/// \param OutRects    �������ַ������
/// \param MHeight     �ο����ָ߶�     
/// \return  �з��Ժ�õ������ָ���
int GetTextLineComptChain::GetHorizontalTextChain(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects,int MHeight)
{
	vector<ImgRECT> TempVec;
	sort(RectVec.begin(),RectVec.end(),LessByX);
	
	DispBoxPositionPic(RectVec,  Pic.Width,  Pic.Height,"step_.bmp");//<---debug

	CombOverlapRectX(RectVec,TempVec, OverLayPix);

	DispBoxPositionPic(TempVec,  Pic.Width,  Pic.Height,"step0.bmp");//<---debug
	
	int cHeight = MHeight; 
	cHeight = GetAignTop(TempVec,cHeight);
	
	GetTextLineHorizontal(cHeight,TempVec,OutRects);
	
	if(DEBUG_OUT)
	{
		printf("Get Obj Num %i\n", OutRects.size());
		FILE* file =fopen("objpos.txt","wt+");
        int i; Loopi(OutRects.size())
			fprintf(file,"%3i,%3i,%3i,%3i\n", OutRects[i].left ,
			OutRects[i].top, OutRects[i].right,OutRects[i].bottom);
		fclose(file);
	}
	
	return OutRects.size();
}

/// \brief �õ����������������з�
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param RectVec     ��ʼ������ͨ�����Ӿ�������
/// \param OutRects    �������ַ������
/// \param MHeight     �ο����ֿ�     
/// \return  �з��Ժ�õ������ָ���
int GetTextLineComptChain::GetVerticalTextChain(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects, int Mwidth)
{
	vector<ImgRECT> TempVec;
	sort(RectVec.begin(),RectVec.end(),LessByY);
	CombOverlapRectY(RectVec,TempVec,OverLayPix);
	int cWidth = Mwidth;
	cWidth = GetMedianWidth(TempVec);
	GetTextLineVertical(cWidth,TempVec,OutRects);
	return OutRects.size();
}



/// \brief �ϲ�X�������ص�����Ӿ���
/// \param  RectVec    ���ϲ�����1
/// \param  OutRects   ���ϲ�����2
/// \param  nOverlap   ��������Ҫ������ص�����
/// \return  ��
int GetTextLineComptChain::CombOverlapRectX(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects,int nOverlap)
{
	int i,j,sz, tmp_left, tmp_right;
	i = j = sz = tmp_left = tmp_right = 0;
	sz = RectVec.size(); //�����ѡ����ĸ���
	
	vector<ImgRECT> RectVecBk;
    RectVecBk = RectVec;
	
	Loopi(sz)  RectVecBk[i].rectValid = true; //��ʼʱ������ÿһ�������Ϊ��Ч����
	
	Loopi(sz-1) 
	{  
		tmp_left  = RectVecBk[i].left;
		tmp_right = RectVecBk[i].right;		
		if(RectVecBk[i].rectValid == false) continue; 
		
		for(j = i+1; j< RectVecBk.size(); j++)
		{        
			int Debug_Val = RectVecBk[j].left;
			if( tmp_right - RectVecBk[j].left >= nOverlap)
			{
				RectVecBk[i] = MergRectBox(RectVecBk[i] , RectVecBk[j]);
				RectVecBk[i].rectValid = true;
				tmp_right = RectVecBk[i].right;
				RectVecBk[j].rectValid = false;
			}
			else
			{  
				i = j -1 ;
				break;
			}
		}
	}
	
	Loopi(RectVecBk.size())
	{
		if(RectVecBk[i].rectValid) 
		{   
			getRectWH(RectVecBk[i]);//���¾��ε�������Ϣ
			OutRects.push_back(RectVecBk[i]);
		}
	}
	
	return 1;
};

/// \brief  �ϲ�Y�������ص�����Ӿ���
/// \param  RectVec    ���ϲ�����1
/// \param  OutRects   ���ϲ�����2
/// \param  nOverlap   ��������Ҫ������ص�����
/// \return  ��
int GetTextLineComptChain::CombOverlapRectY(vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects,int nOverlap)
{
	int i,j; int y1,y2; int sz;
	sz = RectVec.size();
	
	Loopi(sz)  RectVec[i].rectValid = true;
	
	Loopi(sz-1) 
	{  
		y1 = RectVec[i].top;
		y2 = RectVec[i].bottom;
		
		if(RectVec[i].rectValid==false) 
			continue; 
		
		for(j = i+1; j< RectVec.size(); j++)
		{        
			if(y2 - RectVec[j].top >= nOverlap )
			{
				RectVec[i] = MergRectBox(RectVec[i] , RectVec[j]);
				y2 = RectVec[i].bottom;
				RectVec[j].rectValid = false;
			}
			else
			{  
				i = j -1 ;
				break;
			}
		}
	}
	
	Loopi(RectVec.size())
	{
		if(RectVec[i].rectValid) 
		{   
			getRectWH(RectVec[i]);//���¾��ε�������Ϣ
			OutRects.push_back(RectVec[i]);
		}
	}
	
	return 1;
};

//\brief �õ����εĸ߶���ֵ
int GetTextLineComptChain::GetMedianHeight(vector<ImgRECT>&RectVec)
{
	int i;
	vector<int> HeightVec;
	int Mid = RectVec.size()/2;
	
	Loopi(RectVec.size()) HeightVec.push_back(RectVec[i].height);
	
	nth_element(HeightVec.begin(), HeightVec.begin()+Mid, HeightVec.end());
	return HeightVec[Mid];
};

//\brief �õ����εĿ����ֵ
int GetTextLineComptChain::GetMedianWidth(vector<ImgRECT>&RectVec)
{
	int i;
	vector<int> WidthVec;
	int Mid = RectVec.size()/2;
	
	Loopi(RectVec.size()) WidthVec.push_back(RectVec[i].width);
	
	nth_element(WidthVec.begin(), WidthVec.begin()+Mid, WidthVec.end());
	return WidthVec[Mid];
};


//\brief �õ��ʻ�����͸��
bool GetTextLineComptChain::StrokeNumberY(ImgRECT Rect, int StrokeNumTh, int StrokeSumTh, int gap)
{
	int StrokeSum = 0;
	
	for(int x = Rect.left; x <= Rect.right; x++)
	{
		int start = -gap, StrokeNum = 0;
		
		for(int y = Rect.top; y <= Rect.bottom; y++)
		{
			if(*get_pix_color(Pic, x, y) == 0)
			{
				if(y - start >= gap) StrokeNum++;
				start = y;
			}
		}
		if(StrokeNum >= StrokeNumTh)
		{
			StrokeSum++;
			if(StrokeSum >= StrokeSumTh) return true;
		}
	}
	return false;
}

//\brief �õ��ʻ�����͸��
bool GetTextLineComptChain::StrokeNumberX(ImgRECT Rect, int StrokeNumTh, int StrokeSumTh, int gap)
{
	int StrokeSum = 0;
	for(int y = Rect.top; y <= Rect.bottom; y++)
	{
		int start = -gap, StrokeNum = 0;
		for(int x = Rect.left; x <= Rect.right; x++)
		{
			if(*get_pix_color(Pic, x, y) == 0)
			{
				if(y - start >= gap) StrokeNum++;
				start = x;
			}
		}
		if(StrokeNum >= StrokeNumTh)
		{
			StrokeSum++;
			if(StrokeSum >= StrokeSumTh) return true;
		}
	}
	return false;
}

/// \brief  �����з����ϵĵ�һ�κϲ� 
/// \param  ObjVec     �������Ԫ������
/// \param  RectVec    �������Ӿ�������
/// \param  hwRatio    �ַ��ĸ߶ȿ�ȱ���
/// \param  chrHeight  �ַ��Ĳο��߶�
/// \param  RectVec    �����������ӿ�  
/// \param  chObjDone  ��������ӿ��״̬��Ϣ(0��δʶ��,1�Ѿ����Ϊ����)
/// \param  OutRects   ��ѡ������Ӿ��ο�
/// \param  WidthSum   �ϸ��ַ�����ܺ�
/// \param  nChractNum �ϸ��ַ���
/// \return ��
void GetTextLineComptChain::MarkValidRectH1(float hwRatio, int chrHeight,vector<ImgRECT>&RectVec,
											vector<int>&chObjDone,vector<ImgRECT>&OutRects,
											int &WidthSum,int&nChractNum)
{
	int i,Num; Num = chObjDone.size();
	float nW,nH;
	float MinR = CNchr_WH_LowR, MaxR = CNchr_WH_HighR;
	//����MinR, MaxR ����Width/Height������
	if(hwRatio != 1.0f)
	{
		//<--------------hwVet[i] = 1.0f * RectVec[i].height / RectVec[i].width;
		MinR = 1.0 / hwRatio - RatioTolerance1;
		MaxR = 1.0 / hwRatio + RatioTolerance0;
	}
	//single
	//Loopi(Num) 
	for( i = 0; i < Num; i++)
	{
		if(chObjDone[i] != 0)
		{
			continue;
		}
		nW = RectVec[i].width  ; nH = RectVec[i].height ;
		//double v1,v2;v1 =  CHR_BOX_STEP(RectVec[i  ],RectVec[i-1]);v2 = CHR_BOX_STEP(RectVec[i+1],RectVec[i  ]) ;

		if (   (i == 0       || CHR_BOX_STEP(RectVec[i  ],RectVec[i-1]) > chrHeight*(MaxR + RatioTolerance1))
			&& (i == Num - 1 || CHR_BOX_STEP(RectVec[i+1],RectVec[i  ]) > chrHeight*(MaxR + RatioTolerance1))) // may // 1.2???
		{
			if(VAL_IN_RANGE(nW, MinR * nH, MaxR * nH ) && ( nH > chrHeight* MinLowChrHeightR ))
			{
				if(i == 0 && i+1 < Num )
				{
					ImgRECT RectTemp;
					RectTemp = MergRectBox(RectVec[i], RectVec[i+1]);
					if(   fabs(  CHR_WH_RATIO(RectTemp  ) - hwRatio) > 
						  fabs(  CHR_WH_RATIO(RectVec[i]) - hwRatio)
						&& CHR_BOX_GAP(RectVec[i + 1],RectVec[i]) < nH * ChrGapRatio)
					{
						chObjDone[i] = 1;
						WidthSum += nW; nChractNum++;
						OutRects.push_back(RectVec[i]);
					}
				}
				if( i > 0)
				{
					chObjDone[i] = 1;
					WidthSum += nW; nChractNum++;
					OutRects.push_back(RectVec[i]);
				}
			}
			
			if((nW > chrHeight * CNchr_WH_LowR) && (nH * 4 < chrHeight)) // һ 
			{
				chObjDone[i] = 1;
				WidthSum += nW; nChractNum++;
				OutRects.push_back(RectVec[i]);
			}
			
		}		
	}
	
	
}

/// \brief  ����ʺ�Ӣ���ַ���߱ȵ��ַ�
/// \param  ObjVec     �������Ԫ������
/// \param  RectVec    �������Ӿ�������
/// \param  hwRatio    �ַ��ĸ߶ȿ�ȱ���
/// \param  chrHeight  �ַ��Ĳο��߶�
/// \param  RectVec    �����������ӿ�  
/// \param  chObjDone  ��������ӿ��״̬��Ϣ(0��δʶ��,1�Ѿ����Ϊ�ַ�)
/// \param  OutRects   ��ѡ������Ӿ��ο�
/// \param  WidthSum   �ϸ��ַ�����ܺ�
/// \param  nChractNum �ϸ��ַ���
/// \return ��

void GetTextLineComptChain::MarkValidRectEN(float hwRatio,int chrHeight,vector<ImgRECT>&RectVec,
											vector<int>&chObjDone,vector<ImgRECT>&OutRects,
											int &WidthSum,int&nChractNum)
{
	int i,Num; Num = chObjDone.size();
	float nW,nH;
	vector<ImgRECT> TmpRect;
	Loopi(Num) 
	{
		if(chObjDone[i] != 0)
		{
			continue;
		}
		nW = RectVec[i].width  ; nH = RectVec[i].height;
		if( i > 0 && i < Num - 1 &&  CHR_BOX_GAP(RectVec[i+1], RectVec[i-1]) <  chrHeight * 0.6)  // nH * 0.6
		{
			continue;
		}
		
		if(nW *  20 < 3 * nH || nW * 10 < 3 * nH && 
			!((i == 0    || CHR_BOX_GAP(RectVec[i],RectVec[i-1]) >= 3 ) 
			&& (i + 1 == Num || CHR_BOX_GAP(RectVec[i+1],RectVec[i]) >= 3)))
		{
			continue;
		}
		if( i < Num - 1 && CHR_BOX_STEP(RectVec[i+1],RectVec[i]) < 1.1 * RectVec[i+1].height 
			&& RectVec[i + 1].width > 2 * nW)  // nH * 0.6
		{
			continue;
		}
		
		if(nH > chrHeight* 0.6 && nW < CNchr_WH_LowR * chrHeight && nW < nH 
			&& StrokeNumberY(RectVec[i], 4, 5,3) == false
			&& StrokeNumberX(RectVec[i], 3, 5,3) == false)
		{
			chObjDone[i] = 1;
			OutRects.push_back(RectVec[i]);
		}
	}
	
	
	
}

/// \brief  �����з����ϵĵڶ��κϲ� 
/// \param  ObjVec     �������Ԫ������
/// \param  RectVec    �������Ӿ�������
/// \param  hwRatio    �ַ��ĸ߶ȿ�ȱ���
/// \param  chrHeight  �ַ��Ĳο��߶�
/// \param  RectVec    �����������ӿ�  
/// \param  chObjDone  ��������ӿ��״̬��Ϣ(0��δʶ��,1�Ѿ����Ϊ����)
/// \param  OutRects   ��ѡ������Ӿ��ο�
/// \param  WidthSum   �ϸ��ַ�����ܺ�
/// \param  nChractNum �ϸ��ַ���
/// \return ��
void GetTextLineComptChain::MarkValidRectH2(float hwRatio,int chrHeight,vector<ImgRECT>&RectVec,
							 vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum)
{
	int i,Num; Num = chObjDone.size();
	float nW,nH;
	ImgRECT rRect;  
	float MinR = CNchr_WH_LowR, MaxR = CNchr_WH_HighR;
	if(hwRatio > CNchr_WH_HighR)
	{
		MinR = 1.0 / hwRatio - RatioTolerance0;
		MaxR = 1.0 / hwRatio + RatioTolerance1;
	}
	if(hwRatio < CNchr_WH_LowR)
	{
		MinR =  1.0 / hwRatio - RatioTolerance1;
		MaxR =  1.0 / hwRatio + RatioTolerance0;
	}
	Loopi(Num-1){
		if(chObjDone[i] == 0 && chObjDone[i+1] == 0) 
		{		  
			rRect = MergRectBox(RectVec[i],RectVec[i+1]);	
			getRectWH(rRect);
			nW = rRect.width;  nH = rRect.height;
			
			int dis = RectVec[i+1].left - RectVec[i].right;
			
				if(    VAL_IN_RANGE(nW, MinR * nH, MaxR * nH ) 
					&& nW > chrHeight * MinR 
					&& nH > chrHeight * MinLowChrHeightR//0.8 
					//&& dis * 6 <= chrHeight 
					&& (i == 0       ||
					CHR_BOX_STEP(RectVec[i+1],RectVec[i-1])  > chrHeight * MaxR && 
					  (CHR_BOX_GAP(RectVec[i],RectVec[i-1])  > dis || 
					  (CHR_BOX_GAP(RectVec[i],RectVec[i-1])  >= dis - 1 && nH >= chrHeight -  1))
					  )
					&& (i == Num - 2 || 
					CHR_BOX_STEP(RectVec[i+2],RectVec[i  ])  > chrHeight * MaxR &&
					(CHR_BOX_GAP(RectVec[i+2],RectVec[i+1])  > dis || 
					(CHR_BOX_GAP(RectVec[i+2],RectVec[i+1])  >= dis - 1 && nH >= chrHeight - 1))
					)
			   )

			{			   
				chObjDone[i  ] = 2; chObjDone[i+1] = 2;
				WidthSum += rRect.width ; nChractNum++; 
				OutRects.push_back(rRect);
			}
		}
	}
}

/// \brief  �õ��ַ���������ֵ
/// \param  RectVec     �����������ӿ�  
/// \param  hwRatio     �ַ���������ֵ���
/// \return ��
void GetTextLineComptChain::GetHwRatio(vector<ImgRECT>&RectVec, float &hwRatio)
{
	int i;
	if(RectVec.size() < 5)
	{
		return;
	}
	vector<float> hwVet;
	//hwVet.resize(RectVec.size());
	for( i = 0; i < RectVec.size(); i++)
	{
		double Ratio_ = CHR_WH_RATIO(RectVec[i]);
		if(VAL_IN_RANGE(Ratio_, CNchr_WH_LowR,CNchr_WH_HighR))
		{
           hwVet.push_back(Ratio_);
		}
		hwVet.push_back(Ratio_);
	}
	sort(hwVet.begin(), hwVet.end());
	float hwRMidean = hwVet[hwVet.size() / 2];
	if(VAL_IN_RANGE(hwRMidean, CNchr_WH_LowR,CNchr_WH_HighR))
	{
		hwRatio = hwRMidean;
		return;
	}
	int candidate_num = 0;
	for( i = 0; i < RectVec.size(); i++)
	{
		if(VAL_IN_RANGE(  CHR_WH_RATIO(RectVec[i]), hwRMidean - RatioTolerance0, hwRMidean + RatioTolerance0))
		{
			candidate_num ++;
		}
	}
	if(candidate_num > float(RectVec.size()) * 0.6 )
	{
		hwRatio = hwRMidean;
		return;
	}
}


/// \brief ���к������ֵ�ճ������
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param hwRatio     �ο��ַ�����
/// \param MeanWidth   �ο��ַ���ֵ
/// \param RectVec     �������Ԫ������
/// \param chObjDone   ��������ӿ��״̬��Ϣ(0��δʶ��,1�Ѿ����Ϊ����)
/// \param OutRects    ��ѡ������Ӿ��ο�
/// \param WidthSum   �ϸ��ַ�����ܺ�
/// \param nChractNum �ϸ��ַ���
/// \return  �з��Ժ�õ������ָ���
void GetTextLineComptChain::MarkValidRectH_overlap(float hwRatio, int MeanWidth,
					 vector<ImgRECT>&RectVec,vector<int>&chObjDone,vector<ImgRECT>&OutRects,
						  int &WidthSum,int&nChractNum)
{
	int i,j,k,Num;
	Num = chObjDone.size();
	int range;
	range = MeanWidth / 5;
	vector<int> DisVet;
	for(  i = 0; i < Num; i++)
	{
		//if( i > 1  &&  
		//	( (pnDone[i] == 1 && pnDone[i-1] == 1) || (pnDone[i] == 1 && pnDone[i-1] == 2) || (pnDone[i] == 2 && pnDone[i-1] == 1)) 
		//	&& 1) // ���ڵ������ֶ��ѱ�ʶ�� // 1?ͬ�ϣ��Է���������
		if(i > 0 && chObjDone[i] > 0 && chObjDone[i - 1] > 0)
		{
			if( CHR_BOX_GAP(RectVec[i],RectVec[i - 1]) + 1 > 0 )
			{
				DisVet.push_back(CHR_BOX_GAP(RectVec[i], RectVec[i-1]) + 1);
			}
		}
	}
	int nAveW2 = MeanWidth;
	int nAveD = MeanWidth / 5;
	if(DisVet.size() > 0)
	{
		sort(DisVet.begin(),DisVet.end());
		nAveD =DisVet[DisVet.size() / 2];
		nAveW2 = MeanWidth + nAveD;
	}
	range = nAveD + 1;
	for( i = 0; i < Num; i++)
	{
		if(chObjDone[i] == 0)
		{
			if(RectVec[i].width > 1.7 * MeanWidth) // ��ƽ����ȵļ���������
			{
				int no = 0;
				int width = RectVec[i].width;
				//int n = (width + 0.5 * MeanWidth) / (MeanWidth);
				int n = (width + 0.5 * nAveW2) / (nAveW2);
				ImgRECT rRect0 = RectVec[i];
				int post1 = rRect0.left, post2  = rRect0.left;
				for(int j = 0; j < n; j++)
				{
					ImgRECT rRect = rRect0;
					rRect.left = post2;
					if(j != n-1)
					{
						//int post = (j+1) * width / n + RectVec[i].left;
						int post = (j+1) * nAveW2 - nAveD / 2 + RectVec[i].left;
						int mingap = rRect0.bottom - rRect0.top + 1;
						post1 = post2 = post;
						for( int m = post - range; m <= post + range; m++)
						{
							if(m < rRect0.left  ) continue;
							if(m > rRect0.right|| m < rRect.left) break;
						 
							int y1 =  rRect0.bottom, y2 = rRect0.top;
							for(  k = rRect0.top; k <= rRect0.bottom; k++)
							{
								byte *pub = get_pix_color(Pic,m, k);
								if(pub[0] == 0) // ������
								{
									y1 = k;
									break;
								}
							}
							for( k = rRect0.bottom; k >= rRect0.top; k--)
							{
								byte *pub = get_pix_color(Pic,m, k);
								if(pub[0] == 0) // ������
								{
									y2 = k;
									break;
								}
							}
							int gap = y2 - y1 + 1;
							if(y1 > y2)  gap = 0;// û�����ֵ�
							 
							if(gap < mingap)
							{
								mingap = gap;
								post1 = post2 = m;
							}
							else if(mingap == 0 && gap > 0)
							{
								break;
							}
							else if(mingap == 0 && gap == 0)
							{
								post2 = m;
							}
						}
						rRect.right = post1;
					}

					ImgRECT rT;
					rT.left = max(0,rRect.left - 1);
					rT.top = max(0,rRect.top - 1);
					rT.right = min(Pic.Width - 1,rRect.right + 1);
					rT.bottom = min(Pic.Height - 1,rRect.bottom + 1);
					rT.width = rT.right - rT.left + 1;
					rT.height = rT.bottom - rT.top + 1;
					if(j == n - 1)
					{
						if(rT.width >  1.5 * MeanWidth && no == 0)
						{
							no++;
							int wm = 0;
							for(int m = 0; m < j; m++)
							{
								wm += OutRects[OutRects.size() - 1 - m].width;
							}
							if(j > 0)
							{
								wm /= j;
								nAveW2 = nAveW2 + (wm - MeanWidth);
							}
							else
							{
								wm = MeanWidth;
							}
							if(rT.width >  1.7 * wm)
							{
								n++;
								j--;
							}
							else
							{ OutRects.push_back(rT); }
							
						}
						else
						{ OutRects.push_back(rT); 	}
					}
					else
					{ OutRects.push_back(rT); }	
				}
				chObjDone[i] = 4;
			}
		}
	}
	
}


void GetTextLineComptChain::MarkValidRectV_overlap(float hwRatio, int MeanHeight,vector<ImgRECT>&RectVec,
							 vector<int>&chObjDone,vector<ImgRECT>&OutRects,int &WidthSum,int&nChractNum)
{
	int i,j,k, Num;
	Num = chObjDone.size();
	int range;
	range = MeanHeight / 5;
	vector<int> DisVet;
	for(  i = 0; i < Num; i++)
	{
		//if( i > 1  &&  
		//	( (pnDone[i] == 1 && pnDone[i-1] == 1) || (pnDone[i] == 1 && pnDone[i-1] == 2) || (pnDone[i] == 2 && pnDone[i-1] == 1)) 
		//	&& 1) // ���ڵ������ֶ��ѱ�ʶ�� // 1?ͬ�ϣ��Է���������
		if(i > 0 && chObjDone[i] > 0 && chObjDone[i - 1] > 0)
		{
			if( RectVec[i].top - RectVec[i - 1].bottom + 1 > 0 )
			{
				DisVet.push_back(RectVec[i].top - RectVec[i - 1].bottom + 1);
			}
		}
	}
	int nAveH2 = MeanHeight;
	int nAveD = MeanHeight / 5;
	if(DisVet.size() > 0)
	{
		sort(DisVet.begin(),DisVet.end());
		nAveD =DisVet[DisVet.size() / 2];
		nAveH2 = MeanHeight + nAveD;
	}
	range = nAveD + 1;
	int temp = RectVec[0].height;
	for( i = 0; i < Num; i++)
	{
		if(chObjDone[i] == 0)
		{
			if(RectVec[i].height > 1.7 * MeanHeight) // ��ƽ����ȵļ���������
			{
				int no = 0;
				int height = RectVec[i].height;
				int n = (height + 0.5 * nAveH2) / (nAveH2);
				ImgRECT rRect0 = RectVec[i];
				int post1 = rRect0.top, post2  = rRect0.top;
				for(int j = 0; j < n; j++)
				{
					ImgRECT rRect = rRect0;
					rRect.top = post2;
					if(j != n-1)
					{
						//int post = (j+1) * width / n + RectVec[i].left;
						int post = (j+1) * nAveH2 - nAveD / 2 + RectVec[i].top;
						int mingap = rRect0.bottom - rRect0.top + 1;
						post1 = post2 = post;
						for( int m = post - range; m <= post + range; m++)
						{
							if(m < rRect0.top  )
							{
								continue;
							}
							if(m > rRect0.bottom|| m < rRect.top)
							{
								break;
							}
							int x1 =  rRect0.left, x2 = rRect0.right;
							for(k = rRect0.left; k <= rRect0.right; k++)
							{
								byte *pub = get_pix_color(Pic,k, m);
								if(pub[0] == 0) // ������
								{
									x1 = k;
									break;
								}
							}
							for(  k = rRect0.left; k >= rRect0.right; k--)
							{
								byte *pub = get_pix_color(Pic,k, m);
								if(pub[0] == 0) // ������
								{
									x2 = k;
									break;
								}
							}
							int gap = x2 - x1 + 1;
							if(x1 > x2) // û�����ֵ�
							{
								gap = 0;
							}
							if(gap < mingap)
							{
								mingap = gap;
								post1 = post2 = m;
							}
							else if(mingap == 0 && gap > 0)
							{
								break;
							}
							else if(mingap == 0 && gap == 0)
							{
								post2 = m;
							}
						}
						rRect.bottom = post1;
					}
					ImgRECT rT;
					rT.left = max(0,rRect.left - 1);
					rT.top = max(0,rRect.top - 1);
					rT.right = min(Pic.Width - 1,rRect.right + 1);
					rT.bottom = min(Pic.Height - 1,rRect.bottom + 1);
					rT.width = rT.right - rT.left + 1;
					rT.height = rT.bottom - rT.top + 1;
					if(j == n - 1)
					{
						if(rT.height >  1.5 * MeanHeight && no == 0)
						{
							no++;
							int wm = 0;
							for(int m = 0; m < j; m++)
							{
								wm += OutRects[OutRects.size() - 1 - m].height;
							}
							if(j > 0)
							{
								wm /= j;
								nAveH2 = nAveH2 + (wm - MeanHeight);
							}
							else
							{
								wm = MeanHeight;
							}
							if(rT.height >  1.7 * wm)
							{
								n++;
								j--;
							}
							else
							{
								OutRects.push_back(rT);
							}
						}
						else
						{
							OutRects.push_back(rT);
						}
					}
					else
					{
						OutRects.push_back(rT);
					}
				}
				chObjDone[i] = 4;
			}
		}
	}
}


/// \brief �õ�����������ӿ�ĺ����з�
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param RectVec       ��ʼ������ͨ�����Ӿ�������
/// \param OutRects      �������ַ������
/// \param chrHeight     �ο����ָ߶�     
/// \return  �з��Ժ�õ������ָ���
int GetTextLineComptChain::GetTextLineHorizontal(int chrHeight,vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects)
{
	int i;
	vector<int> chObjDone;
	chObjDone.resize(RectVec.size());
	fill(chObjDone.begin(),chObjDone.end(),0);
	int Num = RectVec.size();
	int nChractNum = 0;
	int WidthSum   = 0;
	OutRects.clear();
	
	float nW,nH;
	ImgRECT rRect;  
	float hwRatio = 1.0f;
	GetHwRatio( RectVec, hwRatio);
	
	MarkValidRectH1( hwRatio,chrHeight, RectVec, chObjDone, OutRects, WidthSum, nChractNum);

	DispBoxPositionPic(OutRects,  Pic.Width,  Pic.Height,"step1.bmp");//<---debug

	MarkValidRectH2( hwRatio,chrHeight, RectVec, chObjDone, OutRects, WidthSum, nChractNum);

	DispBoxPositionPic(OutRects,  Pic.Width,  Pic.Height,"step2.bmp");//<---debug
	
	if(ChineseOnly == false)
	{
		MarkValidRectEN( hwRatio,chrHeight, RectVec, chObjDone, OutRects, WidthSum, nChractNum);
	}
	
	float MeadWidth;
	if(nChractNum > 0)
	{
		MeadWidth = float(WidthSum)/float(nChractNum);
	}
	else
	{
		MeadWidth = chrHeight / hwRatio;
	}
	
	Loopi(Num)
	{
		if(chObjDone[i] == 0 )
		{
			nW = RectVec[i].width ; nH = RectVec[i].height;
			
			if( VAL_IN_RANGE(nW, CNchr_WH_LowR * MeadWidth, CNchr_WH_HighR * MeadWidth ))
			{
				OutRects.push_back(RectVec[i]);
				chObjDone[i] = 3;
				WidthSum += nW;
				nChractNum++; 
			}
		}
	}
	if(nChractNum > 0)
	{
		MeadWidth = float(WidthSum)/float(nChractNum);
	}
	else
	{
		MeadWidth = chrHeight / hwRatio;
	}
	
	Loopi(Num-1)
	{
		if(chObjDone[i] == 0 ) // has not been recognized
		{		
			int j = i + 1;
			for(; j < Num; j++)
			{
				if(chObjDone[j] != 0)
				{
					break;
				}
			}
			
			rRect = RectVec[i];
			
			
			for(int k = i + 1; k < j; k++)
			{
				rRect = MergRectBox(rRect,RectVec[k]);	
			}
			
			getRectWH(rRect);
			nW = rRect.width;  nH = rRect.height;
			if( VAL_IN_RANGE(nW, CNchr_WH_LowR * MeadWidth, CNchr_WH_HighR * MeadWidth ))
			{
				for(int k = i; k < j; k++)
				{
					chObjDone[k  ] = 2;;
				}
				
				WidthSum += rRect.width ; nChractNum++; 
				OutRects.push_back(rRect);
			}
		}
	}
    
	MarkValidRectH_overlap( hwRatio,chrHeight, RectVec, chObjDone, OutRects, WidthSum, nChractNum);
	Loopi(Num)
	{
		if(chObjDone[i] == 0)
		{
			OutRects.push_back(RectVec[i]);
		}
	}
	sort(OutRects.begin(),OutRects.end(),LessByX);
	vector<bool> flag;
	flag.resize(OutRects.size());
	Loopi(flag.size())
	{
		flag[i] = true;
	}

	if( OutRects.size() > 1)
	{
		for(int i = 1; i < OutRects.size() - 1; i++)
		{
			if(     OutRects[i].width   < ChrGapRatio * OutRects[i].height 
				&&  OutRects[i+1].width < RatioTolerance1 * OutRects[i+1].height 
				&&  OutRects[i].height  < OutRects[i+1].height
				&&  OutRects[i - 1].width < 0.65 * OutRects[i - 1].height )
			{
				OutRects[i - 1] = MergRectBox(OutRects[i - 1],RectVec[i]);	
				OutRects[i - 1] = MergRectBox(OutRects[i - 1],RectVec[i + 1]);	
				flag[i] = false;
				flag[i + 1] = false;
			}
		}
	}
	if(FirstLineFlag == true)
	{
		int no = 0;
		int i = 0;
		while( no < 9)
		{
			if( OutRects[i].width < CNchr_WH_LowR * OutRects[i].height && i + 1 < OutRects.size() && OutRects[i + 1].width < CNchr_WH_LowR * OutRects[i + 1].height)
			{
				OutRects[i] = MergRectBox(OutRects[i],RectVec[i + 1]);	
				flag[i + 1] = false;
				i++;
			}
			no++;
			i++;
		}
	}
	
	vector<ImgRECT> tmp;
	for( i = 0; i < OutRects.size() ; i++)
	{
		if(flag[i] == true 
		/*&& (OutRects[i].height > MeadWidth * 0.6 || OutRects[i].height < MeadWidth * 0.25)*/)
		{
			tmp.push_back(OutRects[i]);
		}
	}
	
	OutRects = tmp;
	sort(OutRects.begin(),OutRects.end(),LessByX);
	
	if(nChractNum!=0)
		return WidthSum/nChractNum;
	return -1;
};

/// \brief �õ����������������з�
///  ������ʼ������ͨ�����Ӿ�������ϲ��зֺ����Ӿ�������
/// \param RectVec      ��ʼ������ͨ�����Ӿ�������
/// \param OutRects     �������ַ������
/// \param chrWidth     �ο����ֿ��     
/// \return  �з��Ժ�õ������ָ���
int  GetTextLineComptChain::GetTextLineVertical(int chrWidth,vector<ImgRECT>&RectVec,vector<ImgRECT>&OutRects)
{
	int i;
	vector<int> chObjDone;
	chObjDone.resize(RectVec.size());
	fill(chObjDone.begin(),chObjDone.end(),0);
	int Num = RectVec.size();
	int nChractNum = 0;
	int HeightSum   = 0;
	OutRects.clear();
	
	float nW,nH;
	ImgRECT rRect;  
	float hwRatio = 1.0f;
	GetHwRatio( RectVec, hwRatio);
	float MinR = CNchr_WH_LowR, MaxR = CNchr_WH_HighR;
	if(hwRatio != 1.0f)
	{
		MinR = 1.0 / hwRatio - RatioTolerance0;
		MaxR = 1.0 / hwRatio + RatioTolerance0;
	}
	
	Loopi(Num) {
		nW = RectVec[i].width  ; nH = RectVec[i].height ;		
		if (       (i == 0       || (RectVec[i].bottom - RectVec[i-1].top) > chrWidth * hwRatio * 1.6 )
			&& (i == Num - 1 || (RectVec[i+1].bottom - RectVec[i].top) > chrWidth * hwRatio * 1.6 ))
		{
			if(VAL_IN_RANGE(nW, MinR * nH, nH* MaxR ) && ( nH > chrWidth* MinLowChrHeightR ))
			{
				chObjDone[i] = 1;
				HeightSum += nH; nChractNum++;
				OutRects.push_back(RectVec[i]);
				continue;
			}
			
			if((nW > chrWidth * CNchr_WH_LowR) && (nH * 4 < chrWidth)) // һ 
			{
				chObjDone[i] = 1;
				HeightSum += nH; nChractNum++;
				OutRects.push_back(RectVec[i]);
				
			}
			
			
		}		
	}
	
	Loopi(Num-1)
	{
		if(chObjDone[i] == 0 && chObjDone[i+1] == 0) // has not been recognized
		{		  
			rRect = MergRectBox(RectVec[i],RectVec[i+1]);	
			getRectWH(rRect);
			nW = rRect.width;  nH = rRect.height;
			
			int dis = RectVec[i+1].top - RectVec[i].bottom + 1;
			
			if( VAL_IN_RANGE(nW, MinR * nH, MaxR * nH ) && nH > chrWidth * MinCharWidthR && dis * 7 <= chrWidth 
				&& (i == 0       || RectVec[i+1].bottom - RectVec[i-1].top > chrWidth * hwRatio * 1.6 && RectVec[i  ].top - RectVec[i-1].bottom + 1 > dis)
				&& (i == Num - 2 || RectVec[i+2].bottom - RectVec[i  ].top > chrWidth * hwRatio * 1.6 && RectVec[i+2].top - RectVec[i+1].bottom + 1 > dis))
			{
				
				chObjDone[i  ] = 2; chObjDone[i+1] = 2;
				HeightSum += rRect.width ; nChractNum++; 
				OutRects.push_back(rRect);
			}
		}
	}
	
	float MeadHeight;
	if(nChractNum > 0)
	{
		MeadHeight = float(HeightSum)/float(nChractNum);
	}
	else
	{
		MeadHeight = chrWidth * hwRatio;
	}
	
	
	Loopi(Num)
	{
		if(chObjDone[i] == 0 )
		{
			nW = RectVec[i].width ;
			nH = RectVec[i].height;
			
			if( VAL_IN_RANGE(nH, 0.9 * MeadHeight, 1.1 * MeadHeight ))
			{
				OutRects.push_back(RectVec[i]);
				chObjDone[i] = 3;
				HeightSum += nH;
				nChractNum++; 
			}
			
		}		
	}
	if(nChractNum > 0)
	{
		MeadHeight = float(HeightSum)/float(nChractNum);
	}
	else
	{
		MeadHeight = chrWidth * hwRatio;
	}
	Loopi(Num-1)
	{
		if(chObjDone[i] == 0 ) // has not been recognized
		{		
			int j = i + 1;
			for(; j < Num; j++)
			{
				if(chObjDone[j] != 0)
				{
					break;
				}
			}
			rRect = RectVec[i];
			for(int k = i + 1; k < j; k++)
			{
				rRect = MergRectBox(rRect,RectVec[k]);	
			}
			getRectWH(rRect);
			nW = rRect.width;  nH = rRect.height;
			if( VAL_IN_RANGE(nH, CNchr_WH_LowR * MeadHeight, CNchr_WH_HighR * MeadHeight ))
			{
				for(int k = i; k < j; k++)
				{
					chObjDone[k  ] = 2;
				}
				HeightSum += rRect.width ; nChractNum++; 
				OutRects.push_back(rRect);
			}
		}
	}
	MarkValidRectV_overlap( hwRatio,chrWidth, RectVec, chObjDone, OutRects, HeightSum, nChractNum);
	Loopi(Num)
	{
		if(chObjDone[i] == 0)
		{
			OutRects.push_back(RectVec[i]);
		}
	}
	sort(OutRects.begin(),OutRects.end(),LessByY);
	
	if(nChractNum > 0)
	{
		return HeightSum/nChractNum;
	}
	else
	{
		return 0;
	}
	
};
//<========================================================================================================

#endif
