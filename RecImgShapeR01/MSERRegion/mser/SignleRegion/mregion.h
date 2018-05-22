 

#ifndef Mser_Region_Head
#define Mser_Region_Head
#include <vector>
#include <deque>
#include <assert.h>
 
#include "c24bitmap.h"
#define NUM_INFINITY 99999 //infinity�����
//#include "hzmerge.h"
#include <string>

using namespace std;

//\brief ��ɫ����ĺ궨��
#define ColorDist(R1,G1,B1,R2,G2,B2)  ((R1-R2)*(R1-R2) + (G1-G2)*(G1-G2) + (B1-B2)*(B1-B2))
//\brief ��ֵ��Χ�ĺ궨��
#define VAL_IN_RANGE(x, a, b)   ((x) >= (a) && (x) <= (b))
 

//\brief Region ����Ӿ�������
struct ImgRECT  
{
	int left ,  right;
	int top  , bottom;
	int width, height;
	int rectValid;     //�Ƿ���Ч����
    
	int     GeoX,GeoY;
	double  NewX,NewY;
    double  Weight;
 };


//\brief �����ڵĵ�ṹ
struct RPoint
{
  float x,y;
};

 
/// \brief �ж�һ�����Ƿ��ھ�����
///  
/// \param p   ��
/// \param r   ���жϵľ��� 
/// \return    ����ڷ���1, ���򷵻�0
int RgPtInRect(RPoint & p, ImgRECT &r)
{
	return p.x > r.left && p.x < r.right &&
		   p.y > r.top  && p.y < r.bottom ;
};


/*--------------------------------------------
Region��һЩ�������Ի���AnalysisRegionObj������
---------------------------------------------*/ 


/// A Maximally Stable Extremal Region.
//\brief ��������
struct RgFeature
{
	int   perimeter;          // ������ܳ�
	int   InnerHoles;         // �����ڲ��ն���
	float HoleRatio;          // �ն��ı���
	float postProb;           // �����Ƿ���������ĺ������
	float StrokeMean,StrokeStd,StkMin,StkMax,StkMedian;   
	//�ʻ�ƽ��ֵ�������Сֵ�����ֵ���Ȼ���ֵ
	//float StrokeStd;     //
	int   label;              // ����ı��(���)
	float  aspectRatio;       // ����ĳ����
	//vector<float> OccupyVec;//�������񻮷�ռ����Ϣ
	float AngleCnt[4];        // ����ĽǶȷ���ͳ��
	double roundratio;        // �����Բ�Ķ�(�����ж�ϸ������)
};


//\brief MSER�ȶ���������ṹ��
class Region
{
public:
	int level_;   // ���򱻲��ҳ�ʱ��ĻҶȵȼ�
	int r_id;     // �ڽ���MSER����ʱ�������ڶ����е�Ψһ���id
	int pixel_;   // ������������ص�(y * width + x).
	int area_;    // ��������
	int childNum; // �������������
	int RgIndex;  // //�����Index�������ӵ�ʱ���ʾ�������������
	double moments_[5];  // һ�ס�����ת������(x, y, x^2, xy, y^2).
	double variation_;	 // MSER �����variation.(�ȶ��仯��)
	double normal_variation;  // �������MSER variation

//----------------------------------------------------------------------
//-----------------�����һЩ����-------------------------------------
	float x, y, GeoX, GeoY;  //���� x, y ; �������� GeoX, GeoY
	float RR, GG,BB;         //������ɫRGB
	float HH, SS,VV;         //ת�������ӦHSV
	int   left, right, top, bottom;    //boundingbox����������
	int   rwidth, rheight;   //�������Ӿ��ο�͸�
	int   IsValid;           //�����Ƿ���Ч
	float variance_wp;    	 //�����ڻҶ��ϵı�����(variance with penalty)
    RgFeature Feature;       //����ͳ������

	Region(int level = 256, int pixel = 0);
	vector<RPoint> PtVec;              //�����ڵĵ�
	vector<RPoint> ContourPtVec;       //���������������
	vector<RPoint> ConvexHullPtVec;    //�����ڵĵ�
	double ConvexHullArea;
	vector< vector<RPoint> > MContours; //�����ڵ�����Ƕ����Ϣ(����������)
	vector<RPoint> InflectionPts;      //����Ľǵ�(ת�۵�)
	vector<RPoint> CenterPts;          //����ĹǼ����ĵ�
 	vector<RPoint> OutterCurve;        //����������� 
	Region * parent_;      // ����ĸ��ڵ�ָ��
	Region * child_ ;	   // ����ĺ���ָ��
	Region * next_  ;	   // ������ֵ�ָ��
	
	void accumulate(int x, int y);   //��ǰlevel��������
	void merge(Region * child);
	void getvariation(int delta);
	
	friend class MSER;
    void  update_childnum();
	void  SaveMserRegion(vector<Region>&ReginVec, double minVariation);
	void  PrintTree(FILE*file,int start);
	int GRegionMinSize, GRegionMaxSize; //�������
	int PicWidth, PicHeight; //�������
 
	void SaveToFile(vector<Region>&RgVec, char*filename);
	void LoadFromFile(vector<Region>&RgVec, char*filename);
	//void SaveToPic(char*filename);
};


/// \brief �ж�һ�������Ƿ��ھ�����
///  
/// \param p   ����
/// \param r   ���жϵľ��� 
/// \return    ����ڷ���1, ���򷵻�0
int RegionInRect(Region & p, ImgRECT &r)
{
	return p.x > r.left && p.x < r.right &&
		p.y > r.top  && p.y < r.bottom ;
};


//*****************************************************
//**********������ʾRegion��Ϣ��ͼ��*******************
//*****************************************************
//#include "RegionDisp.h"
//*****************************************************
//*****************************************************
//*****************************************************
 

int Region::GRegionMinSize = -1;
int Region::GRegionMaxSize = -1;
int Region::PicWidth  = -1;
int Region::PicHeight = -1;

/// \brief  ��MSER��������ӡ���Ե��ļ���ȥ
/// \param  file   �ļ���
/// \param  start  ����Ĳ����Ϣ
/// \return ��
void Region::PrintTree(FILE*file,int start)
{
	string blankstr;
	blankstr = ""; 
	//if(area_<32) //<-----?���������С������
	//	 return;

	int i;
	for(i=0;i<start;i++) blankstr+="--";
	
	fprintf(file,"%s,id:%i,parent_id:%i;thre:%i, %.6lf ---size:%i , childnum %i, pos:%i,%i,%i,%i\n",
		blankstr.c_str(),
		r_id,
		parent_->r_id,
		level_,
		normal_variation,area_,
		childNum, left, top, right, bottom);
	
	for(Region * child = child_; child; child = child->next_)
		child->PrintTree(file,start+1);
}


/// \brief  �õ�MSER�������ڲ��������ָ��
//  �ݹ�������������������õ�����������ڲ���
//  ���浽(vector<RPoint> PtVec;//�����ڵĵ�)
/// \param  T  ��������ĳ������ 
/// \return �� 
vector<RPoint>* GetSubNode(Region*T)
{	
	int i;
	for(Region * child = T->child_; child; child = child->next_) 
	{
		vector<RPoint>* PtVecPt;
		PtVecPt= GetSubNode(child);
		for(i=0;i< PtVecPt->size();i++)
		{
			T->PtVec.push_back( (*PtVecPt)[i] );
		}
	}
	
	return &T->PtVec;
}


/// \brief  ��������С��minVariation������ݹ鱣�浽ReginVec
//   
/// \param  ReginVec      ��������Region����
/// \param  minVariation  ��С������
/// \return �� 
void Region::SaveMserRegion(vector<Region>&ReginVec,double minVariation)
{
	int i;
	if((this->normal_variation < minVariation) )
	{
		if((this->area_>=GRegionMinSize) &&
			(this->area_<=GRegionMaxSize) )
		{ 
			int RegionWidth,RegionHeight;
			RegionWidth  =  right  -  left +1;
			RegionHeight =  bottom -  top  +1;
			
			//if((RegionWidth *4 < PicWidth ) 
			//	&&
			//	(RegionHeight*3 < PicHeight*2) )

			{Region R;
			R = *(this);
			ReginVec.push_back(R);}
		}
		
	}
	
	for(Region * child = child_; child; child = child->next_) 
		child->SaveMserRegion(ReginVec,minVariation);
	
}



Region::Region(int level, int pixel) : level_(level), pixel_(pixel), area_(0),
variation_(NUM_INFINITY ),  parent_(0), child_(0), next_(0)
{
	fill_n(moments_,5, 0.0);
	
	r_id = 0;
	left =  right =  top =  bottom = -1;
}

inline void Region::accumulate(int x, int y)
{
	++area_;
	moments_[0] += x;
	moments_[1] += y;
	moments_[2] += x * x;
	moments_[3] += x * y;
	moments_[4] += y * y;
	
    
	if(left>=0)
	{
		if(left   > x ) left   = x;
		if(right  < x ) right  = x;
		if(top    > y ) top    = y;
		if(bottom < y ) bottom = y;
	}
	else
	{
		left = right  = x;
		top  = bottom = y;
	}
	
	RPoint Pt; Pt.x = x; Pt.y =y;
	PtVec.push_back(Pt);
	
}



/// \brief  һ���������������(������)�ϲ�
///  
/// \param subNode   ���ϲ�������ָ��
/// \return   
void Region::merge(Region * subNode)  //�ϲ�����
{  
    assert(!subNode->parent_ );
	assert(!subNode->next_   );	
	
	area_ += subNode->area_;
	moments_[0] += subNode->moments_[0];
	moments_[1] += subNode->moments_[1];
	moments_[2] += subNode->moments_[2];
	moments_[3] += subNode->moments_[3];
	moments_[4] += subNode->moments_[4];
	
	subNode->next_ = child_;
	child_ = subNode;
	subNode->parent_ = this;
	
	if(left>=0)
	{
		if(left   > subNode->left    ) left   = subNode->left   ;
		if(right  < subNode->right   ) right  = subNode->right  ;
		if(top    > subNode->top     ) top    = subNode->top    ;
		if(bottom < subNode->bottom  ) bottom = subNode->bottom ;
	}              
	else
	{
		left   = subNode->left   ;
		right  = subNode->right  ;
		top    = subNode->top    ;
		bottom = subNode->bottom ;
	}
}

/// \brief  ���㵱ǰ�Ҷ���ֵV�仯delta����������ı仯
///  
/// \param delta     �Ҷ���ֵ�仯
/// \return   
void Region::getvariation(int delta)
{
	Feature.aspectRatio = (float(right - left + 1) / float(bottom - top + 1));
	
	// Find the last parent with level not higher than level + delta
	const Region * parent = this;
	
	while (parent->parent_ && (parent->parent_->level_ <= (level_ + delta)))
		parent = parent->parent_;
	
	// Calculate variation
	variation_ = static_cast<double>(parent->area_ - area_) / area_;
	
	
	// Process all the children
	childNum = 0;
	for (Region * child = child_; child; child = child->next_) {
		child->getvariation(delta);
		childNum ++;
	}
	
	normal_variation = variation_;//<---------------�д�����
	if(Feature.aspectRatio<0.3)
	{
      normal_variation+= 0.35*(0.3 - Feature.aspectRatio);
	}
	if(Feature.aspectRatio>1.2)
	{
      //normal_variation+= 0.01 * (Feature.aspectRatio-1.2);
	}
}


//---------------------------------------------------------------------------------------------
/// \brief  ����MSER�ڵ���ӽڵ������Ϣ
/// ͨ�������ӽڵ���������ӽڵ�������ݹ����ʵ���������ĸ���
/// \return   
void Region::update_childnum()
{
	childNum = 0;
	for (Region * child = child_; child; child = child->next_) {
		child->update_childnum();
		childNum ++;
	}
	
}


//---------------------------------------------------------------------------------------------
/// \brief  ��һ��MSER��������л��ļ��ж�ȡ����
///  
/// \param RgVec     MSER��������
/// \param filename  ������ļ���
/// \return   
void Region::LoadFromFile(vector<Region>&RgVec, char*filename)
{
	int i,j;
	
	int ObjNum ;
	FILE*file = fopen(filename,"rb");
	
	fread(&rwidth , 4, 1, file);
	fread(&rheight, 4, 1, file);
	fread(&ObjNum , 4, 1, file);
	
	int TotalCnt = rwidth * rheight;
	int * DumpImg =  new int[rwidth * rheight];
	
	fread((void*)DumpImg,rwidth * rheight, 4, file);
	fclose(file);
    
	RgVec.resize(ObjNum);
	
	Loopi(TotalCnt)
	{
		if(DumpImg[i] != -1)
		{int xx,yy;
		xx = i % rwidth;
		yy = i / rwidth;
		RPoint Pt;
		Pt.x = xx;
		Pt.y = yy;
		RgVec[DumpImg[i]].PtVec.push_back(Pt);
		}
	}
	
    delete [] DumpImg;  
}


/// \brief  ��һ��MSER�������л����ļ���
///  
/// \param RgVec     MSER��������
/// \param filename  ������ļ���
/// \return    
void Region::SaveToFile(vector<Region>&RgVec, char*filename)
{
	int i,j;
	int * DumpImg =  new int[rwidth * rheight];
	
	int TotalCnt = rwidth * rheight;
	Loopi(TotalCnt)
		DumpImg[i] = -1;
	
	Loopi(RgVec.size())
	{
		Loopj(RgVec[i].PtVec.size())
		{
			int xx,yy;
			xx = RgVec[i].PtVec[j].x ;
			yy = RgVec[i].PtVec[j].y ;
			int Pos;
			Pos = yy * rwidth + xx;
			DumpImg[Pos] = i;
		}
	}
	
	int ObjSize = RgVec.size();
	FILE*file = fopen(filename,"wb+");
	fwrite(&rwidth , 4, 1, file);
	fwrite(&rheight, 4, 1, file);
	fwrite(&ObjSize, 4, 1, file);
	fwrite((void*)DumpImg,rwidth * rheight, 4, file);
	
	fclose(file);
	
	delete [] DumpImg;
}

/// \brief  ��һ��MSER�������������������ļ�
///  
/// \param RgVec     MSER��������
/// \param filename  ������ļ���
/// \return  
void PrintMSERRegionAxis(vector<Region>&RgVec, char*filename)
{
  int i;
  FILE*file = fopen(filename,"wt+");
  
  Loopi(RgVec.size())
  {
    fprintf(file,"%3i, %i, %i, %i, %i, %i\n", i,
		  RgVec[i].left, RgVec[i].top, RgVec[i].right, RgVec[i].bottom, RgVec[i].PtVec.size());

  }

  fclose(file);
}

void DisPlayRegionCircle(vector<Region>&RgVec, int ImgWidth,int ImgHeight, char*filename)
{
	int i;
	C24BitMap CPic;
	CPic.FormatF(ImgWidth,ImgHeight);
	 
	
	Loopi(RgVec.size())
	{
	 float Radius = sqrt(double(RgVec[i].PtVec.size()))/2.0;

	 CPic.PenColor.R = RgVec[i].RR;
	 CPic.PenColor.G = RgVec[i].GG;
	 CPic.PenColor.B = RgVec[i].BB;

     CPic.DrawCircle(RgVec[i].GeoX, RgVec[i].GeoY, Radius);
		//fprintf(file,"%3i, %i, %i, %i, %i\n", i,
		//	RgVec[i].left, RgVec[i].top, RgVec[i].right, RgVec[i].bottom);
		
	}
	
	CPic.Save(filename);
	 
}

/// \brief     ����һ�����ߵ���������
///            �������
/// \param     Curve                   �����������
/// \param     left,top,right,bottom   ���ߵ���������
/// \param     ww,hh,Area			   ���ߵĿ��ߺ����
/// \return ��
void GetCurveBound(vector<RPoint>&Curve,int&left,int&top,int&right,int&bottom,int&ww,int&hh,int&Area)
{
  int i;
 
  float x,y,xmin,xmax,ymin,ymax;
  int val;
  
  xmin = ymin = 9999;
  xmax = ymax = 0;
  /*-----  integrate results */
  for (i=0; i<Curve.size();i++)
  {  
	  x  = Curve[i].x; y = Curve[i].y;
	  if (xmin > x) xmin = x; if (xmax < x) xmax = x;
	  if (ymin > y) ymin = y; if (ymax < y) ymax = y;
  }   
  
  left   = xmin;  top    = ymin;
  right  = xmax;  bottom = ymax;
  
  ww   = right  - left + 1;
  hh   = bottom - top  + 1;
  Area = ww * hh;
  
}
#endif