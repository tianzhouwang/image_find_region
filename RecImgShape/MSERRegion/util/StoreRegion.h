#ifndef StoreRegion_Head
#define StoreRegion_Head
#include "mregion.h"


//\brief �������л�MSER���м���ʱ�ṹ��
struct RegionInfo
{
	int level_;     ///< Level at which the region is processed.
	int area_;      ///< Area of the region (moment zero).
	int childNum;
	double normal_variation;
	int   left,right,top,bottom;       // boundingbox left,right,top,bottom
	int   rwidth,rheight;
	//-----------------�����һЩ����-------------------------------------
	float  x,y, GeoX, GeoY;
	float     RR,GG,BB;
	RgFeature Feature;
};

/*\brief ��������
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
};*/


int CalcuRegionSize(Region&R)
{
 // int 
	return 1;
}
/*struct RgFeature
{
	int   perimeter;     // perimeter �ܳ�
	int   InnerHoles;    //euler_no;      // euler no  �ڲ��ն���
	float postProb;      // posterior prob
	float StrokeMean,StrokeStd,StkMin,StkMax,StkMedian;    //
	//float StrokeStd;     //
	int   label;           // classified label
	float  aspectRatio;    // aspect ratio
	//vector<float> OccupyVec;
	float AngleCnt[4];
	double roundratio;
};*/


/// \brief  RGBɫ�ʵ�HSVɫ�ʵ�ת��
/// \param  R,G,B  ����ɫ�ʵ�RGB��ֵ
/// \param  H,S,V  �õ�ӳ����HSV��ֵ
/// \return ��
void RGBToHSV(float R,float G, float B, float &H, float &S, float &V)
{
	float maxv = max(max(R,G),B); float minv = min(min(R,G),B);
	
	if( maxv - minv < 1)  {H = 0;} 
	else
	{
		if (R == maxv)  H = (G-B)/(maxv-minv);	
		if (G == maxv)  H = 2 + (B-R)/(maxv-minv);
		if (B == maxv)  H = 4 + (R-G)/(maxv-minv);	 
	}
	
	H = H * 60;  if(H < 0) H = H + 360;
	V = maxv;
    S  =  ( (maxv < 1) ?  0 : (maxv-minv)/maxv) ;
 	
	return;
}


/// \brief  ��ͼ������ȡ�������ɫ��Ϣ  
/// \param  CPic   ԭͼ
/// \param  RgVec  ������Ϣ����
/// \return ��
void  GetObjColor(C24BitMap&CPic,vector<Region>&RgVec)
{
	int i,j;	
	Loopi(RgVec.size())
	{
		RgVec[i].Feature.label = i;
		double RRsum,GGsum,BBsum;
		RRsum = GGsum = BBsum = 0;
		
		Loopj(RgVec[i].PtVec.size())
		{
			C24PixVal Pix;
			Pix = get_pix_color(CPic,RgVec[i].PtVec[j].x,RgVec[i].PtVec[j].y);
			RRsum += *Pix.r;
			GGsum += *Pix.g;
			BBsum += *Pix.b;
		}

		//<----------------����ע�ͼ������Ƿ��������ı߽�ɫ��-------------
		//Loopj(RgVec[i].ContourPtVec.size())
		//{
		//	C24PixVal Pix;
		//	Pix = get_pix_color(CPic,RgVec[i].ContourPtVec[j].x,RgVec[i].ContourPtVec[j].y);
		//	RRsum -= *Pix.r;
		//	GGsum -= *Pix.g;
		//	BBsum -= *Pix.b;
		//}
		//int num = RgVec[i].PtVec.size() - RgVec[i].ContourPtVec.size();
		int num = RgVec[i].PtVec.size();
		RgVec[i].RR = RRsum / float(num); 
		RgVec[i].GG = GGsum / float(num);
		RgVec[i].BB = BBsum / float(num);

		RGBToHSV(RgVec[i].RR, RgVec[i].GG, RgVec[i].BB, RgVec[i].HH, RgVec[i].SS, RgVec[i].VV);
		

		RgVec[i].x = (RgVec[i].left + RgVec[i].right ) /2;
		RgVec[i].y = (RgVec[i].top  + RgVec[i].bottom) /2;
		RgVec[i].rwidth  = RgVec[i].right  - RgVec[i].left + 1;
		RgVec[i].rheight = RgVec[i].bottom - RgVec[i].top  + 1;
	}
	
}




/// \brief  ��������MSER���浽�ļ�
/// \param  PicWidth   ԭͼ�Ŀ��
/// \param  PicHeight  ԭͼ�ĸ߶�
/// \param  ObjVec     ������Ϣ����
/// \param  filename   ��������ļ���Ϣ
/// \return �� 
void SaveMRegionFile(int PicWidth,int PicHeight, vector<Region>&ObjVec ,char*filename)
{
 
	int i,j;
	int * DumpImg =  new int[PicWidth * PicHeight];
				
	int TotalCnt = PicWidth * PicHeight;
	Loopi(TotalCnt) DumpImg[i] = -1;

	RegionInfo*RgPt = new RegionInfo[ObjVec.size()];
	
	Loopi(ObjVec.size())
	{
		Loopj(ObjVec[i].PtVec.size())
		{
			int xx,yy;
			xx = ObjVec[i].PtVec[j].x ;
			yy = ObjVec[i].PtVec[j].y ;
			int Pos;
			Pos = yy * PicWidth + xx;
			DumpImg[Pos] = i;
		}

   RgPt[i].level_           =  ObjVec[i].level_           ;            
   RgPt[i].area_            =  ObjVec[i].area_            ;
   RgPt[i].childNum         =  ObjVec[i].childNum         ;
   RgPt[i].normal_variation =  ObjVec[i].normal_variation ;
   RgPt[i].left             =  ObjVec[i].left             ;
   RgPt[i].right            =  ObjVec[i].right            ;
   RgPt[i].top              =  ObjVec[i].top              ;
   RgPt[i].bottom           =  ObjVec[i].bottom           ;
   RgPt[i].rwidth           =  ObjVec[i].rwidth           ;
   RgPt[i].rheight          =  ObjVec[i].rheight          ;
   RgPt[i].x                =  ObjVec[i].x                ;
   RgPt[i].y                =  ObjVec[i].y                ;
   RgPt[i].GeoX             =  ObjVec[i].GeoX             ;
   RgPt[i].GeoY             =  ObjVec[i].GeoY             ;
   RgPt[i].RR               =  ObjVec[i].RR               ;
   RgPt[i].GG               =  ObjVec[i].GG               ;
   RgPt[i].BB               =  ObjVec[i].BB               ;
   RgPt[i].Feature 	        =  ObjVec[i].Feature 	      ;
	
	}
				
	int ObjSize = ObjVec.size();
	FILE*file = fopen(filename,"wb+");
	fwrite(&PicWidth , 4, 1, file);
	fwrite(&PicHeight, 4, 1, file);
	fwrite(&ObjSize, 4, 1, file);
	fwrite((void*)DumpImg,PicWidth * PicHeight, 4, file);
	fwrite((void*)RgPt, ObjVec.size(),sizeof(RegionInfo),file )	;		
	
	Loopi(TotalCnt) DumpImg[i] = -1;
	
	Loopi(ObjVec.size())
	{
		Loopj(ObjVec[i].ContourPtVec.size())
		{
			int xx,yy;
			xx = ObjVec[i].ContourPtVec[j].x ;
			yy = ObjVec[i].ContourPtVec[j].y ;
			int Pos;
			Pos = yy * PicWidth + xx;
			DumpImg[Pos] = i;
		}
	}

    fwrite((void*)DumpImg,PicWidth * PicHeight, 4, file);
	fclose(file);

	delete [] RgPt;
	delete [] DumpImg;
}


/// \brief  ���ļ��ж��������MSER����
/// \param  PicWidth   ԭͼ�Ŀ��
/// \param  PicHeight  ԭͼ�ĸ߶�
/// \param  ObjVec     ������Ϣ����
/// \param  filename   ��������ļ���Ϣ
/// \return �� 
void LoadMRegionFile(int&PicWidth,int&PicHeight, vector<Region>&ObjVec ,char*filename)
{
	ObjVec.clear();
	int i,j;
	
	int ObjNum ;
	FILE*file = fopen(filename,"rb");
	
	fread(&PicWidth , 4, 1, file);
	fread(&PicHeight, 4, 1, file);
	fread(&ObjNum   , 4, 1, file);
	
	int TotalCnt = PicWidth * PicHeight;
	int * DumpImg =  new int[PicWidth * PicHeight];
	
	fread((void*)DumpImg, PicWidth * PicHeight, 4, file);
	RegionInfo*RgPt = new RegionInfo[ObjNum];
    fread((void*)RgPt, ObjNum, sizeof(RegionInfo),file );
	
    ObjVec.resize(ObjNum);
		  
   Loopi(TotalCnt)
   {
    if(DumpImg[i] != -1)
	{int xx,yy;
	  xx = i % PicWidth;
	  yy = i / PicWidth;
	 RPoint Pt;
	 Pt.x = xx;
     Pt.y = yy;
     ObjVec[DumpImg[i]].PtVec.push_back(Pt);
	}
 }
	

 Loopi(ObjVec.size())
	{
    ObjVec[i].level_           = RgPt[i].level_          ;            
    ObjVec[i].area_            = RgPt[i].area_           ;
    ObjVec[i].childNum         = RgPt[i].childNum        ;
    ObjVec[i].normal_variation = RgPt[i].normal_variation;
    ObjVec[i].left             = RgPt[i].left            ;
    ObjVec[i].right            = RgPt[i].right           ;
    ObjVec[i].top              = RgPt[i].top             ;
    ObjVec[i].bottom           = RgPt[i].bottom          ;
    ObjVec[i].rwidth           = RgPt[i].rwidth          ;
    ObjVec[i].rheight          = RgPt[i].rheight         ;
    ObjVec[i].x                = RgPt[i].x               ;
    ObjVec[i].y                = RgPt[i].y               ;
    ObjVec[i].GeoX             = RgPt[i].GeoX            ;
    ObjVec[i].GeoY             = RgPt[i].GeoY            ;
    ObjVec[i].RR               = RgPt[i].RR              ;
    ObjVec[i].GG               = RgPt[i].GG              ;
    ObjVec[i].BB               = RgPt[i].BB              ;
    ObjVec[i].Feature 	       = RgPt[i].Feature 	     ;
 }	

 fread((void*)DumpImg, PicWidth * PicHeight, 4, file);

 Loopi(TotalCnt)
 {
	 if(DumpImg[i] != -1)
	 {int xx,yy;
	 xx = i % PicWidth;
	 yy = i / PicWidth;
	 RPoint Pt;
	 Pt.x = xx;
     Pt.y = yy;
     ObjVec[DumpImg[i]].ContourPtVec.push_back(Pt);
	 }
 }
 delete [] DumpImg; 
 fclose(file);
}

#endif