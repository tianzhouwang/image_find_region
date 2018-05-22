 
#ifndef Region_Display_Head
#define Region_Display_Head

#define DISPLAY_INNER_OBJECT     1
#define DISPLAY_REGION_RECT      0
#define DISPLAY_REGION_NUMID     1
#define DISPLAY_REGION_CONVEX    0
#define DISPLAY_DBG_REGION_COLOR 1
#include <vector>
#include <algorithm>
//#include "hzmerge.h"

using namespace std;
/// \brief  从文件中读取文件的标记信息
/// \param  flename  图像文件名
/// \param  LabelVec 图像中区域的标记信息
/// \param  width    图像的宽
/// \param  height   图像的高
/// \param  ObjNum   所读取的区域个数
/// \return 无
void ReadRegionsDumpFile(char*filename,vector<int>&LabelVec,
						 int&Width,int&Height,int&ObjNum)
{ 
    int i;
	FILE*file = fopen(filename,"rb");
	
	fread(&Width , 4, 1, file);
	fread(&Height, 4, 1, file);
	fread(&ObjNum , 4, 1, file);
   
    int   TotalCnt = Width * Height;
	int * DumpImg  = new int[TotalCnt];

    fread((void*)DumpImg, TotalCnt, 4, file);
    fclose(file);
    LabelVec.resize(TotalCnt);

	Loopi(TotalCnt)
	{
       LabelVec[i] = DumpImg[i];
	}

	delete []DumpImg;
}

/// \brief  将区域的基本信息输出到每个单独的文件
///         (每个区域将绘制一个自己的文件)
/// \param  ReginVec 一组区域的向量
/// \param  dirNames 文件保存路径
/// \return 无
void PrintfRegionImg(vector<Region>&ReginVec, char*dirNames)
{
	int i,j;	 
	C24BitMap CPic;
	CPic.CleanPic(0);
	char NameBuff[100],regioninfo[100];

    Loopj(ReginVec.size())
	{
	    CPic.FormatF(ReginVec[i].rwidth, ReginVec[i].rheight);
		CPic.CleanPic(0); //  KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
		CPic.RandPenColor();
		//CPic.SetPenColor(BOUND(ReginVec[j].level_*3/2,0,255),0,0);
		//KBPic.RandPenColor();
		sprintf(NameBuff,"%s\\%06i.bmp",dirNames,ReginVec[j].r_id);
		//if(ReginVec[j].PtVec.size()<20)
		//	continue;
		Loopi(ReginVec[j].PtVec.size())
		{
			CPic.SigDot(ReginVec[j].PtVec[i].x + ReginVec[j].left, 
				ReginVec[j].PtVec[i].y + ReginVec[j].top);	 
		} 

		sprintf(regioninfo,"id:%i,area:%i,ratio:%.2lf,thr:%i",ReginVec[j].r_id,
			ReginVec[j].PtVec.size(),ReginVec[j].normal_variation,ReginVec[j].level_);
		//MergeTxtStr(CPic, 10,10,20, regioninfo,255,12,12);
		CPic.Save(NameBuff);
	}
	
 
}

/// \brief  在图像文件显示输出一组MSER区域
///         将绘制区域的轮廓拐点
/// \param  ObjVec    一组MSEG区域 
/// \param  width     图像宽
/// \param  height    图像高
/// \param  filename  输出的图像文件名
/// \return 无
void DispReionObjvec1(vector<Region>&ObjVec,int width,int height,char*filename )
{	
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
	int i,j,t;
	Loopj(ObjVec.size())
	{
		//AnalysisObj(ObjVec[j]);
		//if(ObjVec[j].IfValid==-1)continue;
        char buff[30];
		sprintf(buff,"%i",j);

		KBPic.RandPenColor();
		KBPic.RandPenColor();
        //KBPic.Rectangle(ObjVec[j].left,ObjVec[j].top,ObjVec[j].right,ObjVec[j].bottom);


		Loopi(ObjVec[j].PtVec.size())
		{
		//	KBPic.SigDot(ObjVec[j].PtVec[i].x,ObjVec[j].PtVec[i].y);
		}

		//MergeTxtStr(KBPic, ObjVec[j].GeoX,ObjVec[j].GeoY,20, buff,255,0,0); 

        //KBPic.SetPenColor(255, 255, 0);
        //KBPic.SigDot( (ObjVec[j].left + ObjVec[j].right )/2,
			 //         (ObjVec[j].top  + ObjVec[j].bottom)/2 );

		for(t=0;t<ObjVec[j].MContours.size();t++)
		{
			if(ObjVec[j].MContours[t].size()<20)
				continue;
			KBPic.RandPenColor();
			Loopi(ObjVec[j].MContours[t].size())
			{
				RPoint Pt = ObjVec[j].MContours[t][i];
				
				KBPic.SigDot(Pt.x,Pt.y);

				if(i%12==0)
				{
                  KBPic.SigDot(Pt.x-1,Pt.y);
				  KBPic.SigDot(Pt.x+1,Pt.y);
				  KBPic.SigDot(Pt.x,Pt.y-1);
				  KBPic.SigDot(Pt.x,Pt.y+1);

				  KBPic.SigDot(Pt.x-2,Pt.y);
				  KBPic.SigDot(Pt.x+2,Pt.y);
				  KBPic.SigDot(Pt.x,Pt.y-2);
				  KBPic.SigDot(Pt.x,Pt.y+2);
				}
			}
		}
	
	}

	KBPic.Save(filename);
}



/// \brief  在图像文件显示一个单一Region
/// \param  Obj       待显示的Region
/// \param  width     图像宽
/// \param  height    图像高
/// \param  filename  输出的图像文件名
/// \return 无
void DispReionInPic( Region &Obj ,int width,int height,char*filename)
{	
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
	int i,j,t;
	
	KBPic.RandPenColor();	
	KBPic.SetPenColor(Obj.RR, Obj.GG, Obj.BB);
	
	//KBPic.Rectangle(Obj.left,Obj.top,Obj.right,Obj.bottom);

	if(DISPLAY_INNER_OBJECT)
	Loopi(Obj.PtVec.size())
	{
		KBPic.SigDot(Obj.PtVec[i].x,Obj.PtVec[i].y);
	}
	
	
	KBPic.SetPenColor(255, 255, 0);
	KBPic.SigDot( (Obj.left + Obj.right )/2,
		(Obj.top  + Obj.bottom)/2 );
	
	int Cnt = Obj.MContours.size();
	for(t=0;t<Obj.MContours.size();t++)
	{
		KBPic.RandPenColor();
		Loopi(Obj.MContours[t].size())
		{
			RPoint Pt = Obj.MContours[t][i];
			KBPic.SigDot(Pt.x,Pt.y);
		}
	}
	
	KBPic.Save(filename);
}


//计算区域的八边形轮廓 ----------------------------------
void GetBountpt( Region &Obj, double Dx,double Dy, 
				 double&left, double&top, double&right, double&bottom)
{
	int i;
	
	vector<double> RotX;
	vector<double> RotY;
	
	Loopi( Obj.ContourPtVec.size() )
	{ 
	  double NewX, NewY;
	 
	  NewX = double( Obj.ContourPtVec[i].x - Obj.GeoX) * Dx   + double( Obj.ContourPtVec[i].y - Obj.GeoY) *Dy;
	  NewY = double( Obj.ContourPtVec[i].x - Obj.GeoX) *(-Dy) + double( Obj.ContourPtVec[i].y - Obj.GeoY) *Dx;

	  RotX.push_back(NewX);
	  RotY.push_back(NewY);

	}
	
	left  = * min_element(RotX.begin(),RotX.end());
	right = * max_element(RotX.begin(),RotX.end());
	
	top   = * min_element(RotY.begin(),RotY.end());
	bottom= * max_element(RotY.begin(),RotY.end());
}

/// \brief  显示矩形外接框向量
/// \param  RectVec   一组矩形框向量
/// \param  width     图像宽
/// \param  height    图像高
/// \param  filename  输出的图像文件名
/// \param  selObjCr  所选的显示颜色
/// \return 无
void DispImgRectvec(vector<ImgRECT>& RectVec,int width,int height,char*filename,int selObjCr=0)
{
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
    int i;
	switch (selObjCr)
	{
	case 0: KBPic.PenColor.R = 254; KBPic.PenColor.G =   0; KBPic.PenColor.B = 244;
		break;
	case 1: KBPic.PenColor.R = 150; KBPic.PenColor.G =   0; KBPic.PenColor.B = 253;
		break;
	case 2: KBPic.PenColor.R =   0; KBPic.PenColor.G =   0; KBPic.PenColor.B = 255;
		break;
	case 3: KBPic.PenColor.R =   0; KBPic.PenColor.G = 255; KBPic.PenColor.B = 255;
		break;
	case 4: KBPic.PenColor.R =   0; KBPic.PenColor.G = 255; KBPic.PenColor.B = 128;
		break;
	case 5: KBPic.PenColor.R =   0; KBPic.PenColor.G = 255; KBPic.PenColor.B = 0;
		break;
	case 6: KBPic.PenColor.R = 128; KBPic.PenColor.G = 255; KBPic.PenColor.B = 0;
		break;
	case 7: KBPic.PenColor.R = 255; KBPic.PenColor.G = 255; KBPic.PenColor.B = 0;
		break;
	default:
		KBPic.PenColor.R = 255;KBPic.PenColor.G = 0;KBPic.PenColor.B = 0;
		break;
	}

    Loopi(RectVec.size())
	{
      KBPic.Rectangle(RectVec[i].left, RectVec[i].top, RectVec[i].right,RectVec[i].bottom);
	}

	KBPic.Save(filename);
}


void DrawTiltedRectangle(C24BitMap&CPic,
						 double  CenterX,  double CenterY,
						 double  rot_left , double  rot_top,
						 double  rot_right, double  rot_bottom)
{
  int XX[4];
  int YY[4];

  XX[0] = CenterX + rot_left  * 0.707 + rot_top     * 0.707 ;
  XX[1] = CenterX + rot_left  * 0.707 + rot_bottom  * 0.707 ;
  XX[2] = CenterX + rot_right * 0.707 + rot_bottom  * 0.707 ;
  XX[3] = CenterX + rot_right * 0.707 + rot_top     * 0.707 ;
                                                             
  YY[0] = CenterY + rot_left  * -0.707 + rot_top     * 0.707 ;
  YY[1] = CenterY + rot_left  * -0.707 + rot_bottom  * 0.707 ;
  YY[2] = CenterY + rot_right * -0.707 + rot_bottom  * 0.707 ;
  YY[3] = CenterY + rot_right * -0.707 + rot_top     * 0.707 ;

  CPic.DrawLine(XX[0], YY[0], XX[1],YY[1]);
  CPic.DrawLine(XX[1], YY[1], XX[2],YY[2]);
  CPic.DrawLine(XX[2], YY[2], XX[3],YY[3]);
  CPic.DrawLine(XX[3], YY[3], XX[0],YY[0]);
}

double GetRegionRadius(Region&R)
{
  int i;
  double distance=0;

  Loopi(R.ContourPtVec.size())
  {
	double dx = R.ContourPtVec[i].x - R.GeoX;
	double dy = R.ContourPtVec[i].y - R.GeoY;

    double dis = sqrt(dx*dx + dy*dy);
	if(dis>distance)
		distance = dis;
  }

  return distance;
}

#include "convexhull.h"


void DrawRegionConvexhull(C24BitMap&CPic, Region&R)
{

 CMyConvexHull  ConvexHull;
 int ptstep = 3;
 int size = R.ContourPtVec.size()/ptstep;
 
 ConvexHull.SetNodeNum(size);

 int i;
 for(i=0;i<size;i++)
 {
	ConvexHull.node[i].x = R.ContourPtVec[i*ptstep].x;
	ConvexHull.node[i].y = R.ContourPtVec[i*ptstep].y;
    //CPic.DrawCircle(x,y,2);
 }

 ConvexHull.Solve();

  
	int x1,y1,x2,y2;
	int x0,y0;

	x0  = ConvexHull.node[ConvexHull.V[ConvexHull.back]].x;
	y0  = ConvexHull.node[ConvexHull.V[ConvexHull.back]].y;
   
	int Cnt =0;
	for(i=ConvexHull.back+1;i<ConvexHull.front;i++){
		
		
		x1 = ConvexHull.node[ConvexHull.V[i-1]].x;
		y1 = ConvexHull.node[ConvexHull.V[i-1]].y;
		x2 = ConvexHull.node[ConvexHull.V[i]].x;
		y2 = ConvexHull.node[ConvexHull.V[i]].y;
		
		Cnt++;
		CPic.DrawTkLine(x1,y1,x2,y2,2);
		//printf("%i,%i,%i,%i\n", x1,y1,x2,y2);
	}
   
	printf("%i\n",Cnt);
	//x1 = ConvexHull.node[ConvexHull.V[ConvexHull.back]].x;
	//y1 = ConvexHull.node[ConvexHull.V[ConvexHull.back]].x;

	CPic.DrawTkLine(x2,y2, x0, y0,2);

	//CPic.Save("res.bmp");
}

void DrawRegionConvexhullA(C24BitMap&CPic, Region&R)
{
   int i;

   for(i=0;i<R.ConvexHullPtVec.size()-1;i++)
   {
	   CPic.DrawDashLine(R.ConvexHullPtVec[i  ].x, R.ConvexHullPtVec[i  ].y,
		               R.ConvexHullPtVec[i+1].x, R.ConvexHullPtVec[i+1].y,2);
   }

   if(R.ConvexHullPtVec.size()>2)
   {
     int Last = R.ConvexHullPtVec.size() - 1;
     CPic.DrawDashLine(R.ConvexHullPtVec[0].x,    R.ConvexHullPtVec[0].y,
		             R.ConvexHullPtVec[Last].x, R.ConvexHullPtVec[Last].y,2);
   }
}



/// \brief  在图像文件显示输出一组MSER区域
/// \param  ObjVec    一组MSEG区域 
/// \param  width     图像宽
/// \param  height    图像高
/// \param  filename  输出的图像文件名
/// \param  selObjCr  用随机颜色还是区域自身的颜色绘制区域
/// \return 无


void DispReionObjvec(vector<Region>&ObjVec,int width,int height,char*filename)
{	
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
	int i,j,t;

	Loopj(ObjVec.size())
	{
	  KBPic.RandPenColor();

      double rot_left, rot_top, rot_right, rot_bottom;
 
	  GetBountpt( ObjVec[j], 0.707, -0.707, 
				   rot_left, rot_top, rot_right, rot_bottom);
   
	  //DrawTiltedRectangle(KBPic,
		//  ObjVec[j].GeoX,  ObjVec[j].GeoY,
		//  rot_left, rot_top, rot_right, rot_bottom);

     /*  double c_radius =  GetRegionRadius(ObjVec[j]);

	   KBPic.DrawCircleLine(ObjVec[j].GeoX,ObjVec[j].GeoY,
		                    c_radius);*/
       //if(ObjVec[j].MContours.size()>20) 
		//   continue;

 
	   if(DISPLAY_REGION_CONVEX)
       DrawRegionConvexhullA(KBPic, ObjVec[j]);

		//if(ObjVec[j].rwidth*2>width)continue;
		//if(ObjVec[j].IsValid ==0)continue;
		//AnalysisObj(ObjVec[j]);
		//if(ObjVec[j].IfValid==-1)continue;
        char buff[30];
		sprintf(buff,"%i",j);

		
        if(DISPLAY_DBG_REGION_COLOR)
		   KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
        //if(ObjVec[j].ContourPtVec.size()==0)
		//	continue;
		//if(float(ObjVec[j].InflectionPts.size())/float(ObjVec[j].ContourPtVec.size())>0.4)
         //   KBPic.SetPenColor(0, 255, 255);
		if(DISPLAY_REGION_RECT)
		    KBPic.Rectangle(ObjVec[j].left,ObjVec[j].top,ObjVec[j].right,ObjVec[j].bottom);

		if(DISPLAY_INNER_OBJECT)
		Loopi(ObjVec[j].PtVec.size())
		{
			KBPic.SigDot(ObjVec[j].PtVec[i].x,ObjVec[j].PtVec[i].y);
		}

	

        KBPic.SetPenColor(255, 255, 0);
        KBPic.SigDot( (ObjVec[j].left + ObjVec[j].right )/2,
			          (ObjVec[j].top  + ObjVec[j].bottom)/2 );
				
		int Cnt = ObjVec[j].MContours.size();

        for(t=0;t<ObjVec[j].MContours.size();t++)
		{
			KBPic.RandPenColor();
			Loopi(ObjVec[j].MContours[t].size())
			{
				RPoint Pt = ObjVec[j].MContours[t][i];
				KBPic.SigDot(Pt.x,Pt.y);
			}
		}

		if(DISPLAY_REGION_NUMID)
		MergeTxtStr(KBPic, ObjVec[j].GeoX,ObjVec[j].GeoY,15, buff, 255, 0, 112); 
		/*KBPic.PenColor.R = BOUND(int(KBPic.PenColor.R)*15/10,0,255);
		KBPic.PenColor.G = BOUND(int(KBPic.PenColor.G)*15/10,0,255);
		KBPic.PenColor.B = BOUND(int(KBPic.PenColor.B)*15/10,0,255);

		Loopi(ObjVec[j].ContourPtVec.size())
		{
           KBPic.SigDot(ObjVec[j].ContourPtVec[i].x,ObjVec[j].ContourPtVec[i].y);
		}

		KBPic.PenColor.R = 255; //BOUND(int(KBPic.PenColor.R)*8/10,0,255);
		KBPic.PenColor.G =   0; //BOUND(int(KBPic.PenColor.G)*8/10,0,255);
		KBPic.PenColor.B =   0; //BOUND(int(KBPic.PenColor.B)*8/10,0,255);

		Loopi(ObjVec[j].InflectionPts.size())
		{
           KBPic.SigDot(ObjVec[j].InflectionPts[i].x,ObjVec[j].InflectionPts[i].y);
		}
		*/
	}

	KBPic.Save(filename);
}

/// \brief  在图像显示绘制一组MSER区域
/// \param  KBPic     所画的背景图像文件
/// \param  ObjVec    一组MSEG区域 
/// \param  selObjCr  用随机颜色还是区域自身的颜色绘制区域
/// \return 无
void DispReionObjvec(C24BitMap&KBPic,vector<Region>&ObjVec, int selObjCr=0)
{	
	int i,j,t;
	Loopj(ObjVec.size())
	{
		//AnalysisObj(ObjVec[j]);
		//if(ObjVec[j].IfValid==-1)continue;
		KBPic.RandPenColor();
        if(selObjCr)
			KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
	
		Loopi(ObjVec[j].PtVec.size())
		{
			KBPic.SigDot(ObjVec[j].PtVec[i].x,ObjVec[j].PtVec[i].y);
		}/**/
		
		if(DISPLAY_REGION_RECT)
		KBPic.Rectangle(ObjVec[j].left,ObjVec[j].top,ObjVec[j].right,ObjVec[j].bottom);

		KBPic.PenColor.R = BOUND(int(KBPic.PenColor.R)*15/10,0,255);
		KBPic.PenColor.G = BOUND(int(KBPic.PenColor.G)*15/10,0,255);
		KBPic.PenColor.B = BOUND(int(KBPic.PenColor.B)*15/10,0,255);
		
	/*	Loopi(ObjVec[j].ContourPtVec.size())
		{
			KBPic.SigDot(ObjVec[j].ContourPtVec[i].x,ObjVec[j].ContourPtVec[i].y);
		}
		*/
         
		Loopi(ObjVec[j].MContours.size())
		{ KBPic.RandPenColor();
		  vector<RPoint>  PtVec;
		 PtVec = ObjVec[j].MContours[i];
          //for(t=0;t<ObjVec[j].MContours[i].size();t++)
		  for(t=0;t<PtVec.size();t++)
		  {
           KBPic.SigDot(PtVec[t].x,PtVec[t].y);
		  }
		}

		KBPic.PenColor.R = 255;//BOUND(int(KBPic.PenColor.R)*8/10,0,255);
		KBPic.PenColor.G =   0;//BOUND(int(KBPic.PenColor.G)*8/10,0,255);
		KBPic.PenColor.B =   0;//BOUND(int(KBPic.PenColor.B)*8/10,0,255);
		
		/*Loopi(ObjVec[j].InflectionPts.size())
		{
			KBPic.SigDot(ObjVec[j].InflectionPts[i].x,ObjVec[j].InflectionPts[i].y);
		}
		*/
	}
	
}
	

void PrintfRegionVec(vector<Region>&RVec,char*filename);

#endif