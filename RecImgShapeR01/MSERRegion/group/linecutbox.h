 

#ifndef CUT_TEXTBOX_FROM_LINE_HEAD
#define CUT_TEXTBOX_FROM_LINE_HEAD

#include "GetRectRotAngle.h"

#define LINEBOX_REGION_DEBUG 0
#define GET_MIN_MAX(val, min_v, max_v ) {if(min_v>(val)) min_v = val; if(max_v <(val))max_v = val;}
#define SWAP_VAL(x,y)  { x = (x)+(y); y =(x)-(y); x = (x)-(y); } 


/*
向量结构体
*/
struct vector2D
{
	double Dx,Dy;
};

//对向量进行规范化
void NormalizeVector(vector2D&vec)
{
 double VecLength = vec.Dx * vec.Dx + vec.Dy * vec.Dy;
 VecLength = sqrt(VecLength);
 vec.Dx/= VecLength;
 vec.Dy/= VecLength;
}

//得到对象的连接
void GetObjLink(FastSelCoLinearRegion&RegionCoAngleInfo,
				int Rg1Idx, int Rg2Idx, vector<Region>&RgVec, C24BitMap&CPic)
{
  int i,j;
  vector<int> IdxVec;
  double DisAverage;

  GetInSideRegion(RegionCoAngleInfo,Rg1Idx, Rg2Idx, RgVec, IdxVec, DisAverage);
  printf("test zero %i\n", IdxVec.size());


  CPic.PenColor.R = 0; CPic.PenColor.G = 255; CPic.PenColor.B = 0;  
  
  vector<Region> ObjVec;
  ObjVec.clear();

  Loopi(IdxVec.size())
  { 
	 CPic.Rectangle(RgVec[IdxVec[i]].GeoX-2,RgVec[IdxVec[i]].GeoY-2,RgVec[IdxVec[i]].GeoX+2,RgVec[IdxVec[i]].GeoY+2);
	 CPic.Rectangle(RgVec[IdxVec[i]].left,RgVec[IdxVec[i]].top,RgVec[IdxVec[i]].right,RgVec[IdxVec[i]].bottom);
     ObjVec.push_back(RgVec[IdxVec[i]]);

	//CPic.Rectangle(OutRects[i].left  , OutRects[i].top   , 
	//	             OutRects[i].right , OutRects[i].bottom);
  }

  ImgRECT  RectBound;

  
  GetOutterRect(ObjVec,  RectBound);

  CPic.Rectangle(RectBound.left, RectBound.top, RectBound.right, RectBound.bottom);
  //CPic.Save("link0.bmp");

  ObjVec.clear();
  

  Loopi(RgVec.size())
  {
    if(RegionInRect(RgVec[i], RectBound))
	 {
		if(ColorDist(RgVec[i].RR,RgVec[i].GG,RgVec[i].BB,
			         RgVec[Rg1Idx].RR, RgVec[Rg1Idx].GG,  RgVec[Rg1Idx].BB)<40*40*3)
       ObjVec.push_back(RgVec[i]);
	   printf("~ %i , ",i);
	 }
  } 

  
  GetTextLineComptChain txt_chain;
  vector<ImgRECT> RectVec, OutRects;

  int MHeight;
    
  MHeight = RectBound.bottom - RectBound.top + 1;

  transObjsToRects( ObjVec, RectVec);
   

  txt_chain.Pic.FormatF(CPic.Width,CPic.Height);
  txt_chain.Pic.Clear();
  
  Loopi(RgVec.size())
  {
	  Loopj(RgVec[i].PtVec.size())
      {
		  * get_pix_color(txt_chain.Pic,
			  RgVec[i].PtVec[j].x,
			  RgVec[i].PtVec[j].y) =0;
	  }
  }


  txt_chain.OverLayPix = 0;
  txt_chain.GetHorizontalTextChain( RectVec, OutRects, MHeight);/**/

  //OutRects = RectVec;
  sort(OutRects.begin(),OutRects.end(),LessByX);
   
  /*Loopi(OutRects.size())
  {
	  CPic.Rectangle(OutRects[i].left  , OutRects[i].top   , 
		             OutRects[i].right , OutRects[i].bottom);
  }*/
  
 // CPic.Save("link0.bmp");

  double MeanDistance;
  double lk_uniformity =
	  GetObjLinkUniformity(OutRects, MeanDistance);

 
  printf("\n%.5lf,mean_dis:%.2lf",lk_uniformity,MeanDistance);

  double height_uniformity 
	  = GetObjHeightUniformity(OutRects, MeanDistance, RectBound.top);

  printf("\n%.5lf", height_uniformity);

  //CPic.Save("link.bmp");  
}

//│AXo+BYo+C│/√（A2+B2）
// Y = KX + B ---> KX - Y + B = 0 
float PointToLineDistance(float X,float Y,float K,float B)
{
	double Denominator = sqrt( K * K + 1.0);
	double Numerator   = fabs( K* X - Y +B);
	return Numerator/Denominator;
}

int GetLineFitPosition(vector<int>&Link, vector<Region>&RegionVec,
						double & x1, double &y1,
						double & x2, double &y2, int SepLine=1)
{
	
	if(Link.size()<3)
		return 0;

    int i,j;
    vector<Point2D>  PointVec;
	int LinkMinX, LinkMaxX, LinkMinY, LinkMaxY;
	LinkMinX = LinkMinY =  999;
	LinkMaxX = LinkMaxY = -999;

	int ReverseXY = 0;
	Loopi(Link.size())
	{
		Point2D Pt;
		int Idx= Link[i];
		Pt.x = RegionVec[Idx].x; 
		Pt.y = RegionVec[Idx].y;
		Pt.Weight = RegionVec[Idx].PtVec.size();
        PointVec.push_back(Pt);

		GET_MIN_MAX(Pt.x, LinkMinX, LinkMaxX);
        GET_MIN_MAX(Pt.y, LinkMinY, LinkMaxY);
	}
    
	if((LinkMaxY-LinkMinY) > (LinkMaxX-LinkMinX))
	{
		ReverseXY = 1;
		 
		Loopi(PointVec.size())
		{
			SWAP_VAL(PointVec[i].x,PointVec[i].y);	 
		}
	}

	double k, b, meanx, meany;
	CalculateLineKB( PointVec, k, b, meanx, meany);

	double SumDis =0;
    Loopi(PointVec.size()) 
	{
	  SumDis += PointToLineDistance(PointVec[i].x,PointVec[i].y, k, b);
	}
     SumDis/= double(PointVec.size());

   if(SumDis>5.0 && SepLine)
	 {double Angle  = 
		GetRotAngle(PointVec, 3.0);
	Angle =   Angle / 180 * M_PI ;
	double Dx_,Dy_;
	Dx_ =  cos(Angle); Dy_ = - sin(Angle);
    k = Dy_/Dx_;}

	double normal_len = sqrt( k * k + 1.0 );

	vector<double> ProjLen;

	Loopi(Link.size())
	{
		Point2D Pt; 
		Pt = PointVec[i];
        
		Pt.x -=meanx;
		Pt.y -=meany;
        
		double proj_val = Pt.x * 1.0 + Pt.y * k;
	    ProjLen.push_back(proj_val);
	}

	double MinProjV, MaxProjV;

	MinProjV = *min_element(ProjLen.begin(), ProjLen.end());
	MaxProjV = *max_element(ProjLen.begin(), ProjLen.end());

	x1 = meanx + (MinProjV/normal_len) * 1 / normal_len;
    y1 = meany + (MinProjV/normal_len) * k / normal_len;

	x2 = meanx + (MaxProjV/normal_len) * 1 / normal_len;
    y2 = meany + (MaxProjV/normal_len) * k / normal_len;


	if(ReverseXY)
	{
		SWAP_VAL(x1,y1); 
		SWAP_VAL(x2,y2);
	}

	return 1;
}

#endif