 
#ifndef CalcuStk_Head
#define CalcuStk_Head
#include "mregion.h"
#include "TraceBound.h"
#include "convexhull.h"
//#include "StoreRegion.h"
//#include "StoreRegion.h"


/// \brief  删除内部为空或者边界有问题的区域
/// \param  RVec MSER区域向量 
/// \return 无 处理完以后在原基础上删除
void DeleteNullRegion(vector<Region>&RVec)
{
	int i;
	vector<Region> BkVec;
	Loopi(RVec.size())
	{
		if(RVec[i].PtVec.size()==0||RVec[i].ContourPtVec.size()==0)
		{
			continue;
		}
		//if(RVec[i].PtVec.size()!=0)
		BkVec.push_back(RVec[i]);
	}
	RVec = BkVec;
}


/// \brief  得到Region划好网格后区域在每个网格内的占比
/// \param  R 待计算的区域
/// \return 无 计算完后自动更新区域的参数
void GetOccupyFeature(Region&R)
{
	int i,j;
	
	double rectheight,v1,v2,v3;
	int Cnt[6];
	Loopj(6) Cnt[j] =0;
	rectheight = R.bottom - R.top +1;
	v1 = R.top + rectheight/3;
	v2 = v1 + rectheight/3;
	
	Loopj(R.PtVec.size())
	{
		int y = R.PtVec[j].y;
		
		if(float(R.PtVec[j].y) > v2)
			Cnt[0]+=1;
		else if(float(R.PtVec[j].y) > v1)  
			Cnt[1]+=1;
		else
			Cnt[2]+=1;
		
	}
	
	//Loopj(3)
	//	R.Feature.OccupyVec.push_back(float(Cnt[j])/float(R.PtVec.size()));
	
	Loopj(R.ContourPtVec.size())
	{
		int y = R.ContourPtVec[j].y;
		if(float(R.ContourPtVec[j].y) > v2)
			Cnt[3]+=1;
		else if(float(R.ContourPtVec[j].y) > v1)  
			Cnt[4]+=1;
		else
			Cnt[5]+=1;
		
	}
	
	//Loopj(3)
	//	R.Feature.OccupyVec.push_back(float(Cnt[j+3])/float(R.ContourPtVec.size()));
	
}

/// \brief     将区域所计算的特征输出到文件中
/// \param     RVec  待计算的一组区域
/// \filename  所保存的文件名
/// \return 无  
void  PrintfRegionVec(vector<Region>&RVec,char*filename)
{
	int i;
	FILE*file;
	file = fopen(filename, "wt+");
	
	fprintf(file,"Idx,左,上,右,下,宽,高, centX, centY, 红,绿,蓝, 面积,轮廓面积比,面积凸包比,圆心度,笔画平均,笔画方差,笔画中值,笔画最小,笔画最大,空洞数,空洞比,周长面积比\n");
	
	
	for (i=0;i<RVec.size();i++)
	{
		GetOccupyFeature(RVec[i]);
		if(RVec[i].PtVec.size()==0)
			continue;
		
		float r1,r2;
		int s1,s2,s3; s1 = RVec[i].ContourPtVec.size();s2 = RVec[i].InflectionPts.size();s3 =RVec[i].PtVec.size();
        r1 = float( RVec[i].ContourPtVec.size())  / float(RVec[i].PtVec.size()        );
		//r2 = float( RVec[i].InflectionPts.size()) / float(RVec[i].ContourPtVec.size() );
		r2 =  RVec[i].ConvexHullArea/float(RVec[i].PtVec.size() );
		Region R = RVec[i];
		if(RVec[i].PtVec.size()==0||RVec[i].ContourPtVec.size()==0)continue;
		fprintf(file,"%i,  %3i,%3i,%3i,%3i,   %i,%i,  %.2f,%.2f, %.2f,%.2f,%.2f,  %i,%.3f,%.3f,  %.3lf, %.2f,%.2f,%.2f,%.2f,%.2f,%i,%.2f,%.2f\n",//---%.2f,%.2f,%.2f,%.2f\n",////----%.2f,%.2f,%.2f,%.2f,%.2f,%.2f\n",
			i,
			RVec[i].left,RVec[i].top,RVec[i].right,RVec[i].bottom, //区域的四周位置
            RVec[i].rwidth,RVec[i].rheight,   //区域的宽和高
			RVec[i].x,RVec[i].y,              //区域的中心位置
			RVec[i].RR,RVec[i].GG,RVec[i].BB, //区域的RGB数值
			RVec[i].PtVec.size(),  r1,        //区域的点数
			r2,
			RVec[i].Feature.roundratio, 
			RVec[i].Feature.StrokeMean,RVec[i].Feature.StrokeStd,RVec[i].Feature.StkMedian,RVec[i].Feature.StkMin,RVec[i].Feature.StkMax,
			RVec[i].Feature.InnerHoles,
			RVec[i].Feature.HoleRatio,
			float(RVec[i].ContourPtVec.size())/float(RVec[i].rwidth + RVec[i].rheight));
		//RVec[i].Feature.AngleCnt[0],RVec[i].Feature.AngleCnt[1],RVec[i].Feature.AngleCnt[2],RVec[i].Feature.AngleCnt[3]);
		//RVec[i].Feature.OccupyVec[0],RVec[i].Feature.OccupyVec[1],RVec[i].Feature.OccupyVec[2],
		//RVec[i].Feature.OccupyVec[3],RVec[i].Feature.OccupyVec[4],RVec[i].Feature.OccupyVec[5]
		
		//compact<----边界/总数
	}
	
	fclose(file);
}

/// \brief     计算一个区域内的特征
///            计算区域的上下左右、长宽、圆心度等 
/// \param     Obj  待计算的一组区域
/// \return 无
void AnalysisRegionObj(Region&Obj)
{
	int i;
	float x,y,xmin,xmax,ymin,ymax;
	float xWeight,yWeight,gWeight;
    xWeight = yWeight = gWeight = 0;
	int val;
	xmin = ymin = 9999;
	xmax = ymax = 0;
	/*-----  integrate results */
	for (i=0; i<Obj.PtVec.size();i++)
	{  
		x  = Obj.PtVec[i].x; y = Obj.PtVec[i].y;
		val= 1;//Obj.PtVec[i].pixIntensity;
		
		if (xmin > x) xmin = x; if (xmax < x) xmax = x;
		if (ymin > y) ymin = y; if (ymax < y) ymax = y;
		
		xWeight += x * val; yWeight += y *val; gWeight  +=val;
	}   
	
	/* copy some data to "obj" structure */
	//Mxx=Mxx/rv;Myy=Myy/rv;Mxy=Mxy/rv;
	
	Obj.left  = xmin; Obj.right  = xmax;
	Obj.top   = ymin; Obj.bottom = ymax;
    Obj.GeoX =  ( xmin + xmax ) / 2;
	Obj.GeoY =  ( ymin + ymax ) / 2;
	
    Obj.rwidth  = xmax - xmin +1; 
	Obj.rheight = ymax - ymin +1;
	
	Obj.x   =  (xWeight / gWeight+1.0);
	Obj.y   =  (yWeight / gWeight+1.0);
    
	double Mxx,Myy,Mxy;
    Mxx = Myy = Mxy = 0;
	
	for (i=0; i<Obj.PtVec.size();i++)
	{
		x  = float(Obj.PtVec[i].x)-Obj.x;
		y  = float(Obj.PtVec[i].y)-Obj.y;
		val= 1;//Obj.PtVec[i].Flux;
		Mxx +=  (x * x * val); // / sum (I)
		Myy +=  (y * y * val); // / sum (I)
		Mxy +=  (x * y * val); // / sum (I) 
	}
	
	Obj.Feature.roundratio = sqrt(pow((Mxx - Myy), 2) + pow((2 * Mxy) , 2)) / (Mxx + Myy); 
	
}

//用来记录区域信息的临时结构体
struct TmpRgInfo_ 
{
	int left, top, right, bottom, ww, hh, Area;
};



/// \brief    计算一个区域的内部空洞数
/// R.Feature.InnerHoles  区域内部的空洞数   
/// R.Feature.HoleRatio   区域空洞面积和整个区域的面积比例         
/// \param     R          待计算的区域   
/// \return 无
void GetRegionInnerHolesAndArea( Region &R )
{
  //vector<int> WidthVec;
	vector< TmpRgInfo_ > InsideObjs;
    int i;
	int SumArea,MaxArea;
	SumArea = MaxArea = 0;
	R.Feature.InnerHoles = 0;
    for (i=0;i<R.MContours.size();i++)
    {
	   TmpRgInfo_ tmp_subobj;
       GetCurveBound( R.MContours[i], 
		   tmp_subobj.left, tmp_subobj.top, tmp_subobj.right, tmp_subobj.bottom, 
		   tmp_subobj.ww  , tmp_subobj.hh , tmp_subobj.Area );

	   if(tmp_subobj.Area<50) continue;

	   if(tmp_subobj.Area > MaxArea)
		   MaxArea = tmp_subobj.Area;
	   SumArea += tmp_subobj.Area; 
	   R.Feature.InnerHoles++;
    }

	if(R.Feature.InnerHoles>=1)
	{
		R.Feature.InnerHoles-=1;
		R.Feature.HoleRatio = double(SumArea - MaxArea)/ double(SumArea);
	}
}

/// \brief     计算一组区域的内部颜色
///   得到区域的内部区域颜色    
/// \param   RVec      待计算的区域 
/// \param   CPic      分析的原始图像 
/// \param   LabelVec  图像的联通域标记信息       
/// \return  
void GetObjContourColor(vector<Region> & RVec,C24BitMap&CPic,vector<int>&LabelVec)
{
	//void CMyConvexHull::CalcuConvexhull(vector<RPoint>&ContourPtVec, vector<RPoint>&ConvexhullPt, int ptstep=3)
	CMyConvexHull Cvx;
	TraceRegion(CPic.Width,CPic.Height,RVec,LabelVec);
	DeleteNullRegion(RVec);
    GetObjColor(CPic,RVec);
 	int i;

 	Loopi(RVec.size())
	{
     AnalysisRegionObj(RVec[i]);
	 GetRegionInnerHolesAndArea( RVec[i]);
	 RVec[i].ConvexHullArea = Cvx.CalcuConvexhull(RVec[i].ContourPtVec, RVec[i].ConvexHullPtVec);
	}
	//MarkInvalidRegionByCrossValidate(RVec);
}



#endif

