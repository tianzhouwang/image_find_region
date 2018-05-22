 
#ifndef CalcuStk_Head
#define CalcuStk_Head
#include "mregion.h"
#include "TraceBound.h"
#include "convexhull.h"
//#include "StoreRegion.h"
//#include "StoreRegion.h"


/// \brief  ɾ���ڲ�Ϊ�ջ��߽߱������������
/// \param  RVec MSER�������� 
/// \return �� �������Ժ���ԭ������ɾ��
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


/// \brief  �õ�Region���������������ÿ�������ڵ�ռ��
/// \param  R �����������
/// \return �� ��������Զ���������Ĳ���
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

/// \brief     �����������������������ļ���
/// \param     RVec  �������һ������
/// \filename  ��������ļ���
/// \return ��  
void  PrintfRegionVec(vector<Region>&RVec,char*filename)
{
	int i;
	FILE*file;
	file = fopen(filename, "wt+");
	
	fprintf(file,"Idx,��,��,��,��,��,��, centX, centY, ��,��,��, ���,���������,���͹����,Բ�Ķ�,�ʻ�ƽ��,�ʻ�����,�ʻ���ֵ,�ʻ���С,�ʻ����,�ն���,�ն���,�ܳ������\n");
	
	
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
			RVec[i].left,RVec[i].top,RVec[i].right,RVec[i].bottom, //���������λ��
            RVec[i].rwidth,RVec[i].rheight,   //����Ŀ�͸�
			RVec[i].x,RVec[i].y,              //���������λ��
			RVec[i].RR,RVec[i].GG,RVec[i].BB, //�����RGB��ֵ
			RVec[i].PtVec.size(),  r1,        //����ĵ���
			r2,
			RVec[i].Feature.roundratio, 
			RVec[i].Feature.StrokeMean,RVec[i].Feature.StrokeStd,RVec[i].Feature.StkMedian,RVec[i].Feature.StkMin,RVec[i].Feature.StkMax,
			RVec[i].Feature.InnerHoles,
			RVec[i].Feature.HoleRatio,
			float(RVec[i].ContourPtVec.size())/float(RVec[i].rwidth + RVec[i].rheight));
		//RVec[i].Feature.AngleCnt[0],RVec[i].Feature.AngleCnt[1],RVec[i].Feature.AngleCnt[2],RVec[i].Feature.AngleCnt[3]);
		//RVec[i].Feature.OccupyVec[0],RVec[i].Feature.OccupyVec[1],RVec[i].Feature.OccupyVec[2],
		//RVec[i].Feature.OccupyVec[3],RVec[i].Feature.OccupyVec[4],RVec[i].Feature.OccupyVec[5]
		
		//compact<----�߽�/����
	}
	
	fclose(file);
}

/// \brief     ����һ�������ڵ�����
///            ����������������ҡ�����Բ�Ķȵ� 
/// \param     Obj  �������һ������
/// \return ��
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

//������¼������Ϣ����ʱ�ṹ��
struct TmpRgInfo_ 
{
	int left, top, right, bottom, ww, hh, Area;
};



/// \brief    ����һ��������ڲ��ն���
/// R.Feature.InnerHoles  �����ڲ��Ŀն���   
/// R.Feature.HoleRatio   ����ն����������������������         
/// \param     R          �����������   
/// \return ��
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

/// \brief     ����һ��������ڲ���ɫ
///   �õ�������ڲ�������ɫ    
/// \param   RVec      ����������� 
/// \param   CPic      ������ԭʼͼ�� 
/// \param   LabelVec  ͼ�����ͨ������Ϣ       
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

