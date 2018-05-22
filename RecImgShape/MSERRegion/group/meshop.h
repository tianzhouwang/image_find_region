#ifndef Mesh_Operation_H
#define Mesh_Operation_H
//#include "mregion.h"
//#include "RegionDisp.h"


void GetRegionOutterCurve(Region&R)
{
	if(R.MContours.size()==1)
		R.OutterCurve = R.MContours[0];

    int i;
	int MaxArea =-1;
    int MaxIdx  =-1;
	int left, top, right, bottom, ww, hh, Area;
    for(i=0;i<R.MContours.size();i++)
	{
      GetCurveBound(R.MContours[i],left, top, right, bottom, ww, hh, Area);
	  if(Area>MaxArea)
	  {
		  MaxArea = Area;
		  MaxIdx  = i;
	  }
	}

	if(MaxIdx!=-1)
       R.OutterCurve = R.MContours[MaxIdx];
    else
       R.OutterCurve = R.MContours[0];
}

//////////////////////////////////////////////////////////////////////////////////////
/*!
\brief		computes the area of a 3D planar polygon.
\note		Reference: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#3D%20Polygons.
			Copyright 2000, softSurfer (www.softsurfer.com)
			This code may be freely used and modified for any purpose
			providing that this copyright notice is included with it.
			SoftSurfer makes no warranty for this code, and cannot be held
			liable for any real or imagined damage resulting from its use.
			Users of this code must verify correctness for their application.
\param[in]	n		The number of vertices in the polygon.
\param[out]	V		An array of n vertices in a plane.
\param[out]	N		Unit normal vector of the polygon's plane.
\return		The area of the polygon.
*/
//////////////////////////////////////////////////////////////////////////////////////

struct  Tb3dPoint
{
  float x,y,z;
};

double Cal3DPolygonSignedArea( Tb3dPoint* V, int n,  Tb3dPoint N )
{
	double area = 0;
	double an, ax, ay, az;  // abs value of normal and its coords
	int   coord;           // coord to ignore: 1=x, 2=y, 3=z
	int   i, j, k;         // loop indices

	// select largest abs coordinate to ignore for projection
	// Start modification by Yefeng Zheng 06-14-2007.
	/* Original
	ax = (N.x>0 ? N.x : -N.x);     // abs x-coord
	ay = (N.y>0 ? N.y : -N.y);     // abs y-coord
	az = (N.z>0 ? N.z : -N.z);     // abs z-coord

	coord = 3;                     // ignore z-coord
	if (ax > ay) 
	{
	if (ax > az) coord = 1;    // ignore x-coord
	}
	else if (ay > az) coord = 2;   // ignore y-coord
	*/
	ax = N.x;     // x-coord
	ay = N.y;     // y-coord
	az = N.z;     // z-coord

	coord = 3;                     // ignore z-coord
	if (fabs(ax) > fabs(ay)) 
	{
		if (fabs(ax) > fabs(az)) coord = 1;    // ignore x-coord
	}
	else if (fabs(ay) > fabs(az)) coord = 2;   // ignore y-coord
	// End modification by Yefeng Zheng 06-14-2007.

	// compute area of the 2D projection
	for( k=0; k<n; k++ )
	{
		i = (k+1)%n;
		j = (k+2)%n;
		switch (coord) 
		{
		case 1:
			area += (V[i].y * (V[j].z - V[k].z));
			continue;
		case 2:
			/* Start modification by Yefeng Zheng 06-14-2007.
			/* Original
			area += (V[i].x * (V[j].z - V[k].z));
			*/
			area += (V[i].z * (V[j].x - V[k].x));
			// End modificatin by Yefeng Zheng 06-14-2007.
			continue;
		case 3:
			area += (V[i].x * (V[j].y - V[k].y));
			continue;
		}
	}
	// scale to get area before projection
	an = sqrt( ax*ax + ay*ay + az*az);  // length of normal vector
	switch (coord) {
	case 1:
		area *= (an / (2*ax));
		break;
	case 2:
		area *= (an / (2*ay));
		break;
	case 3:
		area *= (an / (2*az));
	}
	return area;
}


bool IsCounterClockWiseOrder(vector<RPoint>&ContourPts)
{ 
	
	vector<Tb3dPoint>  Pnt;
	
	if( ContourPts.size() == 0 )
		return false;
    int i;
    Pnt.resize(ContourPts.size());

    Loopi(Pnt.size()) 
	{
		Pnt[i].x = ContourPts[i].x;
		Pnt[i].y = ContourPts[i].y;
		Pnt[i].z = 0;
	}	

	Tb3dPoint N;// = tIds3dPoint( Normal[0], Normal[1], Normal[2] );
	N.x = 0;
	N.y = 0;
	N.z = 1;
	return Cal3DPolygonSignedArea( &Pnt[0], (int)Pnt.size(), N ) >= 0;
}

// vector<RPoint> PtVec;              //区域内的点
// vector<RPoint> ContourPtVec;       //区域的所有轮廓点
// vector<RPoint> ConvexHullPtVec;    //区域内的点

//////////////////////////////////////////////////////////////////////////////////////
/*!
\brief		Uniformly resample a curve.
\param[in]	input		Input curve.
\param[out]	output		Curve after resampling.
\param[in]	nPoint		Number of sampling points.
\return		0 for success and others for failure.
*/
//////////////////////////////////////////////////////////////////////////////////////

double diff_two_norm(RPoint&Pt1, RPoint&Pt2) 
{ 
	double  sum = 0;
	double Dx = Pt1.x - Pt2.x;
	double Dy = Pt1.y - Pt2.y;

	 sum  = Dx*Dx + Dy*Dy;
	return  sqrt( sum );
  }


int ResampleCurve(vector<RPoint>& input, vector<RPoint>& output, float nStepSize)
{
	int nPoint;


	if( input.size() == 0 )
		return -1;
	
	int i;
	double full_length  = 0;
	int nInPoints = int(input.size());
	for( i=1; i<nInPoints; i++ )
		full_length += diff_two_norm(input[i-1],input[i]); 
	
    nPoint = full_length/nStepSize;
    if(nPoint==0)
		return 0;

	double crtdelta = full_length / (nPoint-1);
	double cumdelta = 0;
	double cumlength = 0;
	
	output.resize(nPoint);
	int crt_idx = 0;
	output[crt_idx] = input[0];
	++crt_idx;
	cumdelta += crtdelta;
	for( i=1; i<nInPoints; i++ )
	{
		double crtlength = diff_two_norm(input[i-1],input[i]);

		crtlength = max( 1e-5, crtlength );
		cumlength += crtlength;
		while ((cumlength>=cumdelta)&& (crt_idx<(nPoint-1)))
		{
			//add intersection
			double alpha = (cumlength-cumdelta)/crtlength;
			output[crt_idx].x = alpha*input[i-1].x+(1-alpha)*input[i].x;
			output[crt_idx].y = alpha*input[i-1].y+(1-alpha)*input[i].y;
			++crt_idx;
			cumdelta += crtdelta;
		}
	}
	//last point
	output[crt_idx] = input[nInPoints-1];
	return 0;
}

int	ReversePointOrder( vector<RPoint>&Contour )
{
	int nPnt =  Contour.size();

	for(int i=0; i<nPnt/2; i++ )
	{
		RPoint temp = Contour[i];
		Contour[i] = Contour[nPnt-1-i];
		Contour[nPnt-1-i] = temp;
	}
	return 0;
}

 
int EnforceAntiCounterClockWiseOrder(vector<RPoint>&Contour)
{
	if( Contour.size() == 0 )
		return -1;
	
	bool bCounterClockWise = IsCounterClockWiseOrder( Contour);
	if( bCounterClockWise ) // Reverse the order.
		ReversePointOrder( Contour);
	
	return 0;
}

/*int CrossProduct( Tb3dPoint&U, Tb3dPoint&V, Tb3dPoint&W )
{
	W.x = U.y*V.z - U.z*V.y;
	W.y = U.z*V.x - U.x*V.z;
	W.z = U.x*U.y - U.y*V.x;
	return 0;
}*/


Tb3dPoint CrossProduct(Tb3dPoint &u, Tb3dPoint &v)
{
	// cross product u x v
	Tb3dPoint r;
	r.x = u.y*v.z - u.z*v.y;
	r.y = u.z*v.x - u.x*v.z;
	r.z = u.x*v.y - u.y*v.x;
	return r;
}

Tb3dPoint operator-=(Tb3dPoint&m_data,const Tb3dPoint& rhs)
{
	m_data.x -= rhs.x;
	m_data.y -= rhs.y;
	m_data.z -= rhs.z;
	
	return m_data;
}

Tb3dPoint operator-(Tb3dPoint& lhs, const Tb3dPoint& rhs)
{   Tb3dPoint res = lhs;
res -= rhs;
return res;
}
int Normalize( Tb3dPoint&U)
{
	double norm =  sqrt(U.x*U.x+U.y*U.y+U.z*U.z);
	if( norm != 0 )
	{
		U.x /= norm;
		U.y /= norm;
		U.z /= norm;
	}
	return 0;
}


int PointInPolygon(vector<RPoint> PolygonPts, float testx, float testy)
{
	//int nvert, float *vertx, float *verty, 
	int i, j, c = 0;
	int nvert = PolygonPts.size();

	for (i = 0, j = nvert-1; i < nvert; j = i++) {
		if ( ((PolygonPts[i].y>testy) != (PolygonPts[j].y>testy)) &&
			(testx < (PolygonPts[j].x-PolygonPts[i].x) * (testy-PolygonPts[i].y) / (PolygonPts[j].y-PolygonPts[i].y) + PolygonPts[i].x) )
			c = !c;
	}

	return c;
}


int CalculateContourNormal( vector<Tb3dPoint> &MeshPoint, vector<Tb3dPoint> &Norm,vector<Tb3dPoint>&SegNorm )
{
	int i, j;
	
	int nMeshPoint = (int)MeshPoint.size();
	if( nMeshPoint == 0 )
		return -1;
	
	// Calculate the plane normal.
	Tb3dPoint PlaneNorm;
	PlaneNorm.x = 0;
	PlaneNorm.y = 0;
	PlaneNorm.z = 1.0;
	//CalPlaneNormal( MeshPoint, PlaneNorm, true );
	
	// Calculate normal for each segment.
	//vector<Tb3dPoint> SegNorm;
	SegNorm.resize( nMeshPoint );
	for( i=0; i<nMeshPoint; i++ )
	{
		j = (i+1)%nMeshPoint;
		Tb3dPoint u    = MeshPoint[i] - MeshPoint[j];
		Tb3dPoint Norm = CrossProduct(u , PlaneNorm );
		SegNorm[i] = Norm;
		Normalize(SegNorm[i]);//.normalize();
	}
	
	// Calculate the normal for each point.
	Norm.resize( nMeshPoint );
	for( i=0; i<nMeshPoint; i++ )
	{
		j = (i+nMeshPoint-1)%nMeshPoint;

		Norm[i].x =  float(SegNorm[i].x +  SegNorm[j].x)/2.0;
		Norm[i].y =  float(SegNorm[i].y +  SegNorm[j].y)/2.0;
		Norm[i].z =  float(SegNorm[i].z +  SegNorm[j].z)/2.0;

		Normalize(Norm[i]);
	}
	
	return 0;
}


//void GetCurveBound(vector<RPoint>&Curve,int&left,int&top,int&right,int&bottom,int&ww,int&hh,int&Area)
void GetRegionExpContour(Region&R, vector<RPoint>&ContourExp,float explen)
{
  int i,j;
  GetRegionOutterCurve( R);
  vector<RPoint> ReSample;
  
  //int region_scale = max(R.rwidth,R.rheight);
  float stepSize = 5.0;
  ResampleCurve(R.OutterCurve, ReSample, stepSize);
  EnforceAntiCounterClockWiseOrder(ReSample);


  vector<Tb3dPoint>  MeshPoints, ContourNorm, SegmentNorm; 
  MeshPoints.resize(ReSample.size());

  Loopi(ReSample.size())
  {
	  MeshPoints[i].x =  ReSample[i].x; MeshPoints[i].y =  ReSample[i].y;
      MeshPoints[i].z =  0;
  }
  
  ContourNorm.resize(MeshPoints.size());
  ContourExp.resize(MeshPoints.size());
  CalculateContourNormal(MeshPoints, ContourNorm,SegmentNorm);
  
  Loopi(ReSample.size())
  {
	  ContourExp[i].x  = ReSample[i].x + ContourNorm[i].x * explen;
	  ContourExp[i].y  = ReSample[i].y + ContourNorm[i].y * explen;
  }
}

void DrawContour(C24BitMap&CPic,vector<RPoint>&Contour)
{
  int i;
  vector<RPoint> ReSample;
  ResampleCurve(Contour, ReSample, 100);
  
  Contour = ReSample;

  EnforceAntiCounterClockWiseOrder(Contour);

  vector<Tb3dPoint>  MeshPoints,ContourNorm,SegmentNorm; 


  MeshPoints.resize(Contour.size());
  Loopi(Contour.size())
  {
	  MeshPoints[i].x =  Contour[i].x;
	  MeshPoints[i].y =  Contour[i].y;
      MeshPoints[i].z =  0;
  }

  ContourNorm.resize(MeshPoints.size());

  CalculateContourNormal(MeshPoints, ContourNorm,SegmentNorm);
  //ContourNorm = SegmentNorm;

  Loopi(Contour.size())
  {
	 CPic.DrawCircle(Contour[i].x,Contour[i].y,1);
  }

  Loopi(Contour.size()-1)
  {
	  //CPic.DrawCircle(Contour[i].x,Contour[i].y,1);
	  CPic.DrawLine(Contour[i  ].x,Contour[i  ].y,
		            Contour[i+1].x,Contour[i+1].y);
	  CPic.PenColor.R = CPic.PenColor.R * 9/10;
  }



  Loopi(Contour.size())
  {
	  Normalize( ContourNorm[i]);
	  Contour[i].x += ContourNorm[i].x * 15.0;
	  Contour[i].y += ContourNorm[i].y * 15.0;

	   
	  CPic.DrawCircle(Contour[i].x,Contour[i].y,1);
  }
  
  Loopi(Contour.size()-1)
  {
	  //CPic.DrawCircle(Contour[i].x,Contour[i].y,1);
	  CPic.DrawLine(Contour[i  ].x,Contour[i  ].y,
		  Contour[i+1].x,Contour[i+1].y);
	  CPic.PenColor.R = CPic.PenColor.R * 9/10;
  }

  CPic.Save("result\\temp.bmp");
} 

 
void DispReionObjvecExp(vector<Region>&ObjVec, int width, int height, char*filename)
{	
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
	int i,j,t;

	Loopj(ObjVec.size())
	{
	   KBPic.RandPenColor();
 
       if(ObjVec[j].MContours.size()>20) 
		   continue;

        char buff[30];
		sprintf(buff,"%i",j);

		
        //if(DISPLAY_DBG_REGION_COLOR)
		 //  KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
         
		vector<RPoint> ContourExp;
		ContourExp.clear();
	    GetRegionExpContour(ObjVec[j],  ContourExp,5);	 
        
		Loopi(ContourExp.size())
		{
            KBPic.DrawCircle(ContourExp[i].x,ContourExp[i].y,1);
		}

	
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
		   
	}

	KBPic.Save(filename);
}

#endif