#ifndef STW_TYPEs_H
#define STW_TYPEs_H
#include <vector>
#include <math.h>



using namespace std;

struct CtRect
{
 //int x,y,width,height;
 int left,top,right,bottom;
 CtRect(int left_,int top_,int right_,int bottom_)
	{ left  = left_;
      top   = top_;
	  right = right_;
	  bottom= bottom_;}
 int RR,GG,BB;
 int ChildNum;
 CtRect(){};
 vector<CtRect> Childs;
};

struct fSPOINT
{//public:
	//double x,y;
	int ix,iy;
	fSPOINT(){ix=iy=0; };//x=y=0;};
	fSPOINT(int x_,int y_){ix = x_; iy=y_; };//x=x_;y=y_;};
	double Flux;
	double RR,GG,BB;
	float StrokeWidth;
};
 

struct MatType 
{
 double StrokeWidth;
 double RR,GG,BB;
 double Value;
 int CCMapLabel;
};

class Mat
{
public:
	vector<MatType>  data;
	vector<MatType*> Pt;
	Mat(){};//  {data=NULL;Pt=NULL;}
	~Mat() { ;};
	int Width,Height;
	int size(){return Width*Height;}
	MatType At(int x,int y){return Pt[x][y];}
	MatType At(fSPOINT& inp){return Pt[inp.ix][inp.iy];}
	void Set(fSPOINT&inp,MatType V) {Pt[inp.ix][inp.iy] = V;}
	void FormatF(int w,int h);
	void SetImg(C256BitMap&OriPic);
	void SaveFile(char*file);
	void SaveImg(C256BitMap&SvPic);
	
	void FillStw(double Value);
	void FillLbVal(int Value);
	void FillVal(double Value){int i;for (i=0;i<Width*Height;i++)data[i].Value = Value;};
	
	int  GetCCLabel(int i,int j);
	int  GetSkWidth(int i,int j);
	void SetCCLabel(int i,int j,int Val);
	void CopyMat(Mat&StwMap);
	
};

void Mat::CopyMat(Mat&StwMap)
{
this->FormatF(StwMap.Width,StwMap.Height);
int i;
for(i=0;i<Width*Height;i++) data[i] = StwMap.data[i]; 
}

int Mat::GetCCLabel(int i,int j)
{
	return Pt[i][j].CCMapLabel;
}

int Mat::GetSkWidth(int i,int j)
{
   return Pt[i][j].StrokeWidth;
}

void Mat::SetCCLabel(int i,int j,int LbVal)
{
    Pt[i][j].CCMapLabel = LbVal;
}

void Mat::SaveImg(C256BitMap&SvPic)
{
  SvPic.FormatF(Width,Height);
  int i,j;
  for(i=0;i< SvPic.Width;i++)
	 for(j=0;j< SvPic.Height;j++)
	 {      
	  *get_pix_color( SvPic,i,j) = BOUND(Pt[i][j].StrokeWidth*10,0,255);
	 }		
  //Pic.Save(filename);
}

void Mat::FillStw(double Value)
{  int i; for (i=0;i<Width*Height;i++)
 	  data[i].StrokeWidth = Value;
}

void Mat::FillLbVal(int Value)
{ int i;  for (i=0;i<Width*Height;i++)
 	  data[i].CCMapLabel = Value;
}

void Mat::SaveFile(char*filename)
{ C256BitMap Pic;
  Pic.FormatF(Width,Height);
  int i,j;
  for(i=0;i< Pic.Width;i++)
	  for(j=0;j< Pic.Height;j++)
	  {      
		if(Pt[i][j].StrokeWidth!=0)
		{
			int kk = Pt[i][j].StrokeWidth;
		}
	    *get_pix_color( Pic,i,j) = BOUND(Pt[i][j].StrokeWidth*20,0,255);
      }
  Pic.Save(filename);
}

void Mat::SetImg(C256BitMap&OriPic)
{
  FormatF(OriPic.Width,OriPic.Height);
  int i,j;
  for(i=0;i<OriPic.Width;i++)
	  for(j=0;j<OriPic.Height;j++)
	  {      
	     Pt[i][j].Value = *get_pix_color(OriPic,i,j);
	  }	
}

void Mat::FormatF(int w,int h)
{    
  int i;
  Width  = w; Height = h;
  data.resize(w*h); Pt.resize(w);
  for(i=0;i<w;i++) Pt[i] = &data[i*h];
  for(i=0;i<Width*Height;i++) data[i].Value = data[i].StrokeWidth = -1;
}



/*struct ObjInfo {
	double	x,y;	       // 重心 X, Y轴坐标
	double GeoX,GeoY;      // 几何中心   
	double	adu;	       // 象素ADU值
	double  RR,GG,BB;	   // 象素ADU值
	double StdRR,StdGG,StdBB;
	int Width,Height,Radius;      // 宽度、高度
	int pixnum;            // 包含的点数
	double Flux;           // 流量总和
	bool	isMark;  	   // 是否已经标记
	double  max_adu_x,max_adu_y;//最大adu的坐标
	double  peak;          //最大流量
	vector<fSPOINT> PtVec; //所包含的点
	int left,right,top,bottom;//外接矩形
	double  Stdev;         //区域内方差
	double  Average;       //区域内平均比划
	//double ratio;          //圆心度
	int ClassT;            //label标号
	int FindFlag;          //查找标记2
	int IfValid;           //是否有效区域 -1为无效区域
	double StrokeWidth;    //比划的宽度
    double StwStdev;       //比划的方差
	fSPOINT FourCorner[4]; //四个角点
	double AngWidth,AngHeight;//斜角宽，斜角高
	double cxx,cyy,cxy;
	double roundratio;
	int ObjIdx;
};*/


//<===============================================================
/*void  analysis(ObjInfo&Obj)
{
	//
	double    peak,val,powVal;
	float     i,x, y, xmin,xmax, ymin,ymax;
	double xWeight, yWeight, gWeight,MaxAdu,Mxx,Myy,Mxy ; 
	double RRVal,GGVal,BBVal;
	// X, Y, 和ADU的权重
	double v_flux;
	Mxx=0;Myy=0;Mxy=0; RRVal = GGVal = BBVal =0;
 
	//   ~~~~~~ 此处省略  ~~~~
	xWeight = yWeight = gWeight = 0;
	xmin=xmax=Obj.PtVec[0].ix;
    ymin=ymax=Obj.PtVec[0].iy;
	powVal = v_flux =0;
	peak = 0;
 
	for (i=0; i<Obj.PtVec.size();i++)
    {
		x  = Obj.PtVec[i].ix; y  = Obj.PtVec[i].iy;
		val= 1;//Obj.PtVec[i].Flux;
		v_flux+=val;
		powVal += val * val;
		
		RRVal += Obj.PtVec[i].RR;
		GGVal += Obj.PtVec[i].GG;
		BBVal += Obj.PtVec[i].BB;

		if (peak < val) peak = val;
		if (xmin > x)  xmin = x;
		if (xmax < x)  xmax = x;
		if (ymin > y)  ymin = y;
		if (ymax < y)  ymax = y;
		
		//计算物体几何中心/重心
		val = 1.0;//此处注释掉是重心
		xWeight += x * val;
		yWeight += y * val;
		gWeight += val;   
    }  
    
	
	Obj.peak = peak;
	Obj.left = xmin;  Obj.right = xmax;
	Obj.top = ymin;   Obj.bottom = ymax;
	
	Obj.GeoX = (Obj.left + Obj.right )/2;
	Obj.GeoY = (Obj.top  + Obj.bottom)/2;

	Obj.x =  (xWeight/gWeight+0.5);
	Obj.y =  (yWeight/gWeight+0.5);

	//Obj.x = Obj.GeoX; Obj.y = Obj.GeoY;//Obj.bottom; <===调试
	//Obj.x = Obj.left; Obj.y =  Obj.bottom;  

	Obj.Flux    = v_flux;
	Obj.pixnum  = Obj.PtVec.size();
	Obj.RR = RRVal/double(Obj.pixnum); 
	Obj.GG = GGVal/double(Obj.pixnum); 
	Obj.BB = BBVal/double(Obj.pixnum); 
	Obj.Average = v_flux/double(Obj.pixnum); 
	Obj.Stdev    = sqrt(powVal / double(Obj.pixnum) - Obj.Average * Obj.Average);
	Obj.StwStdev = Obj.Stdev;// 笔划的变异率
	Obj.Width  = Obj.right  - Obj.left +1;
	Obj.Height = Obj.bottom - Obj.top  +1;
    Obj.Radius = sqrt(double(Obj.Width * Obj.Width + Obj.Height * Obj.Height));
	//以下是圆心度数计算，检查物体是否是长条形
	for (i=0; i<Obj.PtVec.size();i++)
	{
		x  = float(Obj.PtVec[i].ix)-Obj.x;
		y  = float(Obj.PtVec[i].iy)-Obj.y;
		val= Obj.PtVec[i].Flux;
		
		Mxx +=  (x * x * val);// / sum (I)
		Myy +=  (y * y * val);// / sum (I)
		Mxy +=  (x * y * val);// / sum (I) 
	}
	
	Obj.roundratio =sqrt(pow((Mxx - Myy), 2) + pow((2 * Mxy) , 2)) / (Mxx + Myy); 
	//以下是物体的椭圆拟合参数
	double  xm2,ym2, xym, temp;
	
	xm2 = Mxx; ym2 = Myy; xym = Mxy;
	temp = sqrt((xm2-ym2)*(xm2-ym2)/4 - xym*xym);
	
	Obj.cxx = (float)(ym2/temp);
	Obj.cyy = (float)(xm2/temp);
	Obj.cxy = (float)(-2*xym/temp);
	
} 
*/

#endif
