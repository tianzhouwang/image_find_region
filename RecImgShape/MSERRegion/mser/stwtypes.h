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
	double	x,y;	       // ���� X, Y������
	double GeoX,GeoY;      // ��������   
	double	adu;	       // ����ADUֵ
	double  RR,GG,BB;	   // ����ADUֵ
	double StdRR,StdGG,StdBB;
	int Width,Height,Radius;      // ��ȡ��߶�
	int pixnum;            // �����ĵ���
	double Flux;           // �����ܺ�
	bool	isMark;  	   // �Ƿ��Ѿ����
	double  max_adu_x,max_adu_y;//���adu������
	double  peak;          //�������
	vector<fSPOINT> PtVec; //�������ĵ�
	int left,right,top,bottom;//��Ӿ���
	double  Stdev;         //�����ڷ���
	double  Average;       //������ƽ���Ȼ�
	//double ratio;          //Բ�Ķ�
	int ClassT;            //label���
	int FindFlag;          //���ұ��2
	int IfValid;           //�Ƿ���Ч���� -1Ϊ��Ч����
	double StrokeWidth;    //�Ȼ��Ŀ��
    double StwStdev;       //�Ȼ��ķ���
	fSPOINT FourCorner[4]; //�ĸ��ǵ�
	double AngWidth,AngHeight;//б�ǿ�б�Ǹ�
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
	// X, Y, ��ADU��Ȩ��
	double v_flux;
	Mxx=0;Myy=0;Mxy=0; RRVal = GGVal = BBVal =0;
 
	//   ~~~~~~ �˴�ʡ��  ~~~~
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
		
		//�������弸������/����
		val = 1.0;//�˴�ע�͵�������
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

	//Obj.x = Obj.GeoX; Obj.y = Obj.GeoY;//Obj.bottom; <===����
	//Obj.x = Obj.left; Obj.y =  Obj.bottom;  

	Obj.Flux    = v_flux;
	Obj.pixnum  = Obj.PtVec.size();
	Obj.RR = RRVal/double(Obj.pixnum); 
	Obj.GG = GGVal/double(Obj.pixnum); 
	Obj.BB = BBVal/double(Obj.pixnum); 
	Obj.Average = v_flux/double(Obj.pixnum); 
	Obj.Stdev    = sqrt(powVal / double(Obj.pixnum) - Obj.Average * Obj.Average);
	Obj.StwStdev = Obj.Stdev;// �ʻ��ı�����
	Obj.Width  = Obj.right  - Obj.left +1;
	Obj.Height = Obj.bottom - Obj.top  +1;
    Obj.Radius = sqrt(double(Obj.Width * Obj.Width + Obj.Height * Obj.Height));
	//������Բ�Ķ������㣬��������Ƿ��ǳ�����
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
	//�������������Բ��ϲ���
	double  xm2,ym2, xym, temp;
	
	xm2 = Mxx; ym2 = Myy; xym = Mxy;
	temp = sqrt((xm2-ym2)*(xm2-ym2)/4 - xym*xym);
	
	Obj.cxx = (float)(ym2/temp);
	Obj.cyy = (float)(xm2/temp);
	Obj.cxy = (float)(-2*xym/temp);
	
} 
*/

#endif
