#ifndef _CIMGGEOMETRY_scale_H
#define _CIMGGEOMETRY_scale_H
#include "c256bitmap.h"
#define BOUND(x,a,b) (((x) < (a)) ? (a) : (((x) > (b)) ? (b) : (x)))
typedef unsigned char BYTE;
///lpbySrcXY----����Դ����(x, y)�ĵ�ַ, 
//x, y, ��������任��õ��Ķ�Ӧ��ԭͼ�����������
//nScanWidth,  int nScanHeight, Դɨ���Ⱥ�ɨ��߶�
BYTE interpolate(BYTE* PtPos,int x,int y,
float fu,float fv,int nScanWidth,int nScanHeight)
{
//���ڵ��ĸ����������½ǵ��x, y����ƫ����
 int nx(1),ny(1);
 if((x + 1)>(nScanWidth  - 1)) nx = 0;
 if((y + 1)>(nScanHeight - 1)) ny = 0;
 BYTE abyDot[2][2];//���ص�(x, y)������λ��
 BYTE* pbySrc = PtPos;//��ȡ������ֵ(x, y) = (x, y) + (0, 0)						
 abyDot[0][0] = *pbySrc;
 pbySrc = PtPos + nx ;//(x + 1, y) = (x, y) + (1, 0)
 abyDot[1][0] = *pbySrc;
 pbySrc = (PtPos - ny * nScanWidth);//ָ����һ������
 abyDot[0][1] = *pbySrc;//(x + 1, y + 1) = (x, y) + (1, 1)
 pbySrc = pbySrc+ nx;
 abyDot[1][1] = *pbySrc;
 BYTE Val = (BYTE)(BOUND(((1 - fu) * (1 - fv) * ((float)abyDot[0][0]) + 
				      (1 - fu) * fv * ((float)abyDot[0][1]) + 
					  fu * (1 - fv) * ((float)abyDot[1][0]) +
					  fu * fv * ((float)abyDot[1][1])), 0, 255));
 return Val;
}


void Scale(BYTE* lpbyBitsSrc32,int x,int y,int Width,int Height,
int ScanWidth,int ScanHeight,BYTE* BitsDst,int WidthImgDst,int HeightImgDst)
{
int w(min(Width,ScanWidth- x)),h(min(Height,ScanHeight- y)),i,j;
//fScalex,fScaley����ʾ�����ű�Ϊ��ʵ���űȵĵ���֮������������������,��һ�γ���,�ܱ�һ�γ˷�Ҫ��.
//������ű�
float fScalex((float)w/(float)WidthImgDst),fScaley((float)h/(float)HeightImgDst);
Loopi(HeightImgDst)//��ɱ任
	{float fYInverse = i * fScaley;//����任���õĸ���yֵ
	 int yy = (int)fYInverse;//ȡ��
	 float fv = fYInverse - yy;//�����ֵ
	 yy+=y;//��Ӧ��ԭͼ���y����
     
	 BYTE* pbySrc = lpbyBitsSrc32+(ScanHeight-1-yy)*ScanWidth;
	 BYTE* pbyDst = BitsDst+ 4*(HeightImgDst-1-i)*int((WidthImgDst+3)/4);
     Loopj(WidthImgDst)
	 {float fXInverse=j*fScalex;//����任���õĸ���xֵ
	  int xx = (int)fXInverse;//ȡ��
	  float fu = fXInverse - xx;//�����ֵ
	  xx+=x;//��Ӧ��ԭͼ���x����
	  BYTE* BPtCurrent = pbySrc+xx;//��ȡ����
	  *pbyDst++ = interpolate(BPtCurrent,xx,yy,fu,fv,ScanWidth,ScanHeight);
	 }
	}
}


/// \brief  ���ŻҶ�ͼ��
/// \param  src �����ŵ�ԭͼ
/// \param  x ROI��x ����
/// \param  y ROI��y ����
/// \param  ���ź�ͼ��
/// \destw  ���ź�ͼ����
/// \desth  ���ź�ͼ��߶�
/// \return 
void scale_img(C256BitMap&src,int x,int y,
  int w,int h,C256BitMap&dest,int dstw,int dsth)
{
C256BitMap destP;
destP.FormatF(dstw,dsth);
Scale(src.Buffer,x,y,w,h,src.LineWidth,src.Height,
	destP.Buffer,dstw,dsth);
dest=destP;
}

/// \brief  �������ͼ��
/// \param  SrcPic �����ŵ�ԭͼ
/// \param  DsPic  ���ź��ͼ��
/// \w  ���ź�ͼ����
/// \h  ���ź�ͼ��߶�
/// \return 
void  Scale_BMP(C24BitMap&SrPic,C24BitMap&DsPic,int w,int h)
{
	IJT;DsPic.FormatF(w,h);
	float fScalex=float(SrPic.Width )/(float)w;
	float fScaley=float(SrPic.Height)/(float)h;
	float fYInverse,fXInverse,fu,fv;int xx,yy;
	Loopi(w)Loopj(h)
	{
		fYInverse=j*fScaley;yy=fYInverse;fv=fYInverse-yy;
		fXInverse=i*fScalex;xx=fXInverse;fu=fXInverse-xx;
		C24PixVal Tmp,TmpDs;
		Tmp=C24PtItrp(SrPic,xx,yy,fu,fv);
		TmpDs=get_pix_color(DsPic,i,j);
		*TmpDs.r=Tmp.R;*TmpDs.g=Tmp.G;*TmpDs.b=Tmp.B;
	}
}

#endif
