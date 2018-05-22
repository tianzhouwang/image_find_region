#ifndef Test_Func_Head
#define Test_Func_Head

#define DISPLAY_INNER_OBJECT 1




void TransHsvImg(C24BitMap&CPic,C256BitMap&GPic)
{
 int i,j;
 GPic.FormatF(CPic.Width,CPic.Height);
 Loopj(CPic.Height)
   Loopi(CPic.Width)
 {
   C24PixVal Pix =  get_pix_color(CPic,i,j);
//    float  H,  S,  V;
//    RGBToHSV((*Pix.r),(*Pix.g), (*Pix.b),  H,  S,  V);
//    S = S*128.0;

   //*get_pix_color(GPic,i,j) = 255 - BOUND( V/2+S, 0, 255);
   double Dis = (255 - (*Pix.r))* (255 - (*Pix.r)) + 
	            (*Pix.g) * (*Pix.g)                +
			    (*Pix.b) * (*Pix.b);
   Dis /= 2.5;

   *get_pix_color(GPic,i,j) = BOUND(sqrt(Dis),0,255);

 }

}


/// \brief  在图像文件显示输出一组长宽比像文字的区域 
/// \param  ObjVec    一组MSEG区域 
/// \param  width     图像宽
/// \param  height    图像高
/// \param  filename  输出的图像文件名
/// \param  selObjCr  用随机颜色还是区域自身的颜色绘制区域
/// \return 无

void DispSquareReionObjvec(vector<Region>&ObjVec,int width,int height,char*filename,int selObjCr=0)
{	
	C24BitMap KBPic;
	KBPic.FormatF(width, height);
	KBPic.ClearPic();
	int i,j,t;

	Loopj(ObjVec.size())
	{
	 
        char buff[30];
		sprintf(buff,"%i",j);

		KBPic.RandPenColor();
        if(selObjCr)
		   KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
         
        KBPic.Rectangle(ObjVec[j].left,ObjVec[j].top,ObjVec[j].right,ObjVec[j].bottom);

		if(DISPLAY_INNER_OBJECT)
		Loopi(ObjVec[j].PtVec.size())
		{
			KBPic.SigDot(ObjVec[j].PtVec[i].x,ObjVec[j].PtVec[i].y);
		}

        KBPic.SetPenColor(255, 0, 0);
		if(ObjVec[j].rwidth==0) continue;

        double char_ratio = float(ObjVec[j].rheight) / float(ObjVec[j].rwidth);
		Region R = ObjVec[j];

		KBPic.SetPenColor(255, 255, 0);
		MergeTxtStr(KBPic, ObjVec[j].GeoX,ObjVec[j].GeoY,15, buff,0,255,255); 

	    if(char_ratio >0.85 && char_ratio<1.2)
		{
          KBPic.Rectangle(ObjVec[j].left,ObjVec[j].top,ObjVec[j].right,ObjVec[j].bottom);
		  KBPic.DrawCircle(ObjVec[j].GeoX, ObjVec[j].GeoY, 3);
          
		  KBPic.SetPenColor(0, 255, 0);
		  
		}
        
		
	}

	KBPic.Save(filename);
}

#endif