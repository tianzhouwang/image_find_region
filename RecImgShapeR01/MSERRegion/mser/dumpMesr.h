 
#ifndef Dump_Mser_Info_Head
#define Dump_Mser_Info_Head
//#include "hzmerge.h"

void  CompPictures_A(C256BitMap&GPic, Region*rootRegion,vector<Region>&ReginVec);
void  CompleteRegion_A(FillFindRegion& RgFill,C256BitMap&Pic,Region&Rg);

class DumpMserRg
{
public:
	DumpMserRg();
	void DumpPictures(char*dirNames,vector<Region>&ReginVec);
	void DumpPictures(char*dirNames,Region*RootRg);
	void PrintfRegion(char*filename,Region*TreeRoot);
	void PrintfRegionList(char*filename,Region*TreeRoot);
	void MakParentInfo(Region*TreeRoot);
	C24BitMap  OPic;
	C256BitMap GPic;
};

void DumpMserRg::MakParentInfo(Region*TreeRoot)
{
   
	for(Region * child = TreeRoot->child_; child; child = child->next_)
	{
		child->parent_ = TreeRoot;
		MakParentInfo(child);
		//child->PrintTree(file,start+1);
	}
}

void DumpMserRg::PrintfRegionList(char*filename,Region*TreeRoot)
{
//	vector<Region> ReginVec;
//    TreeRoot->SaveMserRegion(ReginVec, 0.8);
//	PrintfRegionVec(ReginVec,filename);
}

DumpMserRg::DumpMserRg()
{
  InitTextMerge();
}

 void DumpMserRg::PrintfRegion(char*filename,Region*TreeRoot)
 {
	 FILE*file = fopen(filename,"wt+");
	 TreeRoot->parent_ = TreeRoot;
	 MakParentInfo(TreeRoot);
	 TreeRoot->PrintTree(file,0);
	 fclose(file);
 }

void DumpMserRg::DumpPictures(char*dirNames,Region*RootRg)
{
   
   vector<Region> RgVec;
   CompPictures_A( GPic, RootRg, RgVec);
   DumpPictures( dirNames,RgVec);
}

void DumpMserRg::DumpPictures(char*dirNames,vector<Region>&ReginVec)
{
    //IniPic
	int i,j;	 
	C24BitMap GMPic,CPic;
	CPic.FormatF(OPic.Width,OPic.Height);
	GMPic = CPic;
	GMPic.CleanPic(0);
	CPic.CleanPic(0);

	char NameBuff[100],regioninfo[100];
    Loopj(ReginVec.size())
	{
		//AnalysisObj(ObjVec[j]);
		//if(ObjVec[j].IfValid==-1)continue;	
		CPic.CleanPic(0);
		//  KBPic.SetPenColor(ObjVec[j].RR, ObjVec[j].GG, ObjVec[j].BB);
		CPic.RandPenColor();
		GMPic.RandPenColor();
		//CPic.SetPenColor(BOUND(ReginVec[j].level_*3/2,0,255),0,0);
		//KBPic.RandPenColor();
		sprintf(NameBuff,"%s\\%06i.bmp",dirNames,ReginVec[j].r_id);

		//if(ReginVec[j].PtVec.size()<20)
		//	continue;

		Loopi(ReginVec[j].PtVec.size())
		{
			CPic.SigDot(ReginVec[j].PtVec[i].x,ReginVec[j].PtVec[i].y);
			GMPic.SigDot(ReginVec[j].PtVec[i].x,ReginVec[j].PtVec[i].y);
		} 
		
       sprintf(regioninfo,"id:%i,area:%i,ratio:%.2lf,thr:%i",ReginVec[j].r_id,
		   ReginVec[j].PtVec.size(),ReginVec[j].normal_variation,ReginVec[j].level_);

		MergeTxtStr(CPic, 10,10,20, regioninfo,255,12,12);
		CPic.Save(NameBuff);
	}
 
	GMPic.Save("res.bmp");
}


//------------------------------------------------------------------------------

void  CompleteRegion_A(FillFindRegion& RgFill,C256BitMap&Pic,Region&Rg)
//void MserRegion::CompleteRegion(FillFindRegion&Region,C256BitMap&Pic,int thrd,vector<Point> &PtVec)
{
	int i,j;
	RgFill.fill_threshold = Rg.level_;
	
	for (i=Rg.left;i<=Rg.right;i++)
	{
		for(j=Rg.top;j<=Rg.bottom;j++)
			RgFill.Matrix[i][j] = 0;
	}
	
	
	if(Rg.PtVec.size()>0)
	{
		int x, y;
		x = Rg.PtVec[0].x;
		y = Rg.PtVec[0].y;
		
		//RgFill.PtVec.clear();
		//RgFill.Pic.Save("dest1.bmp");
		RgFill.FillSeed(x,y);
		
		
		if(RgFill.PtVec.size()>Rg.PtVec.size())
			Rg.PtVec = RgFill.PtVec;
		
	}
	
} 

void  CompPictures_A(C256BitMap&GPic, Region*rootRegion,vector<Region>&ReginVec)
{   
	int i; 
    FillFindRegion RegionFill(GPic);
	rootRegion->SaveMserRegion(ReginVec, 0.4);
//	MarkInvalidRegionByCrossValidate(ReginVec);
    Loopi(ReginVec.size())
	{
		if(ReginVec[i].IsValid)
			CompleteRegion_A(RegionFill,GPic,ReginVec[i]);
		else
			ReginVec[i].PtVec.clear();
	}
}
#endif