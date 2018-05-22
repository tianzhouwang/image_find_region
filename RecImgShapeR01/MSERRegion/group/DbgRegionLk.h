#ifndef Debug_Region_Link_Head
#define Debug_Region_Link_Head


double GetRegionLink(vector<Region>&RegionVec,int Idx1, int Idx2, int PicWidth, int PicHeight)
{
   double Val = 
	   SkipNotSimilarRegion(RegionVec[Idx1], RegionVec[Idx2]);


   int i,j,t;
   
   //	int xDir[8],yDir[8];
   //	xDir[0] =1; xDir[1] =1; xDir[2] =0; xDir[3] = -1; xDir[4] = -1; xDir[5] =-1; xDir[6] = 0; xDir[7] = 1;
   //  yDir[0] =0; yDir[1] =1; yDir[2] =1; yDir[3] =  1; yDir[4] =  0; yDir[5] =-1; yDir[6] =-1; yDir[7] =-1;
   
   int xDir[12],yDir[12];
   xDir[0] =1; xDir[1] =1; xDir[2] =0; xDir[3] = -1; xDir[4] = -1; xDir[5] =-1; xDir[6] = 0; xDir[7] = 1;
   yDir[0] =0; yDir[1] =1; yDir[2] =1; yDir[3] =  1; yDir[4] =  0; yDir[5] =-1; yDir[6] =-1; yDir[7] =-1;
   
   xDir[8] = 2; xDir[9] =-2; xDir[10] = -2; xDir[11] =  2;
   yDir[8] = 1; yDir[9] = 1; yDir[10] = -1; yDir[11] = -1;
   
   double searchLenNorm[12];
   vector< vector<int> > RegionLink;
   for(i=0 ; i < 12 ; i++)
   {
	   double dx = xDir[i];
	   double dy = yDir[i];
	   searchLenNorm[i] = sqrt(dx * dx + dy * dy);
   }
   
   if(RegionLink.size()!= RegionVec.size())
   {
	   RegionLink.resize(RegionVec.size());
	   for(i=0;i< RegionVec.size();i++)
	   {
		   RegionLink[i].resize(RegionVec.size());
	   }
   }
   
   
   vector<int> LabelVec;
   LabelVec.resize(PicWidth*PicHeight);
   
   Loopi(PicWidth*PicHeight)
	   LabelVec[i]= -1;
   
   for (i=0;i<RegionVec.size();i++)
   {
	   for(j=0;j<RegionVec[i].PtVec.size();j++)
	   {
		   int x,y;
		   x = RegionVec[i].PtVec[j].x;
		   y = RegionVec[i].PtVec[j].y;
		   LabelVec[y*PicWidth+x] = i;
	   }
   }
   
   //for(i=0; i<RegionVec.size(); i++)
   i = Idx1;
   {
	   //if(SkipInvalidRegion(RegionVec[i]))continue;
	   //if(i==86)
	   //{
	   //	int debug_here;debug_here =1;
	   //}
	   int  max_size      = max(RegionVec[i].rwidth, RegionVec[i].rheight);
	   int  stroke_length = RegionVec[i].Feature.StrokeMean * 8;
	   int  size_length   = max_size  * 3/2; 
	   int  search_length = size_length; //min(stroke_length, size_length); 
	   
	   
	   for (j=0;j<12;j++)
	   {
		   int region_seach_len;
		   region_seach_len = search_length / searchLenNorm[j] ;//(xDir[j]*yDir[j]==0) ? search_length : (search_length* 3/2); 
		   int startX,startY;		  
		   
		   startX = RegionVec[i].GeoX + RegionVec[i].rwidth /3 * xDir[j];
		   startY = RegionVec[i].GeoY + RegionVec[i].rheight/3 * yDir[j];
		   
		   
		   for(t=0;t<region_seach_len;t++)
		   {
			   
			   int newX,newY;
			   newX = startX + t * xDir[j];
			   newY = startY + t * yDir[j];
			   
			   if(newY > (PicHeight-1)|| newX >(PicWidth-1) || newX <0 || newY <0 )
				   break;
			   
			   int RgLabel;
			   RgLabel = LabelVec[newY * PicWidth +newX];
			   
			   if(RgLabel==-1) continue;
			   
			   //&& RegionSimilar[i][RgLabel];
			   if(RgLabel!= i )
			   {
				   if(!SkipInvalidRegion(RegionVec[RgLabel])) 
				   {
					   if( RegionXYOverLap(RegionVec[i],RegionVec[RgLabel]))
						   RegionLink[i][RgLabel] = 1;
				   }
				   break;
			   }
		   }
	   }
	   
	}
   return 1;
}


#endif