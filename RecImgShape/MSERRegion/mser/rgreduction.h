 
#ifndef Region_Reduction_Head
#define Region_Reduction_Head
#include "mregion.h"
#include "mser.h"

/// \brief  把树节点的
/// \param  Obj 待查找的元素
/// \return 临近物体的个数, 且临近物体的指针将保存在AdjacentObjs数组中 
void InsertChildNodes(Region*Parent,vector<Region *>&SubNodes)
{
	int i;
	Region *Pt;
	if(SubNodes.size()>0)
	{        
		Parent->child_ = SubNodes[0];
		Parent->child_->next_ = NULL;
		Pt = Parent->child_;
		
		for(i=1; i<SubNodes.size(); i++)
		{
			Pt->next_ = SubNodes[i];
			Pt = SubNodes[i];
		}
		Pt->next_ =NULL;/**/
	}
	else
	{
		Parent->child_   = NULL;
	}
	
}

/* */

/*
procedure LINEAR-REDUCTION(T)
2: if nchildren[T] = 0 then
3:		return T
4: else if nchildren[T] = 1 then
5:		c = LINEAR-REDUCTION(child[T])
6:		if var[T] .<= var[c] then
7:			link-children(T, children[c])
8:			return T
9:		else
10:			return c
11:		end if
12:else // nchildren[T] =. 2
13:		for each c 2 children[T] do
14:			link-children(T, LINEAR-REDUCTION(c))
15:		end for
16:		return T
17:end if
end procedure
*/

Region *linear_reduction(Region *T)
{
 	 
	if (!T->child_ ) {
		// has no child
		return T;

	} else if (T->childNum == 1) {
		 
		Region *c = linear_reduction(T->child_);

		if (c->normal_variation < T->normal_variation ) {
			
			return c; // c is better

		} else {
			// remove c, link T to its new child, "c's child"
			T->child_   = c->child_  ;
			T->childNum = c->childNum;
			
			return T; // T is better
			
		}
	} else {

		//void  changeCild(Region *NodeC, Region *PrevNode,Region *CurNode);
        vector<Region*> TempArr;
        TempArr.clear();

		for(Region * child = T->child_; child; child = child->next_)
		{
		   Region * recNode =  linear_reduction(child);
           //recNode->next_ = NULL;
           TempArr.push_back(recNode);
		}
		
	   InsertChildNodes(T,TempArr);
	 
		return T;
	}

	return T;
}


/*
procedure TREE-ACCUMULATION(T)
2: if nchildren[T] .>= 2 then
3:		C = empty;
4:		for each c in children[T] do
5:			C = C union TREE-ACCUMULATION(c)
6:		end for
7:		if var[T] <=. min-var[C] then
8:			discard-children(T)
9:			return T
10:		else
11:			return C
12:		end if
13: else // nchildren[T] = 0
14:		return T
15: end if
end procedure
*/


void InsertSubNodeVec(vector<Region*>&RegionVec, vector<Region*>&SubNodeVec)
{
  int i;
  for(i=0;i<SubNodeVec.size();i++)
  {
	  RegionVec.push_back(SubNodeVec[i]);
  }
}

vector<Region *> tree_accumulation(Region *T )
{
	vector<Region*> TempArr;
    int i;
	if(1)//T->childNum>=2 && T->PtVec.size()< GRegionMaxSize)
    {
	 
	  /* if(T->r_id == 53329)//9531)//7032 )
		{
			int aa;
			aa = 1;
		}	*/
        //vector<Region*> RegionVec;
		//RegionVec.clear();
        double child_var ;
		child_var = 9999;//T->child_->normal_variation;

		Region *Prev;
		
		TempArr.clear();
		int AddCond;
		AddCond = 1;
        int maxchild;
		maxchild =-1;

		for(Region * child = T->child_; child; child = child->next_)
		{
           //if(child->PtVec.size()<32)continue;
		  
	  	   //  if(child_var > child->normal_variation)
		   //	  child_var = child->normal_variation;   
		   int child_id; 
		   child_id = child->r_id;
		   InsertSubNodeVec(TempArr, tree_accumulation(child));
		   if(child->area_>maxchild)
			   maxchild = child->area_;
		  //----------------------------------------------------------------------
		}

        //********************以下代码处理细线干扰********************
		//************************************************************
		//************************************************************
        double stable_child_area;
		int    stable_child_cnt;
		int    big_child_cnt;

		stable_child_cnt = stable_child_area = big_child_cnt= 0;

		for(i=0;i<TempArr.size();i++)
		{
          if(child_var > TempArr[i]->normal_variation)
			{ // child_var = TempArr[i]->normal_variation;
		//-------------------------------------------------------
		     if( (TempArr[i]->bottom - TempArr[i]->top+1 ) > 2
				  &&
				 (TempArr[i]->right  - TempArr[i]->left+1) > 2 )
			   child_var = TempArr[i]->normal_variation;
			  else
			  {
				 // if( child_var <0.01)
                  //     child_var = 0.01;
			  }/**/
		//-------------------------------------------------------
			}
          
		  if(TempArr[i]->normal_variation< 0.08)//0.008)
		  {
			  stable_child_cnt  ++;
			  stable_child_area += TempArr[i]->area_;
			  int child_id = TempArr[i]->r_id;
			  if(TempArr[i]->area_ > T->area_ /100)
			  {
				  big_child_cnt++;
			  }
		  }
		}

        //************************************************************
		//************************************************************
		//************************************************************
      
		double OcuupyRatio = double(T->area_)/
			                 (double(T->bottom - T->top+1)*double(T->right - T->left+1));
           
        
        if(child_var<0.008 && maxchild >800 )
			AddCond = 0;

		if(OcuupyRatio>0.9 &&T->area_ >20000)
            AddCond = 0;

		int not_picoutter = ( (T->right - T->left+1 ) < T->PicWidth * 4/5 );

		if(T->normal_variation <0.003  &&  not_picoutter)
		{
			if(big_child_cnt<6)
			{ 
				if(T->area_ / stable_child_area > 8)
				{ 
					T->child_    = NULL; TempArr.clear();
					TempArr.push_back(T);
					return TempArr;
				}
			};

			if(big_child_cnt>6 && big_child_cnt< 15)
			{ 
				if(T->area_ / stable_child_area > 20)
				{ 
					T->child_    = NULL; TempArr.clear();
					TempArr.push_back(T);
					return TempArr;
				}
			};
		}

		int NoiseCond =0;
		if(stable_child_area <32 && T->normal_variation <0.01 )//&& T->normal_variation <= child_var)
			 NoiseCond = 1; //子孙是微小噪声

		if((T->normal_variation <= child_var || T->normal_variation<0.003 || NoiseCond) 
			&& T->area_ < T->GRegionMaxSize && AddCond && not_picoutter)//父节点最稳定了
		{   
           
		   T->child_    = NULL;
		   TempArr.clear();
		   TempArr.push_back(T);
		   return TempArr;
		}
		else
		{ 
		  return  TempArr;
		}
	}
    else
	{
		TempArr.clear();
		TempArr.push_back(T);
		return TempArr;
	}

//--------------------------------
}

Region * cut_noiseNode(Region *T,int MinSize)
{
    int childNum,i;
    Region * Pt;
    vector<Region*> TempArr;
    TempArr.clear();
    childNum = 0;

	for(Region * child = T->child_; child; child = child->next_)
	{
		if(child->area_ < MinSize)continue;
		childNum++;
		//----------------------------------------------------------------------
		Region * recNode = cut_noiseNode(child,MinSize);
		TempArr.push_back(recNode); 
	}
    
	InsertChildNodes(T, TempArr);

	return T;
}

#endif