#ifndef Connect_GraphNode_Head
#define Connect_GraphNode_Head
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;


bool GetMaxConnected(const vector< vector<int> >& initialChains, vector< vector<int> >& chains);
void GetBox(vector <KRegion>&ResultVec,vector<int>&vec,int&x1,int&y1,int&x2, int&y2);

class NodeSimilarity 
{
public:
	NodeSimilarity() 
	{id = -1;sameas=-1;
	   group_num = 0;}; 
       //NodeSimilarity(int _id, int _sameas){id=_id;sameas=_sameas};
       NodeSimilarity(int _id)
	   {id=_id;sameas=_id;
	   group_num = 0;};
	   
       int id, sameas, tag;
       int     group_num;
};

#define VertexLabed 123
class ConnectedLinks
{
public:	
	int highest_label;
	ConnectedLinks(){ highest_label=-1;labels.resize(1024);};
	~ConnectedLinks(){};
	vector<NodeSimilarity> labels;	 //标记向量
	bool  is_root_label(int id);
    void  clear();
	int  root_of(int id) ;
	bool  is_equivalent(int id, int as) ;
	bool  merge(int id1, int id2) ;
	int  new_label() ;
	void AddEdge(string item1,string item2);
	void AddEdge(int item1,int item2);
	map<string,int>   LabedMarks;
	map<string,int> VertexTag;
	void GetLabelInfo(vector< vector<int> >&LbChain);
	map<int,int>         IdNumCount;//每个组里面有多少个成员


	vector<string>  VertexList;
	map<string,int> VertexCnt;
    map<string,int> VertexIdx;
};

void ConnectedLinks::AddEdge(int item1,int item2)
{
   char buffer1[50];
   char buffer2[50];
   sprintf(buffer1,"%i",item1);
   sprintf(buffer2,"%i",item2);
   AddEdge(buffer1,buffer2);

}
void ConnectedLinks::GetLabelInfo(vector< vector<int> >&LbChain)
{
	int i;
    LbChain.clear();
	
	int newtag = 0;
	
	for(int id=0; id<this->labels.size(); ++id)
	{
		//fprintf(file,"%i,%i\n",id,labels[id].sameas);
		if(this->is_root_label(id))
		{ this->labels[id].tag = newtag++;
		//printf("%i,%i,tag:%i\n",id,labels[id].sameas,labels[id].tag);
		IdNumCount[this->labels[id].tag] = 0;
		}
	}
	// fclose(file);
	
	
	

	for(i = 0; i<VertexList.size(); i++)		
	{
		//if(i<10)continue;
		string word = VertexList[i];
		int tag = VertexTag[word];
		if(LabedMarks[word]==VertexLabed)
		{
			tag = this->labels[this->root_of(tag)].tag;
			IdNumCount[tag]+=1;
			VertexTag[word] = tag;
			//fprintf(file,"%s,%i\n",word.c_str(),tag);
		}
	} 
	
     LbChain.resize(newtag);
	 
	 for(i = 0; i<VertexList.size(); i++)		
	 {
		 //if(i<10)continue;
		 string word = VertexList[i];
		 int Idx= atoi(word.c_str());
		 int tag = VertexTag[word];
		 //if(IdNumCount[tag]>10)
		 LbChain[tag].push_back(Idx);
	} 

	/*printf("==================================\n");
	
	for(i=0;i<newtag;i++)
	{ 
		
		printf("%i,",IdNumCount[i]);
	}
	
FILE*file = fopen("wordtag.txt","wt+");
	
	for(i = 0; i<VertexList.size(); i++)		
	{
		//if(i<10)continue;
		string word = VertexList[i];
		int tag = VertexTag[word];
		//if(IdNumCount[tag]>10)
		{
			
			fprintf(file,"%s,%i\n",word.c_str(),tag);
		}
	} 
    fclose(file);	*/
}

void ConnectedLinks::AddEdge(string item1,string item2)
{

	   if(VertexCnt[item1]==0)
	   {
		   VertexList.push_back(item1);
		   VertexCnt[item1] =1;
	   }
	   else
		   VertexCnt[item1]++;

	   if(VertexCnt[item2]==0)
	   {
		   VertexList.push_back(item2);
		   VertexCnt[item2] =1;
	   }
	   else
		   VertexCnt[item2]++;

	   bool Lb1,Lb2;
	   Lb1 =  (LabedMarks[item1]==VertexLabed);
	   Lb2 =  (LabedMarks[item2]==VertexLabed);
	   
	   if((!Lb1)&&(!Lb2))
	   {
		   int Lb = this->new_label();
		   VertexTag[item1]=Lb;
		   VertexTag[item2]=Lb;
       }
       if((!Lb1)&&(Lb2))
	   {
		   int Lb = VertexTag[item2];
		   VertexTag[item1]=Lb;
	   }
	   if((Lb1)&&(!Lb2))
	   {
		   int Lb = VertexTag[item1];
		   VertexTag[item2]=Lb;
	   }
	   if((Lb1)&&(Lb2))
	   {
		   int Lb1 = VertexTag[item1];  
		   int Lb2 = VertexTag[item2];
		   this->merge(Lb1, Lb2);
	   }

	    LabedMarks[item1] = VertexLabed;
	    LabedMarks[item2] = VertexLabed;
}


bool ConnectedLinks::is_root_label(int id)
{
	return (labels[id].sameas == id);
};

void ConnectedLinks::clear()
{  
	fill(labels.begin(), labels.end(), NodeSimilarity());
	highest_label = -1;};
	
	
	int ConnectedLinks::root_of(int id) 
	{
		while (!is_root_label(id)) {
			labels[id].sameas = labels[labels[id].sameas].sameas;
			id = labels[id].sameas;
		}
		return id;
	} ;
	
	bool ConnectedLinks::is_equivalent(int id, int as) 
	{
		return (root_of(id) == root_of(as));
	};
	
	bool ConnectedLinks::merge(int id1, int id2) 
	{
		if(!is_equivalent(id1, id2)) 
		{
			labels[root_of(id1)].sameas = root_of(id2);
			return false;
		}
		return true;
	};
	
	int ConnectedLinks::new_label() 
	{
		if(unsigned int(highest_label+2) > labels.size())
			labels.resize(highest_label*2);
		
		highest_label++;
		labels[highest_label] = NodeSimilarity(highest_label);
		
		return highest_label;
}


bool GetMaxConnected(const vector< vector<int> >& initialChains, vector< vector<int> >& chains)
{
    char Buff[40];
	ConnectedLinks Cnnt;
	string str1,str2;
	int i,j;
	int Idx1,Idx2,Idx3;

	for(i=0;i< initialChains.size();i++)
	{
		//if(initialChains[i].size()<1)continue; // 没可能
        
		int Idx1;
		Idx1 = initialChains[i][0];
		sprintf(Buff,"%i",Idx1);str1 = Buff;
		for(j=1;j<initialChains[i].size();j++)
		{
		Idx2 = initialChains[i][j];
		sprintf(Buff,"%i",Idx2);str2 = Buff;
	 	Cnnt.AddEdge(str1,str2);
		}
	}
	
   Cnnt.GetLabelInfo(chains);

	return 1;
}

void  GetBoxRec(vector <Region>&ResultVec,vector<int>&vec,CtRect&Rect)
{ 
	Rect.RR = Rect.GG =Rect.BB =0;
    int x1,y1,x2,y2;
	int i,j,idx;
	idx = vec[0];
	x1 = ResultVec[idx].left;
	y1 = ResultVec[idx].top;
	x2 = ResultVec[idx].right;
	y2 = ResultVec[idx].bottom;
	for(i=0;i<vec.size();i++)
	{    
		Region Rg;
		idx = vec[i];
		Rg = ResultVec[idx];
		
		if(x1>Rg.left)   x1 = Rg.left;
		if(y1>Rg.top)    y1 = Rg.top;
		if(x2<Rg.right)  x2 = Rg.right;
		if(y2<Rg.bottom) y2 = Rg.bottom;
        
		Rect.RR += Rg.RR;
		Rect.GG += Rg.GG;
		Rect.BB += Rg.BB;
	}


	Rect.left    = x1;  Rect.top     = y1;
	Rect.right   = x2;  Rect.bottom  = y2;

	Rect.RR /= float(vec.size());
	Rect.GG /= float(vec.size());
	Rect.BB /= float(vec.size());
	
}
void GetBox(vector <KRegion>&ResultVec,vector<int>&vec,int&x1,int&y1,int&x2, int&y2)
{
	int i,j,idx;
	idx = vec[0];
  	x1 = ResultVec[idx].left;
 	y1 = ResultVec[idx].top;
  	x2 = ResultVec[idx].right;
 	y2 = ResultVec[idx].bottom;
	for(i=0;i<vec.size();i++)
	{
		idx = vec[i];
		if(x1>ResultVec[idx].left)
			x1 = ResultVec[idx].left;
		
		if(y1>ResultVec[idx].top)
			y1 = ResultVec[idx].top;
		
		if(x2<ResultVec[idx].right)
			x2 = ResultVec[idx].right;
		
		if(y2<ResultVec[idx].bottom)
			y2 = ResultVec[idx].bottom;
	}
	
}
// ConnectedLinks

//bool mergePairs()
#endif