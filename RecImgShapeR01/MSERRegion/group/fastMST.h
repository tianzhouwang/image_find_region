  

/*===============DEMO CODE===============
  int i,j;
  vector< vector<int> >  chains;
  MSTConnectedLinks Cnnt;
  Cnnt.AddEdge(1,5); Cnnt.AddEdge(1,6);
  Cnnt.AddEdge(3,5); Cnnt.AddEdge(8,4);
  Cnnt.AddEdge(4,2); Cnnt.AddEdge(7,2);
  Cnnt.GetLabelInfo(chains);

 for(i=0;i<chains.size();i++)
 {  printf("\n================================\n");
	 for (j=0;j<chains[i].size();j++)
	     printf("%i,",chains[i][j]); 
 }*/
#ifndef Connect_GraphNode_Head
#define Connect_GraphNode_Head
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
//#include <string>
using namespace std;


bool GetMaxConnected(const vector< vector<int> >& initialChains, vector< vector<int> >& chains);
//void GetBox(vector <KRegion>&ResultVec,vector<int>&vec,int&x1,int&y1,int&x2, int&y2);

//最大联通分支的节点类
class MSTNodeSimilarity 
{
public:
	MSTNodeSimilarity() 
	{id = -1;sameas=-1;
	   group_num = 0;}; 
       //MSTNodeSimilarity(int _id, int _sameas){id=_id;sameas=_sameas};
       MSTNodeSimilarity(int _id)
	   {id=_id;sameas=_id;
	   group_num = 0;};
	   
       int id, sameas, tag;
       int     group_num;
};

//最大联通分支连接线
#define Int32 int
#define VertexLabed 123
class MSTConnectedLinks
{
public:	
	int highest_label;
	MSTConnectedLinks() { highest_label=-1;labels.resize(1024);};
	~MSTConnectedLinks(){};
	vector<MSTNodeSimilarity> labels;	 //标记向量
	bool  is_root_label(int id);
    void  clear();
	int   root_of(int id) ;
	bool  is_equivalent(int id, int as) ;
	bool  merge(int id1, int id2) ;
	int   new_label() ;
	//void AddEdge(Int32 item1,Int32 item2);
	void AddEdge(int item1,int item2);
	map<Int32,int>   LabedMarks;
	map<Int32,int>   VertexTag;
	void GetLabelInfo(vector< vector<int> >&LbChain);
	map<int,int>     IdNumCount;//每个组里面有多少个成员

	vector<Int32>  VertexList;
	map<Int32,int> VertexCnt;
    map<Int32,int> VertexIdx;
};

 

void MSTConnectedLinks::GetLabelInfo(vector< vector<int> >&LbChain)
{
	int i;
    LbChain.clear();
	
	int newtag = 0;
	
	for(int id=0; id<this->labels.size(); ++id)
	{	 
		if(this->is_root_label(id))
		{ this->labels[id].tag = newtag++;
		  IdNumCount[this->labels[id].tag] = 0;
		}
	}
	 
	for(i = 0; i<VertexList.size(); i++)		
	{
		 
		Int32 word = VertexList[i];
		int tag = VertexTag[word];
		if(LabedMarks[word]==VertexLabed)
		{
			tag = this->labels[this->root_of(tag)].tag;
			IdNumCount[tag]+=1;
			VertexTag[word] = tag;	 
		}
	} 
	
     LbChain.resize(newtag);
	 
	 for(i = 0; i<VertexList.size(); i++)		
	 {	  
		 Int32 word = VertexList[i];
		 int Idx=  word;
		 int tag = VertexTag[word];
		  
		 LbChain[tag].push_back(Idx);
	} 
	 
}

void MSTConnectedLinks::AddEdge(Int32 item1,Int32 item2)
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


bool MSTConnectedLinks::is_root_label(int id)
{
	return (labels[id].sameas == id);
};

void MSTConnectedLinks::clear()
{  
	fill(labels.begin(), labels.end(), MSTNodeSimilarity());
	highest_label = -1;};
		
	int MSTConnectedLinks::root_of(int id) 
	{
		while (!is_root_label(id)) {
			labels[id].sameas = labels[labels[id].sameas].sameas;
			id = labels[id].sameas;
		}
		return id;
	} ;
	
	bool MSTConnectedLinks::is_equivalent(int id, int as) 
	{
		return (root_of(id) == root_of(as));
	};
	
	bool MSTConnectedLinks::merge(int id1, int id2) 
	{
		if(!is_equivalent(id1, id2)) 
		{
			labels[root_of(id1)].sameas = root_of(id2);
			return false;
		}
		return true;
	};
	
	int MSTConnectedLinks::new_label() 
	{
		if( (unsigned int)(highest_label+2) > labels.size())
			labels.resize(highest_label*2);
		
		highest_label++;
		labels[highest_label] = MSTNodeSimilarity(highest_label);
		
		return highest_label;
}


bool GetMaxConnected(const vector< vector<int> >& initialChains, vector< vector<int> >& chains)
{
    char Buff[40];
	MSTConnectedLinks Cnnt;
	Int32 str1,str2;
	int i,j;
	int Idx1,Idx2,Idx3;

	for(i=0;i< initialChains.size();i++)
	{
		//if(initialChains[i].size()<1)continue; // 没可能
		int Idx1;
		Idx1 = initialChains[i][0];
		int str1 =  Idx1; 
		for(j=1;j<initialChains[i].size();j++)
		{
		Idx2 = initialChains[i][j];
		int str2 = Idx2; //sprintf(Buff,"%i",Idx2);str2 = Buff;
	 	Cnnt.AddEdge(str1,str2);
		}
	}
	
   Cnnt.GetLabelInfo(chains);

   return 1;
}

 
#endif
