 

#ifndef LB_CONNECTED_H
#define LB_CONNECTED_H
 
#include <math.h>
#include <vector>
#include <algorithm>
#include "ktypes.h"
#include "kutil.h"


using namespace std;

/*struct ImgPix
{
  int v1,v2,v3;
};*/

int LBCheckBoundPoint(unsigned char*PixPt,int x,int y,int w, int h)
{   
    if(x==0 || x==(w-1)) return 0;  
	if(y==0 || y==(h-1)) return 0;
	
	int val, v1, v2, v3, v4;
	
    val = *PixPt;
    v1 = *(PixPt+1); if(val!=v1) return 1;
    v2 = *(PixPt-1); if(val!=v2) return 1;
    v3 = *(PixPt+w); if(val!=v3) return 1;
    v4 = *(PixPt-w); if(val!=v4) return 1;
	return 0;
}

class TMConnectedComponents
{
public:
    TMConnectedComponents(int soft_maxlabels) : labels(soft_maxlabels) {
		clear();
    }
    void clear() {
		std::fill(labels.begin(), labels.end(), Similarity_());
		highest_label = 0;
    }
     
    int connected(const unsigned char *img, int *out,
				  int width, int height,  int K8_connectivity);



private:
    struct Similarity_ {
		Similarity_() : id(0), sameas(0) {}
		Similarity_(int _id, int _sameas) : id(_id), sameas(_sameas) {}
		Similarity_(int _id) : id(_id), sameas(_id) {}
		int id, sameas, tag;
    };

    bool is_root_label(int id) {
		return (labels[id].sameas == id);
    }
    int root_of(int id) {
	while (!is_root_label(id)) {
	    // link this node to its parent's parent, just to shorten
	    // the tree.
	    labels[id].sameas = labels[labels[id].sameas].sameas;

	    id = labels[id].sameas;
	}
	return id;
    }
    bool is_equivalent(int id, int as) {
		return (root_of(id) == root_of(as));
    }
    bool merge(int id1, int id2) {
	if(!is_equivalent(id1, id2)) {
	    labels[root_of(id1)].sameas = root_of(id2);
	    return false;
	}
	return true;
    }
    int new_label() {
		// Thanks to Pierre Grinspan for suggesting using capacity()
		// instead of size() here.
		if(highest_label+1 > labels.capacity())
			labels.reserve(highest_label*2);
		labels.resize(highest_label+1);
		labels[highest_label] = Similarity_(highest_label);
		return highest_label++;
    }


    
    void label_image(const unsigned char *img, int *out,
					 int width, int height, int K8_connectivity);
   
    int relabel_image(int *out, int width, int height);


    std::vector<Similarity_> labels;
	public:
    int highest_label;
	int SAME(unsigned char c1,unsigned char c2) {if(c1==c2)return 1; return 0;};
	vector<float> lbAngles;
	void computePointsAngle(vector<KRegion>&regionVec,int*Ptlbs);
	int picwidth,picheight;
};

 
void TMConnectedComponents::computePointsAngle(vector<KRegion>&regionVec,int*Ptlbs)
{
  int i,j,t;
  int offsetY[] = {-1, -1, -1,  0, 0,  1, 1, 1};
  int offsetX[] = {-1,  0,  1, -1, 1, -1, 0, 1};
  int pixoffset[8];
  int nNeighbors = 8;
  for (i=0;i<nNeighbors;i++)
  {
    pixoffset[i] = offsetX[i] + picwidth * offsetY[i];
  }
  

  Loopi(regionVec.size())
  {
   Loopj(regionVec[i].contourPts.size())
   {
	   int x,y;
	   x = regionVec[i].contourPts[j].x;
	   y = regionVec[i].contourPts[j].y;
	   int pixel_pos = y*picwidth+x;
	   int pixlabel = Ptlbs[pixel_pos];
       int nbpattern = 0;
	   for(t=0;t<nNeighbors;t++)
	   {  
		   int offset = pixoffset[t];
		   int nblabel = Ptlbs[pixel_pos+offset];
		   if(nblabel != pixlabel )
			   nbpattern+= 1<<t;
	   }

 	  regionVec[i].contourPts[j].PatternIdx = nbpattern;
          regionVec[i].contourPts[j].angle = contour_angle[nbpattern];
	  regionVec[i].contourPts[j].dx    = angle_dx[nbpattern];
	  regionVec[i].contourPts[j].dy    = angle_dy[nbpattern];
   }
  }
   

}

int TMConnectedComponents::connected(const unsigned char *img, int *labelimg,
			       int width, int height, int K8_connectivity)
{
	lbAngles.resize(width*height);
	picwidth  = width;
	picheight = height;

    label_image(img,labelimg, width,height, K8_connectivity);
    return relabel_image(labelimg, width,height);

}



 
void TMConnectedComponents::label_image(const unsigned char *img, int *labelimg,
				 int width, int height, const int K8_CONNECTIVITY)
{
	int c,r;
    const unsigned char *row = img;
    const unsigned char *last_row = 0;
    struct Label_handler {
	Label_handler(const unsigned char *img, int *limg) :
	    piximg(img), labelimg(limg) {}
	int &operator()(const unsigned char *pixp) { return labelimg[pixp-piximg]; }
	const unsigned char *piximg;
	int *labelimg;
    } label(img, labelimg);

    clear();

    label(&row[0]) = new_label();

    // label the first row.
    for(  c=1, r=0; c<width; ++c) {
	if(SAME(row[c], row[c-1]))
	    label(&row[c]) = label(&row[c-1]);
	else
	    label(&row[c]) = new_label();
    }

    // label subsequent rows.
    for( r=1; r<height; ++r)    {
	// label the first pixel on this row.
	last_row = row;
	row = &img[width*r];

	if(SAME(row[0], last_row[0]))
	    label(&row[0]) = label(&last_row[0]);
	else
	    label(&row[0]) = new_label();

	// label subsequent pixels on this row.
	for(int c=1; c<width; ++c)	{
	    int mylab = -1;

	    // inherit label from pixel on the left if we're in the same blob.
	    if(SAME(row[c],row[c-1]))
		mylab = label(&row[c-1]);
	    for(int d=(K8_CONNECTIVITY?-1:0); d<1; ++d) {
		// if we're in the same blob, inherit value from above pixel.
		// if we've already been assigned, merge its label with ours.
		if(SAME(row[c], last_row[c+d])) {
		    if(mylab>=0) merge(mylab, label(&last_row[c+d]));
		    else mylab = label(&last_row[c+d]);
		}
	    }
	    if(mylab>=0) label(&row[c]) = static_cast<int>(mylab);
	    else label(&row[c]) = new_label();

	    if(K8_CONNECTIVITY && SAME(row[c-1], last_row[c]))
		merge(label(&row[c-1]), label(&last_row[c]));
	}
    }
}

 
int TMConnectedComponents::relabel_image(int *labelimg, int width, int height)
{
    int newtag = 0;
    for(int id=0; id<labels.size(); ++id)
	if(is_root_label(id))
	    labels[id].tag = newtag++;

    for(int i = 0; i<width*height; ++i)
	labelimg[i] = labels[root_of(labelimg[i])].tag;

    return newtag;
}


#endif // _CONNECTED_H
