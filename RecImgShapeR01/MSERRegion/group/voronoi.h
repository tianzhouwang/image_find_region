#ifndef  VORONOI_HEAD_H
#define  VORONOI_HEAD_H

#include <stdio.h>
#include <math.h>
#include "c256bitmap.h"
#include <vector>

using namespace std;

#define SelectMinMax(val_,minv_,maxv_) { if(val_>maxv_) maxv_=val_; if(val_<minv_)minv_=val_; }

//--------------------------------------------------------------------
#define	BYTE_WIDTH		8       /* [bit] */
#define	BLACK		1
#define WHITE		0
#define	NAMELEN		500

/* initial size and increments */
#define INITPIXEL	5000000  /* initial size of black pixels */
#define INCPIXEL        100000   /* increments of black pixels */
#define INITNEIGHBOR	50000    /* initial size of the array "neighbor" */
#define INCNEIGHBOR     10000    /* its increments */

#define INITLINE        50000    /* initial size of the array "lineseg" */
#define INCLINE         10000    /* its increments */
#define SITE_BOX        10000

/* MAX valures */
#define UNSIGNED_MAX    0xFFFFFFFF /* max. of unsigned int */
#define NOLABEL		0xFFFF     /* max. of unsigned short */
#define LABELMAX        NOLABEL

#define YES		1
#define NO		0
#define RIGHTANGLE	90
#define M_PI		3.14159265358979323846

/* for sweep2 */
#define LE		0
#define RE		1
#define DELETED		-2

/* An end point of a Voronoi edge is on the border of an image,
   its sitenbr is FRAME (erase.c, output.c) */
#define	FRAME		-1

/* for hash table */
#define M1		10007	/* prime number */
#define M2		3	/* prime number */
#define NIL		0
#define NODATA		-1

/* default values */
#define SAMPLE_RATE	7
#define	NOISE_MAX	20
#define FREQ_RATE	0.5
#define Ta_CONST        40
#define SMWIND          2    /* for smoothing the freq. distribution
			        of distance */
                          
/* voronoi�x�γ������ */
#define OUTPUT		1
#define NO_OUTPUT	0


typedef unsigned short Coordinate;
typedef unsigned short Label;
typedef unsigned int   NumPixel;
typedef unsigned int   Key;
typedef unsigned int   HashVal;
typedef double         Coord;

/* the structure for a vector */
typedef struct{
  int x;
  int y;
} Vector;

/*  the structure for a binary image */ 
typedef struct{
  unsigned char *image;
  Coordinate imax, jmax; /* width and height of the image */
} ImageData;

/*  the structure for a pixel */
typedef struct{
  Label label;     /* label */
  Coordinate xax,yax; /* coordinates */
} BlackPixel;

/* ���B�Y�ɷ֤��������ˤ�, ���������Τοk����L����������� */
/*
typedef struct{
    unsigned short x,y;
    unsigned short dx,dy;
    unsigned int bpn;
} Component;
*/

/* �O�Ӥ����B�Y�ɷ֤��v�S���������
   the structure for a neighboring relation
   between connected components(CC's)

         lab2
       -------
       |  x  |
       -------
    --- /
    | |/ angle
    |x|-----
    | |
    ---
    lab1
*/
typedef struct{
  float dist;       /* min. distance between CC's */
  float angle;      /* angle between CC's */
  Label lab1, lab2; /* labels of CC's */
} Neighborhood;

/* voronoi�x��ʼ��?�K���, �ɤ��O���B�Y�ɷ��g�Ǥ��뤫���������
   the structure for a Voronoi edge */
typedef struct{
  int	sp,ep;   //��ʼ�㣬������
  Coordinate xs,xe,ys,ye;  /* + (xs,ys)
				\
				\  Voronoi edge
				\
				+ (xe,ye)
				*/
  Label lab1,lab2;  /* this Voronoi edge is between
			      CC of a label "lab1" and that of lab2 */
  unsigned short yn;
} LineSegment;

/* voronoi�x�ζ˵�
   the structure for a Voronoi point */
typedef struct node{
  int line;
  struct node *next;
} EndPoint;

/* �O���B�Y�ɷ֤Υ�٥�ˌ�����ϥå����
   the structure for a hash table for
   representing labels of CC's */
typedef struct hash {
  Label lab1;
  Label lab2;
  unsigned int entry;
  struct hash *next;
} HashTable;
/*
typedef struct hash {
    unsigned long id;
    unsigned int entry;
    struct hash *next;
} HashTable;
*/

/* ���Τ��������
   the structure for a rectangle */
typedef struct{
    Coordinate is,ie,js,je;
} RectV;

struct Freenode {
    struct Freenode	*nextfree;
};

struct Freelist {
    struct Freenode	*head;
    int			nodesize;
};

struct Point {
    float		x,y;
};

/* structure used both for sites and for vertices */
struct Site {
    struct Point	coord;
    int			sitenbr;
    int			refcnt;
    unsigned int	label;
};

struct Edge	{
    float		a,b,c;
    struct	Site 	*ep[2];
    struct	Site	*reg[2];
    int			edgenbr;
    unsigned int	lab1,lab2;
    float		dist;
};

struct Halfedge {
    struct Halfedge	*ELleft, *ELright;
    struct Edge		*ELedge;
    int			ELrefcnt;
    char		ELpm;
    struct Site		*vertex;
    float		ystar;
    struct Halfedge	*PQnext;
};
//====================================================================
/* sort sites on y, then x, coord */
int scomp(const void *s1, const void *s2)
{
	float s1x = ((struct Point *)s1)->x;
	float s1y = ((struct Point *)s1)->y;
	float s2x = ((struct Point *)s2)->x;
	float s2y = ((struct Point *)s2)->y;
	
	if(s1y < s2y) return(-1);
	if(s1y > s2y) return(1);
	if(s1x < s2x) return(-1);
	if(s1x > s2x) return(1);
	return(0);
}
//LABELnbr = ++ln; /* �B�Y�ɷ֤��� */

//====================================================================
/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
deltax, deltay (can all be estimates).
Performance suffers if they are wrong; better to make nsites,
deltax, and deltay too big than too small.  (?) */
class MkVoronoi
{
public:
	float	xmin, xmax, ymin, ymax, deltax, deltay;
	void voronoi( Coordinate imax, Coordinate jmax);
	
	void geominit();
	struct Edge *bisect( struct Site *, struct Site * );
	struct Site *intersect( struct Halfedge *, struct Halfedge * );
	int right_of( struct Halfedge *, struct Point * );
	void endpoint( struct Edge *, int, struct Site *, Coordinate, Coordinate );
	float dist( struct Site *, struct Site * );
	void makevertex( struct Site * );
	void deref( struct Site * );
	void ref( struct Site * );
	struct Site * nextsite();
	//==========================================
	void CleanSites();
	
	//============================================
	int			nsites;
	int			siteidx;
	int			sqrt_nsites;
	int			nvertices;
	struct 	Freelist 	sfl;
	struct	Site		*bottomsite;
	
	int 			nedges;
	struct	Freelist 	efl;
	
	struct  Freelist	hfl;
	struct	Halfedge	*ELleftend, *ELrightend;
	int 			ELhashsize;
	struct	Halfedge	**ELhash;
	
	int 			PQhashsize;
	struct	Halfedge 	*PQhash;
	int 			PQcount;
	int 			PQmin;
	//=============================================
	Neighborhood	*neighbor;	/* �O���B�Y�ɷ��g���؏��� */
	LineSegment	*lineseg;	    /* ʼ��?�K������ˤȥ�٥� */
	struct	Site		*sites;
	unsigned int	point_edge;	/* ��voronoi�x�α��� */
	unsigned int	edge_nbr;	/* ��ȥ��Υ��ꥢvoronoi�x�� */
	/* ���֤α��� */
	
	unsigned int	NEIGHnbr;	/* �O�Ӥ����B�Y�ɷ֤νM���� */
	unsigned int	LINEnbr;	/* ��ȥǰ�Υ��ꥢvoronoi�x�� ���֤��� */
	unsigned int	Enbr;		/* voronoi�x����ȥ������B�Y �ɷ֤νM���� */
	long		SiteMax;	/* voronoi��η��Ť���� */
	//-------------------------------------------
	int ntry, totalsearch;
	
	/* edgelist.c */
	void ELinitialize();
	struct Halfedge *HEcreate( struct Edge *, int );
	void ELinsert( struct Halfedge *, struct Halfedge * );
	struct Halfedge *ELgethash( int );
	struct Halfedge *ELleftbnd( struct Point * );
	void ELdelete( struct Halfedge * );
	struct Halfedge *ELright( struct Halfedge * );
	struct Halfedge *ELleft( struct Halfedge * );
	struct Site *leftreg( struct Halfedge * );
	struct Site *rightreg( struct Halfedge * );
	/* memory.c */
	void freeinit( struct Freelist *, int );
	char *getfree( struct Freelist * );
	void makefree( struct Freenode *, struct Freelist * );
	char *myalloc( unsigned );
	char *myrealloc( void *, unsigned, unsigned int, size_t);
	
	/* output.c */
	void in_frame( float *, float *, float, struct Edge *, int,
		Coordinate, Coordinate );
	void s_in_frame( float *, float *, float *, float *,
		struct Edge *, Coordinate );
	void e_in_frame( float *, float *, float *, float *,
		struct Edge *, Coordinate, Coordinate );
	void frameout( float *, float *, float *, float *,
		int *, int *, struct Edge *, Coordinate, Coordinate );
	void out_ep2( struct Edge *, struct Site *, Coordinate, Coordinate );
	
	/* heap.c */
	void PQinsert( struct Halfedge *, struct Site *, float );
	void PQdelete( struct Halfedge * );
	int PQbucket( struct Halfedge * );
	int PQempty();
	struct Point PQ_min();
	struct Halfedge *PQextractmin();
	void PQinitialize();
	
	/* hash.c */
	HashVal hash1( Key );
	HashVal hash2( Key );
	void init_hash();
	Key key( Label, Label );
	int search( Label, Label );
	void enter( Label, Label, unsigned int );
	HashTable	*hashtable[M1+M2];
	/*
	free(area);
	free(sites);
	free(lineseg);
	*/
	//--------------------------------------
	unsigned int total_alloc; 
	
	MkVoronoi(int Width,int Height);
	~MkVoronoi();
};

MkVoronoi::MkVoronoi(int Width,int Height)
{
	total_alloc;// = 0;��Ҫ��ʼ��
    NEIGHnbr = LINEnbr = Enbr = SiteMax = 0;
	/* �O���B�Y�ɷ��g���؏���neighbor ���I��_�� */
	neighbor = (Neighborhood *)myalloc(sizeof(Neighborhood)* INITNEIGHBOR);
	
	/* ����lineseg ���I��_�� */
	lineseg = (LineSegment *)myalloc(sizeof(LineSegment)* INITLINE);
    sites=(struct Site *)malloc( Width* Height*sizeof(Site));
    
	ymax = xmax = 0;
	ymin = xmin = 9999; 
	init_hash();

}
MkVoronoi::~MkVoronoi()
{

}
void MkVoronoi::CleanSites()
{   
	int index = 0;
	int i;
	qsort(sites, nsites, sizeof *sites, scomp);
	xmin=sites[0].coord.x; 
	xmax=sites[0].coord.x;
	for(i=1; i<nsites; i+=1) {
		if(sites[i].coord.x < xmin) xmin = sites[i].coord.x;
		if(sites[i].coord.x > xmax) xmax = sites[i].coord.x;
	}
	ymin = sites[0].coord.y;
	ymax = sites[nsites-1].coord.y;
	for(i=0;i<nsites;i++){
		sites[i].sitenbr = i;
	}
	
	/* ���}���Ƥ��� sites ���ȥ���� */
	for(i=0;i<nsites;i++){
		if(!(sites[i].coord.x==sites[i+1].coord.x&&
			sites[i].coord.y==sites[i+1].coord.y)){
			sites[index].coord.x=sites[i].coord.x;
			sites[index].coord.y=sites[i].coord.y;
			sites[index].label=sites[i].label;
			sites[index].sitenbr=index;
			index++;
		}
	}
	nsites = index;
}

struct Site * MkVoronoi::nextsite()
{
	struct Site *s;
	if(siteidx < nsites) {
		s = &sites[siteidx];
		siteidx += 1;
		return(s);
	}
	else return( (struct Site *)NULL);
}

void MkVoronoi::voronoi( Coordinate imax, Coordinate jmax)
{
  struct Site *newsite, *bot, *top, *temp, *p;
  struct Site *v;
  struct Point newintstar;
  int pm;
  struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
  struct Edge *e;

  freeinit(&sfl, sizeof *sites);
  siteidx = 0;
  geominit();
  PQinitialize();
  bottomsite = nextsite();//(*nextsite)();
  ELinitialize();

  newsite = nextsite();//(*nextsite)();

  while(1) {
    if(!PQempty()) newintstar = PQ_min();

    if (newsite != (struct Site *)NULL 
	&& (PQempty() 
	    || newsite -> coord.y < newintstar.y
	    || (newsite->coord.y == newintstar.y 
		&& newsite->coord.x < newintstar.x))) {
      /* new site is smallest */
      lbnd = ELleftbnd(&(newsite->coord));
      rbnd = ELright(lbnd);
      bot = rightreg(lbnd);
      e = bisect(bot, newsite);
      bisector = HEcreate(e, LE);
      ELinsert(lbnd, bisector);
      if ((p = intersect(lbnd, bisector)) != (struct Site *) NULL) {
	PQdelete(lbnd);
	PQinsert(lbnd, p, dist(p,newsite));
      }
      lbnd = bisector;
      bisector = HEcreate(e, RE);
      ELinsert(lbnd, bisector);
      if ((p = intersect(bisector, rbnd)) != (struct Site *) NULL) {
	PQinsert(bisector, p, dist(p,newsite));	
      }
      newsite = nextsite();//(*nextsite)();	
    }
    else if (!PQempty()) {
      /* intersection is smallest */
      lbnd = PQextractmin();
      llbnd = ELleft(lbnd);
      rbnd = ELright(lbnd);
      rrbnd = ELright(rbnd);
      bot = leftreg(lbnd);
      top = rightreg(rbnd);
	    
      v = lbnd->vertex;
      makevertex(v);
      endpoint(lbnd->ELedge,lbnd->ELpm,v,imax,jmax);
      endpoint(rbnd->ELedge,rbnd->ELpm,v,imax,jmax);
      ELdelete(lbnd); 
      PQdelete(rbnd);
      ELdelete(rbnd); 
      pm = LE;
      if (bot->coord.y > top->coord.y) {
	temp = bot;
	bot = top;
	top = temp;
	pm = RE;
      }
      e = bisect(bot, top);
      bisector = HEcreate(e, pm);
      ELinsert(llbnd, bisector);
      endpoint(e, RE-pm, v, imax, jmax);
      deref(v);
      if((p = intersect(llbnd, bisector)) != (struct Site *) NULL){
	PQdelete(llbnd);
	PQinsert(llbnd, p, dist(p,bot));
      }
      if ((p = intersect(bisector, rrbnd)) != (struct Site *) NULL){
	PQinsert(bisector, p, dist(p,bot));
      }
    }
    else break;
  }

  for(lbnd=ELright(ELleftend); lbnd != ELrightend; lbnd=ELright(lbnd)) {
    e = lbnd -> ELedge;
    out_ep2(e,v,imax,jmax); /* Voronoi �x������ */
  }
}

//=============================================================================
void MkVoronoi::geominit()
{
  struct Edge e;
  float sn;

  freeinit(&efl, sizeof e);
  nvertices = 0;
  nedges = 0;
  sn = nsites+4;
  sqrt_nsites = sqrt(sn);
  deltay = ymax - ymin;
  deltax = xmax - xmin;
}


struct Edge *MkVoronoi::bisect(struct Site *s1, struct Site *s2)
{
  float dx,dy,adx,ady;
  struct Edge *newedge;
    
  newedge = (struct Edge *) getfree(&efl);

  newedge -> reg[0] = s1;
  newedge -> reg[1] = s2;
  ref(s1); 
  ref(s2);
  newedge -> ep[0] = (struct Site *) NULL;
  newedge -> ep[1] = (struct Site *) NULL;

  dx = s2->coord.x - s1->coord.x;
  dy = s2->coord.y - s1->coord.y;
  adx = dx>0 ? dx : -dx;
  ady = dy>0 ? dy : -dy;
  newedge -> c = s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5;
  if (adx>ady) {
    newedge -> a = 1.0;
    newedge -> b = dy/dx;
    newedge -> c /= dx;
  }
  else {
    newedge -> b = 1.0;
    newedge -> a = dx/dy;
    newedge -> c /= dy;
  }
  newedge -> edgenbr = nedges;
  nedges += 1;

  /* voronoi�x�ˤ��΁I�Ȥ��B�Y�ɷ֤Υ�٥�򴢴椹�� */
  if(s1->label < s2 -> label) {
    newedge -> lab1 = s1 -> label;
    newedge -> lab2 = s2 -> label;
  }
  else {
    newedge -> lab1 = s2 -> label;
    newedge -> lab2 = s1 -> label;
  }
    
  /* ĸ���g�ξ��x�򴢴椹�� */
  newedge -> dist =
    sqrt(pow((double)dx,(double)2) + pow((double)dy,(double)2));
                
  return(newedge);
}

struct Site *MkVoronoi::intersect(struct Halfedge *el1, struct Halfedge *el2)
{
  struct Edge *e1,*e2, *e;
  struct Halfedge *el;
  float d, xint, yint;
  int right_of_site;
  struct Site *v;

  e1 = el1 -> ELedge;
  e2 = el2 -> ELedge;
  if(e1 == (struct Edge*)NULL || e2 == (struct Edge*)NULL) 
    return ((struct Site *) NULL);
  if (e1->reg[1] == e2->reg[1]) return ((struct Site *) NULL);

  d = e1->a * e2->b - e1->b * e2->a;
  if (-1.0e-10<d && d<1.0e-10) return ((struct Site *) NULL);

  xint = (e1->c*e2->b - e2->c*e1->b)/d;
  yint = (e2->c*e1->a - e1->c*e2->a)/d;

  if( (e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
      (e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
       e1->reg[1]->coord.x < e2->reg[1]->coord.x) ) {
    el = el1;
    e = e1;
  }
  else {
    el = el2;
    e = e2;
  }
  right_of_site = xint >= e -> reg[1] -> coord.x;
  if ((right_of_site && el -> ELpm == LE) ||
      (!right_of_site && el -> ELpm == RE)) return ((struct Site *) NULL);

  v = (struct Site *) getfree(&sfl);
  v -> refcnt = 0;
  v -> coord.x = xint;
  v -> coord.y = yint;
  return(v);
}

/* returns 1 if p is to right of halfedge e */
int MkVoronoi::right_of(struct Halfedge *el, struct Point *p)
{
  struct Edge *e;
  struct Site *topsite;
  int right_of_site, above, fast;
  float dxp, dyp, dxs, t1, t2, t3, yl;

  e = el -> ELedge;
  topsite = e -> reg[1];
  right_of_site = p -> x > topsite -> coord.x;
  if(right_of_site && el -> ELpm == LE) return(1);
  if(!right_of_site && el -> ELpm == RE) return (0);

  if (e->a == 1.0) {
    dyp = p->y - topsite->coord.y;
    dxp = p->x - topsite->coord.x;
    fast = 0;
    if ((!right_of_site &e->b<0.0) | (right_of_site&e->b>=0.0) ) {
      above = dyp>= e->b*dxp;	
      fast = above;
    }
    else {
      above = p->x + p->y*e->b > e-> c;
      if(e->b<0.0) above = !above;
      if (!above) fast = 1;
    }
    if (!fast) {
      dxs = topsite->coord.x - (e->reg[0])->coord.x;
      above = e->b * (dxp*dxp - dyp*dyp) <
	dxs*dyp*(1.0+2.0*dxp/dxs + e->b*e->b);
      if(e->b<0.0) above = !above;
    }
  }
  else { /*e->b==1.0 */
    yl = e->c - e->a*p->x;
    t1 = p->y - yl;
    t2 = p->x - topsite->coord.x;
    t3 = yl - topsite->coord.y;
    above = t1*t1 > t2*t2 + t3*t3;
  }
  return (el->ELpm==LE ? above : !above);
}

void MkVoronoi::endpoint(struct Edge *e, int lr, struct Site *s,
	      Coordinate imax, Coordinate jmax)
{
  e -> ep[lr] = s;
  ref(s);
  if(e -> ep[RE-lr]== (struct Site *) NULL) return;
  out_ep2(e,s,imax,jmax);	/* voronoi�x���� */

  deref(e->reg[LE]);
  deref(e->reg[RE]);
  makefree((Freenode *)e, &efl); //modified 
}


float MkVoronoi::dist(struct Site *s, struct Site *t)
{
  float dx,dy;

  dx = s->coord.x - t->coord.x;
  dy = s->coord.y - t->coord.y;
  return(sqrt(dx*dx + dy*dy));
}

void MkVoronoi::makevertex(struct Site *v)
{
  v -> sitenbr = nvertices;
  nvertices += 1;
}

void MkVoronoi::deref(struct Site *v)
{
  v -> refcnt -= 1;
  if (v -> refcnt == 0 ) makefree((Freenode *)v, &sfl);//modified 
}

void MkVoronoi::ref(struct Site *v)
{
  v -> refcnt += 1;
}

//================================================
void MkVoronoi::ELinitialize()
{
  int i;
  freeinit(&hfl, sizeof **ELhash);
  ELhashsize = 2 * sqrt_nsites;
  ELhash = (struct Halfedge **) myalloc ( sizeof *ELhash * ELhashsize);
  for(i=0; i<ELhashsize; i +=1) ELhash[i] = (struct Halfedge *)NULL;
  ELleftend = HEcreate( (struct Edge *)NULL, 0);
  ELrightend = HEcreate( (struct Edge *)NULL, 0);
  ELleftend -> ELleft = (struct Halfedge *)NULL;
  ELleftend -> ELright = ELrightend;
  ELrightend -> ELleft = ELleftend;
  ELrightend -> ELright = (struct Halfedge *)NULL;
  ELhash[0] = ELleftend;
  ELhash[ELhashsize-1] = ELrightend;
}


struct Halfedge *MkVoronoi::HEcreate(struct Edge *e, int pm)
{
  struct Halfedge *answer;
  answer = (struct Halfedge *) getfree(&hfl);
  answer -> ELedge = e;
  answer -> ELpm = pm;
  answer -> PQnext = (struct Halfedge *) NULL;
  answer -> vertex = (struct Site *) NULL;
  answer -> ELrefcnt = 0;
  return(answer);
}


void MkVoronoi::ELinsert(struct Halfedge *lb, struct Halfedge *newlb)
{
  newlb -> ELleft = lb;
  newlb -> ELright = lb -> ELright;
  (lb -> ELright) -> ELleft = newlb;
  lb -> ELright = newlb;
}

/* Get entry from hash table, pruning any deleted nodes */
struct Halfedge *MkVoronoi::ELgethash(int b)
{
  struct Halfedge *he;

  if(b<0 || b>=ELhashsize) return((struct Halfedge *) NULL);
  he = ELhash[b]; 
  if (he == (struct Halfedge *) NULL || 
      he -> ELedge != (struct Edge *) DELETED ) return (he);

  /* Hash table points to deleted half edge.  Patch as necessary. */
  ELhash[b] = (struct Halfedge *) NULL;
  if ((he -> ELrefcnt -= 1) == 0) makefree((Freenode *)he, &hfl); // error C2664: 'makefree' : cannot convert parameter 1 from 'struct Halfedge * NOTE WTZ
  return ((struct Halfedge *) NULL);
}	

struct Halfedge *MkVoronoi::ELleftbnd(struct Point *p)
{
  int i, bucket;
  struct Halfedge *he;

  /* Use hash table to get close to desired halfedge */
  bucket = (p->x - xmin)/deltax * ELhashsize;
  if(bucket<0) bucket =0;
  if(bucket>=ELhashsize) bucket = ELhashsize - 1;
  he = ELgethash(bucket);
  if(he == (struct Halfedge *) NULL)
    {   for(i=1; 1 ; i += 1)
      {	if ((he=ELgethash(bucket-i)) != (struct Halfedge *) NULL) break;
      if ((he=ELgethash(bucket+i)) != (struct Halfedge *) NULL) break;
      };
    totalsearch += i;
    };
  ntry += 1;
  /* Now search linear list of halfedges for the corect one */
  if (he==ELleftend  || (he != ELrightend && right_of(he,p)))
    {do {he = he -> ELright;} while (he!=ELrightend && right_of(he,p));
    he = he -> ELleft;
    }
  else 
    do {he = he -> ELleft;} while (he!=ELleftend && !right_of(he,p));

  /* Update hash table and reference counts */
  if(bucket > 0 && bucket <ELhashsize-1)
    {	if(ELhash[bucket] != (struct Halfedge *) NULL) 
      ELhash[bucket] -> ELrefcnt -= 1;
    ELhash[bucket] = he;
    ELhash[bucket] -> ELrefcnt += 1;
    };
  return (he);
}

	
/* This delete routine can't reclaim node, since pointers from hash
   table may be present.   */
void MkVoronoi::ELdelete(struct Halfedge *he)
{
  (he -> ELleft) -> ELright = he -> ELright;
  (he -> ELright) -> ELleft = he -> ELleft;
  he -> ELedge = (struct Edge *)DELETED;
}


struct Halfedge *MkVoronoi::ELright(struct Halfedge *he)
{
  return (he -> ELright);
}

struct Halfedge *MkVoronoi::ELleft(struct Halfedge *he)
{
  return (he -> ELleft);
}


struct Site *MkVoronoi::leftreg(struct Halfedge *he)
{
  if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
  return( he -> ELpm == LE ? 
	  he -> ELedge -> reg[LE] : he -> ELedge -> reg[RE]);
}

struct Site *MkVoronoi::rightreg(struct Halfedge *he)
{
  if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
  return( he -> ELpm == LE ? 
	  he -> ELedge -> reg[RE] : he -> ELedge -> reg[LE]);
}

void MkVoronoi::freeinit(struct Freelist *fl, int size)
{
	fl -> head = (struct Freenode *) NULL;
	fl -> nodesize = size;
}

char *MkVoronoi::getfree(struct Freelist *fl)
{
	int i; struct Freenode *t;
	if(fl->head == (struct Freenode *) NULL) {
		t =  (struct Freenode *) myalloc(sqrt_nsites * fl->nodesize);
		for(i=0; i<sqrt_nsites; i+=1) 	
			makefree((struct Freenode *)((char *)t+i*fl->nodesize), fl);
	}
	t = fl -> head;
	fl -> head = (fl -> head) -> nextfree;
	return((char *)t);
}

void MkVoronoi::makefree(struct Freenode *curr, struct Freelist *fl)
{
	curr -> nextfree = fl -> head;
	fl -> head = curr;
}

char *MkVoronoi::myalloc(unsigned n)
{
	char *t;
	if ((t= (char*)malloc((size_t) n)) == (char *) 0) {
		fprintf(stderr,
			"Insufficient memory (%d bytes in use)\n",
			total_alloc);
		exit(0);
	}
	total_alloc += n;
	return(t);
}

char *MkVoronoi::myrealloc(void *ptr, unsigned current, unsigned inc, size_t unit)
{
	char *t;
	if ((t= (char*)realloc(ptr,(current+inc)*unit)) == (char *) 0) {
		fprintf(stderr,
			"Insufficient memory (%d bytes in use)\n",
			total_alloc);
		exit(0);
	}
	total_alloc += inc;
	return(t);
}


//--------------------------------------------------------------
/* float pxmin, pxmax, pymin, pymax, cradius; */

/*
 * (x,y)�����ˤ�����Ǥ����, ���ڤ����ˤˉ�Q�����v��
 * �¤���voronoi��(x,y)�η��Ť�FRAME(-1)�ˤ���.
 */
void MkVoronoi::in_frame(float *x, float *y, float d,
	      struct Edge *e, int lr,
	      Coordinate max_x, Coordinate max_y)
{
  /* x �����ˤ�0 ���¤Έ��� */
  if(*x < 0){
    *y += (-*x)*d;
    *x = 0;
    e->ep[lr]->sitenbr = FRAME;
    if(*y < 0){
      *x += (-*y)/d;
      *y = 0;
    }
    else if(*y > max_y){
      *x += ((float)max_y - *y)/d;
      *y = (float)max_y;
    }
  }
  /* x �����ˤ�max_x ���ϤΈ��� */
  else if(*x > max_x){
    *y += ((float)max_x - *x)*d;
    *x = (float)max_x;
    e->ep[lr]->sitenbr = FRAME;
    if(*y < 0){
      *x += (-*y)/d;
      *y = 0;
    }
    else if(*y > max_y){
      *x += ((float)max_y - *y)/d;
      *y = (float)max_y;
    }
  }
  /* y �����ˤ�0 ���¤Έ��� */
  else if(*y < 0){
    *x += (-*y)/d;
    *y = 0;
    e->ep[lr]->sitenbr = FRAME;
  }
  /* y �����ˤ�max_y ���ϤΈ��� */
  else if(*y > max_y){
    *x += ((float)max_y - *y)/d;
    *y = (float)max_y;
    e->ep[lr]->sitenbr = FRAME;
  }
}

/* �o���h��Ǥ���ʼ����ڤ����ˤˉ�Q�����v�� */
void MkVoronoi::s_in_frame(float *xsf, float *xef, float *ysf, float *yef,
		struct Edge *e, Coordinate max_y)
{
  float d;
    
  if((e->b) == 0) {	/* �A�����ޤΈ��� */
    *xsf = *xef;
    *ysf = (float)max_y;
  }
  else {
    d = -(e->a)/(e->b); /* �A�� */
	    
    if(d == 0){	/* �A����0 �Έ��� */
      *xsf = 0;
      *ysf = *yef;
    }
    else if(d > 0){	/* �A�������Έ��� */
      if(*yef > *xef*d){
	*xsf = 0;
	*ysf = *yef - (*xef*d);
      }
      else {
	*xsf = *xef - (*yef/d);
	*ysf = 0;
      }
    }
    else {		/* �A����ؓ�Έ��� */
      if(((float)max_y - *yef) > ((-*xef) * d)){
	*xsf = 0;
	*ysf = *yef + ((-*xef) * d);
      }
      else {
	*xsf = *xef - (-((float)max_y - *yef)) / d;
	*ysf = (float)max_y;
      }
    }
  }
}

/* �o���h��Ǥ���K����ڤ����ˤˉ�Q�����v�� */
void MkVoronoi::e_in_frame(float *xsf, float *xef, float *ysf, float *yef,
		struct Edge *e, Coordinate max_x, Coordinate max_y)
{
  float d;

  if((e->b) ==0) {	/* �A�����ޤΈ��� */
    *xef = *xsf;
    *yef = 0;
  }
  else {
    d = -(e->a)/(e->b); /* �A�� */

    if(d == 0){	/* �A����0 �Έ��� */
      *xef = (float)max_x;
      *yef = *ysf;
    }
    else if(d > 0){	/* �A�������Έ��� */
      if(((float)max_y - *ysf) > (((float)max_x - *xsf) * d)){
	*xef = (float)max_x;
	*yef = *ysf + ((float)max_x - *xsf) * d;
      }
      else {
	*xef = *xsf + ((float)max_y - *ysf) / d;
	*yef = (float)max_y;
      }
    }
    else {		/* �A����ؓ�Έ��� */
      if(*ysf > (- ((float)max_x - *xsf) * d)){
	*xef = (float)max_x;
	*yef = *ysf - (- ((float)max_x - *xsf) * d);
      }
      else {
	*xef = *xsf + (-*ysf) / d;
	*yef =0;
      }
    }
  }
}

/*
 * ʼ��, �K������ˤ��o���h��, �ޤ���, ����Έ���,
 * ���ˤ�������, ����voronoi��η��Ť�FRAME(-1) �ˤ����v��
 */
void MkVoronoi::frameout(float *xsf, float *xef, float *ysf, float *yef,
	      int *sp, int *ep, struct Edge *e,
	      Coordinate max_x, Coordinate max_y)
{
  float d;

  /* ʼ��, �K��Ȥ�˟o���h��Ǥʤ����� */
  if((e->ep[LE] != (struct Site *)NULL) &&
     (e->ep[RE] != (struct Site *)NULL)){
    *xsf = e->ep[LE]->coord.x; /* ʼ���x ���� */
    *xef = e->ep[RE]->coord.x; /* �K���x ���� */
    *ysf = e->ep[LE]->coord.y; /* ʼ���y ���� */
    *yef = e->ep[RE]->coord.y; /* �K���y ���� */

    /* �A������ �ΤȤ� */
    if(*xsf == *xef) {
      if(*ysf < 0) {
	*ysf =0;
	e->ep[LE]->sitenbr = FRAME;
      }
      else if(*ysf > max_y) {
	*ysf = (float)max_y;
	e->ep[LE]->sitenbr = FRAME;
      }
      if(*yef < 0) {
	*yef = 0;
	e->ep[RE]->sitenbr = FRAME;
      }
      else if(*yef > max_y) {
	*yef = (float)max_y;
	e->ep[RE]->sitenbr = FRAME;
      }
    }

    /* �A�������ޤΤȤ� */
    else {
      d = (*yef - *ysf)/(*xef - *xsf); /* �A�� */

      /* ʼ�㤬����Ζ���Έ��� */
      in_frame(xsf,ysf,d,e,LE,max_x,max_y);
      /* �K�㤬����Ζ���Έ��� */	    
      in_frame(xef,yef,d,e,RE,max_x,max_y);
    }
	
    *sp = e->ep[LE]->sitenbr; /* ʼ��η��Ť���� */
    *ep = e->ep[RE]->sitenbr; /* �K��η��Ť���� */
	
  }
    
  /* ʼ�㤬�o���h��Έ��� */
  else if((e->ep[LE] == (struct Site *)NULL)) {
    *xef = e->ep[RE]->coord.x; /* �K���x ���� */
    *yef = e->ep[RE]->coord.y; /* �K���y ���� */
    *sp = FRAME;	/* ʼ��η���FRAME ����� */
    *ep = e->ep[RE]->sitenbr; /* �K��η��Ť���� */

    /* ʼ����ڤˤ��� */
    s_in_frame(xsf,xef,ysf,yef,e,max_y);
  }
        
  /* �K�㤬�o���h��Έ��� */
  else if((e->ep[RE] == (struct Site *)NULL)){
    *xsf = e->ep[LE]->coord.x; /* ʼ���x ���� */
    *ysf = e->ep[LE]->coord.y; /* ʼ���y ���� */
    *sp = e->ep[LE]->sitenbr; /* ʼ��η��Ť���� */
    *ep = FRAME;	/* �K��η���FRAME ����� */

    /* �K����ڤˤ��� */
    e_in_frame(xsf,xef,ysf,yef,e,max_x,max_y);
  }
}

/*
 * �B�Y�ɷ��g��voronoi�x�Τߤ�lineseg �˴��椷,
 * �B�Y�ɷ��g���v�Sneighbor ��Ĥ����v��.
 */
void MkVoronoi::out_ep2(struct Edge *e, struct Site *v,
	     Coordinate imax, Coordinate jmax)
{
  int i,sp,ep;
  float xsf,xef,ysf,yef;
//  float si,sj,ei,ej;
  Coordinate max_x=imax-1;
  Coordinate max_y=jmax-1;

  static unsigned int current_neighbor_size = INITNEIGHBOR;
  static unsigned int current_lineseg_size  = INITLINE;

  /* double i1,j1,i2,j2; */



  /* ��٥뤬ͬ�����Ϥϳ������ʤ� */
  //if(output_pvor == NO && e->lab1 == e->lab2) 
  if(  e->lab1 == e->lab2) 
  {
    point_edge++;
    return;
  }
  else {
    /*
     * voronoi�㤬�o���h��, �ޤ���, ���⤫�ɤ����ж���,
     * �⤷, �����ʤ閘�ڤ����ˤ�����
     */
    frameout(&xsf,&xef,&ysf,&yef,&sp,&ep,e,imax-1,jmax-1);
	
    /* �Ϥ΄I��ˤ�äƤ�ʼ��, �K���һ��������Έ��ϳ������ʤ� */
    if((xsf < 0.0) || (xsf > (float)max_x) ||
       (xef < 0.0) || (xef > (float)max_x) ||
       (ysf < 0.0) || (ysf > (float)max_y) ||
       (yef < 0.0) || (yef > (float)max_y))
      return;
  }

// ���ꥸ�ʥ���`���Ȥ���
// һ�ȥ����ȥ����Ȥ��Ƥߤ�
//
//  if((e->ep[LE] != (struct Site *)NULL) &&
//     (e->ep[RE] != (struct Site *)NULL)){
//    /* voronoi�x��ʼ��?�K��Ȥ�˟o���h��Ǥʤ����� */
//		
//    si = e->ep[LE]->coord.x; /* ʼ���x ���� */
//    ei = e->ep[RE]->coord.x; /* �K���x ���� */
//    sj = e->ep[LE]->coord.y; /* ʼ���y ���� */
//    ej = e->ep[RE]->coord.y; /* �K���y ���� */
//		
//    if(!((si < (float)imax)&&(ei < (float)imax)&&
//	 (sj < (float)jmax)&&(ej < (float)jmax)
//	 &&(si > 0)&&(ei > 0)&&(sj > 0)&&(ej > 0))){
//      /* �⤷ʼ��?�K��Ȥ�˻�����ˤ��ä��� */
//      return;
//      /* neighbor �Υ��åȤ��Ф�ʤ� */
//    }
//  }
//  else{
//    /* voronoi�x��ʼ��?�K��Τɤ��餫���o���h��Ǥ������ */
//    return;
//    /* neighbor �Υ��åȤ��Ф�ʤ� */
//  }
	
  /*
   * �B�Y�ɷ֤��O���v�Sneighbor �˼Ȥ˵��h���Ƥ����,
   * ���h����Ƥ�����x��, ��ؤ�voronoi�x�򘋳ɤ���ĸ���g��
   * ���x�δ󤭤�����^��,�B�Y�ɷ��g�ξ��x���ٶ��x����
   */

  i = search(e->lab1,e->lab2); /* �ϥå����˵��h����Ƥ��뤫�{�٤� */

  /* δ���h�Έ��� */

  if(i == NODATA){
    enter(e->lab1,e->lab2,NEIGHnbr); /* �ϥå����˵��h */
    /*
    i1 = (double)(component[e->lab1].xc);
    j1 = (double)(component[e->lab1].yc);
    i2 = (double)(component[e->lab2].xc);
    j2 = (double)(component[e->lab2].yc);
    */

    /* ĸ���g�ξ��x����� */
    neighbor[NEIGHnbr].dist = e->dist; /* �B�Y�ɷ��g����С���x */
    neighbor[NEIGHnbr].lab1 = e->lab1;
    neighbor[NEIGHnbr].lab2 = e->lab2;
	
    /* �����g�νǶȤ���� */
    /*	if(i1 == i2) {
	neighbor[NEIGHnbr].angle = -RIGHTANGLE;
	}
	else {
	angle = atan2((j2-j1),(i2-i1))*2*RIGHTANGLE/M_PI;
			
	if(angle > RIGHTANGLE){
	neighbor[NEIGHnbr].angle = (float)(angle-2*RIGHTANGLE);
	}
	else if(angle <= -RIGHTANGLE){
	neighbor[NEIGHnbr].angle = (float)(angle+2*RIGHTANGLE);
	}
	else {
	neighbor[NEIGHnbr].angle = (float)angle;
	}
	}
	*/
    /*	neighbor[NEIGHnbr].dist = 
	(float)sqrt((component[e->lab1].xc - component[e->lab2].xc) *
	            (component[e->lab1].xc - component[e->lab2].xc) +
		    (component[e->lab1].yc - component[e->lab2].yc) *
		    (component[e->lab1].yc - component[e->lab2].yc));*/
    /* �B�Y�ɷ֤������g�ξ��x */
    /*	dx = (double)(component[e->lab2].xc - component[e->lab1].xc);
	dy = (double)(component[e->lab2].yc - component[e->lab1].yc);
	neighbor[NEIGHnbr].angle = (float)atan2(dy,dx);*/
	
    NEIGHnbr++;
    if(NEIGHnbr >= current_neighbor_size) {
      neighbor=(Neighborhood *)myrealloc(neighbor,
					 current_neighbor_size,
					 INCNEIGHBOR,
					 sizeof(Neighborhood));
      current_neighbor_size+=INCNEIGHBOR;
    }
  }

  /* ���h�g�ߤΈ��� */
  else {
    if(neighbor[i].dist > e->dist) /* ���x����^���ƶ̤������B�Y */
      neighbor[i].dist = e->dist; /* �ɷ��g�ξ��x�Ȥ��� */
  }

  if(sp > SiteMax) SiteMax = sp;
  if(ep > SiteMax) SiteMax = ep;

  /* voronoi�x�����򴢴椹�� */
  lineseg[LINEnbr].sp = sp;
  lineseg[LINEnbr].ep = ep;
  lineseg[LINEnbr].xs = (unsigned int)(xsf+0.5);
  lineseg[LINEnbr].xe = (unsigned int)(xef+0.5);
  lineseg[LINEnbr].ys = (unsigned int)(ysf+0.5);
  lineseg[LINEnbr].ye = (unsigned int)(yef+0.5);
  lineseg[LINEnbr].lab1 = e->lab1;
  lineseg[LINEnbr].lab2 = e->lab2;
  lineseg[LINEnbr].yn = OUTPUT;
  LINEnbr++;
  point_edge++;
  if(LINEnbr >= current_lineseg_size) {
      lineseg=(LineSegment *)myrealloc(lineseg,
				       current_lineseg_size,
				       INCLINE,
				       sizeof(LineSegment));
      current_lineseg_size+=INCLINE;
  }
}
//-----------------------------------------------------
void MkVoronoi::PQinsert(struct Halfedge *he, struct Site *v, float offset)
{
	struct Halfedge *last, *next;
	
	he -> vertex = v;
	ref(v);
	he -> ystar = v -> coord.y + offset;
	last = &PQhash[PQbucket(he)];
	while ((next = last -> PQnext) != (struct Halfedge *) NULL &&
		(he -> ystar  > next -> ystar  ||
		(he -> ystar == next -> ystar && v -> coord.x >
		next->vertex->coord.x))) {
		last = next;};
		he -> PQnext = last -> PQnext; 
		last -> PQnext = he;
		PQcount += 1;
}

void MkVoronoi::PQdelete(struct Halfedge *he)
{
	struct Halfedge *last;
	
	if(he ->  vertex != (struct Site *) NULL) {
		last = &PQhash[PQbucket(he)];
		while (last -> PQnext != he) last = last -> PQnext;
		last -> PQnext = he -> PQnext;
		PQcount -= 1;
		deref(he -> vertex);
		he -> vertex = (struct Site *) NULL;
	}
}

int MkVoronoi::PQbucket(struct Halfedge *he)
{
	int bucket;
	
	bucket = (he->ystar - ymin)/deltay * PQhashsize;
	if (bucket<0) bucket = 0;
	if (bucket>=PQhashsize) bucket = PQhashsize-1 ;
	if (bucket < PQmin) PQmin = bucket;
	return(bucket);
}

int MkVoronoi::PQempty()
{
	return(PQcount==0);
}

struct Point MkVoronoi::PQ_min()
{
	struct Point answer;
	
	while(PQhash[PQmin].PQnext == (struct Halfedge *)NULL) {PQmin += 1;};
	answer.x = PQhash[PQmin].PQnext -> vertex -> coord.x;
	answer.y = PQhash[PQmin].PQnext -> ystar;
	return (answer);
}

struct Halfedge *MkVoronoi::PQextractmin()
{
	struct Halfedge *curr;
	
	curr = PQhash[PQmin].PQnext;
	PQhash[PQmin].PQnext = curr -> PQnext;
	PQcount -= 1;
	return(curr);
}

void MkVoronoi::PQinitialize()
{
	int i;
	
	PQcount = 0;
	PQmin = 0;
	PQhashsize = 4 * sqrt_nsites;
	PQhash = (struct Halfedge *) myalloc(PQhashsize * sizeof *PQhash);
	for(i=0; i<PQhashsize; i+=1) PQhash[i].PQnext = (struct Halfedge *)NULL;
}
//========================================

/* �ϥå����v�� 1 */
HashVal MkVoronoi::hash1(Key key)
{
  return(key % M1);
}

/* �ϥå����v�� 2 */
HashVal MkVoronoi::hash2(Key key)
{
  return(key % M2);
}

/* ���ڻ��v�� */
void MkVoronoi::init_hash()
{
  HashVal i;

  for(i=0;i<M1+M2;i++)
    hashtable[i]=NIL;
}

/*
 * �ϥå����v���ˌ����륭�`��Ĥ����v��
 * 2�ĤΥ�٥���ܤ��Ȥ�, ����ˌ�����id �򷵤�
 */
Key MkVoronoi::key(Label lab1, Label lab2)
{
  unsigned long key,tmp;

  key = lab1;
  key = key << 2*BYTE_WIDTH;
  tmp = lab2;
  key = key | tmp;
  return key;
}

/* id ���ϥå����˵��h����Ƥ��뤫���{�٤��v�� */ 
int MkVoronoi::search(Label lab1, Label lab2)
{
  Key id;
  HashVal x;
  HashTable *p;

  id = key(lab1,lab2);
  x = hash1(id)+hash2(id);	/* �ϥå��傎��Ӌ�� */
  p = hashtable[x];

  while(p != NIL) {
    if((lab1 == p->lab1) && (lab2 == p->lab2)) /* ���h����Ƥ���� */
      return(p->entry);	                       /* ����entry �΂��򷵤� */
    p = p->next;
  }
    
  return NODATA;	/* ���h����Ƥ��ʤ��Ȥ���NODATA �΂��򷵤� */
}

/*
 * ���h����Ƥ��ʤ�id �Ȥ���ˌ�����entry �΂���
 * �ϥå����˵��h�����v��
 */ 
void MkVoronoi::enter(Label lab1, Label lab2, unsigned int entry)
{
  Key id;
  HashVal x;
  HashTable *p;

  id = key(lab1,lab2);
  x = hash1(id)+hash2(id);	/* �ϥå��傎��Ӌ�� */
    
  /* ���h���뤿����I���_������ */
  p = (HashTable *)myalloc(sizeof(HashTable)* 1);

  /* �_�������I����뤷, ������h���� */
  p->next = hashtable[x];
  hashtable[x] = p;
  hashtable[x]->lab1 = lab1;
  hashtable[x]->lab2 = lab2;
  hashtable[x]->entry = entry;
}

#endif