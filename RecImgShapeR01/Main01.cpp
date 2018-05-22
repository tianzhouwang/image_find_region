#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "c24bitmap.h"
#include "c256bitmap.h"
 
int main(int argc, char*argv[])
{
 int i,j;
 C256BitMap GPic;
 //C24BitMap CPic;

 
	   GPic.Load("shapeRec.bmp");
	   GPic.Save("tmpA.bmp");
 

 //CPic.Save("dest.bmp");
 return 1;
}
