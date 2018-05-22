#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <math.h>
#include "c24bitmap.h"
#include "c256bitmap.h"
#include "mregion.h"
#include "PicMser.h"

void main(int argc,char* argv[])
{
  //Region R;
  C24BitMap CPic;
  CPic.Load(argv[1]);
  vector<Region> RegionB, RegionW, RegionR;
  GetMserRegion( CPic, RegionB, RegionW, RegionR);
}