////////////////////////////////////////////////////////////////////////////////////////////////
// Description: This is the file for comuting time consuming of a procedure
// Usage: 
/*	   C256BitMap InPic, OutPic;
 	   InPic.Load("input.bmp");
 	   TimeCounter timer;
 	   ImgErode(InPic,OutPic);
 	   timer.gettime();
 	   OutPic.Save("output.bmp");*/

#ifndef TIME_COUNT_HEAD
#define TIME_COUNT_HEAD
#include <time.h>
#include <string.h>
 

class TimeCounter
{
 public:
 int tcnt;
 char txtbuffer[100];
 void settime(char*txt);
 void gettime();
 int  addTime();
 int  counters[200];
 int  countersTmp[200];
 string commentstr[200];
 TimeCounter();
 void AddCounterIdx(int idx);
};

//TimeCounter TCnt;

void TimeCounter::AddCounterIdx(int idx)
{
  counters[idx] += addTime();
}

TimeCounter::TimeCounter()
{
  tcnt = clock();
  strcpy(txtbuffer,"process time:");
  int i;
  for(i=0;i<200;i++)
  {
    counters[i]=0;
  }
}

void TimeCounter::settime(char*txt=NULL)
{
  tcnt = clock();
 if(txt!=NULL)
   strcpy(txtbuffer,txt);
 else
   strcpy(txtbuffer,"process time:");
}

void TimeCounter::gettime()
{
  // printf("%s %i\n",txtbuffer,clock()-tcnt);
   tcnt = clock();
}

int TimeCounter::addTime()
{
  int retVal = clock()-tcnt;
    tcnt = clock();
  return  retVal;
}

#endif
