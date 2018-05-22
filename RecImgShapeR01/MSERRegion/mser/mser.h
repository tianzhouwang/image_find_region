 

#ifndef MSER_BASIC_PROCESS_H
#define MSER_BASIC_PROCESS_H
#include <algorithm>
#include <cassert>
#include <limits>

using namespace std;

#include "mregion.h"
#include <deque>
#include <vector>

using namespace std;

//================================================================
#define M_PI   3.1415926535897932384626433832795 
typedef signed char       int8_t;
typedef signed short      int16_t;
typedef signed int        int32_t;
typedef unsigned char     uint8_t;
typedef unsigned short    uint16_t;
typedef unsigned int      uint32_t;
//================================================================

//#include <stdint.h>

/// The MSER class extracts maximally stable extremal regions from a grayscale (8 bits) image.
/// @note The MSER class is not reentrant, so if you want to extract regions in parallel, each
/// thread needs to have its own MSER class instance.


class MSER
{
public:
	/// Constructor.
	/// @param[in] eight Use 8-connected pixels instead of 4-connected.
	/// @param[in] delta DELTA parameter of the MSER algorithm. Roughly speaking, the stability of a
	/// region is the relative variation of the region area when the intensity is changed by delta.
	/// @param[in] minArea Minimum area of any stable region relative to the image domain area.
	/// @param[in] maxArea Maximum area of any stable region relative to the image domain area.
	/// @param[in] maxVariation Maximum variation (absolute stability score) of the regions.
	/// @param[in] minDiversity Minimum diversity of the regions. When the relative area of two
	/// nested regions is below this threshold, then only the most stable one is selected.
	MSER(bool eight = false, int delta = 2, double minArea = 0.0001, double maxArea = 0.5,
	  	 double maxVariation = 0.5, double minDiversity = 0.33);
	
	/// Extracts maximally stable extremal regions from a grayscale (8 bits) image.
	/// @param[in] bits Pointer to the first scanline of the image.
	/// @param[in] width Width of the image.
	/// @param[in] height Height of the image.
	/// @param[out] regions Detected MSER.
	void GenerateMSER(const uint8_t * bits, int width, int height);//, std::vector<Region> & regions);
	Region*TopRegion;
	// Implementation details (could be moved outside this header file)
//private:
	// Helper method
	void processStack(int newPixelGreyLevel, int pixel, std::vector<Region *> & regionStack);
	
	// Parameters
	bool eight_;
	int delta_;
	double minArea_;
	double maxArea_;
	double maxVariation_;
	double minDiversity_;
	
	// Memory pool of regions for faster allocation
	std::deque<Region> pool_;        //存孩子
	unsigned int	poolIndex_;      //队列顶端的元素Index数值
	
};


MSER::MSER(bool eight, int delta, double minArea, double maxArea, double maxVariation,
		   double minDiversity) : eight_(eight), delta_(delta), minArea_(minArea),
 maxArea_(maxArea), maxVariation_(maxVariation), minDiversity_(minDiversity), pool_(256), poolIndex_(0)
{
	// Parameter check
	assert(delta   > 0   );
	assert(minArea >= 0.0);
	assert(maxArea <= 1.0);
	assert(minArea < maxArea);
	assert(maxVariation > 0.0);
	assert(minDiversity >= 0.0);
}

void MSER::GenerateMSER(const uint8_t * bits, int width, int height)//, vector<Region> & regions)
{
	Region *R;
	if (!bits || (width <= 0) || (height <= 0))
		return;
	
	// 1. Clear the accessible pixel mask, the heap of boundary pixels and the component stack. Push
	// a dummy-component onto the stack, with grey-level higher than any allowed in the image.
	vector<bool> accessible(width * height);
	vector<int>  boundaryPixels[256];
	int priority = 256;
	vector<Region *> regionStack;
	
	//Region *R =  new(&pool_[poolIndex_]) Region;
	
	regionStack.push_back( &pool_[poolIndex_++]);
	
	// 2. Make the source pixel (with its first edge) the current pixel, mark it as accessible and
	// store the grey-level of it in the variable current level.
	int curPixel = 0;
	int curEdge  = 0;
	int curLevel = bits[0];
	accessible[0] = true;
	
	// 3. Push an empty component with current level onto the component stack.
step_3:
    R         = &pool_[poolIndex_++];
	//R->r_id   = poolIndex_ - 1;
	R->level_ = curLevel; 
	R->pixel_ = curPixel;
	
	regionStack.push_back(R);
	
	if (poolIndex_ == pool_.size())
		pool_.resize(pool_.size() + 256);
	
	// 4. Explore the remaining edges to the neighbors of the current pixel, in order, as follows:
	// For each neighbor, check if the neighbor is already accessible. If it is not, mark it as
	// accessible and retrieve its grey-level. If the grey-level is not lower than the current one,
	// push it onto the heap of boundary pixels. If on the other hand the grey-level is lower than
	// the current one, enter the current pixel back into the queue of boundary pixels for later
	// processing (with the next edge number), consider the new pixel and its grey-level and go to 3.
	for (;;) {
		const int x = curPixel % width;
		const int y = curPixel / width;
		
		for (; curEdge < (eight_ ? 8 : 4); ++curEdge) {
			int neighborPixel = curPixel;
			
			if (eight_) {
				switch (curEdge) {
				case 0:
					if (x < width - 1)
						neighborPixel = curPixel + 1;
					break;
				case 1:
					if ((x < width - 1) && (y > 0))
						neighborPixel = curPixel - width + 1;
					break;
				case 2:
					if (y > 0)
						neighborPixel = curPixel - width;
					break;
				case 3:
					if ((x > 0) && (y > 0))
						neighborPixel = curPixel - width - 1;
					break;
				case 4:
					if (x > 0)
						neighborPixel = curPixel - 1;
					break;
				case 5:
					if ((x > 0) && (y < height - 1))
						neighborPixel = curPixel + width - 1;
					break;
				case 6:
					if (y < height - 1)
						neighborPixel = curPixel + width;
					break;
				default:
					if ((x < width - 1) && (y < height - 1))
						neighborPixel = curPixel + width + 1;
					break;
				}
			}
			else {
				switch (curEdge) {
				case 0:
					if (x < width - 1)
						neighborPixel = curPixel + 1;
					break;
				case 1:
					if (y < height - 1)
						neighborPixel = curPixel + width;
					break;
				case 2:
					if (x > 0)
						neighborPixel = curPixel - 1;
					break;
				default:
					if (y > 0)
						neighborPixel = curPixel - width;
					break;
				}
			}
			
			if (neighborPixel != curPixel && !accessible[neighborPixel]) {
				const int neighborLevel = bits[neighborPixel];
				accessible[neighborPixel] = true;
				
				if (neighborLevel >= curLevel) {//将高的像素设置成boundaryPixels
					boundaryPixels[neighborLevel].push_back(neighborPixel << 4);
					
					if (neighborLevel < priority)
						priority = neighborLevel;
				}
				else {
					boundaryPixels[curLevel].push_back((curPixel << 4) | (curEdge + 1));
					
					if (curLevel < priority)
						priority = curLevel;
					
					curPixel = neighborPixel;
					curEdge = 0;
					curLevel = neighborLevel;
					
					goto step_3;
				}
			}
		}
		
		// 5. Accumulate the current pixel to the component at the top of the stack (water
		// saturates the current pixel).
		regionStack.back()->accumulate(x, y);
		
		// 6. Pop the heap of boundary pixels. If the heap is empty, we are done. If the returned
		// pixel is at the same grey-level as the previous, go to 4.
		if (priority == 256) {
			TopRegion = regionStack.back();
			/*regionStack.back()->detect(delta_, minArea_ * width * height,
				maxArea_ * width * height, maxVariation_, minDiversity_,
				regions);
            printf("region stack size:%i\n",regions.size());
			printf("Idx: %i\n",poolIndex_);
			*/
			poolIndex_ = 0;
			return;
		}
		
		curPixel = boundaryPixels[priority].back() >> 4;
		curEdge  = boundaryPixels[priority].back() & 15;
		
		boundaryPixels[priority].pop_back();
		
		while (boundaryPixels[priority].empty() && (priority < 256))
			++priority;
		
		const int newPixelGreyLevel = bits[curPixel];
		
		if (newPixelGreyLevel != curLevel) {
			curLevel = newPixelGreyLevel;
			
			// 7. The returned pixel is at a higher grey-level, so we must now process
			// all components on the component stack until we reach the higher
			// grey-level. This is done with the processStack sub-routine, see below.
			// Then go to 4.
			processStack(newPixelGreyLevel, curPixel, regionStack);
			// printf("%i,%i\n",x,y);
		}
	}
	//printf("Idx%i\n",poolIndex_);
}

void MSER::processStack(int newPixelGreyLevel, int pixel, vector<Region *> & regionStack)
{
	// 1. Process component on the top of the stack. The next grey-level is the minimum of
	// newPixelGreyLevel and the grey-level for the second component on the stack.
	do {
		Region * top = regionStack.back();
		
		regionStack.pop_back();
		
		// 2. If newPixelGreyLevel is smaller than the grey-level on the second component on the
		// stack, set the top of stack grey-level to newPixelGreyLevel and return from sub-routine
		// (This occurs when the new pixel is at a grey-level for which there is not yet a component
		// instantiated, so we let the top of stack be that level by just changing its grey-level.
		if (newPixelGreyLevel < regionStack.back()->level_) {
			Region*R;
			R = &pool_[poolIndex_++];
			R->level_ = newPixelGreyLevel;R->pixel_ = pixel;
			regionStack.push_back(R);
			
			if (poolIndex_ == pool_.size())
				pool_.resize(pool_.size() + 256);
			
			regionStack.back()->merge(top);			
			return;
		}
		
		// 3. Remove the top of stack and merge it into the second component on stack as follows:
		// Add the first and second moment accumulators together and/or join the pixel lists.
		// Either merge the histories of the components, or take the history from the winner. Note
		// here that the top of stack should be considered one ’time-step’ back, so its current
		// size is part of the history. Therefore the top of stack would be the winner if its
		// current size is larger than the previous size of second on stack.
		regionStack.back()->merge(top);
	}
	// 4. If(newPixelGreyLevel>top of stack grey-level) go to 1.
	while (newPixelGreyLevel > regionStack.back()->level_);
}



#endif



