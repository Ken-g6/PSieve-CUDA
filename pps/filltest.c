// Test the filling of the bitskip array, to multiply by 2^-b at once.
#include <stdio.h>
#include <stdint.h>
#define BITSATATIME 5
#define BITSMAX (1<<BITSATATIME)
#define MYP 10007
int bitsatatime = BITSATATIME;

// Fill the bitskip array, to multiply by 2^-b at once.
// O(2^bitsatatime) performance, only 2^bitsatatime writes.
void fillbitskip(uint64_t *bitskip, uint64_t p) {
	int len = 1<<bitsatatime;
	int halflen=len/2; 
	int j, k; 

	// Initialize the first two entries.
	bitskip[halflen] = (p+1)/2;	// Needed first.
	bitskip[0] = 0;			// Ignored through the end.

	// Fill in the intervening spaces, two numbers at a time.
	for(j=halflen; j > 1; j >>= 1) {
		for(k=j/2; k < halflen; k+=j) {
			register uint64_t bl = bitskip[2*k];
			//printf("Filling k=%d from bitskip=%lu\n", k, bl);
			bitskip[k] = (bl+((bl&1)?p:(uint64_t)0))/2;
			//printf("Filling k=%d\n", k+halflen);
			bitskip[k+halflen] = (bl+1+((bl&1)?(uint64_t)0:p))/2;
		}
	}
}
int main(void) {
	uint64_t bitskip[BITSMAX];
	int i;

	fillbitskip(bitskip, MYP);
	for(i=0; i < BITSMAX; i++) {
		printf("%d: %lu\n", i, bitskip[i]);
	}
	return 0;
}
