#pragma once

#include "miscStdHeaders.h"

#define DRAW_WIDTH 4
#define DEFAULT_PRNG_SEED0 0x2545f4914f6cdd1d
#define DEFAULT_PRNG_SEED1 (0x2545f4914f6cdd1d + 1595721457336359713)
#define DEFAULT_PRNG_SEED2 (0x2545f4914f6cdd1d*3)
#define DEFAULT_PRNG_SEED3 1595721457336359713
#define DEFAULT_TICK_COUNT 0
#define DEFAULT_TOTAL_MULTIPLIER 0

#pragma warning(push) //pop at end of file
#pragma warning(disable : 4244) //uint64_t to uint32_t conversion is intended
namespace AZ{
    //PERFORMANCE: on a phenon II x4 820, 2,8 Ghz, 8gb DDR2 1333Mhz ram (for large n):
    //x86 debug: ~29 - 39 (AZ test) nanos per PRN
    //x64 debug: ~13,5 - 20 (AZ test) nanos per PRN
    //x86 release: ~11,5 - 12,5 (with load) nanos per PRN
    //x64 release: ~2,5 - 2,8 (with load) nanos per PRN
    //(times to draw numbers, multiply by a float and assign to an array)
    //(about 3-3,5 times then draw1spcg32 on release x64 for test system)

	inline void draw4spcg32s(uint64_t* s0, uint64_t* s1, uint64_t* s2, uint64_t* s3,
                                                   uint32_t* dest0, uint32_t* dest1, 
                                                   uint32_t* dest2, uint32_t* dest3) { 
        uint64_t m = 0x9b60933458e17d7d;
        uint64_t a = 0xd737232eeccdf7ed;

        *s0 = *s0 * m + a;
        *s1 = *s1 * m + a;
        *s2 = *s2 * m + a;
        *s3 = *s3 * m + a;
        
        int shift0 = 29 - (*s0 >> 61);
        int shift1 = 29 - (*s1 >> 61);
        int shift2 = 29 - (*s2 >> 61);
        int shift3 = 29 - (*s3 >> 61); 

        *dest0 = *s0 >> shift0;
        *dest1 = *s1 >> shift1;
        *dest2 = *s2 >> shift2;
        *dest3 = *s3 >> shift3;
    }


    //from https://nullprogram.com/blog/2017/09/21/
    //on test system:
    //debug x64: ~16,5 - 26.5 nanos per prn
    //release x86: ~53 nanos per prn
    //debug x86: ~16,5 nanos per prn
    //release x64: ~9 nanos per prn
	inline uint32_t draw1spcg32(uint64_t* seed) {
        uint64_t m = 0x9b60933458e17d7d;
        uint64_t a = 0xd737232eeccdf7ed;
        *seed = *seed * m + a;
        int shift = 29 - (*seed >> 61);
        return *seed >> shift;
    }
}
#pragma warning(pop)


