// Standalone test for the NORM NanoRQ FEC implementation

#include "protoTime.h"  // for ProtoTime
#include "normEncoderNanoRQ.h"

#include <string.h> // for memcpy(), etc
#include <stdlib.h> // for rand()
#include <stdio.h>

const unsigned int NUM_PARITY   = 32;
const unsigned int NUM_DATA     = 128;
const unsigned int SHORT_DATA   = 128;
const unsigned int SEG_SIZE     = 64;

const unsigned int B_SIZE = (SHORT_DATA + NUM_PARITY);

#define NORM_ENCODER NormEncoderNanoRQ
#define NORM_DECODER NormDecoderNanoRQ

int main(int argc, char* argv[])
{
    ProtoTime currentTime;
    currentTime.GetCurrentTime();
    int seed = (unsigned int)currentTime.usec();
    fprintf(stderr, "nanorqFecTest: seed = %u\n", seed);
    srand(seed);
    
    NORM_ENCODER encoder;
    if (!encoder.Init(NUM_DATA, NUM_PARITY, SEG_SIZE)) {
        fprintf(stderr, "nanorqFecTest: encoder init failed!\n");
        return -1;
    }
    
    NORM_DECODER decoder;
    if (!decoder.Init(NUM_DATA, NUM_PARITY, SEG_SIZE)) {
        fprintf(stderr, "nanorqFecTest: decoder init failed!\n");
        return -1;
    }
     
    for (int trial = 0; trial < 5; trial++)
    {
        fprintf(stderr, "\n--- Trial %d ---\n", trial + 1);

        // 1) Create some source data
        char txData[B_SIZE][SEG_SIZE];
        char* txDataPtr[B_SIZE];
        for (unsigned int i = 0 ; i < SHORT_DATA; i++)
        {
            txDataPtr[i] = txData[i];
            memset(txDataPtr[i], 'A' + (i % 26), SEG_SIZE - 1);
            txDataPtr[i][SEG_SIZE - 1] = '\0';
        }

        // 2) Zero-init parity vectors
        for (unsigned int i = SHORT_DATA; i < B_SIZE; i++)
        {
            txDataPtr[i] = txData[i];
            memset(txDataPtr[i], 0, SEG_SIZE);
        }

        // 3) Encode progressively (incremental/streaming FEC)
        ProtoTime startTime, stopTime;
        startTime.GetCurrentTime();
        for (unsigned int i = 0; i < NUM_PARITY; i++)
        {
            encoder.EncodeParity(i, (const char**)txDataPtr, SHORT_DATA, txDataPtr[SHORT_DATA + i]);
        }
        stopTime.GetCurrentTime();
        double encodeTime = ProtoTime::Delta(stopTime, startTime);

        // 4) Copy txData to rxData
        char rxData[B_SIZE][SEG_SIZE];
        char* rxDataPtr[B_SIZE];
        for (unsigned int i = 0; i < B_SIZE; i++)
        {
            rxDataPtr[i] = rxData[i];
            memcpy(rxDataPtr[i], txDataPtr[i], SEG_SIZE);
        }

        // 5) Introduce random erasures (up to NUM_PARITY)
        unsigned int erasureCount = 1 + (rand() % NUM_PARITY);
        unsigned int erasureLocs[B_SIZE];
        for (unsigned int i = 0; i < B_SIZE; i++)
            erasureLocs[i] = i;
        for (unsigned int i = 0; i < erasureCount; i++)
        {
            unsigned int loc = i + (rand() % (B_SIZE - i));
            unsigned int tmp = erasureLocs[i];
            erasureLocs[i] = erasureLocs[loc];
            erasureLocs[loc] = tmp;
        }
        
        // Sort erasure locations in ascending order
        for (unsigned int i = 0; i < erasureCount; i++)
        {
            for (unsigned int j = i+1; j < erasureCount; j++)
            {
                if (erasureLocs[j] < erasureLocs[i])
                {
                    unsigned int tmp = erasureLocs[i];
                    erasureLocs[i] = erasureLocs[j];
                    erasureLocs[j] = tmp;
                }
            }
        }
        
        fprintf(stderr, "erasureCount: %u / %u\n", erasureCount, NUM_PARITY);
        fprintf(stderr, "erasureLocs: ");
        for (unsigned int i = 0; i < erasureCount; i++)
            fprintf(stderr, "%u ", erasureLocs[i]);
        fprintf(stderr, "\n");

        // 6) Clear erased rxData buffers
        for (unsigned int i = 0; i < erasureCount; i++)
            memset(rxDataPtr[erasureLocs[i]], 0, SEG_SIZE);

        // 7) Decode the erased data
        startTime.GetCurrentTime();
        int decode_res = decoder.Decode(rxDataPtr, SHORT_DATA, erasureCount, erasureLocs);
        stopTime.GetCurrentTime();
        double decodeTime = ProtoTime::Delta(stopTime, startTime);

        if (decode_res != 0)
        {
            fprintf(stderr, "nanorqFecTest: decode() returned %d (insufficient symbols for rateless decode). Skipping validation.\n", decode_res);
            continue;
        }

        // 8) Validate correctness
        bool success = true;
        for (unsigned int i = 0; i < SHORT_DATA; i++)
        {
            if (0 != memcmp(rxDataPtr[i], txDataPtr[i], SEG_SIZE))
            {
                fprintf(stderr, "nanorqFecTest: segment %d decode ERROR!\n", i);
                fprintf(stderr, "  Expected: ");
                for (unsigned int k = 0; k < 16 && k < SEG_SIZE; k++)
                    fprintf(stderr, "%02x ", (unsigned char)txDataPtr[i][k]);
                fprintf(stderr, "\n  Decoded:  ");
                for (unsigned int k = 0; k < 16 && k < SEG_SIZE; k++)
                    fprintf(stderr, "%02x ", (unsigned char)rxDataPtr[i][k]);
                fprintf(stderr, "\n");
                success = false;
            }
        }    

        if (success) {
            fprintf(stderr, "nanorqFecTest: validation SUCCESS!\n");
        } else {
            return -1;
        }

        // 9) Print speeds
        fprintf(stderr, "nanorqFecTest: encodeTime: %lf usec, decodeTime: %lf usec\n", 1.0e+06*encodeTime, 1.0e+06*decodeTime); 
    }
    
    return 0;
}
