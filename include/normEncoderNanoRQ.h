#ifndef _NORM_ENCODER_NANORQ
#define _NORM_ENCODER_NANORQ

#include "normEncoder.h"
#include "protoDefs.h"  // for UINT16
#ifdef __cplusplus
extern "C" {
#endif
#include "nanorq_core.h"
#include "nanorq_ops.h"
#ifdef __cplusplus
}
#endif

class NormEncoderNanoRQ : public NormEncoder
{
    public:
        NormEncoderNanoRQ();
        virtual ~NormEncoderNanoRQ();
        
        virtual bool Init(unsigned int numData, unsigned int numParity, UINT16 vectorSize);
        virtual void Destroy();
        virtual void Encode(unsigned int segmentId, const char* dataVector, char** parityVectorList);
        virtual void EncodeParity(unsigned int parityId, const char** sourceVectorList, unsigned int numData, char* parityVector);
        
        virtual bool IsRateless() const { return true; }
        
        unsigned int GetNumData() {return ndata;}
        unsigned int GetNumParity() {return npar;}
        unsigned int GetVectorSize() {return vector_size;}
	
    private:
        unsigned int    ndata;        // max data pkts per block (k)
        unsigned int    npar;         // No. of parity packets (n-k)
        unsigned int    vector_size;  // Size of biggest vector to encode
        nanorq_core     rq;
        bool            rq_initialized;
        uint8_t*        prep_mem;
        uint8_t*        work_mem;
        uint8_t*        D;
        schedule        S;
        uint32_t        stride;
        unsigned int    encoded_count;
};


class NormDecoderNanoRQ : public NormDecoder
{
    public:
        NormDecoderNanoRQ();
        virtual ~NormDecoderNanoRQ();
        virtual bool Init(unsigned int numData, unsigned int numParity, UINT16 vectorSize);
        virtual void Destroy();
        virtual int Decode(char** vectorList, unsigned int numData,  unsigned int erasureCount, unsigned int* erasureLocs);
        
        unsigned int GetNumParity() {return npar;}
        unsigned int GetVectorSize() {return vector_size;}
        
    private:
        unsigned int    ndata;        // max data pkts per block (k)
        unsigned int    npar;         // No. of parity packets (n-k)
        UINT16          vector_size;  // Size of biggest vector to encode
        nanorq_core     rq;
        bool            rq_initialized;
};

#ifdef __cplusplus
extern "C" {
#endif

NormEncoder* NormCreateNanoRQEncoder();
NormDecoder* NormCreateNanoRQDecoder();

#ifdef __cplusplus
}
#endif

#endif // _NORM_ENCODER_NANORQ
