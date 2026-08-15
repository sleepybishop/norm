#include <stdlib.h>
#include <string.h>
#include "normEncoderNanoRQ.h"
#include "protoDebug.h"
extern "C" {
#include "obl/oblas_lite.h"
}
#include "util.h" // for nanorq_oblas

NormEncoderNanoRQ::NormEncoderNanoRQ()
  : ndata(0), npar(0), vector_size(0), rq_initialized(false),
    prep_mem(NULL), work_mem(NULL), D(NULL), stride(0), encoded_count(0)
{
    S.ops.a = NULL;
}

NormEncoderNanoRQ::~NormEncoderNanoRQ()
{
    Destroy();
}

bool NormEncoderNanoRQ::Init(unsigned int numData, unsigned int numParity, UINT16 vecSizeMax)
{
    Destroy();
    
    ndata = numData;
    npar = numParity;
    vector_size = vecSizeMax;
    
    if (!nanorq_core_encoder_new(ndata, 0, &rq))
    {
        PLOG(PL_FATAL, "NormEncoderNanoRQ::Init() error: failed to create nanorq_core instance\n");
        return false;
    }
    rq_initialized = true;
    encoded_count = 0;
    
    size_t prep_len = nanorq_core_calculate_prepare_memory(&rq);
    prep_mem = new uint8_t[prep_len];
    if (!nanorq_core_prepare(&rq, prep_mem, prep_len))
    {
        PLOG(PL_FATAL, "NormEncoderNanoRQ::Init() error: nanorq_core_prepare failed\n");
        return false;
    }
    
    size_t work_len = nanorq_core_calculate_work_memory(&rq);
    work_mem = new uint8_t[work_len];
    
    size_t sched_bytes = ops_estimate_schedule_bytes(ndata);
    void* sched_buf = malloc(sched_bytes);
    if (!schedule_init(&S, sched_buf, sched_bytes))
    {
        PLOG(PL_FATAL, "NormEncoderNanoRQ::Init() error: schedule_init failed\n");
        free(sched_buf);
        return false;
    }
    
    nanorq_core_set_op_callback(&rq, &S, ops_push);
    if (!nanorq_core_precalculate(&rq, work_mem, work_len))
    {
        PLOG(PL_FATAL, "NormEncoderNanoRQ::Init() error: nanorq_core_precalculate failed\n");
        return false;
    }
    
    u32 rows = nanorq_core_get_pc_rows(&rq);
    stride = nanorq_core_recommended_stride(vector_size);
    D = (uint8_t*)obl_alloc(rows, stride, nanorq_oblas.align_size);
    if (!D)
    {
        PLOG(PL_FATAL, "NormEncoderNanoRQ::Init() error: D allocation failed\n");
        return false;
    }
    nanorq_core_init_matrix(&rq, D, stride);
    
    return true;
}

void NormEncoderNanoRQ::Destroy()
{
    if (D)
    {
        obl_free(D);
        D = NULL;
    }
    if (work_mem)
    {
        delete[] work_mem;
        work_mem = NULL;
    }
    if (prep_mem)
    {
        delete[] prep_mem;
        prep_mem = NULL;
    }
    if (S.ops.a)
    {
        free(S.ops.a);
        S.ops.a = NULL;
    }
    rq_initialized = false;
}

void NormEncoderNanoRQ::Encode(unsigned int segmentId, const char* dataVector, char** parityVectorList)
{
    ASSERT(rq_initialized);
    // NORM calls Encode() sequentially as source symbols arrive, but since we synthesize
    // repair symbols on demand in EncodeParity(), we actually do nothing here.
    // The fountain code only needs the source vectors at repair time.
}

void NormEncoderNanoRQ::EncodeParity(unsigned int parityId, const char** sourceVectorList, unsigned int numData, char* parityVector)
{
    ASSERT(rq_initialized);
    // A rateless encoder must not rely on internal cached state from Encode() because
    // NORM can ask it to generate parity for *any* block at any time.
    // 1. Initialize the NanoRQ matrix for the provided `numData`.
    nanorq_core_init_matrix(&rq, D, stride);
    
    // 2. Place all the source symbols into the matrix
    for (unsigned int i = 0; i < numData; i++)
    {
        if (sourceVectorList[i])
        {
            nanorq_core_place_symbol(&rq, D, stride, i, (const uint8_t*)sourceVectorList[i], vector_size);
        }
    }
    
    // 3. Precalculate the schedule and run it
    ops_run(&rq, D, stride, &S);
    
    // 4. Mix the requested parity symbol!
    // The parityId in NORM is 0-indexed for repair symbols, so for NanoRQ the ESI is numData + parityId
    uint8_t* tmp_buf = new uint8_t[stride];
    ops_mix(&rq, D, stride, numData + parityId, tmp_buf);
    memcpy(parityVector, tmp_buf, vector_size);
    delete[] tmp_buf;
}

NormDecoderNanoRQ::NormDecoderNanoRQ()
  : ndata(0), npar(0), vector_size(0), rq_initialized(false)
{
}

NormDecoderNanoRQ::~NormDecoderNanoRQ()
{
    Destroy();
}

bool NormDecoderNanoRQ::Init(unsigned int numData, unsigned int numParity, UINT16 vecSizeMax)
{
    Destroy();
    ndata = numData;
    npar = numParity;
    vector_size = vecSizeMax;
    rq_initialized = true;
    return true;
}

void NormDecoderNanoRQ::Destroy()
{
    rq_initialized = false;
}

int NormDecoderNanoRQ::Decode(char** vectorList, unsigned int numData, unsigned int erasureCount, unsigned int* erasureLocs)
{
    ASSERT(rq_initialized);
    
    unsigned int missing_source_count = 0;
    unsigned int received_repair_count = 0;
    
    unsigned int* missing_source_esis = new unsigned int[numData];
    unsigned int* received_repair_esis = new unsigned int[npar];
    
    bool* is_erasure = new bool[numData + npar];
    memset(is_erasure, 0, sizeof(bool) * (numData + npar));
    for (unsigned int i = 0; i < erasureCount; i++)
    {
        if (erasureLocs[i] < (numData + npar))
        {
            is_erasure[erasureLocs[i]] = true;
        }
    }
    
    for (unsigned int i = 0; i < numData; i++)
    {
        if (is_erasure[i])
        {
            missing_source_esis[missing_source_count++] = i;
        }
    }
    
    if (missing_source_count == 0)
    {
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    
    for (unsigned int i = numData; i < numData + npar; i++)
    {
        if (vectorList[i] != NULL && !is_erasure[i])
        {
            received_repair_esis[received_repair_count++] = i;
        }
        else
        {
            is_erasure[i] = true;
        }
    }
    
    if (received_repair_count < missing_source_count)
    {
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        PLOG(PL_ERROR, "NormDecoderNanoRQ::Decode() error: not enough received packets to decode\n");
        return erasureCount;
    }
    
    unsigned int overhead = received_repair_count - missing_source_count;
    
    nanorq_core dec_rq;
    if (!nanorq_core_encoder_new(numData, overhead, &dec_rq))
    {
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    
    size_t prep_len = nanorq_core_calculate_prepare_memory(&dec_rq);
    uint8_t* prep_mem = new uint8_t[prep_len];
    if (!nanorq_core_prepare(&dec_rq, prep_mem, prep_len))
    {
        delete[] prep_mem;
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    
    for (unsigned int i = 0; i < missing_source_count; i++)
    {
        nanorq_core_replace_symbol(&dec_rq, missing_source_esis[i], received_repair_esis[i]);
    }
    for (unsigned int i = 0; i < overhead; i++)
    {
        nanorq_core_replace_symbol(&dec_rq, dec_rq.P.Kprime + i, received_repair_esis[missing_source_count + i]);
    }
    
    nanorq_core_patch_matrix(&dec_rq);
    
    size_t work_len = nanorq_core_calculate_work_memory(&dec_rq);
    uint8_t* work_mem = new uint8_t[work_len];
    schedule S_dec;
    size_t sched_bytes = ops_estimate_schedule_bytes(numData + overhead);
    void* sched_buf = malloc(sched_bytes);
    if (!schedule_init(&S_dec, sched_buf, sched_bytes))
    {
        free(sched_buf);
        delete[] work_mem;
        delete[] prep_mem;
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    
    nanorq_core_set_op_callback(&dec_rq, &S_dec, ops_push);
    if (!nanorq_core_precalculate(&dec_rq, work_mem, work_len))
    {
        free(sched_buf);
        delete[] work_mem;
        delete[] prep_mem;
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    
    u32 rows = nanorq_core_get_pc_rows(&dec_rq);
    u32 stride = nanorq_core_recommended_stride(vector_size);
    uint8_t* D = (uint8_t*)obl_alloc(rows, stride, nanorq_oblas.align_size);
    if (!D)
    {
        free(sched_buf);
        delete[] work_mem;
        delete[] prep_mem;
        delete[] missing_source_esis;
        delete[] received_repair_esis;
        delete[] is_erasure;
        return erasureCount;
    }
    nanorq_core_init_matrix(&dec_rq, D, stride);
    
    for (unsigned int i = 0; i < numData; i++)
    {
        if (!is_erasure[i])
        {
            nanorq_core_place_symbol(&dec_rq, D, stride, i, (const uint8_t*)vectorList[i], vector_size);
        }
    }
    for (unsigned int i = 0; i < missing_source_count; i++)
    {
        nanorq_core_place_symbol(&dec_rq, D, stride, missing_source_esis[i], (const uint8_t*)vectorList[received_repair_esis[i]], vector_size);
    }
    for (unsigned int i = 0; i < overhead; i++)
    {
        nanorq_core_place_symbol(&dec_rq, D, stride, dec_rq.P.Kprime + i, (const uint8_t*)vectorList[received_repair_esis[missing_source_count + i]], vector_size);
    }
    
    ops_run(&dec_rq, D, stride, &S_dec);
    
    uint8_t* tmp_buf = new uint8_t[stride];
    for (unsigned int i = 0; i < missing_source_count; i++)
    {
        unsigned int missing_esi = missing_source_esis[i];
        ops_mix(&dec_rq, D, stride, missing_esi, tmp_buf);
        memcpy(vectorList[missing_esi], tmp_buf, vector_size);
    }
    delete[] tmp_buf;
    
    obl_free(D);
    free(sched_buf);
    delete[] work_mem;
    delete[] prep_mem;
    delete[] missing_source_esis;
    delete[] received_repair_esis;
    delete[] is_erasure;
    
    return 0;
}

NormEncoder* NormCreateNanoRQEncoder()
{
    return new NormEncoderNanoRQ();
}

NormDecoder* NormCreateNanoRQDecoder()
{
    return new NormDecoderNanoRQ();
}
