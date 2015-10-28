#ifndef DATATRAN_H
#define DATATRAN_H
#include <Partition.h>
void dataTransport(float *h_U, const Partition& pt, int tag);
void copy_UVW_buf(float *UVW, float *buf, const Partition& pt, int transportlen_side, int flag);
#endif // DATATRAN_H
