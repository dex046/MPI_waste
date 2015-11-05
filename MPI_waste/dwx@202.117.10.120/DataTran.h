#ifndef DATATRAN_H
#define DATATRAN_H
#include "Partition.h"
#include <mpi.h>
#include "RWsgy.h"

#define RIGHT_TO_LEFT -2
#define LEFT_TO_RIGHT 2
#define BOTTOM_TO_TOP -1
#define TOP_TO_BOTTOM 1

#define STEP_U 0
#define STEP_VW 1
#define STEP_VP 2

//#define UpdateVpPML 3
void dataTransport(float *data, const Partition& pt, int tag, MPI_Request *request);
void dataTransport_Vp(float *data, const Partition& pt, int tag, const AFDPU2D &Pa, MPI_Request *request);

void dataGather(float *data, const Partition& pt, int tag, MPI_Request *request);

void copydatatobuf(float *data, float *buf, const Partition& pt, uint transportlen_side, int tag, int flag);
void copydatatobuf_Vp(float *data, float *buf, const Partition& pt, uint transportlen_side, int flag, const AFDPU2D &Pa);
void copybuftodata(float *buf, float *data, const Partition& pt, uint transportlen_side, int tag, int flag);
#endif // DATATRAN_H
