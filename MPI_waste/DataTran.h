#ifndef DATATRAN_H
#define DATATRAN_H
#include "Partition.h"
#include <mpi.h>
#include "RWsgy.h"

//#define U_TOP_TO_BOTTOM 0
//#define U_LEFT_TO_RIGHT 1
//#define U_BOTTOM_TO_TOP 2
//#define U_RIGHT_TO_LEFT 3

//#define V_TOP_TO_BOTTOM 4
//#define V_LEFT_TO_RIGHT 5
//#define V_BOTTOM_TO_TOP 6
//#define V_RIGHT_TO_LEFT 7

//#define W_TOP_TO_BOTTOM 8
//#define W_LEFT_TO_RIGHT 9
//#define W_BOTTOM_TO_TOP 10
//#define W_RIGHT_TO_LEFT 11

//#define VP_TOP_TO_BOTTOM 12
//#define VP_LEFT_TO_RIGHT 13
//#define VP_BOTTOM_TO_TOP 14
//#define VP_RIGHT_TO_LEFT 15

#define TOP_TO_BOTTOM 0
#define LEFT_TO_RIGHT 1
#define BOTTOM_TO_TOP 2
#define RIGHT_TO_LEFT 3

#define STEP_U 0
#define STEP_V 4
#define STEP_W 8
#define STEP_VP 12

//#define UpdateVpPML 3
void dataTransport(float *data, const Partition& pt, int tag, MPI_Request *request);
void dataTransport_Vp(float *data, const Partition& pt, int tag, const AFDPU2D &Pa, MPI_Request *request);

void dataGather(float *data, const Partition& pt, int tag);

void copydatatobuf(float *data, float *buf, const Partition& pt, uint transportlen_side, int tag, int flag);
void copydatatobuf_Vp(float *data, float *buf, const Partition& pt, uint transportlen_side, int flag, const AFDPU2D &Pa);
void copybuftodata(float *buf, float *data, const Partition& pt, uint transportlen_side, int tag, int flag);
#endif // DATATRAN_H
