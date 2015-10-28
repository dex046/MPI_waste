/******************************************
 * author:dwx
 * 2015.10.15
 ******************************************/
#include <mpi.h>
#include <DataTran.h>
#include <TDFWI.h>
#include <RWsgy.h>
#define RIGHT_TO_LEFT -2
#define LEFT_TO_RIGHT 2
#define BOTTOM_TO_TOP -1
#define TOP_TO_BOTTOM 1

#define STEP_U 0
#define STEP_VW 1

#define UpdateVpPML 3
void copy_UVW_buf(float *UVW, float *buf, const Partition& pt, int transportlen_side, int flag)
{
    int length_x = pt.getblockLength_x();
    int length_z = pt.getblockLength_z();
    if(flag == RIGHT_TO_LEFT)
    {
        for(int i = 0; i < length_z; ++i)
        {
            for(int j = 0; j < transportlen_side; ++j)
            {
                buf[i * transportlen_side + j] = UVW[i * length_x + j];
            }
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        for(int i = 0; i < length_z; ++i)
        {
            for(int j = length_x - transportlen_side; j < length_x; ++j)
            {
                buf[i * transportlen_side + j] = UVW[i * length_x + j];
            }
        }
    }
    if(flag == TOP_TO_BOTTOM)
    {
        for(int i = length_z - transportlen_side; i < length_z; ++i)
        {
            for(int j = 0; j < length_x; ++j)
            {
                buf[i - length_z + transportlen_side + j] = UVW[i * length_x + j];
            }
        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        for(int i = 0; i < transportlen_side; ++i)
        {
            for(int j = 0; j < length_x; ++j)
            {
                buf[i * transportlen_side + j] = UVW[i * length_x + j];
            }
        }
    }
}

void dataTransport(float *UVW, const Partition& pt, int tag)
{
    if(!pt.isfirstblock_x())
    {

        int transportlength_x;
        switch(tag)
        {
        case STEP_U: transportlength_x = 4;break;
        case STEP_VW: transportlength_x = 3;break;
        };

        int length = transportlength_x * pt.getblockLength_z();
        float *buf = new float[length];
        copy_UVW_buf(UVW, buf, pt, transportlength_x, RIGHT_TO_LEFT);
        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - 1, RIGHT_TO_LEFT, MPI_COMM_WORLD);
    }
    if(!pt.islastblock_x())
    {
        int transportlength_x;
        switch(tag)
        {
        case STEP_U: transportlength_x = 3;break;
        case STEP_VW: transportlength_x = 4;break;
        };

        int length = transportlength_x * pt.getblockLength_z();
        float *buf = new float[length];
        copy_UVW_buf(UVW, buf, pt, transportlength_x, LEFT_TO_RIGHT);
        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + 1, LEFT_TO_RIGHT, MPI_COMM_WORLD);
    }
    if(!pt.isfirstblock_z())
    {
        int transportlength_z;
        switch(tag)
        {
        case STEP_U: transportlength_z = 4;break;
        case STEP_VW: transportlength_z = 3;break;
        };

        int length = transportlength_z * pt.getblockLength_x();
        float *buf = new float[length];
        copy_UVW_buf(UVW, buf, pt, transportlength_z, BOTTOM_TO_TOP);
        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - pt.getsumBlock_x(), BOTTOM_TO_TOP, MPI_COMM_WORLD);
    }
    if(!pt.islastblock_z())
    {
        int transportlength_z;
        switch(tag)
        {
        case STEP_U: transportlength_z = 3;break;
        case STEP_VW: transportlength_z = 4;break;
        };

        int length = transportlength_z * pt.getblockLength_x();
        float *buf = new float[length];
        copy_UVW_buf(UVW, buf, pt, transportlength_z, TOP_TO_BOTTOM);
        MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + pt.getsumBlock_x(), TOP_TO_BOTTOM, MPI_COMM_WORLD);
    }
}

void dataTransport_Vp(float *data, AFDPU2D &Pa, Partition &pt, int tag)
{
    if(pt.getindexmin_x() <= Pa.PMLx && pt.getindexmax_x() >= Pa.PMLx)
    {
        if(!pt.isfirstblock_x())
        {
            int length = pt.getblockLength_z();
            float *buf = new float[length];
            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - 1, UpdateVpPML, MPI_COMM_WORLD);
        }
    }
    else if(pt.getindexmin_x() <= Pa.PMLx + Pa.Nx && pt.getindexmax_x() >= Pa.PMLx + Pa.Nx)
    {
        if(!pt.islastblock_x())
        {
            int length = pt.getblockLength_z();
            float *buf = new float[length];
            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + 1, UpdateVpPML, MPI_COMM_WORLD);
        }
    }
    else if(pt.getindexmin_z() <= Pa.PMLz && pt.getindexmax_z() >= Pa.PMLz)
    {
        if(!pt.isfirstblock_z())
        {
            int length = pt.getblockLength_x();
            float *buf = new float[length];
            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - pt.getsumBlock_x(), UpdateVpPML, MPI_COMM_WORLD);
        }
    }
    else if(pt.getindexmin_z() <= Pa.PMLz + Pa.Nz && pt.getindexmax_z() >= Pa.PMLz + Pa.Nz)
    {
        if(!pt.isfirstblock_z())
        {
            int length = pt.getblockLength_x();
            float *buf = new float[length];
            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + pt.getsumBlock_x(), UpdateVpPML, MPI_COMM_WORLD);
        }
    }
}
