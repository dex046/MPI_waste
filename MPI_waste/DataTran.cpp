/******************************************
 * author:dwx
 * 2015.10.15
 ******************************************/

#include "DataTran.h"
//#include <testTDFWI.h>


void copydatatobuf(float *data, float *buf, const Partition& pt, uint transportlen_side, int tag, int flag)
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint topborder = 0;
    uint leftborder = 0;
    uint bottomborder = 0;
    uint rightborder = 0;
    uint length_x = block_x;
    uint length_z = block_z;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    if(flag == RIGHT_TO_LEFT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            memcpy(buf + i * transportlen_side, data + (i + topborder) * length_x + leftborder, sizeof(float) * transportlen_side);
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            memcpy(buf + i * transportlen_side, data + (i + 1 + topborder) * length_x - rightborder - transportlen_side, sizeof(float) * transportlen_side);
        }
    }
    if(flag == TOP_TO_BOTTOM)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            memcpy(buf + i * block_x, data + (block_z + topborder + i - transportlen_side) * length_x + leftborder, sizeof(float) * block_x);
        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            memcpy(buf + i * block_x, data + (topborder + i) * length_x + leftborder, sizeof(float) * block_x);
        }
    }
}
void copydatatobuf_Vp(float *data, float *buf, const Partition& pt, uint transportlen_side, int flag, const AFDPU2D &Pa)
{
    int gap_x = 0;
    int gap_z = 0;

    H_Border temph_Vp = pt.geth_Vp();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    //int length_x = temph_Vp.length_x;
    //int length_z = temph_Vp.length_z;

    uint topborder = temph_Vp.topborder;
    uint leftborder = temph_Vp.leftborder;
    //int bottomborder = temph_Vp.bottomborder;
    uint rightborder = temph_Vp.rightborder;

    if(flag == RIGHT_TO_LEFT)
    {
        gap_x = Pa.PMLx - indexmin_x;
        for(int i = 0; i < block_z; ++i)
        {
            memcpy(buf + i * transportlen_side, data + (i + topborder) * temph_Vp.length_x + leftborder + gap_x, sizeof(float) * transportlen_side);
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        gap_x = Pa.PMLx + Pa.Nx - 1 - indexmin_x;
        for(int i = 0; i < block_z; ++i)
        {
            memcpy(buf + i * transportlen_side, data + (i + 1 + topborder) * temph_Vp.length_x - rightborder - transportlen_side + gap_x, sizeof(float) * transportlen_side);
        }
    }
    if(flag == TOP_TO_BOTTOM)
    {
        gap_z = Pa.PMLz + Pa.Nz - 1 - indexmin_z;
        for(int i = 0; i < transportlen_side; ++i)
        {
            memcpy(buf + i * block_x, data + (topborder + i - transportlen_side + 1 + gap_z) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        gap_z = Pa.PMLz - indexmin_z;
        for(int i = 0; i < transportlen_side; ++i)
        {
            memcpy(buf + i * block_x, data + (topborder + i + gap_z) * temph_Vp.length_x + leftborder, sizeof(float) * block_x);
        }
    }
}



void copybuftodata(float *buf, float *data, const Partition& pt, uint transportlen_side, int tag, int flag)
{
    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint topborder = 0;
    uint leftborder = 0;
    uint bottomborder = 0;
    uint rightborder = 0;
    uint length_x = block_x;
    uint length_z = block_z;
    transportlen_side = 0;
    //uint length = 0;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;

        length_x = leftborder + block_x + rightborder;
        length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    if(flag == TOP_TO_BOTTOM)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            memcpy(data + i * length_x + leftborder, buf + i * block_x, sizeof(float) * block_x);
        }
    }
    if(flag == LEFT_TO_RIGHT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            memcpy(data + (i + topborder) * length_x, buf + i * transportlen_side, sizeof(float) * transportlen_side);
        }
    }
    if(flag == BOTTOM_TO_TOP)
    {
        for(uint i = 0; i < transportlen_side; ++i)
        {
            memcpy(data + (i + TOP_BORDER + block_z) * length_x, buf + i * block_x, sizeof(float) * block_x);
        }
    }
    if(flag == RIGHT_TO_LEFT)
    {
        for(uint i = 0; i < block_z; ++i)
        {
            memcpy(data + (i + topborder) * length_x + leftborder + block_x, buf + i * transportlen_side, sizeof(float) * transportlen_side);
        }
    }
}
void dataGather(float *data, const Partition& pt, int tag)
{
    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint topborder = 0;
    uint leftborder = 0;
    uint bottomborder = 0;
    uint rightborder = 0;
//    uint length_x = block_x;
//    uint length_z = block_z;
    uint transportlen_side = 0;
    uint length = 0;

    if(tag == STEP_U)
    {
        topborder = temph_U.topborder;
        leftborder = temph_U.leftborder;
        bottomborder = temph_U.bottomborder;
        rightborder = temph_U.rightborder;

//        length_x = leftborder + block_x + rightborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_V || tag == STEP_W)
    {
        topborder = temph_VW.topborder;
        leftborder = temph_VW.leftborder;
        bottomborder = temph_VW.bottomborder;
        rightborder = temph_VW.rightborder;

//        length_x = leftborder + block_x + rightborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else if(tag == STEP_VP)
    {
        topborder = temph_Vp.topborder;
        leftborder = temph_Vp.leftborder;
        bottomborder = temph_Vp.bottomborder;
        rightborder = temph_Vp.rightborder;

//        length_x = leftborder + block_x + rightborder;
//        length_z = topborder + block_z + bottomborder;
    }
    else
    {

    }

    MPI_Status status;
    MPI_Request request;


    if(topborder && !pt.isfirstblock_z())
    {
        length = topborder * block_x;
        float *buf = new float[length];
        MPI_Irecv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + TOP_TO_BOTTOM, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
//        for(int i = 0; i < length; ++i)
//            cout << *(buf++) << " ";
//        cout << status.MPI_SOURCE << endl;
//        cout << pt.getrank() << endl;
        transportlen_side = topborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, TOP_TO_BOTTOM);
        //cout << pt.getrank() << " top" << endl;
    }

    if(leftborder && !pt.isfirstblock_x())
    {

        length = leftborder * block_z;
        float *buf = new float[length];
        MPI_Irecv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + LEFT_TO_RIGHT, MPI_COMM_WORLD, &request);
//cout << pt.getrank() << "=================" << endl;
        MPI_Wait(&request, &status);


        transportlen_side = leftborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, LEFT_TO_RIGHT);
        //cout << pt.getrank() << " left" << endl;
    }

    if(bottomborder && !pt.islastblock_z())
    {
        length = bottomborder * block_x;
        float *buf = new float[length];
        MPI_Irecv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + BOTTOM_TO_TOP, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        transportlen_side = bottomborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, BOTTOM_TO_TOP);
        //cout << pt.getrank() << " bottom" << endl;
    }

    if(rightborder && !pt.islastblock_x())
    {
        length = rightborder * block_z;
        float *buf = new float[length];
        MPI_Irecv(buf, length, MPI_FLOAT, MPI_ANY_SOURCE, tag + RIGHT_TO_LEFT, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, &status);
        transportlen_side = rightborder;
        copybuftodata(buf, data, pt, transportlen_side, tag, RIGHT_TO_LEFT);
        //cout << pt.getrank() << " right"  << tag << endl;
        //cout << pt.getrank() << "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW" << endl;
    }
}


void dataTransport(float *data, const Partition& pt, int tag, MPI_Request *request)
{
    H_Border temph_U = pt.geth_U();
    H_Border temph_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();
    if(!pt.isfirstblock_x())
    {

        uint transportlength_x;
        switch(tag)
        {
        case STEP_U: transportlength_x = temph_U.rightborder;break;
        case STEP_V: transportlength_x = temph_VW.rightborder;break;
        case STEP_W: transportlength_x = temph_VW.rightborder;break;
        //case STEP_VP: transportlength_x = 1;break;
        };

        uint length = transportlength_x * pt.getblockLength_z();
        float *buf = new float[length];
        copydatatobuf(data, buf, pt, transportlength_x, tag, RIGHT_TO_LEFT);
        MPI_Isend(buf, length, MPI_FLOAT, pt.getrank() - 1, tag + RIGHT_TO_LEFT, MPI_COMM_WORLD, request);
    }
    if(!pt.islastblock_x())
    {
        uint transportlength_x;
        switch(tag)
        {
        case STEP_U: transportlength_x = temph_U.leftborder;break;
        case STEP_V: transportlength_x = temph_VW.leftborder;break;
        case STEP_W: transportlength_x = temph_VW.leftborder;break;
        //case STEP_VP: transportlength_x = 1;break;
        };

        uint length = transportlength_x * pt.getblockLength_z();
        float *buf = new float[length];
        copydatatobuf(data, buf, pt, transportlength_x, tag, LEFT_TO_RIGHT);
        MPI_Isend(buf, length, MPI_FLOAT, pt.getrank() + 1, tag + LEFT_TO_RIGHT, MPI_COMM_WORLD, request);
    }
    if(!pt.isfirstblock_z())
    {
        uint transportlength_z;
        switch(tag)
        {
        case STEP_U: transportlength_z = temph_U.bottomborder;break;
        case STEP_V: transportlength_z = temph_VW.bottomborder;break;
        case STEP_W: transportlength_z = temph_VW.rightborder;break;
        //case STEP_VP: transportlength_x = 1;break;
        };

        uint length = transportlength_z * pt.getblockLength_x();
        float *buf = new float[length];
        copydatatobuf(data, buf, pt, transportlength_z, tag, BOTTOM_TO_TOP);
        MPI_Isend(buf, length, MPI_FLOAT, pt.getrank() - pt.getsumBlock_x(), tag + BOTTOM_TO_TOP, MPI_COMM_WORLD, request);
    }
    if(!pt.islastblock_z())
    {
        uint transportlength_z;
        switch(tag)
        {
        case STEP_U: transportlength_z = temph_U.topborder;break;
        case STEP_V: transportlength_z = temph_VW.topborder;break;
        case STEP_W: transportlength_z = temph_VW.rightborder;break;
        //case STEP_VP: transportlength_x = 1;break;
        };

        uint length = transportlength_z * pt.getblockLength_x();
        float *buf = new float[length];
        copydatatobuf(data, buf, pt, transportlength_z, tag, TOP_TO_BOTTOM);
        MPI_Isend(buf, length, MPI_FLOAT, pt.getrank() + pt.getsumBlock_x(), tag + TOP_TO_BOTTOM, MPI_COMM_WORLD, request);
    }
}


void dataTransport_Vp(float *data, const Partition& pt, int tag, const AFDPU2D &Pa, MPI_Request *request)
{
    uint rank = pt.getrank();
    uint blockPosition_x = pt.getblockPosition_x();
    uint blockPosition_z = pt.getblockPosition_z();
    uint sumBlock_x = pt.getsumBlock_x();
    uint sumBlock_z = pt.getsumBlock_z();
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int transportlength_z = 1;
    if(tag == RIGHT_TO_LEFT && !pt.isfirstblock_x())
    {
        for(uint i = 1; i < blockPosition_x; ++i)
        {
            uint length = 1 * block_z;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_z, RIGHT_TO_LEFT, Pa);
            MPI_Isend(buf, length, MPI_FLOAT, rank - i, STEP_VP + RIGHT_TO_LEFT, MPI_COMM_WORLD, request);
        }
    }
    if(tag == LEFT_TO_RIGHT && !pt.islastblock_x())
    {
        for(uint i = 1; i <= sumBlock_x - blockPosition_x; ++i)
        {
            uint length = 1 * block_z;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_z, LEFT_TO_RIGHT, Pa);
            MPI_Isend(buf, length, MPI_FLOAT, rank + i, STEP_VP + LEFT_TO_RIGHT, MPI_COMM_WORLD, request);
        }
    }
    if(tag == BOTTOM_TO_TOP && !pt.isfirstblock_z())
    {
        for(uint i = 1; i < blockPosition_z ; ++i)
        {
            uint length = 1 * block_x;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_z, BOTTOM_TO_TOP, Pa);
            MPI_Isend(buf, length, MPI_FLOAT, rank - (i * sumBlock_x), STEP_VP + BOTTOM_TO_TOP, MPI_COMM_WORLD, request);
        }
    }
    if(tag == TOP_TO_BOTTOM && !pt.islastblock_z())
    {
        for(uint i = 1; i <= sumBlock_z - blockPosition_z; ++i)
        {
            uint length = 1 * block_x;
            float *buf = new float[length];
            copydatatobuf_Vp(data, buf, pt, transportlength_z, TOP_TO_BOTTOM, Pa);
            MPI_Isend(buf, length, MPI_FLOAT, rank + (i * sumBlock_x), STEP_VP + TOP_TO_BOTTOM, MPI_COMM_WORLD, request);
        }
    }
}

//void dataTransport_Vp(float *data, AFDPU2D &Pa, Partition &pt, int tag)
//{
//    if(pt.getindexmin_x() <= Pa.PMLx && pt.getindexmax_x() >= Pa.PMLx)
//    {
//        if(!pt.isfirstblock_x())
//        {
//            int length = pt.getblockLength_z();
//            float *buf = new float[length];
//            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
//            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - 1, UpdateVpPML, MPI_COMM_WORLD);
//        }
//    }
//    else if(pt.getindexmin_x() <= Pa.PMLx + Pa.Nx && pt.getindexmax_x() >= Pa.PMLx + Pa.Nx)
//    {
//        if(!pt.islastblock_x())
//        {
//            int length = pt.getblockLength_z();
//            float *buf = new float[length];
//            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
//            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + 1, UpdateVpPML, MPI_COMM_WORLD);
//        }
//    }
//    else if(pt.getindexmin_z() <= Pa.PMLz && pt.getindexmax_z() >= Pa.PMLz)
//    {
//        if(!pt.isfirstblock_z())
//        {
//            int length = pt.getblockLength_x();
//            float *buf = new float[length];
//            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
//            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() - pt.getsumBlock_x(), UpdateVpPML, MPI_COMM_WORLD);
//        }
//    }
//    else if(pt.getindexmin_z() <= Pa.PMLz + Pa.Nz && pt.getindexmax_z() >= Pa.PMLz + Pa.Nz)
//    {
//        if(!pt.isfirstblock_z())
//        {
//            int length = pt.getblockLength_x();
//            float *buf = new float[length];
//            copy_UVW_buf(data, buf, pt, 1, UpdateVpPML);
//            MPI_Send(buf, length, MPI_FLOAT, pt.getrank() + pt.getsumBlock_x(), UpdateVpPML, MPI_COMM_WORLD);
//        }
//    }
//    else
//    {

//    }
//}
