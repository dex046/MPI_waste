/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#include    "mpi.h"
#include	"testTDFWI.h"
#define     TOP     0
#define     LEFT    1
#define     BOTTOM  2
#define     RIGHT   3

#define     ADD         1
#define     SUBSTRACT   0

using namespace std;

/*------------------------------------------------------------------------
function: Ricker

Information:
返回t时刻Ricker子波的值

Inputs:
f0: Ricker子波的主频
t: 时刻

Outputs:
t时刻Ricker子波的值
------------------------------------------------------------------------*/
float Ricker(float f0,
            float t)
{
    float RV = 0.0f, temp = 0.0f;

    temp = PI * f0 * (t - 1.0f / f0);
    RV = (1 - 2 * temp * temp) * exp(-1.0f * temp * temp);

    return RV;
}

/*------------------------------------------------------------------------
function: ReadData

Information:
read data from Sgy file

Inputs:
FileName: name of the Sgy file
Data: output data
------------------------------------------------------------------------*/
void ReadData(char FileName[],
            float *Data,
            const Partition &pt,
            usht flag)//flag ?
{
    FILE *f1;
    unsigned char f3200[3200];
    //Struct_reelb400 FileHeader;
    REEL *Reel;
    unsigned short *TraceNum, *SampleNum, *SampleInt;
    short *DFormat;
    bool BReel, BIBM;
    uint n,m;

    int rank = pt.getrank();

    TraceNum = new unsigned short[1];
    memset((void *)TraceNum, 0, sizeof(unsigned short));
    SampleNum = new unsigned short[1];
    memset((void *)SampleNum, 0, sizeof(unsigned short));
    SampleInt = new unsigned short[1];
    memset((void *)SampleInt, 0, sizeof(unsigned short));
    DFormat = new short[1];
    memset((void *)DFormat, 0, sizeof(short));
    Reel = new REEL[1];
    memset((void *)Reel, 0, sizeof(REEL));
//cout << "ww" << endl;
    f1 = fopen(FileName, "rb");
    if (f1 == NULL)
    {
        cout << "File Open error!" << "\n" << endl;
    }
    if(rank == ROOT_ID)
    {
        // read the first 3200 bytes
        fread(f3200, 3200, 1, f1);
    }

    // read the information of Sgy file
    InfoOfSgy(FileName, *Reel, TraceNum,
        SampleNum, SampleInt, DFormat,
        &BReel, &BIBM);
//coiut << 1231421 << endl;
//    cout << BReel << endl;
//    cout << *TraceNum << endl;
//    cout << *SampleNum << endl;
//    cout << *SampleInt << endl;
//    cout << *DFormat << endl;
//    cout << BIBM << endl;

    int length_x = pt.getblockLength_x();
    int length_z = pt.getblockLength_z();

//cout << "rank" << rank << "rank" << length_x << " " << length_z << " " << endl;
    Trace *trace;
    trace = new Trace[length_x];
    memset((void *)trace, 0, sizeof(Trace) * length_x);
    for (n = 0; n < length_x; n++)
    {
        trace[n].data = new float[length_z];
        memset((void *)trace[n].data, 0, sizeof(float) * length_z);
    }

    // read the data
    ReadSgyData(FileName, trace, *Reel, SampleNum,
        DFormat, &BReel, &BIBM, pt);

    // write the trace data to Data
    for(n = 0; n < length_x; ++n)
    {
        for(m = 0; m < length_z; ++m)
        {
            if(flag == 0)
            {
                Data[m * length_x + n] = trace[n].data[m];
                //cout << trace[n].data[m];
            }
            else
            {
                Data[n * length_z + m] = trace[n].data[m];
                //cout << trace[n].data[m];
            }
        }
    }

    // free memory
    for (n = 0; n < length_x; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;
    delete []Reel;
    delete []DFormat;
    delete []SampleInt;
    delete []TraceNum;
    delete []SampleNum;

}
void write_sgs_t_Data(const char * const FileName, usht SampleNum, usht TraceNum, usht SampleInt, float *data, const Partition &pt, const AFDPU2D &Pa, usht flag)
{
    unsigned char f3200[3200];
    memset((void*)f3200, 0, 3200);
    uint m, n;

    uint RL_num = pt.getRL_num();
    Trace *trace;
    trace = new Trace[RL_num];
    memset((void *)trace, 0, sizeof(Trace) * RL_num);

    int begin_num = pt.getRL_beginnum();
    int end_num = pt.getRL_endnum();

    int rank = pt.getrank();

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, FileName, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_Offset offset = 0, filesize = 0;
    MPI_Status status;


    filesize = 3600 + (SampleNum * sizeof(float) + 240) * TraceNum;

    MPI_File_set_size(fh, filesize);

    for (int n = 0; n < RL_num; n++)
    {
        trace[n].data = new float[SampleNum];
        memset((void *)trace[n].data, 0, sizeof(float) * SampleNum);
        for (m = 0; m < SampleNum; m++)
        {
            if (flag == 0)
            {
                trace[n].data[m] = *(data + m * RL_num + n);
            }
            else
            {
                trace[n].data[m] = *(data + n * SampleNum + m);
            }
        }
    }

    // trace header
    for (n = 0; n < RL_num; n++)
    {
        for (int m = 0; m < 60; m++)
        {
            trace[n].head.h4[m] = 0;
        }
    }

    if(rank == ROOT_ID)
    {
        // reel header
        for (n = 0; n < 3200; n++)
        {
            f3200[n] = 1;
        }
    }

    //cout << rank << endl;

    if(rank == ROOT_ID)
    {
        MPI_File_write_at(fh, offset, &f3200[0], 3200, MPI_BYTE, &status);

        // 写卷头中400个字节
        REEL reel;
        reel.reelstruct.hns = SampleNum;
        reel.reelstruct.hdt = SampleInt;
        reel.reelstruct.format = 5; // IEEE float
        reel.reelstruct.mfeet = 1;
        //fwrite(&reel.reelstruct, 400, 1, fp);
        MPI_File_write_at(fh, offset + 3200, &reel.reelstruct, 400, MPI_BYTE, &status);
    }

    offset += 3600 + begin_num * (240 + SampleNum * sizeof(float));

    for (int i = 0; i < RL_num; i++)
    {
        // 写道头
        trace[i].head.headstruct.cdp = i + begin_num;//
        trace[i].head.headstruct.ns = SampleNum;
        trace[i].head.headstruct.dt = SampleInt;
        trace[i].head.headstruct.sx = 100000;
        trace[i].head.headstruct.sy = 1000000 + (i + begin_num) * 4;//

        MPI_File_write_at(fh, offset, &trace[i].head.headstruct, 240, MPI_BYTE, &status);

        offset += 240;
        MPI_File_write_at(fh, offset, &trace[i].data[0], SampleNum, MPI_FLOAT, &status);
        offset += SampleNum * sizeof(float);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&fh);
    MPI_Barrier(MPI_COMM_WORLD);

    // free memory
    for (n = 0; n < RL_num; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;
}

/*------------------------------------------------------------------------
function: WriteData

Information:
write data to Sgy file

Inputs:
FileName: name of the Sgy file
SampleNum: number of samples in each trace
TraceNum: number of traces
SampleInt: sample interval
data: output data
------------------------------------------------------------------------*/
void WriteData(const char * const FileName,
            usht SampleNum, //nnz 234
            usht TraceNum, //nnx 610
            usht SampleInt,
            float *data,
            const Partition& pt,
            const AFDPU2D Pa,
            usht flag,
            usht tag)
{
    unsigned char f3200[3200];
    memset((void*)f3200, 0, 3200);
    uint m, n;

    MPI_Offset filesize = (SampleNum * sizeof(float) + 240) * TraceNum + 3600;

    int indexmin_x = pt.getindexmin_x();
    int indexmin_z = pt.getindexmin_z();

    uint block_x = 0;
    uint block_z = 0;

    int rank = pt.getrank();

    if(tag == WRITE_INTER)
    {
        block_x = pt.getinteriorLength_x();
        block_z = pt.getinteriorLength_z();

//        if(!block_x || !block_z)
//            return;
    }
    else if(tag == WRITE_ALL)
    {
        block_x = pt.getblockLength_x();
        block_z = pt.getblockLength_z();
    }
    else
    {

    }

    Trace *trace;
    trace = new Trace[block_x];
    memset((void *)trace, 0, sizeof(Trace) * block_x);

    for (int n = 0; n < block_x; n++)
    {
        trace[n].data = new float[block_z];
        memset((void *)trace[n].data, 0, sizeof(float) * block_z);
        for (m = 0; m < block_z; m++)
        {
            if (flag == 0)
            {
                trace[n].data[m] = *(data + m * block_x + n);
            }
            else
            {
                trace[n].data[m] = *(data + n * block_z + m);
            }
        }
    }

    // trace header
    for (n = 0; n < block_x; n++)
    {
        for (int m = 0; m < 60; m++)
        {
            trace[n].head.h4[m] = 0;
        }
    }

    if(rank == ROOT_ID)
    {
        // reel header
        for (n = 0; n < 3200; n++)
        {
            f3200[n] = 1;
        }
    }

    WriteSgy(FileName, &f3200[0], trace,
        TraceNum, SampleNum, SampleInt, pt, Pa, filesize, tag);///hai

    // free memory
    for (n = 0; n < block_x; n++)
    {
        delete []trace[n].data;
    }
    delete []trace;

}


/*------------------------------------------------------------------------
function: MallocVariables

Information:
在内存上为变量申请空间

Inputs:
Pa: 正演参数
ip: 反演参数
plan: 全局变量
------------------------------------------------------------------------*/
void MallocVariables(AFDPU2D Pa,
                    IP *ip,
                    CPUVs *plan, const Partition &pt)
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint RL_num = pt.getRL_num();

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_z = pt.getinteriorLength_z();

    H_Border h_U = pt.geth_U();
    H_Border h_Vp = pt.geth_Vp();
    H_Border h_VW = pt.geth_VW();

    uint sum_h_Coord = pt.geth_Coord_length();

    try
    {
        // 申请变量// 与NPML有关的变量
        plan->h_dx = new float[block_x];
        plan->h_dz = new float[block_z];
        plan->h_Bx = new float[block_x];
        plan->h_Bz = new float[block_z];
        plan->h_PHIx_V_x = new float[block_z * block_x];
        plan->h_PHIz_W_z = new float[block_z * block_x];
        plan->h_PHIx_U_x = new float[block_z * block_x];//h_U right 4, left 3, seft 1
        plan->h_PHIz_U_z = new float[block_z * block_x];//h_U top 3, bottom 4, seft 1
        plan->h_PHIx_V_x_r = new float[block_z * block_x];//h_V left 4, right 3, seft 1
        plan->h_PHIz_W_z_r = new float[block_z * block_x];//h_W top 4, bottom 3, seft 1
        plan->h_PHIx_U_x_r = new float[block_z * block_x];
        plan->h_PHIz_U_z_r = new float[block_z * block_x];
        // 波场相关变量
        plan->h_U_past = new float[block_z * block_x];//h_V left 4, right 3, h_W top 4, bottom 3
        plan->h_U_now = new float[h_U.length_z * h_U.length_x];
        plan->h_U_next = new float[block_z * block_x];//h_V left 4, right 3, h_W top 4, bottom 3
        plan->h_U_past_r = new float[block_z * block_x];
        plan->h_U_now_r = new float[h_U.length_z * h_U.length_x];
        plan->h_U_next_r = new float[block_z * block_x];
        plan->h_V = new float[h_VW.length_z * h_VW.length_x];//h_U left 3, right 4
        plan->h_W = new float[h_VW.length_z * h_VW.length_x];//h_U top 3, bottom 4
        plan->h_V_r = new float[h_VW.length_z * h_VW.length_x];
        plan->h_W_r = new float[h_VW.length_z * h_VW.length_x];
        plan->h_Vp = new float[h_Vp.length_z * h_Vp.length_x];



        if(RL_num)
        {
            plan->h_TrueWF = new float[RL_num * Pa.Nt];
            plan->h_CurrWF = new float[RL_num * Pa.Nt];
            plan->h_ResWF = new float[RL_num * Pa.Nt];

            plan->h_Obj = new float[RL_num * Pa.Nt];
            // 计算步长的相关变量
            plan->h_TrailWF = new float[RL_num * Pa.Nt];
            plan->h_SumResTrial = new float[RL_num * Pa.Nt];
            plan->h_SumResCurr = new float[RL_num * Pa.Nt];
        }


        if(interiorlength_x && interiorlength_z)
        {
            plan->h_Grad = new float[interiorlength_z * interiorlength_x];//Pa.Nz * Pa.Nx
            plan->h_U_Der = new float[interiorlength_z * interiorlength_x];//Pa.Nz * Pa.Nx
        }


        if(sum_h_Coord)
        {//if(pt.getrank() == 39)cout << pt.getrank() << " " << sum_h_Coord << endl;
            plan->h_Coord = new uint[sum_h_Coord];
            plan->h_RU = new float[sum_h_Coord * (Pa.Nt - 1)];
        }
    }
    catch(const std::bad_alloc& e)
    {
        cout << "2Allocation failed: www   " << pt.getrank() << endl;
    }

    // 初始化空间
    memset((void *)plan->h_dx,			0,	sizeof(float) * block_x);
    memset((void *)plan->h_dz,			0,	sizeof(float) * block_z);
    memset((void *)plan->h_Bx,			0,	sizeof(float) * block_x);
    memset((void *)plan->h_Bz,			0,	sizeof(float) * block_z);
    memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_U_now,		0,	sizeof(float) * h_U.length_z * h_U.length_x);
    memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_U_past_r,	0,	sizeof(float) * block_z * block_x);
    memset((void *)plan->h_U_now_r,		0,	sizeof(float) * h_U.length_z * h_U.length_x);
    memset((void *)plan->h_U_next_r,	0,	sizeof(float) * block_z * block_x);

    memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
    memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
    memset((void *)plan->h_V_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
    memset((void *)plan->h_W_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
    memset((void *)plan->h_Vp,			0,	sizeof(float) * h_Vp.length_z * h_Vp.length_x);

    if(RL_num)
    {
        memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);

        memset((void *)plan->h_Obj,			0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_TrailWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_SumResTrial,	0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_SumResCurr,	0,	sizeof(float) * RL_num * Pa.Nt);
    }

    if(interiorlength_x && interiorlength_z)
    {
        memset((void *)plan->h_Grad,		0,	sizeof(float) * interiorlength_z * interiorlength_x);
        memset((void *)plan->h_U_Der,		0,	sizeof(float) * interiorlength_z * interiorlength_x);
    }



    if(sum_h_Coord)
    {
        //cout << "AAAAA" << pt.getrank() << " " << sum_h_Coord << "AAAAA" << endl;
        memset((void *)plan->h_Coord,		0,	sizeof(uint) * sum_h_Coord);

        memset((void *)plan->h_RU,			0,	sizeof(float) * sum_h_Coord * (Pa.Nt - 1));
    }

        //cout << pt.getrank() << endl;
}


/*------------------------------------------------------------------------
function: GenerateNPML

Information:
生成NPML系数

Inputs:
Pa: 正演参数
plan: 全局变量
------------------------------------------------------------------------*/
void GenerateNPML(AFDPU2D Pa,
                CPUVs *plan,
                 const Partition &pt)
{
    uint length_z = pt.getblockLength_z();
    uint length_x = pt.getblockLength_x();

    // 对NPML参数进行赋值
    float d0 = -3.0f * 2500.0f * logf(1.0e-6f)
        / (2.0f * powf((float)(Pa.PMLz * Pa.dz), 3.0f));
    float d1 = -3.0f * 2500.0f * log(1.0e-6f)
        / (2.0f * powf((float)(Pa.PMLx * Pa.dx), 3.0f));

//    uint n = pt.getinsidenum();
//    Inside *ins = pt.getInside();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    for (uint iz = 0; iz < length_z; iz++)
    {
        int actual_iz = iz + min_z;
        if (actual_iz < Pa.PMLz)
        {
            plan->h_dz[iz] = d0 * (Pa.PMLz - 1 - actual_iz + 0.5f)
                * (Pa.PMLz - 1 - actual_iz + 0.5f) * Pa.dz * Pa.dz;// NPML纵向参数
        }
        else if (actual_iz > Pa.PMLz + Pa.Nz - 1)
        {
            plan->h_dz[iz] = d0 * (actual_iz - Pa.PMLz - Pa.Nz + 0.5f)
                * (actual_iz - Pa.PMLz - Pa.Nz + 0.5f) * Pa.dz * Pa.dz;// NPML纵向参数
        }
        else
        {
            plan->h_dz[iz] = 0.0f;
        }
        plan->h_Bz[iz] = expf(-1.0f * plan->h_dz[iz] * Pa.dt);// NPML纵向参数
    }


    for (uint ix = 0; ix < length_x; ix++)
    {
        uint actual_ix = ix + min_x;
        if (actual_ix < Pa.PMLx)
        {
            plan->h_dx[ix] = d1 * (Pa.PMLx - 1 - actual_ix + 0.5f) *
                (Pa.PMLx - 1 - actual_ix + 0.5f) * Pa.dx * Pa.dx;// NPML横向参数
        }
        else if (actual_ix > Pa.PMLx + Pa.Nx - 1)
        {
            plan->h_dx[ix] = d1 * (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) *
                (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) * Pa.dx * Pa.dx;// NPML横向参数
        }
        else
        {
            plan->h_dx[ix] = 0.0f;// NPML横向参数
        }
        plan->h_Bx[ix] = expf(-1.0f * plan->h_dx[ix] * Pa.dt);// NPML横向参数
    }
}

/*------------------------------------------------------------------------
function: MatMax

Information:
找出矩阵中的最大值

Inputs:
Mat: 矩阵
Len：矩阵的维数
------------------------------------------------------------------------*/
float MatMax(float *Mat,
            uint Len)
{
    float MaxValue = fabs(Mat[0]);

    for (uint n = 0; n < Len; n++)
    {
        if (MaxValue < fabs(Mat[n]))
        {
            MaxValue = fabs(Mat[n]);
        }
    }

    return MaxValue;
}


/*------------------------------------------------------------------------
function: AddSource

Information:
正演加震源

Inputs:
Pa：正演参数
h_U：当前时刻的波场
s: 震源位置
Wave: 当前时刻波场的值
h_Vp: 当前模型速度值的平方
------------------------------------------------------------------------*/
void AddSource(AFDPU2D Pa,
                float *h_U,//h_U_next
                float Wave,
                float *h_Vp,
                SL s,
               const Partition &pt)
{
    H_Border temph_Vp = pt.geth_Vp();
//    for(int i = 0; i < pt.getblockLength_z(); ++i)
//    {
//        for(int j = 0; j < pt.getblockLength_x(); ++j)
//        {
//            cout << *(h_Vp + j + temph_Vp.leftborder + (i + temph_Vp.topborder) * temph_Vp.length_x);
//        }
//    }
    uint block_x = pt.getblockLength_x();

//    vector<pair<uint, uint>> vec = pt.getShot();
//    pair<uint, uint> temp;

//    for(vector<pair<uint, uint>>::iterator iter = vec.begin(); iter != vec.end(); ++iter)
//    {
//        temp = *(iter + is);
//        h_U[temp.second * block_x + temp.first] += Wave * (Pa.dt * Pa.dt * h_Vp[temph_Vp.leftborder + (temp.second + temph_Vp.topborder) * temph_Vp.length_x + temp.first]);
//    //cout << h_Vp[temph_Vp.leftborder + (temp.second + temph_Vp.topborder) * temph_Vp.length_x + temp.first];
//    }

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    uint ix = s.Sx - indexmin_x;
    uint iz = s.Sz - indexmin_z;

    h_U[iz * block_x + ix] += Wave * (Pa.dt * Pa.dt * h_Vp[temph_Vp.leftborder + (iz + temph_Vp.topborder) * temph_Vp.length_x + ix]);
}

/*------------------------------------------------------------------------
function: StepPHIU

Information:
一步更新波场U的卷积项

Inputs:
Pa：正演参数
h_U：当前时刻的波场
h_PHIx_U_x: 横向卷积项
h_PHIz_U_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
------------------------------------------------------------------------*/
void StepPHIU(AFDPU2D Pa,
            float *h_U,
            float *h_PHIx_U_x,
            float *h_PHIz_U_z,
            float *h_Bx,
            float *h_Bz,
            const Partition &pt, int it)///h_PHI_U h_U
{
    H_Border temph_U = pt.geth_U();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    float dUx = 0.0f;
    float dUz = 0.0f;

    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    // 交错网格有限差分
    //yuan cheng xu li cong 4 kai shi,zhe li ye shi cong 4 kai shi
    for (uint iz = 0 + gap_min_z + temph_U.topborder; iz < block_z + temph_U.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + temph_U.leftborder; ix < block_x + temph_U.leftborder + gap_max_x; ix++)
        {
            dUx = C1_4 * (h_U[iz * temph_U.length_x + ix + 1] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[iz * temph_U.length_x + ix + 2] - h_U[iz * temph_U.length_x + ix - 1])
                + C3_4 * (h_U[iz * temph_U.length_x + ix + 3] - h_U[iz * temph_U.length_x + ix - 2])
                + C4_4 * (h_U[iz * temph_U.length_x + ix + 4] - h_U[iz * temph_U.length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * temph_U.length_x + ix] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[(iz + 2) * temph_U.length_x + ix] - h_U[(iz - 1) * temph_U.length_x + ix])
                + C3_4 * (h_U[(iz + 3) * temph_U.length_x + ix] - h_U[(iz - 2) * temph_U.length_x + ix])
                + C4_4 * (h_U[(iz + 4) * temph_U.length_x + ix] - h_U[(iz - 3) * temph_U.length_x + ix]);

            h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] =
                h_Bx[ix - temph_U.leftborder] * (h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] + dUx / Pa.dx) - dUx / Pa.dx;
            h_PHIz_U_z[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] =
                h_Bz[iz - temph_U.topborder] * (h_PHIz_U_z[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] + dUz / Pa.dz) - dUz / Pa.dz;

            if(it == 45 && pt.getrank() == 1 && iz == temph_U.topborder + 53 && ix == 1 + temph_U.leftborder)
            {
                //cout << h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] << endl;
                //cout << ix - temph_U.leftborder + pt.getindexmin_x() << " " << h_Bx[ix - temph_U.leftborder] << endl;
                //cout << dUx << " " << dUz << endl;
//                cout << h_U[iz * temph_U.length_x + ix + 1]
//                     << " " << h_U[iz * temph_U.length_x + ix]
//                     << " " << h_U[iz * temph_U.length_x + ix + 2]
//                     << " " << h_U[iz * temph_U.length_x + ix - 1]
//                     << " " << h_U[iz * temph_U.length_x + ix + 3]
//                     << " " << h_U[iz * temph_U.length_x + ix - 2]
//                     << " " << h_U[iz * temph_U.length_x + ix + 4]
//                     << " " << h_U[iz * temph_U.length_x + ix - 3]
//                     << endl;
            }
        }
    }
}

/*------------------------------------------------------------------------
function: StepPHIVW

Information:
一步更新波场V和W的卷积项

Inputs:
Pa：正演参数
h_V：当前时刻的波场V
h_W：当前时刻的波场W
h_PHIx_V_x: 横向卷积项
h_PHIz_W_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
------------------------------------------------------------------------*/
void StepPHIVW(AFDPU2D Pa,
                float *h_V,
                float *h_W,
                float *h_PHIx_V_x,
                float *h_PHIz_W_z,
                float *h_Bx,
                float *h_Bz,
                const Partition &pt)
{
    H_Border h_VW = pt.geth_VW();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    float dVx = 0.0f;
    float dWz = 0.0f;

    // 交错网格有限差分
    for (uint iz = 0 + gap_min_z + h_VW.topborder; iz < block_z + h_VW.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + h_VW.leftborder; ix < block_x + h_VW.leftborder + gap_max_x; ix++)
        {
            dVx = C1_4 * (h_V[iz * h_VW.length_x + ix] - h_V[iz * h_VW.length_x + ix - 1])
                + C2_4 * (h_V[iz * h_VW.length_x + ix + 1] - h_V[iz * h_VW.length_x + ix - 2])
                + C3_4 * (h_V[iz * h_VW.length_x + ix + 2] - h_V[iz * h_VW.length_x + ix - 3])
                + C4_4 * (h_V[iz * h_VW.length_x + ix + 3] - h_V[iz * h_VW.length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * h_VW.length_x + ix] - h_W[(iz - 1) * h_VW.length_x + ix])
                + C2_4 * (h_W[(iz + 1) * h_VW.length_x + ix] - h_W[(iz - 2) * h_VW.length_x + ix])
                + C3_4 * (h_W[(iz + 2) * h_VW.length_x + ix] - h_W[(iz - 3) * h_VW.length_x + ix])
                + C4_4 * (h_W[(iz + 3) * h_VW.length_x + ix] - h_W[(iz - 4) * h_VW.length_x + ix]);

            h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] =
                h_Bx[ix - h_VW.leftborder] * (h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + dVx / Pa.dx) - dVx / Pa.dx;
            h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] =
                h_Bz[iz - h_VW.topborder] * (h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + dWz / Pa.dz) - dWz / Pa.dz;
        }
    }
}

/*------------------------------------------------------------------------
function: StepU

Information:
一步更新波场U

Inputs:
Pa：正演参数
h_U_next： 下一时刻的波场U
h_U_now： 当前时刻的波场U
h_U_past： 前一时刻的波场U
h_W：当前时刻的波场W
h_V：当前时刻的波场V
h_PHIx_V_x: 横向卷积项
h_PHIz_W_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void StepU(AFDPU2D Pa,
            float *h_U_next,
            float *h_U_now,//
            float *h_U_past,//
            float *h_V,
            float *h_W,
            float *h_PHIx_V_x,//
            float *h_PHIz_W_z,//
            float *h_Bx,
            float *h_Bz,
            float *h_Vp,
            const Partition &pt, int it)//
{
    H_Border h_VW = pt.geth_VW();
    H_Border temph_U = pt.geth_U();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

//        MPI_Barrier(MPI_COMM_WORLD);
//        cout << "rank = " << pt.getrank() << endl;
//        cout << gap_min_x + h_VW.leftborder << endl;
//        cout << block_x + h_VW.leftborder + gap_max_x << endl;
//        MPI_Barrier(MPI_COMM_WORLD);

    float dVx = 0.0f;
    float dWz = 0.0f;



    // 交错网格有限差分
    for (uint iz = 0 + gap_min_z + h_VW.topborder; iz < block_z + h_VW.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + h_VW.leftborder; ix < block_x + h_VW.leftborder + gap_max_x; ix++)
        {
            dVx = C1_4 * (h_V[iz * h_VW.length_x + ix] - h_V[iz * h_VW.length_x + ix - 1])
                + C2_4 * (h_V[iz * h_VW.length_x + ix + 1] - h_V[iz * h_VW.length_x + ix - 2])
                + C3_4 * (h_V[iz * h_VW.length_x + ix + 2] - h_V[iz * h_VW.length_x + ix - 3])
                + C4_4 * (h_V[iz * h_VW.length_x + ix + 3] - h_V[iz * h_VW.length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * h_VW.length_x + ix] - h_W[(iz - 1) * h_VW.length_x + ix])
                + C2_4 * (h_W[(iz + 1) * h_VW.length_x + ix] - h_W[(iz - 2) * h_VW.length_x + ix])
                + C3_4 * (h_W[(iz + 2) * h_VW.length_x + ix] - h_W[(iz - 3) * h_VW.length_x + ix])
                + C4_4 * (h_W[(iz + 3) * h_VW.length_x + ix] - h_W[(iz - 4) * h_VW.length_x + ix]);

            if(pt.getrank() == 1)
            if(it == 45 && iz == 53 + h_VW.topborder && ix == h_VW.leftborder)
            {
                //cout << dWz << " "<< 8.40779e-45 / Pa.dz << endl;

//                cout << h_W[iz * h_VW.length_x + ix]
//                     << " " << h_W[(iz - 1) * h_VW.length_x + ix]
//                     << " " << h_W[(iz + 1) * h_VW.length_x + ix]
//                     << " " << h_W[(iz - 2) * h_VW.length_x + ix]
//                     << " " << h_W[(iz + 2) * h_VW.length_x + ix]
//                     << " " << h_W[(iz - 3) * h_VW.length_x + ix]
//                     << " " << h_W[(iz + 3) * h_VW.length_x + ix]
//                     << " " << h_W[(iz - 4) * h_VW.length_x + ix]
//                     << " " << dWz << endl;

                //cout << dVx / Pa.dx << " " << h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << " " << dWz / Pa.dz << h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << endl;
//                cout << (Pa.dt * Pa.dt * h_Vp[(iz - h_VW.topborder + temph_Vp.topborder) * temph_Vp.length_x + ix - h_VW.leftborder + temph_Vp.leftborder]) << " "
//                     << (dVx / Pa.dx + h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + dWz / Pa.dz + h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder]) << " "
//                     << h_U_past[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << " "
//                     << h_U_now[(iz - h_VW.topborder + temph_U.topborder) * temph_U.length_x + ix - h_VW.leftborder + temph_U.leftborder] << endl;

                //cout << dVx << " " << h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << " " << dWz << " " << h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << endl;
//                cout << h_V[iz * h_VW.length_x + ix]
//                     << " " << h_V[iz * h_VW.length_x + ix - 1]
//                     << " " << h_V[iz * h_VW.length_x + ix + 1]//
//                     << " " << h_V[iz * h_VW.length_x + ix - 2]
//                     << " " << h_V[iz * h_VW.length_x + ix + 2]
//                     << " " << h_V[iz * h_VW.length_x + ix - 3]
//                     << " " << h_V[iz * h_VW.length_x + ix + 3]
//                     << " " << h_V[iz * h_VW.length_x + ix - 4] << endl << dVx << endl;
            }


            h_U_next[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] = (Pa.dt * Pa.dt * h_Vp[(iz - h_VW.topborder + temph_Vp.topborder) * temph_Vp.length_x + ix - h_VW.leftborder + temph_Vp.leftborder])
                * (dVx / Pa.dx + h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + dWz / Pa.dz + h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder])
                - h_U_past[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + 2.0f * h_U_now[(iz - h_VW.topborder + temph_U.topborder) * temph_U.length_x + ix - h_VW.leftborder + temph_U.leftborder];
//            if(pt.getrank() == 1)
//            if(it == 45 && iz == 53 + h_VW.topborder && ix == h_VW.leftborder)
//                cout << h_U_next[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] << endl;
        }
    }


//    if(pt.getrank() == 2)
//    if(it == 52)
//    for(int i = 0; i < block_z; ++i)
//        for(int j = 0; j < block_x; ++j)
//            if(*(h_U_next + i * block_x + j) )
//                cout << i << " " << j << " " << *(h_U_next + i * block_x + j) << endl;
}

/*------------------------------------------------------------------------
function: StepVW

Information:
一步更新波场V和W

Inputs:
Pa：正演参数
h_U： 当前时刻的波场U
h_W：当前时刻的波场W
h_V：当前时刻的波场V
h_PHIx_U_x: 横向卷积项
h_PHIz_U_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
------------------------------------------------------------------------*/
void StepVW(AFDPU2D Pa,
            float *h_U,
            float *h_V,
            float *h_W,
            float *h_PHIx_U_x,
            float *h_PHIz_U_z,
            float *h_Bx,
            float *h_Bz,
            const Partition &pt,int it)
{
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }


    float dUx = 0.0f;
    float dUz = 0.0f;

//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "rank = " << pt.getrank() << endl;
//    cout << gap_min_x + temph_U.leftborder << endl;
//    cout << block_x + temph_U.leftborder + gap_max_x << endl;
//    MPI_Barrier(MPI_COMM_WORLD);

//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "rank = " << pt.getrank() << endl;
//    cout << min_x << max_x << endl;
//    cout << min_z << max_z << endl;
//    MPI_Barrier(MPI_COMM_WORLD);

    // 交错网格有限差分
    for (uint iz = 0 + gap_min_z + temph_U.topborder; iz < block_z + temph_U.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + temph_U.leftborder; ix < block_x + temph_U.leftborder + gap_max_x; ix++)
        {
            dUx = C1_4 * (h_U[iz * temph_U.length_x + ix + 1] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[iz * temph_U.length_x + ix + 2] - h_U[iz * temph_U.length_x + ix - 1])
                + C3_4 * (h_U[iz * temph_U.length_x + ix + 3] - h_U[iz * temph_U.length_x + ix - 2])
                + C4_4 * (h_U[iz * temph_U.length_x + ix + 4] - h_U[iz * temph_U.length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * temph_U.length_x + ix] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[(iz + 2) * temph_U.length_x + ix] - h_U[(iz - 1) * temph_U.length_x + ix])
                + C3_4 * (h_U[(iz + 3) * temph_U.length_x + ix] - h_U[(iz - 2) * temph_U.length_x + ix])
                + C4_4 * (h_U[(iz + 4) * temph_U.length_x + ix] - h_U[(iz - 3) * temph_U.length_x + ix]);


            h_V[(iz - temph_U.topborder + h_VW.topborder) * h_VW.length_x + ix - temph_U.leftborder + h_VW.leftborder] = dUx / Pa.dx + h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder];
            h_W[(iz - temph_U.topborder + h_VW.topborder) * h_VW.length_x + ix - temph_U.leftborder + h_VW.leftborder] = dUz / Pa.dz + h_PHIz_U_z[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder];
            if(pt.getrank() == 1 && it == 45 && iz == 53 + temph_U.topborder && ix == 0 + temph_U.leftborder)
            {
//                cout << h_U[(iz + 1) * temph_U.length_x + ix] << " "
//                     << h_U[iz * temph_U.length_x + ix] << " "
//                     << h_U[(iz + 2) * temph_U.length_x + ix] << " "
//                     << h_U[(iz - 1) * temph_U.length_x + ix] << " "
//                     << h_U[(iz + 3) * temph_U.length_x + ix] << " "
//                     << h_U[(iz - 2) * temph_U.length_x + ix] << " "
//                     << h_U[(iz + 4) * temph_U.length_x + ix] << " "
//                     << h_U[(iz - 3) * temph_U.length_x + ix] << " " << dUz << endl;
            }

            if(pt.getrank() == 1 && it == 45 && iz == 53 + temph_U.topborder && ix == 1 + temph_U.leftborder)
            {
//                cout << h_U[iz * temph_U.length_x + ix + 1]
//                     << " " << h_U[iz * temph_U.length_x + ix]
//                     << " " << h_U[iz * temph_U.length_x + ix + 2]
//                     << " " << h_U[iz * temph_U.length_x + ix - 1]
//                     << " " << h_U[iz * temph_U.length_x + ix + 3]
//                     << " " << h_U[iz * temph_U.length_x + ix - 2]
//                     << " " << h_U[iz * temph_U.length_x + ix + 4]
//                     << " " << h_U[iz * temph_U.length_x + ix - 3] << endl << dUx << endl;
//                cout << h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder] << endl;
            }
        }
    }
}

/*------------------------------------------------------------------------
function: StepShotGather

Information:
一步记录检波器位置的波场值

Inputs:
Pa：正演参数
h_U：当前时刻的波场U
h_SG: 炮集
h_re：检波器位置
nstep：循环到第几个时刻
rn：检波器个数
------------------------------------------------------------------------*/
void StepShotGather(AFDPU2D Pa,
                    float *h_U,
                    float *h_SG,
                    uint nstep,
                    const Partition &pt)
{
    uint id = 0;

    H_Border temph_U = pt.geth_U();

    uint rn = pt.getRL_num();
    vector<pair<uint, uint>> vec = pt.getRL();
    vector<pair<uint, uint>>::iterator begin = vec.begin();
    pair<uint, uint> temp;
    for(uint ir = 0; ir < rn; ++ir)
    {
        temp = *(begin + ir);
        id = (temp.second + temph_U.topborder) * temph_U.length_x + temp.first + temph_U.leftborder;
        h_SG[ir * Pa.Nt + nstep] = h_U[id];
        //cout << h_U[id];
    }
}

/*------------------------------------------------------------------------
function: StepRecordU

Information:
一步记录有效边界存储策略中需要记录的波场值

Inputs:
Pa：正演参数
h_Coord：需要记录的波场位置的坐标
RecNum: 要记录的波场的个数
h_U：波场
h_RU: 记录的波场
------------------------------------------------------------------------*/
void StepRecordU(AFDPU2D Pa,
                uint *h_Coord,
                uint RecNum,
                float *h_U,
                float *h_RU, const Partition &pt)
{
    uint n = pt.get_h_Coord_num();

    H_Coord *coord = pt.get_h_Coord();
    int before = 0;
    for(int i = 0; i < n; ++i)
    {
        int length_x = coord[i].getlength_x();
        int length_z = coord[i].getlength_z();

        for (uint iz = 0; iz < length_z; iz++)
        {
            for (uint ix = 0; ix < length_x; ix++)
            {
                h_RU[before + iz * length_x + ix] = h_U[h_Coord[before + iz * length_x + ix]];
            }
        }

        before += length_x * length_z;
    }
}

/*------------------------------------------------------------------------
function: StepReplaceU

Information:
一步替换有效边界内的波场

Inputs:
Pa：正演参数
h_Coord：需要记录的波场位置的坐标
RecNum: 要记录的波场的个数
h_U：波场
h_RU: 记录的波场
------------------------------------------------------------------------*/
void StepReplaceU(AFDPU2D Pa,
                uint *h_Coord,
                uint RecNum,
                float *h_U,
                float *h_RU, const Partition& pt)
{
    uint n = pt.get_h_Coord_num();

    H_Coord *coord = pt.get_h_Coord();
    int before = 0;

    for(int i = 0; i < n; ++i)
    {
        int length_x = coord[i].getlength_x();
        int length_z = coord[i].getlength_z();

        for (uint iz = 0; iz < length_z; iz++)
        {
            for (uint ix = 0; ix < length_x; ix++)
            {
                h_U[h_Coord[before + iz * length_x + ix]] = h_RU[before + iz * length_x + ix];
            }
        }

        before += length_x * length_z;
    }
}

/*------------------------------------------------------------------------
function: StepRtVW

Information:
一步逆时间更新波场V和W

Inputs:
Pa：正演参数
h_U： 当前时刻的波场U
h_W：当前时刻的波场W
h_V：当前时刻的波场V
h_PHIx_U_x: 横向卷积项
h_PHIz_U_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
------------------------------------------------------------------------*/
void StepRtVW(AFDPU2D Pa,
                float *h_U,
                float *h_V,
                float *h_W,
                float *h_PHIx_U_x,
                float *h_PHIz_U_z,
                const Partition &pt)
{
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();


    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;

    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    float dUx = 0.0f;
    float dUz = 0.0f;

    // 交错网格有限差分
    for (uint iz = 0 + gap_min_z + temph_U.topborder; iz < block_z + temph_U.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + temph_U.leftborder; ix < block_x + temph_U.leftborder + gap_max_x; ix++)
        {
            dUx = C1_4 * (h_U[iz * temph_U.length_x + ix + 1] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[iz * temph_U.length_x + ix + 2] - h_U[iz * temph_U.length_x + ix - 1])
                + C3_4 * (h_U[iz * temph_U.length_x + ix + 3] - h_U[iz * temph_U.length_x + ix - 2])
                + C4_4 * (h_U[iz * temph_U.length_x + ix + 4] - h_U[iz * temph_U.length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * temph_U.length_x + ix] - h_U[iz * temph_U.length_x + ix])
                + C2_4 * (h_U[(iz + 2) * temph_U.length_x + ix] - h_U[(iz - 1) * temph_U.length_x + ix])
                + C3_4 * (h_U[(iz + 3) * temph_U.length_x + ix] - h_U[(iz - 2) * temph_U.length_x + ix])
                + C4_4 * (h_U[(iz + 4) * temph_U.length_x + ix] - h_U[(iz - 3) * temph_U.length_x + ix]);

            h_V[(iz - temph_U.topborder + h_VW.topborder) * h_VW.length_x + ix - temph_U.leftborder + h_VW.leftborder] = dUx / Pa.dx + h_PHIx_U_x[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder];
            h_W[(iz - temph_U.topborder + h_VW.topborder) * h_VW.length_x + ix - temph_U.leftborder + h_VW.leftborder] = dUz / Pa.dz + h_PHIz_U_z[(iz - temph_U.topborder) * block_x + ix - temph_U.leftborder];
        }
    }
}

/*------------------------------------------------------------------------
function: StepRtU

Information:
一步逆时间更新波场U

Inputs:
Pa：正演参数
h_U_next： 下一时刻的波场U
h_U_now： 当前时刻的波场U
h_U_past： 前一时刻的波场U
h_W：当前时刻的波场W
h_V：当前时刻的波场V
h_PHIx_V_x: 横向卷积项
h_PHIz_W_z: 纵向卷积项
h_Bx: 横向NPML参数
h_Bz：纵向NPML参数
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void StepRtU(AFDPU2D Pa,
            float *h_U_next,
            float *h_U_now,
            float *h_U_past,
            float *h_V,
            float *h_W,
            float *h_PHIx_V_x,
            float *h_PHIz_W_z,
            float *h_Vp,
             const Partition &pt)
{
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();
    H_Border temph_Vp = pt.geth_Vp();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    int gap_min_x = 0, gap_min_z = 0, gap_max_x = 0, gap_max_z = 0;
    if(min_x < 4)
    {
        gap_min_x = 4 - min_x;
    }
    if(min_z < 4)
    {
        gap_min_z = 4 - min_z;
    }
    if(max_x > totallength_x - 5)
    {
        gap_max_x = totallength_x - 5 - max_x;
    }
    if(max_z > totallength_z - 5)
    {
        gap_max_z = totallength_z - 5 - max_z;
    }

    float dVx = 0.0f;
    float dWz = 0.0f;

    // 交错网格有限差分
    for (uint iz = 0 + gap_min_z + h_VW.topborder; iz < block_z + h_VW.topborder + gap_max_z; ++iz)
    {
        for (uint ix = 0 + gap_min_x + h_VW.leftborder; ix < block_x + h_VW.leftborder + gap_max_x; ix++)
        {
            dVx = C1_4 * (h_V[iz * h_VW.length_x + ix] - h_V[iz * h_VW.length_x + ix - 1])
                + C2_4 * (h_V[iz * h_VW.length_x + ix + 1] - h_V[iz * h_VW.length_x + ix - 2])
                + C3_4 * (h_V[iz * h_VW.length_x + ix + 2] - h_V[iz * h_VW.length_x + ix - 3])
                + C4_4 * (h_V[iz * h_VW.length_x + ix + 3] - h_V[iz * h_VW.length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * h_VW.length_x + ix] - h_W[(iz - 1) * h_VW.length_x + ix])
                + C2_4 * (h_W[(iz + 1) * h_VW.length_x + ix] - h_W[(iz - 2) * h_VW.length_x + ix])
                + C3_4 * (h_W[(iz + 2) * h_VW.length_x + ix] - h_W[(iz - 3) * h_VW.length_x + ix])
                + C4_4 * (h_W[(iz + 3) * h_VW.length_x + ix] - h_W[(iz - 4) * h_VW.length_x + ix]);

            h_U_past[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] = (Pa.dt * Pa.dt * h_Vp[(iz - h_VW.topborder + temph_Vp.topborder) * temph_Vp.length_x + ix - h_VW.leftborder + temph_Vp.leftborder])
                * (dVx / Pa.dx + h_PHIx_V_x[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + dWz / Pa.dz + h_PHIz_W_z[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder])
                - h_U_next[(iz - h_VW.topborder) * block_x + ix - h_VW.leftborder] + 2.0f * h_U_now[(iz - h_VW.topborder + temph_U.topborder) * temph_U.length_x + ix - h_VW.leftborder + temph_U.leftborder];
        }
    }
}

/*------------------------------------------------------------------------
function: StepCal2Der

Information:
一步计算正传波场对时间的2阶导

Inputs:
Pa：正演参数
h_U_next： 下一时刻的波场U
h_U_now： 当前时刻的波场U
h_U_past： 前一时刻的波场U
h_U_Der：波场对时间的2阶导
------------------------------------------------------------------------*/
void StepCal2Der(AFDPU2D Pa,
                float *h_U_past,
                float *h_U_now,
                float *h_U_next,
                float *h_U_Der,
                const Partition &pt)
{
    uint id1 = 0, id2 = 0;

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_z = pt.getinteriorLength_z();
    uint length_x = pt.getblockLength_x();

    H_Border temph_U = pt.geth_U();

    uint gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    uint gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    for (uint iz = 0; iz < interiorlength_z; iz++)
    {
        for (uint ix = 0; ix < interiorlength_x; ix++)
        {
            id1 = iz * interiorlength_x + ix;
            id2 = (iz + gap_z) * length_x + ix + gap_x;
            h_U_Der[id1] =
                (h_U_next[id2] - 2.0f * h_U_now[(temph_U.topborder + iz + gap_z) * temph_U.length_x + temph_U.leftborder + ix + gap_x] + h_U_past[id2])
                / (Pa.dt * Pa.dt);
        }
    }
}

/*------------------------------------------------------------------------
function: StepCalGrad

Information:
一步求取梯度

Inputs:
Pa：正演参数
h_Grad: 梯度
h_U_r： 反传波场
h_U_Der：正传波场对时间的2阶导
------------------------------------------------------------------------*/
void StepCalGrad(AFDPU2D Pa,
                float *h_Grad,
                float *h_U_Der,
                float *h_U_r,
                const Partition &pt)
{
    uint id1 = 0, id2 = 0;

    H_Border temph_U = pt.geth_U();

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_z = pt.getinteriorLength_z();
    uint length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();


//    ofstream fout("h_U_Der.txt");
//    ofstream fout1("h_U_r.txt");
    //ofstream fout2("")

    for (uint iz = 0; iz < interiorlength_z; iz++)
    {
        for (uint ix = 0; ix < interiorlength_x; ix++)
        {
            id1 = iz * interiorlength_x + ix;
            //id2 = (iz + gap_z) * length_x + ix + gap_x;
            h_Grad[id1] += -1.0f * h_U_Der[id1] * h_U_r[(iz + gap_z + temph_U.topborder) * temph_U.length_x + ix + gap_x + temph_U.leftborder];


//            if(pt.getrank() == 1)
//            {
//                fout << h_U_Der[id1] << " ";
//                fout1 << h_U_r[(iz + gap_z + temph_U.topborder) * temph_U.length_x + ix + gap_x + temph_U.leftborder] << " ";
//            }
        }
    }
}

/*------------------------------------------------------------------------
function: StepResidual

Information:
求取残差波场

Inputs:
Pa：正演参数
h_TrueWF: 观测波场
h_CurrWF：正演波场
h_ResWF：残差波场
Rn: 检波器个数
------------------------------------------------------------------------*/
void StepResidual(AFDPU2D Pa,
                float *h_TrueWF,
                float *h_CurrWF,
                float *h_ResWF,
                int RL_num)//yi ding de fan wei nei zhi you yi ding geshu de jianboqi
{
    uint id = 0;
    for (uint n = 0; n < RL_num; n++)
    {
        for (uint m = 0; m < Pa.Nt; m++)
        {
            id = n * Pa.Nt + m;
            h_ResWF[id] = h_TrueWF[id] - h_CurrWF[id];
            //cout << h_TrueWF[id] << "  " << h_CurrWF[id];
        }
    }
}

/*------------------------------------------------------------------------
function: AddResidual

Information:
残差反传加震源

Inputs:
Pa：正演参数
h_ResWF：残差波场
h_U_r: 反传波场
h_re: 检波器位置
nstep: 时刻
Rn: 检波器个数
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void AddResidual(AFDPU2D Pa,
                float *h_ResWf,
                float *h_U_r,
                uint nstep,
                float *h_Vp,
                const Partition &pt)
{
    H_Border temph_Vp =pt.geth_Vp();

    uint block_x = pt.getblockLength_x();

    vector<pair<uint, uint>> vec = pt.getRL();
    pair<uint, uint> temp;
    vector<pair<uint, uint>>::iterator begin = vec.begin();
    for(uint ir = 0; ir < vec.size(); ++ir)
    {
        temp = *(begin + ir);
        h_U_r[temp.second * block_x + temp.first] += -1.0f * h_ResWf[ir * Pa.Nt + nstep] * (Pa.dt * Pa.dt * h_Vp[(temp.second + temph_Vp.topborder) * temph_Vp.length_x + temp.first + temph_Vp.leftborder]);
    }
//    for (uint ir = 0; ir < vec.size(); ir++)
//    {
//        temp = vec.back();
//        vec.pop_back();
//        h_U_r[temp.second * block_x + temp.first] += -1.0f * h_ResWf[ir * Pa.Nt + nstep] * (Pa.dt * Pa.dt * h_Vp[(temp.second + temph_Vp.topborder) * temph_Vp.length_x + temp.first + temph_Vp.leftborder]);
//    }
}

/*------------------------------------------------------------------------
function: PostProcessGrad

Information:
梯度后处理

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
------------------------------------------------------------------------*/
void PostProcessGrad(AFDPU2D Pa,
                    float *h_Grad,
                    float *h_Vp,
                    const Partition& pt)
{
    uint id1 = 0, id2 = 0;

    H_Border temph_Vp =pt.geth_Vp();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_z = pt.getinteriorLength_z();
    uint length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    uint interiormin_z = pt.getinteriormin_z();
    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint ix = 0; ix < interiorLength_x; ix++)
        {
            id1 = iz * interiorLength_x + ix;
            id2 = (iz + gap_z + temph_Vp.topborder) * temph_Vp.length_x + ix + gap_x + temph_Vp.leftborder;

            h_Grad[id1] = (powf((float)(interiormin_z + iz - Pa.PMLz), PowerGrad)
                / (h_Vp[id2] * sqrtf(h_Vp[id2]))) * h_Grad[id1];

            if (interiormin_z + iz < Pa.PMLz + 5)
            {
                h_Grad[id1] = 0.0f;
            }
        }
    }
}

/*------------------------------------------------------------------------
function: UpdateVp

Information:
更新内部网格内的速度

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
e: 步长
------------------------------------------------------------------------*/
void UpdateVp(AFDPU2D Pa,
            float *h_Vp,
            float *h_Grad,
            float e,
             const Partition &pt)
{
    uint id1 = 0, id2 = 0;

    H_Border temph_Vp =pt.geth_Vp();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_z = pt.getinteriorLength_z();
    uint length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint ix = 0; ix < interiorLength_x; ix++)
        {
            id1 = iz * interiorLength_x + ix;
            id2 = (iz + gap_z + temph_Vp.topborder) * temph_Vp.length_x + ix + gap_x + temph_Vp.leftborder;

            h_Vp[id2] = powf((e * h_Grad[id1] + sqrtf(fabs(h_Vp[id2]))), 2.0f);
        }
    }
}

/*------------------------------------------------------------------------
function: UpdateVpPML

Information:
更新PML网格内的速度

Inputs:
Pa：正演参数
h_Grad：梯度
h_Vp: 当前速度的平方
e: 步长
------------------------------------------------------------------------*/
void UpdateVpPML(AFDPU2D Pa,
                float *h_Vp,
                float *h_Grad,
                float e,
                const Partition &pt)
{
    uint id = 0;
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    H_Border temph_Vp = pt.geth_Vp();

    uint min_x = pt.getindexmin_x();
    uint min_z = pt.getindexmin_z();
    uint max_x = pt.getindexmax_x();
    uint max_z = pt.getindexmax_z();

    for (uint iz = 0; iz < block_z; iz++)
    {
        for (uint ix = 0; ix < block_x; ix++)
        {

            id = (iz + temph_Vp.topborder) * temph_Vp.length_x + ix + temph_Vp.leftborder;
            int temp_z = iz + min_z;
            int temp_x = ix + min_x;
            // 纵向
            if (temp_z < Pa.PMLz)
            {
                if(!temph_Vp.bottomborder)
                    h_Vp[id] = h_Vp[(Pa.PMLz - min_z + temph_Vp.topborder) * temph_Vp.length_x + temph_Vp.leftborder + ix];
                else
                {
                    h_Vp[id] = h_Vp[(block_z + temph_Vp.topborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];///tong yi ge fang xiang shang zhi neng youyige border
                }
            }
            if (temp_z > Pa.Nz + Pa.PMLz - 1)
            {
                if(!temph_Vp.topborder)
                    h_Vp[id] = h_Vp[(Pa.Nz + Pa.PMLz - min_z - 1 + temph_Vp.topborder) * temph_Vp.length_x + ix + temph_Vp.leftborder];
                else
                {
                    h_Vp[id] = h_Vp[ix + temph_Vp.leftborder];
                }
            }
            // 横向
            if (temp_x < Pa.PMLx)
            {
                if(!temph_Vp.rightborder)
                    h_Vp[id] = h_Vp[temph_Vp.length_x * (iz + temph_Vp.topborder) + Pa.PMLx - min_x + temph_Vp.leftborder];
                else
                    h_Vp[id] = h_Vp[temph_Vp.length_x * (iz + temph_Vp.topborder) + temph_Vp.leftborder + block_x];////?
            }
            if (temp_x > Pa.Nx + Pa.PMLx - 1)
            {
                if(!temph_Vp.leftborder)
                    h_Vp[id] = h_Vp[(iz + temph_Vp.topborder) * temph_Vp.length_x + Pa.Nx + Pa.PMLx - 1 - min_x + temph_Vp.leftborder];
                else
                    h_Vp[id] = h_Vp[temph_Vp.length_x * (iz + temph_Vp.topborder)];
            }
        }
    }
}

/*------------------------------------------------------------------------
function: SetCoord

Information:
给有效边界存储策略需要存储的波场坐标赋值

Inputs:
Pa：正演参数
h_Coord: 坐标
------------------------------------------------------------------------*/
void SetCoord(AFDPU2D *Pa,
            uint *h_Coord,
              const Partition &pt)
{
    uint nnx = Pa->Nx + 2 * Pa->PMLx;
    H_Border temph_U = pt.geth_U();

    uint n = pt.get_h_Coord_num();
    H_Coord *coord = pt.get_h_Coord();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    // set Coord
    int before = 0;
    for(int i = 0; i < n; ++i)
    {
        uint length_x = coord[i].getlength_x();
        uint length_z = coord[i].getlength_z();

        uint min_x = coord[i].getindexmin_x();
        uint min_z = coord[i].getindexmin_z();

        for (uint iz = 0; iz < length_z; iz++)
        {
            for (uint ix = 0; ix < length_x; ix++)
            {

                h_Coord[before + iz * length_x + ix] = (iz + min_z - indexmin_z + temph_U.topborder) * temph_U.length_x + temph_U.leftborder + ix + min_x - indexmin_x;////?
            }
        }
        before += length_x * length_z;
    }
}

/*------------------------------------------------------------------------
function: MatAdd

Information:
矩阵加法

Inputs:
Mat: 输出矩阵
Mat1：输入矩阵1
Mat2：输入矩阵2
row：行数
col：列数
coeff1：矩阵1的系数
coeff2：矩阵2的系数
flag：加法或减法
------------------------------------------------------------------------*/
void MatAdd(float *Mat,
            float *Mat1,
            float *Mat2,
            uint row,
            uint col,
            float coeff1,
            float coeff2,
            uint flag)
{
    uint id = 0;
    for (uint iz = 0; iz < row; iz++)
    {
        for (uint ix = 0; ix < col; ix++)
        {
            id = iz * col + ix;
            if (flag == 1)
            {
                Mat[id] = Mat1[id] * coeff1
                    + Mat2[id] * coeff2;
            }
            if (flag == 0)
            {
                Mat[id] = Mat1[id] * coeff1
                    - Mat2[id] * coeff2;
            }
        }
    }
}

/*------------------------------------------------------------------------
function: MatMul

Information:
矩阵点对点乘除法

Inputs:
Mat: 输出矩阵
Mat1：输入矩阵1
Mat2：输入矩阵2
row：行数
col：列数
coeff1：矩阵1的系数
coeff2：矩阵2的系数
flag：加法或减法
------------------------------------------------------------------------*/
void MatMul(float *Mat,
            float *Mat1,
            float *Mat2,
            uint row,
            uint col,
            float coeff1,
            float coeff2,
            uint flag)
{
    uint id = 0;
    for (uint iz = 0; iz < row; iz++)
    {
        for (uint ix = 0; ix < col; ix++)
        {
            id = iz * col + ix;
            if (flag == 1)
            {
                Mat[id] = Mat1[id] * coeff1
                    * Mat2[id] * coeff2;
            }
            if (flag == 0)
            {
                Mat[id] = (Mat1[id] * coeff1)
                    / (Mat2[id] * coeff2);
            }
        }
    }
}

/*------------------------------------------------------------------------
function: CalTrueWF

Information:
求取观测波场

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
------------------------------------------------------------------------*/
void CalTrueWF(AFDPU2D Pa,
            IP *ip,
            CPUVs *plan,
            float *sgs_t,
            const Partition &pt)//plan  sgs_t
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_z = pt.getindexmax_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    uint pos_x = pt.getblockPosition_x();
    uint pos_z = pt.getblockPosition_z();

    float Wavelet = 0.0f;

    int rank = pt.getrank();

    // 给速度赋值
    memset((void *)plan->h_Vp, 0, sizeof(float) * temph_Vp.length_z * temph_Vp.length_x);// h_Vp 正演中使用的速度


    for(int i = 0; i < block_z; ++i)
    {
        memcpy(plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_x + temph_Vp.leftborder, ip->TrueVp + i * block_x, block_x * sizeof(float));
    }
    //cout << *(ip->TrueVp + (52 + temph_Vp.topborder)* temph_Vp.length_x +  260 + temph_Vp.leftborder);
//    if(pt.getrank() == 1)
//        for(int i = 0; i < block_z * block_x; ++i)
//        cout << *(plan->h_Vp + i);
    uint RL_num = pt.getRL_num();


    for (uint is = 0; is < ip->ShotN; is++)// 反演中的炮数
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == ROOT_ID)
            cout << is << " CalTrueWF" << endl;

        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        //
        if(RL_num)
        {
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            // 给检波器的位置信息赋值
        }

//        MPI_Request request_send_U, request_send_V, request_send_W, request_recv_U, request_recv_V, request_recv_W;
//        MPI_Status status_send, status_recv;

        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)// Pa.Nt正演的时间步数
        {
//            if(pt.getrank() == 1)
//                cout <<it << endl;
//            if(it == 32)
//            {
//                //cout << pt.getindexmin_z() << endl;
//                //cout << *(plan->h_U_now + temph_U.topborder * temph_U.length_x + temph_U.leftborder + 247) << endl;
////                for(int i = 0; i < temph_U.length_z * temph_U.length_x; ++i)
////                    cout << *(plan->h_U_now + i);
//                //MPI_Barrier(MPI_COMM_WORLD);
//                for(int i = temph_U.topborder; i < temph_U.topborder + block_z; ++i)
//                    for(int j = temph_U.leftborder; j < temph_U.leftborder + block_x; ++j)
//                        if(*(plan->h_U_now + i * temph_U.length_x + j) != 0)
//                            cout << i - temph_U.topborder + pt.getindexmin_z() << " " << j - temph_U.leftborder + pt.getindexmin_x() << " " << *(plan->h_U_now + i * temph_U.length_x + j) << endl;
//                        //cout << *(plan->h_U_now + i * temph_U.length_x + j) << " ";
//                //MPI_Barrier(MPI_COMM_WORLD);
//            }


            // 波场时刻转换
            for(int i = 0; i < block_z; ++i)
            {
                memcpy(plan->h_U_past + i * block_x, plan->h_U_now + temph_U.leftborder + (i + temph_U.topborder) * temph_U.length_x, block_x * sizeof(float));
                memcpy(plan->h_U_now + temph_U.leftborder + (i + temph_U.topborder) * temph_U.length_x, plan->h_U_next + i * block_x, block_x * sizeof(float));
            }
//if(rank == ROOT_ID)
//    cout << it << endl;
            dataTransport(plan->h_U_now, pt, STEP_U, it);//send data
//cout << pt.getrank() << endl;
//            if(pos_z % 2)
//            {
//                dataTransport(plan->h_U_now, pt, STEP_U);//send data

//                dataGather(plan->h_U_now, pt, STEP_U);//recv data
//            }
//            else
//            {
//                dataGather(plan->h_U_now, pt, STEP_U);//recv data

//                dataTransport(plan->h_U_now, pt, STEP_U);//send data
//            }


           MPI_Barrier(MPI_COMM_WORLD);
//            if(it == 50 && rank == 1)
//            {
//                //MPI_Barrier(MPI_COMM_WORLD);
//                for(int i = 0; i < temph_U.length_z * temph_U.length_x; ++i)
//                    cout << *(plan->h_U_now + i) << " ";
//                //MPI_Barrier(MPI_COMM_WORLD);
//            }


//              if(it == 50 && pt.getrank() == 0)
//              {//cout << temph_U.length_x << endl;
//                  cout << *(plan->h_U_now + (temph_U.topborder + 117) * temph_U.length_x + 247 + temph_U.leftborder) << endl;
//              }
//            if(pt.getrank() == 1)
//            if(it == Pa.Nt - 1)
//            {
//                //cout << pt.getindexmin_z() << endl;
//                //cout << *(plan->h_U_now + temph_U.topborder * temph_U.length_x + temph_U.leftborder + 247) << endl;
////                for(int i = 0; i < temph_U.length_z * temph_U.length_x; ++i)
////                    cout << *(plan->h_U_now + i);
//                //MPI_Barrier(MPI_COMM_WORLD);
//                for(int i = temph_U.topborder; i < temph_U.topborder + block_z; ++i)
//                    for(int j = temph_U.leftborder; j < temph_U.leftborder + block_x; ++j)
//                        if(*(plan->h_U_now + i * temph_U.length_x + j) != 0)
//                            cout << i - temph_U.topborder + pt.getindexmin_z() << " " << j - temph_U.leftborder + pt.getindexmin_x() << " " << *(plan->h_U_now + i * temph_U.length_x + j) << endl;
//                        //cout << *(plan->h_U_now + i * temph_U.length_x + j) << " ";
//                //MPI_Barrier(MPI_COMM_WORLD);
//            }


//            if(it == 0)
//            {
//                cout << RL_num << endl;
//            }

            if(RL_num)
            {
                // 一步记录炮集
                StepShotGather(Pa, plan->h_U_now, plan->h_TrueWF, it, pt);//h_TrueWF
            }

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);//h_PHIx_U_x   h_PHIz_U_z ... h_U_now

//            if(pt.getrank() == 2)
//            if(it == 52)
//            for(int i = 0; i < block_z; ++i)
//                for(int j = 0; j < block_x; ++j)
//                    if(*(plan->h_PHIx_U_x + i * block_x + j) )
//                        cout << i << " " << j << " " << *(plan->h_PHIx_U_x + i * block_x + j) << endl;


//            if(pt.getrank() == 1)
//            {
//                for(int i = 0; i < block_z; ++i)
//                {
//                    for(int j = 0; j < block_x; ++j)
//                    {
//                        cout << plan->h_PHIx_U_x[i * block_x + j] << " ";
//                    }
//                }
//            }
            //cout << pt.getrank() << endl;

//                        if(it == 0)
//                        {
//                            MPI_Barrier(MPI_COMM_WORLD);
//                            for(int i = 0; i < temph_U.length_z * temph_U.length_x; ++i)
//                                cout << *(plan->h_U_now + i);
//                            MPI_Barrier(MPI_COMM_WORLD);
//                        }

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);//h_V  h_W ... h_U   h_PHIx_U_x   h_PHIz_U_z


//            if(pt.getrank() == 0)
//            if(it == 50)
//            {
//                //cout << block_x << endl;
//                for(int i = h_VW.topborder; i < block_z + h_VW.topborder; ++i)
//                    for(int j = h_VW.leftborder; j < block_x + h_VW.leftborder; ++j)
//                        if(*(plan->h_W + i * h_VW.length_x + j) )
//                            cout << i - h_VW.topborder << " " << j - h_VW.leftborder << " " << *(plan->h_W + i * h_VW.length_x + j) << endl;
//            }

            dataTransport(plan->h_V, pt, STEP_V, it);//send data
            dataTransport(plan->h_W, pt, STEP_W, it);//send data

//            if(rank % 2)
//            {
//                dataTransport(plan->h_V, pt, STEP_V);//send data
//                dataTransport(plan->h_W, pt, STEP_W);//send data
//                dataGather(plan->h_V, pt, STEP_V);//recv_v data
//                dataGather(plan->h_W, pt, STEP_W);//recv_w data
//            }
//            else
//            {
//                dataGather(plan->h_V, pt, STEP_V);//recv_v data
//                dataGather(plan->h_W, pt, STEP_W);//recv_w data
//                dataTransport(plan->h_V, pt, STEP_V);//send data
//                dataTransport(plan->h_W, pt, STEP_W);//send data
//            }

            MPI_Barrier(MPI_COMM_WORLD);
            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);//h_PHIx_V_x    h_PHIz_W_z ... h_V, h_W





            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt, it);//h_U_next ...  h_V h_W

//            for(int i = 0; i < pt.getblockLength_z(); ++i)
//            {
//                for(int j = 0; j < pt.getblockLength_x(); ++j)
//                {
//                    cout << *(plan->h_Vp + j + temph_Vp.leftborder + (i + temph_Vp.topborder) * temph_Vp.length_x);
//                }
//            }

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            //cout << "ppppppp" << rank << " " << pt.getShot_num() << endl;
            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, Wavelet, plan->h_Vp, ip->St[is].s, pt);//h_U_next
            }
//            if(it == Pa.Nt - 1 && pt.getrank() == 0)
//            {
//                //if(pt.getrank() == 0)
//                for(int i = 0; i < block_x * block_z; ++i)
//                    cout << *(plan->h_U_next + i);
//                //cout << Wavelet << endl;
//            }
//            if(pt.getrank() == 1)
//            if(it == 45)
//            for(int i = 0; i < block_z; ++i)
//                for(int j = 0; j < block_x; ++j)
//                    if(*(plan->h_U_next + i * block_x + j) )
//                        cout << i << " " << j + pt.getindexmin_x() << " " << *(plan->h_U_next + i * block_x + j) << endl;
        }

        if(RL_num)
        {
            // 输出炮集
            memcpy(sgs_t + is * (Pa.Nt * RL_num), plan->h_TrueWF, Pa.Nt * RL_num * sizeof(float));/////?/\?
            if(rank == ROOT_ID)
                cout << is << " size " << is * (Pa.Nt * RL_num) << endl;
        }

//        if(is == 0)
//        {



//            cout << "lailai" << endl;
//            if(pt.getrank() == 0)
//            {
//                ofstream fout("sgs_t2.txt");
//                for(int is = 0; is < ip->ShotN; ++is)
//                for(int i = 0; i < RL_num; ++i)
//                {
//                    for(int j = 0; j < Pa.Nt; ++j)
//                        fout << *(sgs_t + is * RL_num * Pa.Nt + i * Pa.Nt + j) << " ";
//                }
//                fout << endl;
//            }
//            if(pt.getrank() == 1)
//            {
//                ofstream fout1("sgs_t3.txt");
//                for(int is = 0; is < ip->ShotN; ++is)
//                for(int i = 0; i < RL_num; ++i)
//                {
//                    for(int j = 0; j < Pa.Nt; ++j)
//                        fout1 << *(sgs_t + is * RL_num * Pa.Nt + i * Pa.Nt + j) << " ";
//                }
//                fout1 << endl;
//            }
//        }

    }

//    if(pt.getrank() == 0)
//    {
//        ofstream fout("sgs_t0.txt");
//        for(int is = 0; is < ip->ShotN; ++is)
//        for(int i = 0; i < RL_num; ++i)
//        {
//            for(int j = 0; j < Pa.Nt; ++j)
//                fout << *(sgs_t + is * RL_num * Pa.Nt + i * Pa.Nt + j) << " ";
//            //
//        }
//        fout << endl;
//        fout.flush();
//        fout.close();
//    }
//    if(pt.getrank() == 1)
//    {
//        ofstream fout("sgs_t1.txt");
//        for(int is = 0; is < ip->ShotN; ++is)
//        for(int i = 0; i < RL_num; ++i)
//        {
//            for(int j = 0; j < Pa.Nt; ++j)
//                fout << *(sgs_t + is * RL_num * Pa.Nt + i * Pa.Nt + j) << " ";
//        }
//        fout << endl;
//        fout.flush();
//        fout.close();
//    }
}

/*------------------------------------------------------------------------
function: CalGrad

Information:
求取梯度

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
sgs_c: 当前正演波场
sgs_r: 残差波场
It：第几次迭代
------------------------------------------------------------------------*/
void CalGrad(AFDPU2D Pa,
            IP *ip,
            CPUVs *plan,
            float *sgs_t,
            float *sgs_c,
            float *sgs_r,
            uint It, const Partition &pt)
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_z = pt.getindexmax_z();

    uint interior_min_z = pt.getinteriormin_z();
    uint interior_min_x = pt.getinteriormin_x();
    
    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_z = pt.getinteriorLength_z();

    float Wavelet = 0.0f;

    int rank = pt.getrank();
    int rank_size = pt.getsize();

    uint RecNum = pt.geth_Coord_length();

    if(interiorLength_z && interiorLength_x)
    {
        // 梯度变量空间清零
        memset((void *)plan->h_Grad, 0, sizeof(float) * interiorLength_z * interiorLength_x);
    }

    // 给速度赋值
    memset((void *)plan->h_Vp, 0, sizeof(float) * temph_Vp.length_z * temph_Vp.length_x);
    for(int i = 0; i < block_z; ++i)
    {
        memcpy(plan->h_Vp + (temph_Vp.topborder + i) * temph_Vp.length_x + temph_Vp.leftborder,	ip->CurrVp + i * block_x, sizeof(float) * block_x);
    }

    // 对炮进行循环
    for (uint is = 0; is < ip->ShotN; is++)
    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_past_r,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_now_r,		0,	sizeof(float) * temph_U.length_z * temph_U.length_x);
        memset((void *)plan->h_U_next_r,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_V_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        memset((void *)plan->h_W_r,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);

        int RL_num = pt.getRL_num();//jian bo qi ge shu
        if(RL_num)
        {
            //memset((void *)plan->h_re,			0,	sizeof(RL) * RL_num);
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_Obj,			0,	sizeof(float) * RL_num * Pa.Nt);
        }

        if(RecNum)
        {
            memset((void *)plan->h_RU,			0,	sizeof(uint) * RecNum * (Pa.Nt - 1));
        }


        if(interiorLength_x && interiorLength_z)
            memset((void *)plan->h_U_Der,		0,	sizeof(float) * interiorLength_z * interiorLength_x);

//        MPI_Request request_send_U, request_send_V, request_send_W, request_send, request_recv_OBJ, request_recv, request_recv_min_RL, request_recv_max_RL;
        MPI_Status status_send, status_recv;

        // 计算正演波场，对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            // 波场时刻转换
            for(int iz = 0; iz < block_z; ++iz)
            {
                memcpy(plan->h_U_past + iz * block_x, plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, block_x * sizeof(float));
                memcpy(plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next + iz * block_x, block_x * sizeof(float));
            }

            dataTransport(plan->h_U_now, pt, STEP_U, it);//send data
            //dataGather(plan->h_U_now, pt, STEP_U);//recv data

            if(RL_num)
            {
                // 一步记录炮集
                StepShotGather(Pa, plan->h_U_now, plan->h_CurrWF, it, pt);//plan->h_CurrWF
            }

            if (it < Pa.Nt - 1 && RecNum)
            {
                StepRecordU(Pa, plan->h_Coord, RecNum, plan->h_U_now,
                    &plan->h_RU[it * RecNum], pt);//plan->h_RU[it * RecNum]////
            }

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);//plan->h_PHIx_U_x, plan->h_PHIz_U_z,

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);//plan->h_V, plan->h_W

            dataTransport(plan->h_V, pt, STEP_V, it);//send data
            dataTransport(plan->h_W, pt, STEP_W, it);//send data
            //dataGather(plan->h_V, pt, STEP_V);//recv_v data
            //dataGather(plan->h_W, pt, STEP_W);//recv_w data

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);//h_PHIx_V_x, h_PHIz_W_z

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt, it);//plan->h_U_next

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, Wavelet, plan->h_Vp, ip->St[is].s, pt);//h_U_next
            }
        }


//        MPI_Barrier(MPI_COMM_WORLD);
//        cout << "xxxxxxxxxx" << pt.getrank() << endl;
        //cout << "ppppppp" << pt.getrank() << endl;

        // 输出炮集
        if(RL_num)
        {
            memcpy(sgs_c + is * (Pa.Nt * RL_num), plan->h_CurrWF, Pa.Nt * RL_num * sizeof(float));
            // 重新读入观测数据
            memcpy(plan->h_TrueWF, sgs_t + is * (Pa.Nt * RL_num), Pa.Nt * RL_num * sizeof(float));

            // 计算残差波场
            StepResidual(Pa, plan->h_TrueWF, plan->h_CurrWF,
                plan->h_ResWF, RL_num);//h_ResWF

            //  输出残差波场
            memcpy(sgs_r + is * (Pa.Nt * RL_num),
                plan->h_ResWF,
                Pa.Nt * RL_num * sizeof(float));
        }



//        if(rank == 0)
//        {
//            ofstream fout("h_CurrWF0.txt");
//            for(uint i = 0; i < RL_num; ++i)
//            {
//                for(uint j = 0; j < Pa.Nt; ++j)
//                {
//                    fout << plan->h_CurrWF[i * Pa.Nt + j] << " ";
//                }
//            }
//            fout << endl;
//            fout.flush();
//            fout.close();
//        }



//        if(rank == 1)
//        {
//            ofstream fout1("h_CurrWF1.txt");
//            cout << "rank=" << rank << " " << RL_num << endl;
//            for(uint i = 0; i < RL_num; ++i)
//            {
//                for(uint j = 0; j < Pa.Nt; ++j)
//                {

//                    fout1 << plan->h_CurrWF[i * Pa.Nt + j] << " ";
//                }
//            }
//            fout1 << endl;
//            fout1.flush();
//            fout1.close();
//        }



        //  残差反传以及正传波场逆时间反推
        for (uint it = Pa.Nt - 1; it > 0; it--)
        {
            // 正传波场逆时间反推
            if (it == Pa.Nt - 1)
            {
                // 消除震源
                Wavelet = Ricker(Pa.f0, it * Pa.dt);
                if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
                {
                    AddSource(Pa, plan->h_U_next, -1.0f * Wavelet, plan->h_Vp, ip->St[is].s, pt);
                }

                if(interiorLength_x && interiorLength_z)
                {
                    // 求取波场对时间的2阶导数
                    StepCal2Der(Pa, plan->h_U_past, plan->h_U_now,
                        plan->h_U_next, plan->h_U_Der, pt);

                }

                // 波场转换
                for(int iz = 0; iz < block_z; ++iz)
                {
                    memcpy(plan->h_U_next + iz * block_x, plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, block_x * sizeof(float));
                    memcpy(plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_past + iz * block_x, block_x * sizeof(float));
                }

            }
            else
            {
                // 消除震源
                Wavelet = Ricker(Pa.f0, it * Pa.dt);
                if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
                {
                    AddSource(Pa, plan->h_U_next, -1.0f * Wavelet, plan->h_Vp, ip->St[is].s, pt);
                }

                if(RecNum)
                {
                    // 替换有效边界内的波场
                    StepReplaceU(Pa, plan->h_Coord, RecNum, plan->h_U_now,
                        &plan->h_RU[it * RecNum], pt);

                }

                dataTransport(plan->h_U_now, pt, STEP_U, it);//send data
                //dataGather(plan->h_U_now, pt, STEP_U);//recv data

                // 逆时间更新V和W
                StepRtVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                    plan->h_PHIx_U_x, plan->h_PHIz_U_z, pt);


                dataTransport(plan->h_V, pt, STEP_V, it);//send data
                dataTransport(plan->h_W, pt, STEP_W, it);//send data
                //MPI_Wait(&request_send_V, &status_send);
                //MPI_Wait(&request_send_W, &status_send);
                //dataGather(plan->h_W, pt, STEP_W);//recv_w data
                //dataGather(plan->h_V, pt, STEP_V);//recv_v data

                // 逆时间更新U
                StepRtU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                    plan->h_V, plan->h_W,plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                    plan->h_Vp, pt);

                if(interiorLength_x && interiorLength_z)
                {
                    // 求取波场对时间的2阶导数
                    StepCal2Der(Pa, plan->h_U_past, plan->h_U_now,
                        plan->h_U_next, plan->h_U_Der, pt);
                }


                // 波场转换
                for(int iz = 0; iz < block_z; ++iz)
                {
                    memcpy(plan->h_U_next + iz * block_x, plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, block_x * sizeof(float));
                    memcpy(plan->h_U_now + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_past + iz * block_x, block_x * sizeof(float));
                }
            }

            //  残差反传
            // 波场时刻转换
            for(int iz = 0; iz < block_z; ++iz)
            {
                memcpy(plan->h_U_past_r + iz * block_x, plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, block_x * sizeof(float));
                memcpy(plan->h_U_now_r + (iz + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next_r + iz * block_x, block_x * sizeof(float));
            }

            dataTransport(plan->h_U_now_r, pt, STEP_U, it);//send data
            MPI_Barrier(MPI_COMM_WORLD);
            //dataGather(plan->h_U_now, pt, STEP_U);//recv data



            if(interiorLength_x && interiorLength_z)
            {
                // 一步求取梯度
                StepCalGrad(Pa, plan->h_Grad, plan->h_U_Der, plan->h_U_now_r, pt);//plan->h_Grad

            }

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now_r, plan->h_PHIx_U_x_r,
                plan->h_PHIz_U_z_r, plan->h_Bx, plan->h_Bz, pt, it);

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now_r, plan->h_V_r, plan->h_W_r,
                plan->h_PHIx_U_x_r, plan->h_PHIz_U_z_r, plan->h_Bx, plan->h_Bz, pt, it);

            dataTransport(plan->h_V_r, pt, STEP_V, it);//send data
            dataTransport(plan->h_W_r, pt, STEP_W, it);//send data
            //dataGather(plan->h_V, pt, STEP_V);//recv_v data
            //dataGather(plan->h_W, pt, STEP_W);//recv_w data

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V_r, plan->h_W_r, plan->h_PHIx_V_x_r,
                plan->h_PHIz_W_z_r, plan->h_Bx, plan->h_Bz, pt);

            MPI_Barrier(MPI_COMM_WORLD);
//            if(rank == ROOT_ID && it == 1)
//                cout << "qqqqqqqq" << it << endl;

            // 一步更新波场U
            StepU(Pa, plan->h_U_next_r, plan->h_U_now_r, plan->h_U_past_r,
                plan->h_V_r, plan->h_W_r, plan->h_PHIx_V_x_r, plan->h_PHIz_W_z_r,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt, it);

            // 加震源
            if(RL_num)
            {
                AddResidual(Pa, plan->h_ResWF, plan->h_U_next_r, it, plan->h_Vp, pt);
            }
        }

//        MPI_Barrier(MPI_COMM_WORLD);
//        //cout << "xxxxxxxxxx" << pt.getrank() << endl;
//        cout << "ppppppp" << pt.getrank() << endl;

        //for(uint i = 0; i < RL)

        uint indexmin_x = pt.getindexmin_x();
        uint indexmax_x = pt.getindexmax_x();

        int max_rank_RL = rank, min_rank_RL = rank;

        if(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx)
        {//cout << rank << endl;
            min_rank_RL = rank;
            if(min_rank_RL != ROOT_ID)
                MPI_Send(&min_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MIN_RL, MPI_COMM_WORLD);
            //cout << "AAAA" << min_rank_RL << " " << max_rank_RL << "AAAA" << endl;
        }
        if(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1)
        {
            max_rank_RL = rank;
            if(max_rank_RL != ROOT_ID)
                MPI_Send(&max_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MAX_RL, MPI_COMM_WORLD);
        }

        if(rank == ROOT_ID)
        {
            if(!(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx))
                MPI_Recv(&min_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MIN_RL, MPI_COMM_WORLD, &status_recv);
            if(!(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1))
                MPI_Recv(&max_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MAX_RL, MPI_COMM_WORLD, &status_recv);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&max_rank_RL, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
        MPI_Bcast(&min_rank_RL, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);

//cout << max_rank_RL << endl;
        float *buf_obj = new float[1];
        // 求取目标函数
        if(RL_num)
        {
            if(rank != min_rank_RL)
            {
                MPI_Recv(buf_obj, 1, MPI_FLOAT, rank - 1, STEP_OBJ, MPI_COMM_WORLD, &status_recv);
                ip->ObjIter[It] = *buf_obj;
                for (uint m = 0; m < RL_num * Pa.Nt; m++)
                {
                    //cout << *ip->ObjIter << endl;
                    ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);//ci shi It = 0// zhe li suan de zhi shi jubu de han shu
                }
                if(rank != max_rank_RL)
                {
                    MPI_Send(&ip->ObjIter[It], 1, MPI_FLOAT, rank + 1, STEP_OBJ, MPI_COMM_WORLD);
                }
            }
            else
            {
                for (uint m = 0; m < RL_num * Pa.Nt; m++)
                {
                    //cout << *ip->ObjIter << endl;
                    ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);//ci shi It = 0// zhe li suan de zhi shi jubu de han shu
                }
                MPI_Send(&ip->ObjIter[It], 1, MPI_FLOAT, rank + 1, STEP_OBJ, MPI_COMM_WORLD);
            }
        }

        delete buf_obj;

//        if(rank >= min_rank_RL && rank <= max_rank_RL)
//cout << *ip->ObjIter << " " << rank << endl;

//        if(rank != ROOT_ID)
//        {
//            MPI_Isend(&ip->ObjIter[It], 1, MPI_FLOAT, ROOT_ID, STEP_OBJ, MPI_COMM_WORLD, &request_send);
//        }
//        else
//        {
//            for(uint i = 0; i < rank_size - 1; ++i)
//            {
//                MPI_Irecv(buf_obj, 1, MPI_FLOAT, MPI_ANY_SOURCE, STEP_OBJ, MPI_COMM_WORLD, &request_recv_OBJ);
//                MPI_Wait(&request_recv_OBJ, &status_recv);
//                //cout << *buf_obj << endl;
//                ip->ObjIter[It] += *buf_obj;
//            }
//        }

        MPI_Barrier(MPI_COMM_WORLD);
//        if(rank == max_rank_RL)
//        {
//            cout << ip->ObjIter[It] << " "  << rank << endl;
//        }
        MPI_Bcast(&ip->ObjIter[It], 1, MPI_FLOAT, max_rank_RL, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == ROOT_ID)
        {
            cout << is << " " << ip->ObjIter[It] << endl;
        }
        //cout << "ppppppp" << pt.getrank() << endl;
    }

//    ofstream fout("GradVp0.txt");
//    if(rank == 0)
//    for(int i = 0; i < interiorLength_z; ++i)
//    {
//        for(int j = 0; j < interiorLength_x; ++j)
//        {
//            fout << plan->h_Grad[i * interiorLength_x + j] << " ";
//        }
//    }
//    fout << endl;
//    fout.flush();
//    fout.close();

//    ofstream fout1("GradVp1.txt");
//    if(rank == 1)
//    for(int i = 0; i < interiorLength_z; ++i)
//    {
//        for(int j = 0; j < interiorLength_x; ++j)
//        {
//            fout1 << plan->h_Grad[i * interiorLength_x + j] << " ";
//        }
//    }
//    fout1 << endl;
//    fout1.flush();
//    fout1.close();

//    ofstream fout2("GradVp2.txt");
//    if(rank == 0)
//    for(int i = 0; i < interiorLength_z; ++i)
//    {
//        for(int j = 0; j < interiorLength_x; ++j)
//        {
//            fout2 << plan->h_Grad[i * interiorLength_x + j] << " ";
//        }
//    }
//    fout2 << endl;
//    fout2.flush();
//    fout2.close();

//    ofstream fout3("GradVp3.txt");
//    if(rank == 1)
//    for(int i = 0; i < interiorLength_z; ++i)
//    {
//        for(int j = 0; j < interiorLength_x; ++j)
//        {
//            fout3 << plan->h_Grad[i * interiorLength_x + j] << " ";
//        }
//    }
//    fout3 << endl;
//    fout3.flush();
//    fout3.close();

    // 输出梯度
    for(uint i = 0; i < interiorLength_z; ++i)
    {
        memcpy(ip->GradVp + i * interiorLength_x, plan->h_Grad + i * interiorLength_x, interiorLength_x * sizeof(float));
    }
}

/*------------------------------------------------------------------------
function: CalStepLength

Information:
求取步长

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
sgs_t: 观测波场
sgs_c: 当前正演波场
e: 试探步长
------------------------------------------------------------------------*/
void CalStepLength(AFDPU2D Pa,
                IP *ip,
                CPUVs *plan,
                float *sgs_t,
                float *sgs_c,
                float e,
                  const Partition &pt)
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint interiorLength_x = pt.getinteriorLength_x();
    uint interiorLength_z = pt.getinteriorLength_z();

    uint interior_min_x = pt.getinteriormin_x();
    uint interior_min_z = pt.getinteriormin_z();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmin_z = pt.getindexmin_z();

    uint indexmax_x = pt.getindexmax_x();
    uint indexmax_z = pt.getindexmax_z();

    H_Border temph_Vp = pt.geth_Vp();
    H_Border temph_U = pt.geth_U();
    H_Border h_VW = pt.geth_VW();

    MPI_Request request_send_U, request_send_V, request_send_W, request_recv_U, request_recv_V, request_recv_W, request_send, request_recv_MAX;
    MPI_Status status_send, status_recv;

    int rank = pt.getrank();

    uint sumblock_x = pt.getsumBlock_x();
    //uint sumblock_z = pt.getsumBlock_z();
    uint rank_size = pt.getsize();

    uint RL_num = pt.getRL_num();
    if(RL_num)
    {
        memset((void *)plan->h_SumResTrial,	0,	sizeof(float) * RL_num * Pa.Nt);
        memset((void *)plan->h_SumResCurr,	0,	sizeof(float) * RL_num * Pa.Nt);
    }

    float Wavelet = 0.0f;
    float fenmu = 0.0f, fenzi = 0.0f;
    float *max_buf = new float[rank_size];
    float MaxValue = 0;

    if(interiorLength_z && interiorLength_x)
    {
        // 归一化梯度
        MaxValue = MatMax(ip->GradVp, interiorLength_z * interiorLength_x);
        //cout << MaxValue << endl;
    }

    if(rank != ROOT_ID)
    {
        MPI_Send(&MaxValue, 1, MPI_FLOAT, ROOT_ID, STEP_MAX, MPI_COMM_WORLD);
    }
    else
    {
        max_buf[rank_size - 1] = MaxValue;
        for(uint i = 0; i < rank_size - 1; ++i)
        {
            MPI_Recv(max_buf + i, 1, MPI_FLOAT, MPI_ANY_SOURCE, STEP_MAX, MPI_COMM_WORLD, &status_recv);
        }
        MaxValue = MatMax(max_buf, rank_size);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MaxValue, 1, MPI_FLOAT, ROOT_ID, MPI_COMM_WORLD);

//    if(rank == ROOT_ID)
//    cout << MaxValue << endl;
    if (MaxValue > 1.0e-10f)
    {
        for (uint m = 0; m < interiorLength_z * interiorLength_x; m++)
        {
            ip->GradVp[m] /= MaxValue;
        }
    }

    //cout << "rank" << rank << "MaxValue=" << MaxValue << endl;
    MPI_Barrier(MPI_COMM_WORLD);
//    if(rank == 0)
//    {
//        cout << "rank" << rank << "MaxValue=" << MaxValue << endl;
//    }

    // 生成试探个体
    memset((void *)plan->h_Grad, 0, sizeof(float) * interiorLength_z * interiorLength_x);
    for(uint i = 0; i < interiorLength_z; ++i)
    {
        memcpy(plan->h_Grad + i * interiorLength_x, ip->GradVp + i * interiorLength_x, interiorLength_x * sizeof(float));
    }

    for(uint i = 0; i < block_z; ++i)
    {
        memcpy(plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_x + temph_Vp.leftborder, ip->CurrVp + i * block_x, sizeof(float) * block_x);
    }

    UpdateVp(Pa, plan->h_Vp, plan->h_Grad, e, pt);

    vector<uint> trans_h_Vp = pt.gettrans_h_Vp();
    auto begin = trans_h_Vp.begin();
    //cout << "pppppppp" << rank << endl;
//    if(rank == 15)
//    {
//        cout << *(begin + TOP) << endl;
//        cout << *(begin + LEFT) << endl;
//        cout << *(begin + BOTTOM) << endl;
//        cout << *(begin + RIGHT) << endl;
//    }
    if(*(begin + TOP))
    {
        dataTransport_Vp(plan->h_Vp, pt, BOTTOM_TO_TOP, Pa);
    }
    if(*(begin + LEFT))
    {
        dataTransport_Vp(plan->h_Vp, pt, RIGHT_TO_LEFT, Pa);
    }
    if(*(begin + BOTTOM))
    {
        dataTransport_Vp(plan->h_Vp, pt, TOP_TO_BOTTOM, Pa);
    }
    if(*(begin + RIGHT))
    {
        dataTransport_Vp(plan->h_Vp, pt, LEFT_TO_RIGHT, Pa);
    }

    //cout << rank << endl;

    dataGather(plan->h_Vp, pt, STEP_VP);
    MPI_Barrier(MPI_COMM_WORLD);
//cout << "xxxxxxxxx" << pt.getrank() << endl;
//    if(temph_Vp.topborder)
//    {
//        dataGather(plan->h_Vp, pt, STEP_VP);
//    }
//    if(temph_Vp.leftborder)
//    {
//        dataGather(plan->h_Vp, pt, STEP_VP);
//    }
//    if(temph_Vp.bottomborder)
//    {
//        dataGather(plan->h_Vp, pt, STEP_VP);
//    }
//    if(temph_Vp.rightborder)
//    {
//        dataGather(plan->h_Vp, pt, STEP_VP);
//    }
//cout << "PPPPPPP" << rank << endl;
    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, e, pt);
//cout << "qqqqqqqqq" << rank << endl;
    // 对炮进行循环
    for (uint is = 0; is < ip->ShotN; is = is + 2)
    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * temph_U.length_z * temph_U.length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);
        memset((void *)plan->h_W,			0,	sizeof(float) * h_VW.length_z * h_VW.length_x);

        if(RL_num)
        {
            memset((void *)plan->h_TrueWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_TrailWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_CurrWF,		0,	sizeof(float) * RL_num * Pa.Nt);
            memset((void *)plan->h_ResWF,		0,	sizeof(float) * RL_num * Pa.Nt);
        }

        // 给检波器的位置信息赋值
        //memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));

        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            // 波场时刻转换
            for(uint i = 0; i < block_z; ++i)
            {
                memcpy(plan->h_U_past + i * block_x, plan->h_U_now + (i + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, block_x * sizeof(float));
                memcpy(plan->h_U_now + (i + temph_U.topborder) * temph_U.length_x + temph_U.leftborder, plan->h_U_next + i * block_x, block_x * sizeof(float));
            }
//cout << "xxxxxxxxx" << pt.getrank() << endl;
            dataTransport(plan->h_U_now, pt, STEP_U, it);//send data
            //dataGather(plan->h_U_now, pt, STEP_U);//recv data

            // 一步记录炮集
            StepShotGather(Pa, plan->h_U_now, plan->h_TrailWF, it, pt);

            //    ofstream fout("h_TrailWF.txt");
            //    ofstream fout1("h_ResWF.txt");
            //    if(pt.getrank() == 0)
            //    for(int i = 0; i < RL_num; ++i)
            //    {
            //        for(int j = 0; j < Pa.Nt; ++j)
            //        {
            //            fout << plan->h_TrailWF[i * Pa.Nt + j] << " ";
            //            fout1 << plan->h_ResWF[i * Pa.Nt + j] << " ";
            //        }
            //    }

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt, it);

            dataTransport(plan->h_V, pt, STEP_V, it);//send data
            dataTransport(plan->h_W, pt, STEP_W, it);//send data
            //dataGather(plan->h_V, pt, STEP_V);//recv_v data
            //dataGather(plan->h_W, pt, STEP_W);//recv_w data

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt, it);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            if(ip->St[is].s.Sx >= indexmin_x && ip->St[is].s.Sx <= indexmax_x && ip->St[is].s.Sz >= indexmin_z && ip->St[is].s.Sz <= indexmax_z)
            {
                AddSource(Pa, plan->h_U_next, Wavelet, plan->h_Vp, ip->St[is].s, pt);
            }
        }


        //MPI_Barrier(MPI_COMM_WORLD);
        //cout << "nnnnnnnnnnn" << pt.getrank() << endl;

        // 读取当前和真实波场

        if(RL_num)
        {
            memcpy(plan->h_CurrWF, sgs_c + is * Pa.Nt * RL_num,
                sizeof(float) * Pa.Nt * RL_num);
            memcpy(plan->h_TrueWF, sgs_t + is * Pa.Nt * RL_num,
                sizeof(float) * Pa.Nt * RL_num);


//            ofstream fout("h_TrailWF.txt");
//            ofstream fout1("h_ResWF.txt");
//            if(pt.getrank() == 0)
//            for(int i = 0; i < RL_num; ++i)
//            {
//                for(int j = 0; j < Pa.Nt; ++j)
//                {
//                    fout << plan->h_TrailWF[i * Pa.Nt + j] << " ";
//                    fout1 << plan->h_ResWF[i * Pa.Nt + j] << " ";
//                }
//            }


            // 求取上述波场和当前波场的残差
            MatAdd(plan->h_TrailWF, plan->h_TrailWF, plan->h_CurrWF,
                RL_num, Pa.Nt, 1.0f, 1.0f, 0);

            // 求取残差波场
            StepResidual(Pa, plan->h_TrueWF, plan->h_CurrWF,
                plan->h_ResWF, RL_num);

            // 将上述两种残差分别累加起来
            MatAdd(plan->h_SumResTrial, plan->h_SumResTrial,
                plan->h_TrailWF, RL_num, Pa.Nt, 1.0f, 1.0f, 1);
            MatAdd(plan->h_SumResCurr, plan->h_SumResCurr,
                plan->h_ResWF, RL_num, Pa.Nt, 1.0f, 1.0f, 1);
        }
    }




//    ofstream fout("h_TrailWF.txt");
//    ofstream fout1("h_ResWF.txt");
//    if(pt.getrank() == 0)
//    for(int i = 0; i < RL_num; ++i)
//    {
//        for(int j = 0; j < Pa.Nt; ++j)
//        {
//            fout << plan->h_TrailWF[i * Pa.Nt + j] << " ";
//            fout1 << plan->h_ResWF[i * Pa.Nt + j] << " ";
//        }
//    }


//    if(pt.getrank() == 0)
//    {
//        ofstream fout("h_CurrWF.txt");
//        ofstream fout1("h_TrueWF.txt");

//        for(int i = 0; i < RL_num; ++i)
//        {
//            for(int j = 0; j < Pa.Nt; ++j)
//            {
//                fout << plan->h_CurrWF[i * Pa.Nt + j] << " ";
//                fout1 << plan->h_TrueWF[i * Pa.Nt + j] << " ";
//            }
//        }
//    }




//    ofstream fout("h_SumResTrial.txt");
//    ofstream fout1("h_SumResCurr.txt");
//    if(pt.getrank() == 0)
//    for(int i = 0; i < RL_num; ++i)
//    {
//        for(int j = 0; j < Pa.Nt; ++j)
//        {
//            fout << plan->h_SumResTrial[i * Pa.Nt + j] << " ";
//            fout1 << plan->h_SumResCurr[i * Pa.Nt + j] << " ";
//        }
//    }


    int *ranks = new int[pt.getsumBlock_x()];
    int cpu_index_z = (Pa.PMLz + 2) / pt.getblockLength_z();
    int max_rank_RL = rank, min_rank_RL = rank;
    float buf[2];// = {fenzi, fenmu};
    MPI_Request request, request_recv_fen, request_recv_max_RL, request_recv_min_RL;
    MPI_Status status;
    // 求取最终的步长

    MPI_Barrier(MPI_COMM_WORLD);
//cout << "mmmmmmmm" << pt.getrank() << endl;


    //cout << "oooo" << rank << endl;
    if(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx)
    {
        min_rank_RL = rank;
        if(min_rank_RL != ROOT_ID)
            MPI_Send(&min_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MIN_RL, MPI_COMM_WORLD);
    }
    if(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1)
    {
        max_rank_RL = rank;
        if(max_rank_RL != ROOT_ID)
            MPI_Send(&max_rank_RL, 1, MPI_INT, ROOT_ID, STEP_MAX_RL, MPI_COMM_WORLD);
    }

    if(rank == ROOT_ID)
    {
        if(!(RL_num && indexmin_x <= Pa.PMLx + ip->St[0].rn - 1 && indexmax_x >= Pa.PMLx + ip->St[0].rn - 1))
            MPI_Recv(&max_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MAX_RL, MPI_COMM_WORLD, &status_recv);
        if(!(RL_num && indexmin_x <= Pa.PMLx && indexmax_x >= Pa.PMLx))
            MPI_Recv(&min_rank_RL, 1, MPI_INT, MPI_ANY_SOURCE, STEP_MIN_RL, MPI_COMM_WORLD, &status_recv);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&max_rank_RL, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);
    MPI_Bcast(&min_rank_RL, 1, MPI_INT, ROOT_ID, MPI_COMM_WORLD);

//    if(rank == ROOT_ID)
//        cout << max_rank_RL << " " << min_rank_RL << endl;

    if(RL_num)
    {
        if(rank != min_rank_RL)
        {
            MPI_Recv(buf, 2, MPI_FLOAT, rank - 1, STEP_FEN, MPI_COMM_WORLD, &status_recv);
            fenzi = buf[0];
            fenmu = buf[1];
            for (uint m = 0; m < Pa.Nt * RL_num; m++)
            {
                fenzi += plan->h_SumResCurr[m] * plan->h_SumResTrial[m];
                fenmu += plan->h_SumResTrial[m] * plan->h_SumResTrial[m];
            }
            buf[0] = fenzi;
            buf[1] = fenmu;

            if(rank != max_rank_RL)
            {
                MPI_Send(buf, 2, MPI_FLOAT, rank + 1, STEP_FEN, MPI_COMM_WORLD);
            }
            else
            {
                ip->Alpha = (fenzi / fenmu) * e;
            }
        }
        else
        {
            for (uint m = 0; m < Pa.Nt * RL_num; m++)
            {
                fenzi += plan->h_SumResCurr[m] * plan->h_SumResTrial[m];
                fenmu += plan->h_SumResTrial[m] * plan->h_SumResTrial[m];
            }
            buf[0] = fenzi;
            buf[1] = fenmu;
            MPI_Send(buf, 2, MPI_FLOAT, rank + 1, STEP_FEN, MPI_COMM_WORLD);
        }
    }



//    uint cpu_rank = cpu_index_z * sumblock_x;

//    MPI_Request request, request_recv_fen;
//    MPI_Status status;
//    if(rank > cpu_rank && rank < cpu_rank + sumblock_x)
//        MPI_Isend(buf, 2, MPI_FLOAT, cpu_rank, STEP_FEN, MPI_COMM_WORLD, &request);

//    if(rank == cpu_rank)
//    {//cout << "wwwww" << endl;
//        for(int i = 0; i < sumblock_x - 1; ++i)
//        {
//            MPI_Irecv(buf, 2, MPI_FLOAT, MPI_ANY_SOURCE, STEP_FEN, MPI_COMM_WORLD, &request_recv_fen);
//            MPI_Wait(&request_recv_fen, &status);
//            fenzi += buf[0];
//            fenmu += buf[1];
//            //cout << "fenzi = " << fenzi << "fenmu" << fenmu << endl;
//        }
//        ip->Alpha = (fenzi / fenmu) * e;
//    }

    //cout << ip->Alpha << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&ip->Alpha, 1, MPI_FLOAT, max_rank_RL, MPI_COMM_WORLD);


//    delete &buf;
//    delete ranks;
//    delete max_buf;
}


/*------------------------------------------------------------------------
function: PreProcess

Information:
下次迭代前的预处理

Inputs:
Pa: 正演参数
ip：反演参数
plan：全局变量
----------------------------------------------------------*/
void PreProcess(AFDPU2D Pa,
                IP *ip,
                CPUVs *plan,
                const Partition &pt)
{
    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    int interiorLength_x = pt.getinteriorLength_x();
    int interiorLength_z = pt.getinteriorLength_z();

    if(interiorLength_x == 0 || interiorLength_z == 0)
        return;

    int interior_min_x = pt.getinteriormin_x();
    int interior_min_z = pt.getinteriormin_z();

    H_Border temph_Vp = pt.geth_Vp();

    for(uint i = 0; i < interiorLength_z; ++i)
    {
        memcpy(plan->h_Grad + i * interiorLength_x, ip->GradVp + i * interiorLength_x, interiorLength_x * sizeof(float));
    }
    for(uint i = 0; i < block_z; ++i)
    {
        memcpy(plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_x + temph_Vp.leftborder, ip->CurrVp + i * block_x, block_x * sizeof(float));
    }
//cout << pt.getrank() << endl;
    // 更新速度
    UpdateVp(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);

    vector<uint> trans_h_Vp = pt.gettrans_h_Vp();
    auto begin = trans_h_Vp.begin();
    if(*(begin + TOP))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, BOTTOM_TO_TOP, Pa);
    }
    if(*(begin + LEFT))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, RIGHT_TO_LEFT, Pa);
    }
    if(*(begin + BOTTOM))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, TOP_TO_BOTTOM, Pa);
    }
    if(*(begin + RIGHT))
    {//cout << "what" << endl;
        dataTransport_Vp(plan->h_Vp, pt, LEFT_TO_RIGHT, Pa);
    }

    dataGather(plan->h_Vp, pt, STEP_VP);

//    if(temph_Vp.topborder)
//    {//cout << "eee" << endl;
//        dataGather(plan->h_Vp, pt, BOTTOM_TO_TOP);
//    }
//    if(temph_Vp.leftborder)
//    {//cout << "eee" << endl;
//        dataGather(plan->h_Vp, pt, RIGHT_TO_LEFT);
//    }
//    if(temph_Vp.bottomborder)
//    {//cout << "eee" << endl;
//        dataGather(plan->h_Vp, pt, TOP_TO_BOTTOM);
//    }
//    if(temph_Vp.rightborder)
//    {//cout << "eee" << endl;
//        dataGather(plan->h_Vp, pt, LEFT_TO_RIGHT);
//    }
//cout << "pppppp" << pt.getrank() << endl;
    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);
//cout << "qqqqqq" << pt.getrank() << endl;





    for(uint i = 0; i < block_z; ++i)
    {
        memcpy(ip->CurrVp + i * block_x, plan->h_Vp + (i + temph_Vp.topborder) * temph_Vp.length_x + temph_Vp.leftborder, block_x * sizeof(float));
    }

    ip->Alpha = 0.0f;
}

