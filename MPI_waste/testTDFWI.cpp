
/**************************************************************************
* 2维时间域全波形反演（Time domain Full Waveform Inversion）

* 使用非分裂完全匹配曾（NPML）技术处理吸收边界

* 使用2阶位移运动方程

* 正演使用空间8阶时间2阶精度的交错网格有限差分技术

* 作者：高照奇
**************************************************************************/

#include	"testTDFWI.h"
#define     TOP     0
#define     LEFT    1
#define     BOTTOM  2
#define     RIGHT   3
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

    f1 = fopen(FileName, "rb");
    if (f1 == NULL)
    {
        cout << "File Open error!" << "\n" << endl;
    }
    // read the first 3200 bytes
    fread(f3200, 3200, 1, f1);

    // read the information of Sgy file
    InfoOfSgy(FileName, *Reel, TraceNum,
        SampleNum, SampleInt, DFormat,
        &BReel, &BIBM);

    int length_x = pt.getblockLength_x();
    int length_z = pt.getblockLength_z();
    int indexmin_x = pt.getindexmin_x();
    int indexmin_z = pt.getindexmin_z();

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
    for(n = 0; n < length_z; ++n)
    {
        for(m = 0; m < length_x; ++m)
        {
            if(flag == 0)
            {
                Data[n * length_x + m] = trace[m].data[n];
            }
            else
            {
                Data[m * length_x + n] = trace[m].data[n];
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
void WriteData(char FileName[],
            usht SampleNum,
            usht TraceNum,
            usht SampleInt,
            float *data,
            usht flag)
{
    unsigned char f3200[3200];
    uint m, n;

    Trace *trace;
    trace = new Trace[TraceNum];
    memset((void *)trace, 0, sizeof(Trace) * TraceNum);
    for (int n = 0; n < TraceNum; n++)
    {
        trace[n].data = new float[SampleNum];
        memset((void *)trace[n].data, 0, sizeof(float) * SampleNum);
        for (m = 0; m < SampleNum; m++)
        {
            if (flag == 0)
            {
                trace[n].data[m] = *(data + m * TraceNum + n);
            }
            else
            {
                trace[n].data[m] = *(data + n * SampleNum + m);
            }
        }
    }

    // trace header
    for (n = 0; n < TraceNum; n++)
    {
        for (int m = 0; m < 60; m++)
        {
            trace[n].head.h4[m] = 0;
        }
    }

    // reel header
    for (n = 0; n < 3200; n++)
    {
        f3200[n] = 1;
    }

    WriteSgy(FileName, &f3200[0], trace,
        TraceNum, SampleNum, SampleInt);///hai

    // free memory
    for (n = 0; n < TraceNum; n++)
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
    //uint nnz = Pa.Nz + 2 * Pa.PMLz;// 纵向网格数 + 2*纵向NPML内的网格数
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    //uint RecNum = 2 * (Pa.Nz + Pa.Nx) * 8;

    //uint nnx = pt.getblockLength_x() + 2 * pt.getborderlength_x();
    //uint nnz = pt.getblockLength_z() + 2 * pt.getborderlength_z();
    uint length_x = pt.getblockLength_x();
    uint length_z = pt.getblockLength_z();

//    uint PMLlength_x = 0;
//    uint PMLlength_z = 0;

//    int insidenum = pt.getinsidenum();
//    Inside *ins = pt.getInside();
//    for(int i = 0; i < insidenum; ++i)
//    {
//        PMLlength_x += ins[i].getlength_x();
//        PMLlength_z += ins[i].getlength_z();
//    }

    // 申请变量// 与NPML有关的变量
    plan->h_dx = new float[length_x];
    plan->h_dz = new float[length_z];
    plan->h_Bx = new float[length_x];
    plan->h_Bz = new float[length_z];
    plan->h_PHIx_V_x = new float[length_z * length_x];
    plan->h_PHIz_W_z = new float[length_z * length_x];
    plan->h_PHIx_U_x = new float[length_z * length_x];//h_U right 4, left 3, seft 1
    plan->h_PHIz_U_z = new float[length_z * length_x];//h_U top 3, bottom 4, seft 1
    plan->h_PHIx_V_x_r = new float[length_z * length_x];//h_V left 4, right 3, seft 1
    plan->h_PHIz_W_z_r = new float[length_z * length_x];//h_W top 4, bottom 3, seft 1
    plan->h_PHIx_U_x_r = new float[length_z * length_x];
    plan->h_PHIz_U_z_r = new float[length_z * length_x];
    // 波场相关变量
    plan->h_U_past = new float[length_z * length_x];//h_V left 4, right 3, h_W top 4, bottom 3
    plan->h_U_now = new float[(length_z + pt.geth_U_topborder() + pt.geth_U_bottomborder()) * (length_x + pt.geth_U_leftborder() + pt.geth_U_rightborder())];
    plan->h_U_next = new float[length_z * length_x];//h_V left 4, right 3, h_W top 4, bottom 3
    plan->h_U_past_r = new float[length_z * length_x];
    plan->h_U_now_r = new float[(length_z + pt.geth_U_topborder() + pt.geth_U_bottomborder()) * (length_x + pt.geth_U_leftborder() + pt.geth_U_rightborder())];
    plan->h_U_next_r = new float[length_z * length_x];
    plan->h_V = new float[length_z * (length_x + pt.geth_V_leftborder() + pt.geth_V_rightborder())];//h_U left 3, right 4
    plan->h_W = new float[(length_z + pt.geth_W_topborder() + pt.geth_W_bottomborder()) * length_x];//h_U top 3, bottom 4
    plan->h_V_r = new float[length_z * (length_x + pt.geth_V_leftborder() + pt.geth_V_rightborder())];
    plan->h_W_r = new float[(length_z + pt.geth_W_topborder() + pt.geth_W_bottomborder()) * length_x];
    plan->h_Vp = new float[(length_z + pt.geth_Vp_border()[0] + pt.geth_Vp_border()[2]) * (length_x + pt.geth_Vp_border()[1] + pt.geth_Vp_border()[3])];
    plan->h_re = new RL[ip->St[0].rn];

    plan->h_TrueWF = new float[ip->St[0].rn * Pa.Nt];
    plan->h_CurrWF = new float[ip->St[0].rn * Pa.Nt];
    plan->h_ResWF = new float[ip->St[0].rn * Pa.Nt];

    plan->h_Grad = new float[pt.getinteriorLength_z() * pt.getinteriorLength_x()];//Pa.Nz * Pa.Nx
    plan->h_U_Der = new float[pt.getinteriorLength_z() * pt.getinteriorLength_x()];//Pa.Nz * Pa.Nx

    int sum_h_Coord = 0;
    int n = pt.get_h_Coord_num();
    H_Coord *coord = pt.get_h_Coord_num();
    for(int i = 0; i < n; ++i)
    {
        sum_h_Coord += coord->getlength_x() * coord->getlength_z();
    }
    plan->h_Coord = new uint[sum_h_Coord];
    plan->h_RU = new float[sum_h_Coord * (Pa.Nt - 1)];
    plan->h_Obj = new float[ip->St[0].rn * Pa.Nt];
    // 计算步长的相关变量
    plan->h_TrailWF = new float[ip->St[0].rn * Pa.Nt];
    plan->h_SumResTrial = new float[ip->St[0].rn * Pa.Nt];
    plan->h_SumResCurr = new float[ip->St[0].rn * Pa.Nt];

    // 初始化空间
    memset((void *)plan->h_dx,			0,	sizeof(float) * PMLlength_x);
    memset((void *)plan->h_dz,			0,	sizeof(float) * PMLlength_z);
    memset((void *)plan->h_Bx,			0,	sizeof(float) * length_x);
    memset((void *)plan->h_Bz,			0,	sizeof(float) * length_z);
    memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_U_past,		0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_U_now,		0,	sizeof(float) * (length_z + pt.geth_U_topborder() + pt.geth_U_bottomborder()) * (length_x + pt.geth_U_leftborder() + pt.geth_U_rightborder()));
    memset((void *)plan->h_U_next,		0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_U_past_r,	0,	sizeof(float) * length_z * length_x);
    memset((void *)plan->h_U_now_r,		0,	sizeof(float) * (length_z + pt.geth_U_topborder() + pt.geth_U_bottomborder()) * (length_x + pt.geth_U_leftborder() + pt.geth_U_rightborder()));
    memset((void *)plan->h_U_next_r,	0,	sizeof(float) * length_z * length_x);

    memset((void *)plan->h_V,			0,	sizeof(float) * length_z * (length_x + pt.geth_V_leftborder() + pt.geth_V_rightborder()));
    memset((void *)plan->h_W,			0,	sizeof(float) * (length_z + pt.geth_W_topborder() + pt.geth_W_bottomborder()) * length_x);
    memset((void *)plan->h_V_r,			0,	sizeof(float) * length_z * (length_x + pt.geth_V_leftborder() + pt.geth_V_rightborder()));
    memset((void *)plan->h_W_r,			0,	sizeof(float) * (length_z + pt.geth_W_topborder() + pt.geth_W_bottomborder()) * length_x);
    memset((void *)plan->h_Vp,			0,	sizeof(float) * (length_z + pt.geth_Vp_border()[0] + pt.geth_Vp_border()[2]) * (length_x + pt.geth_Vp_border()[1] + pt.geth_Vp_border()[3]));
    memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);

    memset((void *)plan->h_TrueWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_CurrWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_ResWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);

    memset((void *)plan->h_Grad,		0,	sizeof(float) * pt.getinteriorLength_z() * pt.getinteriorLength_x());
    memset((void *)plan->h_U_Der,		0,	sizeof(float) * pt.getinteriorLength_z() * pt.getinteriorLength_x());
    memset((void *)plan->h_Coord,		0,	sizeof(uint) * sum_h_Coord);
    memset((void *)plan->h_RU,			0,	sizeof(uint) * sum_h_Coord * (Pa.Nt - 1));
    memset((void *)plan->h_Obj,			0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_TrailWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_SumResTrial,	0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_SumResCurr,	0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
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
    //uint nnz = Pa.Nz + 2 * Pa.PMLz;
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;

    uint length_z = pt.getblockLength_z();
    uint length_x = pt.getblockLength_x();

    // 对NPML参数进行赋值
    float d0 = -3.0f * 2500.0f * logf(1.0e-6f)
        / (2.0f * powf((float)(Pa.PMLz * Pa.dz), 3.0f));
    float d1 = -3.0f * 2500.0f * log(1.0e-6f)
        / (2.0f * powf((float)(Pa.PMLx * Pa.dx), 3.0f));

    int n = pt.getinsidenum();
    Inside *ins = pt.getInside();

    for(int i = 0; i < n; ++i)
    {
        int min_x = ins[i].getindexmin_x();
        int min_z = ins[i].getindexmin_z();
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
            int actual_ix = ix + min_x;
            if (actual_ix < Pa.PMLx)
            {
                plan->h_dx[ix] = d0 * (Pa.PMLx - 1 - actual_ix + 0.5f) *
                    (Pa.PMLx - 1 - actual_ix + 0.5f) * Pa.dx * Pa.dx;// NPML横向参数
            }
            else if (actual_ix > Pa.PMLx + Pa.Nx - 1)
            {
                plan->h_dx[ix] = d0 * (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) *
                    (actual_ix - Pa.PMLx - Pa.Nx + 0.5f) * Pa.dx * Pa.dx;// NPML横向参数
            }
            else
            {
                plan->h_dx[ix] = 0.0f;// NPML横向参数
            }
            plan->h_Bx[ix] = expf(-1.0f * plan->h_dx[ix] * Pa.dt);// NPML横向参数
        }
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
                SL s,
                float Wave,
                float *h_Vp,
               const Partition &pt)
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint indexmin_x = pt.getblockLength_x();
    uint block_x = pt.getblockLength_x();
    uint length_x = pt.geth_Vp_border()[LEFT] + pt.geth_Vp_border()[RIGHT];
    vector<pair<int, int>> vec = pt.getShot();
    pair<int, int> temp;
    for(int i = 0; i < vec.size(); ++i)
    {
        temp = vec.back();
        h_U[temp.second * length_x + temp.first] += Wave * (Pa.dt * Pa.dt * h_Vp[id]);
    }
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
            const Partition &pt)
{
    int top_border = pt.geth_U_topborder();
    int bottom_border = pt.geth_U_bottomborder();
    int left_border = pt.geth_U_leftborder();
    int right_border = pt.geth_U_rightborder();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dUx = 0.0f;
    float dUz = 0.0f;

    // 交错网格有限差分
    //yuan cheng xu li cong 4 kai shi,zhe li ye shi cong 4 kai shi
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dUx = C1_4 * (h_U[iz * length_x + ix + 1] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[iz * length_x + ix + 2] - h_U[iz * length_x + ix - 1])
                + C3_4 * (h_U[iz * length_x + ix + 3] - h_U[iz * length_x + ix - 2])
                + C4_4 * (h_U[iz * length_x + ix + 4] - h_U[iz * length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * length_x + ix] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[(iz + 2) * length_x + ix] - h_U[(iz - 1) * length_x + ix])
                + C3_4 * (h_U[(iz + 3) * length_x + ix] - h_U[(iz - 2) * length_x + ix])
                + C4_4 * (h_U[(iz + 4) * length_x + ix] - h_U[(iz - 3) * length_x + ix]);

            h_PHIx_U_x[iz * length_x + ix] =
                h_Bx[ix] * (h_PHIx_U_x[iz * length_x + ix] + dUx / Pa.dx) - dUx / Pa.dx;
            h_PHIz_U_z[iz * length_x + ix] =
                h_Bz[iz] * (h_PHIz_U_z[iz * length_x + ix] + dUz / Pa.dz) - dUz / Pa.dz;
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
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    int top_border = pt.geth_W_topborder();
    int bottom_border = pt.geth_W_bottomborder();
    int left_border = pt.geth_V_leftborder();
    int right_border = pt.geth_V_rightborder();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;
    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dVx = 0.0f;
    float dWz = 0.0f;

    // 交错网格有限差分
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dVx = C1_4 * (h_V[iz * length_x + ix] - h_V[iz * length_x + ix - 1])
                + C2_4 * (h_V[iz * length_x + ix + 1] - h_V[iz * length_x + ix - 2])
                + C3_4 * (h_V[iz * length_x + ix + 2] - h_V[iz * length_x + ix - 3])
                + C4_4 * (h_V[iz * length_x + ix + 3] - h_V[iz * length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * length_x + ix] - h_W[(iz - 1) * length_x + ix])
                + C2_4 * (h_W[(iz + 1) * length_x + ix] - h_W[(iz - 2) * length_x + ix])
                + C3_4 * (h_W[(iz + 2) * length_x + ix] - h_W[(iz - 3) * length_x + ix])
                + C4_4 * (h_W[(iz + 3) * length_x + ix] - h_W[(iz - 4) * length_x + ix]);

            h_PHIx_V_x[iz * length_x + ix] =
                h_Bx[ix] * (h_PHIx_V_x[iz * length_x + ix] + dVx / Pa.dx) - dVx / Pa.dx;
            h_PHIz_W_z[iz * length_x + ix] =
                h_Bz[iz] * (h_PHIz_W_z[iz * length_x + ix] + dWz / Pa.dz) - dWz / Pa.dz;
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
            const Partition &pt)//
{
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    int top_border = pt.geth_W_topborder();
    int bottom_border = pt.geth_W_bottomborder();
    int left_border = pt.geth_V_leftborder();
    int right_border = pt.geth_V_rightborder();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;
    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dVx = 0.0f;
    float dWz = 0.0f;

    // 交错网格有限差分
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dVx = C1_4 * (h_V[iz * length_x + ix] - h_V[iz * length_x + ix - 1])
                + C2_4 * (h_V[iz * length_x + ix + 1] - h_V[iz * length_x + ix - 2])
                + C3_4 * (h_V[iz * length_x + ix + 2] - h_V[iz * length_x + ix - 3])
                + C4_4 * (h_V[iz * length_x + ix + 3] - h_V[iz * length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * length_x + ix] - h_W[(iz - 1) * length_x + ix])
                + C2_4 * (h_W[(iz + 1) * length_x + ix] - h_W[(iz - 2) * length_x + ix])
                + C3_4 * (h_W[(iz + 2) * length_x + ix] - h_W[(iz - 3) * length_x + ix])
                + C4_4 * (h_W[(iz + 3) * length_x + ix] - h_W[(iz - 4) * length_x + ix]);

            h_U_next[iz * length_x + ix] = (Pa.dt * Pa.dt * h_Vp[iz * length_x + ix])
                * (dVx / Pa.dx + h_PHIx_V_x[iz * length_x + ix] + dWz / Pa.dz + h_PHIz_W_z[iz * length_x + ix])
                - h_U_past[iz * length_x + ix] + 2.0f * h_U_now[iz * length_x + ix];
        }
    }
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
            const Partition &pt)
{
    //uint nnz = Pa.Nz + 2 * Pa.PMLz;
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    int top_border = pt.geth_U_topborder();
    int bottom_border = pt.geth_U_bottomborder();
    int left_border = pt.geth_U_leftborder();
    int right_border = pt.geth_U_rightborder();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dUx = 0.0f;
    float dUz = 0.0f;

    // 交错网格有限差分
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dUx = C1_4 * (h_U[iz * length_x + ix + 1] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[iz * length_x + ix + 2] - h_U[iz * length_x + ix - 1])
                + C3_4 * (h_U[iz * length_x + ix + 3] - h_U[iz * length_x + ix - 2])
                + C4_4 * (h_U[iz * length_x + ix + 4] - h_U[iz * length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * length_x + ix] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[(iz + 2) * length_x + ix] - h_U[(iz - 1) * length_x + ix])
                + C3_4 * (h_U[(iz + 3) * length_x + ix] - h_U[(iz - 2) * length_x + ix])
                + C4_4 * (h_U[(iz + 4) * length_x + ix] - h_U[(iz - 3) * length_x + ix]);

            h_V[iz * length_x + ix] = dUx / Pa.dx + h_PHIx_U_x[iz * length_x + ix];
            h_W[iz * length_x + ix] = dUz / Pa.dz + h_PHIz_U_z[iz * length_x + ix];
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
                    RL *h_re,
                    uint nstep,
                    uint rn,
                    const Partition &pt)
{
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id = 0;

    int indexmin_x = pt.getindexmin_x();
    int indexmax_x = pt.getindexmax_x();
    int indexmin_z = pt.getindexmin_z();
    int indexmax_z = pt.getindexmax_z();
    int length_x = pt.getblockLength_x();

    for (uint ir = 0; ir < rn; ir++)// 检波器个数
    {
        if(indexmin_z <= h_re[ir].Rz && h_re[ir].Rz <= indexmax_z && h_re[ir].Rx > indexmin_x && h_re[ir].Rx < indexmax_x)
        {
            id = (h_re[ir].Rz - indexmin_z) * length_x + h_re[ir].Rx - indexmin_x;
            h_SG[ir * Pa.Nt + nstep] = h_U[id];
        }
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
                float *h_RU)//zhe bu fen ye yao gai
{
    for (uint ix = 0; ix < RecNum; ix++)
    {
        h_RU[ix] = h_U[h_Coord[ix]];
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
                float *h_RU)
{
    for (uint ix = 0; ix < RecNum; ix++)
    {
        h_U[h_Coord[ix]] = h_RU[ix];
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
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    int top_border = pt.geth_U_topborder();
    int bottom_border = pt.geth_U_bottomborder();
    int left_border = pt.geth_U_leftborder();
    int right_border = pt.geth_U_rightborder();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dUx = 0.0f;
    float dUz = 0.0f;

    // 交错网格有限差分
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dUx = C1_4 * (h_U[iz * length_x + ix + 1] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[iz * length_x + ix + 2] - h_U[iz * length_x + ix - 1])
                + C3_4 * (h_U[iz * length_x + ix + 3] - h_U[iz * length_x + ix - 2])
                + C4_4 * (h_U[iz * length_x + ix + 4] - h_U[iz * length_x + ix - 3]);
            dUz = C1_4 * (h_U[(iz + 1) * length_x + ix] - h_U[iz * length_x + ix])
                + C2_4 * (h_U[(iz + 2) * length_x + ix] - h_U[(iz - 1) * length_x + ix])
                + C3_4 * (h_U[(iz + 3) * length_x + ix] - h_U[(iz - 2) * length_x + ix])
                + C4_4 * (h_U[(iz + 4) * length_x + ix] - h_U[(iz - 3) * length_x + ix]);

            h_V[iz * length_x + ix] = dUx / Pa.dx + h_PHIx_U_x[iz * length_x + ix];
            h_W[iz * length_x + ix] = dUz / Pa.dz + h_PHIz_U_z[iz * length_x + ix];
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
//    uint nnx = Pa.Nx + 2 * Pa.PMLx;
//    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    int top_border = pt.geth_W_topborder();
    int bottom_border = pt.geth_W_bottomborder();
    int left_border = pt.geth_V_leftborder();
    int right_border = pt.geth_V_rightborder();

    uint totallength_x = pt.gettotallength_x();
    uint totallength_z = pt.gettotallength_z();

    uint length_x = pt.getblockLength_x() + left_border + right_border;
    //uint length_z = pt.getblockLength_z() + top_border + bottom_border;
    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();

    float dVx = 0.0f;
    float dWz = 0.0f;

    // 交错网格有限差分
    for (uint iz = (top_border > 4 - min_z + top_border ? top_border : 4 - min_z + top_border); iz < totallength_z - bottom_border - min_z; iz++)
    {
        for (uint ix = (left_border > 4 - min_x +left_border ? left_border : 4 - min_x + left_border); ix < totallength_x - right_border - min_z; ix++)
        {
            dVx = C1_4 * (h_V[iz * length_x + ix] - h_V[iz * length_x + ix - 1])
                + C2_4 * (h_V[iz * length_x + ix + 1] - h_V[iz * length_x + ix - 2])
                + C3_4 * (h_V[iz * length_x + ix + 2] - h_V[iz * length_x + ix - 3])
                + C4_4 * (h_V[iz * length_x + ix + 3] - h_V[iz * length_x + ix - 4]);
            dWz = C1_4 * (h_W[iz * length_x + ix] - h_W[(iz - 1) * length_x + ix])
                + C2_4 * (h_W[(iz + 1) * length_x + ix] - h_W[(iz - 2) * length_x + ix])
                + C3_4 * (h_W[(iz + 2) * length_x + ix] - h_W[(iz - 3) * length_x + ix])
                + C4_4 * (h_W[(iz + 3) * length_x + ix] - h_W[(iz - 4) * length_x + ix]);

            h_U_past[iz * length_x + ix] = (Pa.dt * Pa.dt * h_Vp[iz * length_x + ix])
                * (dVx / Pa.dx + h_PHIx_V_x[iz * length_x + ix] + dWz / Pa.dz + h_PHIz_W_z[iz * length_x + ix])
                - h_U_next[iz * length_x + ix] + 2.0f * h_U_now[iz * length_x + ix];
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
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id1 = 0, id2 = 0;

    int interiorlength_x = pt.getinteriorLength_x();
    int interiorlength_z = pt.getinteriorlength_z();
    int length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    for (uint iz = 0; iz < interiorlength_z; iz++)
    {
        for (uint ix = 0; ix < interiorlength_x; ix++)
        {
            id1 = iz * interiorlength_x + ix;
            id2 = (iz + gap_z) * length_x + ix + gap_x;
            h_U_Der[id1] =
                (h_U_next[id2] - 2.0f * h_U_now[id2] + h_U_past[id2])
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

    int interiorlength_x = pt.getinteriorLength_x();
    int interiorlength_z = pt.getinteriorlength_z();
    int length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    for (uint iz = 0; iz < interiorlength_z; iz++)
    {
        for (uint ix = 0; ix < interiorlength_x; ix++)
        {
            id1 = iz * interiorlength_x + ix;
            id2 = (iz + gap_z) * length_x + ix + gap_x;
            h_Grad[id1] += -1.0f * h_U_Der[id1] * h_U_r[id2];
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
                uint Rn)//yi ding de fan wei nei zhi you yi ding geshu de jianboqi
{
    uint id = 0;
    for (uint n = 0; n < Rn; n++)
    {
        for (uint m = 0; m < Pa.Nt; m++)
        {
            id = n * Pa.Nt + m;
            h_ResWF[id] = h_TrueWF[id] - h_CurrWF[id];
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
                RL *h_re,
                uint nstep,
                uint Rn,
                float *h_Vp,
                const Partition &pt)
{
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id = 0;

    int indexmin_x = pt.getindexmin_x();
    int indexmin_z = pt.getindexmin_z();
    int indexmax_x = pt.getindexmax_x();
    int indexmax_z = pt.getindexmax_z();
    int length_x = pt.getblockLength_x();
    //int length_z = pt.getblockLength_z();

    for (uint ir = 0; ir < Rn; ir++)
    {
        if(h_re[ir].Rz >= indexmin_z && h_re[ir].Rx <= indexmax_x && h_re[ir].Rz >= indexmin_z && h_re[ir].Rz <= indexmax_z)
        {
            id = (h_re[ir].Rz - indexmin_z) * length_x + h_re[ir].Rx - indexmin_x;
            h_U_r[id] += -1.0f * h_ResWf[ir * Pa.Nt + nstep] * (Pa.dt * Pa.dt * h_Vp[id]);
        }
    }
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
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id1 = 0, id2 = 0;


    int interiorLength_x = pt.getinteriorLength_x();
    int interiorLength_z = pt.getinteriorLength_z();
    int length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    int interiormin_z = pt.getinteriormin_z();
    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint ix = 0; ix < interiorLength_x; ix++)
        {
            id1 = iz * interiorLength_x + ix;
            id2 = (iz + gap_z) * length_x + ix + gap_x;

            h_Grad[id1] = (powf((float)(interiormin_z + iz), PowerGrad)
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
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint id1 = 0, id2 = 0;

    int interiorLength_x = pt.getinteriorLength_x();
    int interiorLength_z = pt.getinteriorLength_z();
    int length_x = pt.getblockLength_x();

    int gap_x = pt.getinteriormin_x() - pt.getindexmin_x();
    int gap_z = pt.getinteriormin_z() - pt.getindexmin_z();
    for (uint iz = 0; iz < interiorLength_z; iz++)
    {
        for (uint ix = 0; ix < interiorLength_x; ix++)
        {
            id1 = iz * interiorLength_x + ix;
            id2 = (iz + gap_z) * length_x + ix + gap_x;

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
    //uint nnx = Pa.Nx + 2 * Pa.PMLx;
    //uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint id = 0;
    int block_x = pt.getblockLength_x();
    int block_z = pt.getblockLength_z();
    int length_x = block_x + pt.geth_Vp_border()[1] + pt.geth_Vp_border()[3];
    int length_z = block_z + pt.geth_Vp_border()[0] + pt.geth_Vp_border()[2];

    int min_x = pt.getindexmin_x();
    int min_z = pt.getindexmin_z();
    int max_x = pt.getindexmax_x();
    int max_z = pt.getindexmax_z();

    for (uint iz = 0; iz < block_z; iz++)
    {
        for (uint ix = 0; ix < block_x; ix++)
        {

            id = iz * length_x + ix;
            int temp_z = iz + min_z;
            int temp_x = ix + min_x;
            // 纵向
            if (temp_z < Pa.PMLz)
            {
                if(max_z >= Pa.PMLz)
                    h_Vp[id] = h_Vp[(Pa.PMLz - min_z) * length_x + ix];
                else
                {
                    h_Vp[id] = h_Vp[(length_z - 1) * length_x + ix];///tong yi ge fang xiang shang zhi neng youyige border
                }
            }
            if (temp_z > Pa.Nz + Pa.PMLz - 1)
            {
                if(min_z <= Pa.Nz + Pa.PMLz - 1)
                    h_Vp[id] = h_Vp[(Pa.Nz + Pa.PMLz - min_z - 1) * length_x + ix];
                else
                {
                    h_Vp[id] = h_Vp[(length_z - 1) * length_x + ix];
                }
            }
            // 横向
            if (temp_x < Pa.PMLx)
            {
                if(max_x >= Pa.PMLx)
                    h_Vp[id] = h_Vp[length_x * iz + Pa.PMLx - min_x];
                else
                    h_Vp[id] = h_Vp[length_x * (iz + 1) - 1];////?
            }
            if (temp_x > Pa.Nx + Pa.PMLx - 1)
            {
                if(min_x <= Pa.Nx + Pa.PMLx - 1)
                    h_Vp[id] = h_Vp[iz * length_x + Pa.Nx + Pa.PMLx - 1 - min_x];
                else
                    h_Vp[id] = h_Vp[length_x * (iz + 1) - 1];
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

    int n = pt.get_h_Coord_num();
    H_Coord *coord = pt.get_h_Coord();
    // set Coord
    int before = 0;
    for(int i = 0; i < n; ++i)
    {
        int indexmin_x = coord[i].getindexmin_x();
        int indexmin_z = coord[i].getindexmin_z();
        int length_x = coord[i].getlength_x();
        int length_z = coord[i].getlength_z();

        if(i > 0)
        {
            before += coord[i-1].getlength_x() * coord[i-1].getlength_z();
        }
        for (uint iz = 0; iz < length_z; iz++)
        {
            for (uint ix = 0; ix < length_x; ix++)
            {

                h_Coord[before + iz * length_x + ix] = (indexmin_z + iz) * nnx + ix + indexmin_x;////?
            }
        }
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
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    uint nnz = Pa.Nz + 2 * Pa.PMLz;

    uint block_x = pt.getblockLength_x();
    uint block_z = pt.getblockLength_z();

    uint length_x = block_x + pt.geth_Vp_border()[1] + pt.geth_Vp_border()[3];
    uint length_z = block_z + pt.geth_Vp_border()[0] + pt.geth_Vp_border()[2];

    float Wavelet = 0.0f;

    // 给速度赋值
    memset((void *)plan->h_Vp, 0, sizeof(float) * length_z * length_x);// h_Vp 正演中使用的速度
    memcpy(plan->h_Vp,	ip->TrueVp,	 nnz * nnx * sizeof(float));

    for(int i = 0; i < length_z; ++i)
    {
        memcpy(plan->h_Vp[i * length_x], ip->TrueVp[i * block_x], block_x * sizeof(float));
    }
    int left_border = pt.geth_U_leftborder();
    int top_border = pt.geth_U_topborder();
    length_x = block_x + left_border + pt.geth_U_rightborder();
    length_z = block_z + top_border + pt.geth_U_bottomborder();
    for (uint is = 0; is < pt.getShot_num(); is++)// 反演中的炮数
    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * length_z * length_x);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * block_z * block_x);
        memset((void *)plan->h_V,			0,	sizeof(float) * block_z * (block_x + pt.geth_V_leftborder() + pt.geth_V_rightborder()));
        memset((void *)plan->h_W,			0,	sizeof(float) * (block_z + pt.geth_W_topborder() + pt.geth_W_bottomborder()) * block_x);
        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);
        memset((void *)plan->h_TrueWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);

        // 给检波器的位置信息赋值
        memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));// 检波器位置数组

        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)// Pa.Nt正演的时间步数
        {
            // 波场时刻转换
            for(int i = top_border; i < top_border + length_z; ++i)
            {
                memcpy(plan->h_U_past[i * block_x], plan->h_U_now[left_border + i * length_x], block_x * sizeof(float));
                memcpy(plan->h_U_now[left_border + i * length_x], plan->h_U_next[i * block_x], block_x * sizeof(float));
            }

            // 一步记录炮集
            StepShotGather(Pa, plan->h_U_now, plan->h_TrueWF,
                plan->h_re, it, ip->St[is].rn);//h_TrueWF

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);//h_PHIx_U_x   h_PHIz_U_z

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);//h_V  h_W

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);//h_PHIx_V_x    h_PHIz_W_z

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt);//h_U_next

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            if(pt.getShot_num())
            {
                AddSource(Pa, plan->h_U_next, ip->St[is].s,
                          Wavelet, plan->h_Vp);//h_U_next
            }
        }

        // 输出炮集
        memcpy(sgs_t + is * (Pa.Nt * ip->St[is].rn),
            plan->h_TrueWF,
            Pa.Nt * ip->St[is].rn * sizeof(float));
    }
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
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;
    float Wavelet = 0.0f;
    uint RecNum = (2 * Pa.Nx + 2 * Pa.Nz) * 8;

    // 梯度变量空间清零
    memset((void *)plan->h_Grad, 0, sizeof(float) * Pa.Nz * Pa.Nx);

    // 给速度赋值
    memset((void *)plan->h_Vp,			0,	sizeof(float) * nnz * nnx);
    memcpy(plan->h_Vp,	ip->CurrVp,	 nnz * nnx * sizeof(float));

    // 对炮进行循环
    for (uint is = 0; is < ip->ShotN; is++)
    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_V,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_W,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIx_U_x_r,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_U_z_r,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIx_V_x_r,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_W_z_r,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_past_r,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_now_r,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_next_r,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_V_r,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_W_r,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);
        memset((void *)plan->h_TrueWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_CurrWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_ResWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_U_Der,		0,	sizeof(float) * Pa.Nz * Pa.Nx);
        memset((void *)plan->h_RU,			0,	sizeof(uint) * RecNum * (Pa.Nt - 1));
        memset((void *)plan->h_Obj,			0,	sizeof(float) * ip->St[0].rn * Pa.Nt);

        // 给检波器位置赋值
        memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));

        // 计算正演波场，对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            // 波场时刻转换
            memcpy(plan->h_U_past, plan->h_U_now, nnz * nnx * sizeof(float));
            memcpy(plan->h_U_now, plan->h_U_next, nnz * nnx * sizeof(float));

            // 一步记录炮集
            StepShotGather(Pa, plan->h_U_now, plan->h_CurrWF,
                plan->h_re, it, ip->St[is].rn);

            // 一步记录有效边界内的波场
            if (it < Pa.Nt - 1)
            {
                StepRecordU(Pa, plan->h_Coord, RecNum, plan->h_U_now,
                    &plan->h_RU[it * RecNum]);
            }

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            AddSource(Pa, plan->h_U_next, ip->St[is].s,
                Wavelet, plan->h_Vp);//h_U_next
        }

        // 输出炮集
        memcpy(sgs_c + is * (Pa.Nt * ip->St[is].rn),
            plan->h_CurrWF,
            Pa.Nt * ip->St[is].rn * sizeof(float));

        // 重新读入观测数据
        memcpy(plan->h_TrueWF,
            sgs_t + is * (Pa.Nt * ip->St[is].rn),
            Pa.Nt * ip->St[is].rn * sizeof(float));

        // 计算残差波场
        StepResidual(Pa, plan->h_TrueWF, plan->h_CurrWF,
            plan->h_ResWF, ip->St[is].rn);

        //  输出残差波场
        memcpy(sgs_r + is * (Pa.Nt * ip->St[is].rn),
            plan->h_ResWF,
            Pa.Nt * ip->St[is].rn * sizeof(float));

        //  残差反传以及正传波场逆时间反推
        for (uint it = Pa.Nt - 1; it > 0; it--)
        {
            // 正传波场逆时间反推
            if (it == Pa.Nt - 1)
            {
                // 消除震源
                Wavelet = Ricker(Pa.f0, it * Pa.dt);
                AddSource(Pa, plan->h_U_next, ip->St[is].s,
                    -1.0f * Wavelet, plan->h_Vp);

                // 求取波场对时间的2阶导数
                StepCal2Der(Pa, plan->h_U_past, plan->h_U_now,
                    plan->h_U_next, plan->h_U_Der, pt);

                // 波场转换
                memcpy(plan->h_U_next, plan->h_U_now, nnz * nnx * sizeof(float));
                memcpy(plan->h_U_now, plan->h_U_past, nnz * nnx * sizeof(float));
            }
            else
            {
                // 消除震源
                Wavelet = Ricker(Pa.f0, it * Pa.dt);
                AddSource(Pa, plan->h_U_next, ip->St[is].s,
                    -1.0f * Wavelet, plan->h_Vp);

                // 替换有效边界内的波场
                StepReplaceU(Pa, plan->h_Coord, RecNum, plan->h_U_now,
                    &plan->h_RU[it * RecNum]);

                // 逆时间更新V和W
                StepRtVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                    plan->h_PHIx_U_x, plan->h_PHIz_U_z, pt);

                // 逆时间更新U
                StepRtU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                    plan->h_V, plan->h_W,plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                    plan->h_Vp, pt);

                // 求取波场对时间的2阶导数
                StepCal2Der(Pa, plan->h_U_past, plan->h_U_now,
                    plan->h_U_next, plan->h_U_Der, pt);

                // 波场转换
                memcpy(plan->h_U_next, plan->h_U_now, nnz * nnx * sizeof(float));
                memcpy(plan->h_U_now, plan->h_U_past, nnz * nnx * sizeof(float));
            }

            //  残差反传
            // 波场时刻转换
            memcpy(plan->h_U_past_r, plan->h_U_now_r, nnz * nnx * sizeof(float));
            memcpy(plan->h_U_now_r, plan->h_U_next_r, nnz * nnx * sizeof(float));

            // 一步求取梯度
            StepCalGrad(Pa, plan->h_Grad, plan->h_U_Der, plan->h_U_now_r, pt);

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now_r, plan->h_PHIx_U_x_r,
                plan->h_PHIz_U_z_r, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now_r, plan->h_V_r, plan->h_W_r,
                plan->h_PHIx_U_x_r, plan->h_PHIz_U_z_r, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V_r, plan->h_W_r, plan->h_PHIx_V_x_r,
                plan->h_PHIz_W_z_r, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next_r, plan->h_U_now_r, plan->h_U_past_r,
                plan->h_V_r, plan->h_W_r, plan->h_PHIx_V_x_r, plan->h_PHIz_W_z_r,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt);

            // 加震源
            AddResidual(Pa, plan->h_ResWF, plan->h_U_next_r, plan->h_re,
                it, ip->St[is].rn, plan->h_Vp);
        }

        // 求取目标函数
        for (uint m = 0; m < ip->St[is].rn * Pa.Nt; m++)
        {
            ip->ObjIter[It] += 0.5f * powf(plan->h_ResWF[m], 2.0f);
        }
    }
    // 输出梯度
    memcpy(ip->GradVp, plan->h_Grad, Pa.Nz * Pa.Nx * sizeof(float));
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
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    memset((void *)plan->h_SumResTrial,	0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
    memset((void *)plan->h_SumResCurr,	0,	sizeof(float) * ip->St[0].rn * Pa.Nt);

    float Wavelet = 0.0f;
    float fenmu = 0.0f, fenzi = 0.0f;

    // 归一化梯度
    float MaxValue = MatMax(ip->GradVp, Pa.Nz * Pa.Nx);
    if (MaxValue > 1.0e-10f)
    {
        for (uint m = 0; m < Pa.Nz * Pa.Nx; m++)
        {
            ip->GradVp[m] /= MaxValue;
        }
    }

    // 生成试探个体
    memset((void *)plan->h_Grad, 0, sizeof(float) * Pa.Nz * Pa.Nx);
    memcpy(plan->h_Grad, ip->GradVp, sizeof(float) * Pa.Nz * Pa.Nx);
    memcpy(plan->h_Vp, ip->CurrVp, sizeof(float) * nnx * nnz);
    UpdateVp(Pa, plan->h_Vp, plan->h_Grad, e, pt);
    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, e, pt);

    // 对炮进行循环
    for (uint is = 0; is < ip->ShotN; is = is + 2)
    {
        memset((void *)plan->h_PHIx_U_x,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_U_z,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIx_V_x,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_PHIz_W_z,	0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_past,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_now,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_U_next,		0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_V,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_W,			0,	sizeof(float) * nnz * nnx);
        memset((void *)plan->h_re,			0,	sizeof(RL) * ip->St[0].rn);
        memset((void *)plan->h_TrueWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_TrailWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_CurrWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);
        memset((void *)plan->h_ResWF,		0,	sizeof(float) * ip->St[0].rn * Pa.Nt);

        // 给检波器的位置信息赋值
        memcpy(plan->h_re,	ip->St[is].re, ip->St[is].rn * sizeof(RL));

        // 对时间进行循环
        for (uint it = 0; it < Pa.Nt; it++)
        {
            // 波场时刻转换
            memcpy(plan->h_U_past, plan->h_U_now, nnz * nnx * sizeof(float));
            memcpy(plan->h_U_now, plan->h_U_next, nnz * nnx * sizeof(float));

            // 一步记录炮集
            StepShotGather(Pa, plan->h_U_now, plan->h_TrailWF,
                plan->h_re, it, ip->St[is].rn);

            // 一步更新波场U的卷积项
            StepPHIU(Pa, plan->h_U_now, plan->h_PHIx_U_x,
                plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场V和W
            StepVW(Pa, plan->h_U_now, plan->h_V, plan->h_W,
                plan->h_PHIx_U_x, plan->h_PHIz_U_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新V和W的卷积项
            StepPHIVW(Pa, plan->h_V, plan->h_W, plan->h_PHIx_V_x,
                plan->h_PHIz_W_z, plan->h_Bx, plan->h_Bz, pt);

            // 一步更新波场U
            StepU(Pa, plan->h_U_next, plan->h_U_now, plan->h_U_past,
                plan->h_V, plan->h_W, plan->h_PHIx_V_x, plan->h_PHIz_W_z,
                plan->h_Bx, plan->h_Bz, plan->h_Vp, pt);

            // 加震源
            Wavelet = Ricker(Pa.f0, it * Pa.dt);

            AddSource(Pa, plan->h_U_next, ip->St[is].s,
                Wavelet, plan->h_Vp);
        }

        // 读取当前和真实波场
        memcpy(plan->h_CurrWF, sgs_c + is * Pa.Nt * ip->St[0].rn,
            sizeof(float) * Pa.Nt * ip->St[0].rn);
        memcpy(plan->h_TrueWF, sgs_t + is * Pa.Nt * ip->St[0].rn,
            sizeof(float) * Pa.Nt * ip->St[0].rn);

        // 求取上述波场和当前波场的残差
        MatAdd(plan->h_TrailWF, plan->h_TrailWF, plan->h_CurrWF,
            ip->St[0].rn, Pa.Nt, 1.0f, 1.0f, 0);

        // 求取残差波场
        StepResidual(Pa, plan->h_TrueWF, plan->h_CurrWF,
            plan->h_ResWF, ip->St[is].rn);

        // 将上述两种残差分别累加起来
        MatAdd(plan->h_SumResTrial, plan->h_SumResTrial,
            plan->h_TrailWF, ip->St[0].rn, Pa.Nt, 1.0f, 1.0f, 1);
        MatAdd(plan->h_SumResCurr, plan->h_SumResCurr,
            plan->h_ResWF, ip->St[0].rn, Pa.Nt, 1.0f, 1.0f, 1);

    }

    // 求取最终的步长
    for (uint m = 0; m < Pa.Nt * ip->St[0].rn; m++)
    {
        fenzi += plan->h_SumResCurr[m] * plan->h_SumResTrial[m];
        fenmu += plan->h_SumResTrial[m] * plan->h_SumResTrial[m];
    }

    ip->Alpha = (fenzi / fenmu) * e;
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
    uint nnz = Pa.Nz + 2 * Pa.PMLz;
    uint nnx = Pa.Nx + 2 * Pa.PMLx;

    memcpy(plan->h_Grad, ip->GradVp, sizeof(float) * Pa.Nz * Pa.Nx);
    memcpy(plan->h_Vp, ip->CurrVp, sizeof(float) * nnz * nnx);

    // 更新速度
    UpdateVp(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);
    UpdateVpPML(Pa, plan->h_Vp, plan->h_Grad, ip->Alpha, pt);

    memcpy(ip->CurrVp, plan->h_Vp, sizeof(float) * nnz * nnx);

    ip->Alpha = 0.0f;
}

