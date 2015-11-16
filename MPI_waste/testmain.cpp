/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#include "mpi.h"
#include "testTDFWI.h"
#include "RWsgy.h"
#include "Partition.h"


using namespace std;

int main(int argc, char ** argv)
{
    int rank, p_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p_size);

    //cout << 111 << endl;
    // 有限差分正演参数
    AFDPU2D *Pa;
    Pa = new AFDPU2D[1];
    memset((void *)Pa, 0, sizeof(AFDPU2D));
    Pa->dt = 0.001f;// 采样间隔
    Pa->dx = 10.0f;// 横向网格间距
    Pa->dz = 10.0f;// 纵向网格间距
    Pa->f0 = 10.0f;// 正演子波的主频
    Pa->Nt = 3000;// 正演的时间步数
    Pa->Nx = 510;// 横向的网格数
    Pa->Nz = 134;// 纵向网格数
    Pa->PMLx = 50;// 横向NPML内的网格数
    Pa->PMLz = 50;// 纵向NPML内的网格数
    uint nnx = Pa->Nx + 2 * Pa->PMLx;
    uint nnz = Pa->Nz + 2 * Pa->PMLz;


    H_Border temph_U(3, 3, 4, 4);//top left bottom right
    H_Border temph_VW(4, 4, 3, 3);

    // 反演参数
    IP *ip;
    ip = new IP[1];
    memset((void *)ip, 0, sizeof(IP));
    ip->ShotN = 1;
    ip->IterN = 1;
    ip->Alpha = 0.0f;

    // 炮信息
    ip->St = new Shot[ip->ShotN];
    memset((void *)ip->St, 0, ip->ShotN * sizeof(Shot));
    for(uint i = 0; i < ip->ShotN; ++i)
    {
        ip->St[i].rn = 510;
    }

    uint cpu_x = 2, cpu_z = 2;
    Partition pt(Pa, ip, nnx, nnz, cpu_x, cpu_z, temph_U, temph_VW, 8, rank, p_size);//borderlength = 4


    uint length_x = pt.getblockLength_x();
    uint length_z = pt.getblockLength_z();

    //cout << "rank=" << rank << "aaaa" << length_x << " " << length_z << endl;

    uint interiorlength_x = pt.getinteriorLength_x();
    uint interiorlength_z = pt.getinteriorLength_z();

    //MPI_Barrier(MPI_COMM_WORLD);
    if(rank == 0)
    cout << rank << " " << length_x << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "rank=" << rank << " " << interiorlength_z << endl;

    uint RL_num = pt.getRL_num();

    //cout << "rank=" << rank << " " << pt.getRL_beginnum() << "wwww" << endl;
    uint indexmin_x = pt.getindexmin_x();
    uint indexmax_x = pt.getindexmax_x();
    uint indexmin_z = pt.getindexmin_z();
    uint indexmax_z = pt.getindexmax_z();

    try
    {
        ip->TrueVp = new float[length_x * length_z];
        memset((void *)ip->TrueVp, 0, sizeof(float) * length_x * length_z);
        ip->CurrVp = new float[length_x * length_z];
        memset((void *)ip->CurrVp, 0, sizeof(float) * length_x * length_z);

        ip->GradVp = new float[interiorlength_z * interiorlength_x];
        memset((void *)ip->GradVp, 0, sizeof(float) * interiorlength_z * interiorlength_x);
        ip->ObjIter = new float[ip->IterN];
        memset((void *)ip->ObjIter, 0, sizeof(float) * ip->IterN);
    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: " << e.what() << endl;
    }
    //cout << 111 << endl;
    char TrueVp[] = "TrueVp-510-134-marmousi.sgy";
    char InitVp[] = "InitVp-510-134-marmousi.sgy";

    // 读取速度信息
    ReadData(TrueVp, ip->TrueVp, pt, 0);
    ReadData(InitVp, ip->CurrVp, pt, 0);
//    for(int i = 0; i < length_x * length_z; ++i)
//        cout << ip->TrueVp[i];

    // 反演中使用到的全局变量
    CPUVs *plan;
    plan = new CPUVs[1];
    memset((void *)plan, 0, sizeof(CPUVs));

    float *sgs_t = NULL, *sgs_c = NULL, *sgs_r = NULL;

    try
    {
        sgs_t = new float[ip->ShotN * Pa->Nt * RL_num];// 反演中的炮数 * 正演的时间步数 * 检波器个数
        memset((void *)sgs_t, 0,
            sizeof(float) * ip->ShotN * Pa->Nt * RL_num);
        sgs_c = new float[ip->ShotN * Pa->Nt * RL_num];
        memset((void *)sgs_c, 0,
            sizeof(float) * ip->ShotN * Pa->Nt * RL_num);
        sgs_r = new float[ip->ShotN * Pa->Nt * RL_num];
        memset((void *)sgs_r, 0,
            sizeof(float) * ip->ShotN * Pa->Nt * RL_num);
    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: " << e.what() << endl;
    }

    // 给全局变量开辟空间
    MallocVariables(*Pa, ip, plan, pt);

    // 求取NPML的参数
    GenerateNPML(*Pa, plan, pt);

    // 给有效边界的坐标赋值
    SetCoord(Pa, plan->h_Coord, pt);// h_Coord 有效边界存储策略需要存储的点的坐标

    if(rank == ROOT_ID)
    {
        cout << "********************************************************************************************" << endl;
        cout << "Doing Time domain Full Waveform Inversion ..." << endl;
        cout << "Parameters of Inversion are as follows:" << endl;
        cout << "\tNx = " << Pa->Nx << endl;
        cout << "\tNz = " << Pa->Nz << endl;
        cout << "\tdx = " << Pa->dx << "m" << endl;
        cout << "\tdz = " << Pa->dz << "m" << endl;
        cout << "\tNt = " << Pa->Nt << endl;
        cout << "\tdt = " << Pa->dt << "s" << endl;
        cout << "\tNpml = " << Pa->PMLx << endl;
        cout << "\tf0 = " << Pa->f0 << endl;
        cout << "\tNshot = " << ip->ShotN << endl;
        cout << "\tIteration number = " << ip->IterN << endl;
        cout << "********************************************************************************************" << endl;
    }

    clock_t begin, duration;

    // 求取观测波场
    if(rank == ROOT_ID)
        cout << "\tCalculating the observed data..." << endl;

    begin = clock();
//cout << rank << endl;
    CalTrueWF(*Pa, ip, plan, sgs_t, pt);//求取观测波场 plan  sgs_t
    duration = clock() - begin;

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == ROOT_ID)
        cout << "\tCalculating the observed data used:\t" << duration / CLOCKS_PER_SEC << "s" << endl;

    // 迭代过程中为了求取步长使用的试探步长
    float e = 24.0f;


    for (uint It = 0; It < ip->IterN; It++)
    {
        // 求取梯度
        if(rank == ROOT_ID)
        {
            cout << "\n\tDoing the " << It << "th iteration" << endl;
            cout << "\tCalculating the Gradient..." << endl;
        }

        begin = clock();
        CalGrad(*Pa, ip, plan, sgs_t, sgs_c, sgs_r, It, pt);

        MPI_Barrier(MPI_COMM_WORLD);

        //cout << ip->ObjIter[It] << endl;
        // 梯度后处理
        PostProcessGrad(*Pa, ip->GradVp, plan->h_Vp, pt);

        // 求取步长
        CalStepLength(*Pa, ip, plan, sgs_t, sgs_c, e, pt);
        duration = clock() - begin;

        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == ROOT_ID)
        {
            cout << "\tObjective function value:\t" << ip->ObjIter[It] << endl;
            cout << "\tStep length:\t" << ip->Alpha << endl;
            cout << "\tThe " << It << "th iteration used " << duration / CLOCKS_PER_SEC << "s" << endl;
        }

        // 下一步迭代预处理
        //PreProcess(*Pa, ip, plan, pt);
    }

    if(rank == ROOT_ID)
    {
        cout << "\tWriting data to .sgy1" << endl;
    }

    const char * const TrueSg = "TrueSG.sgy";
    const char * const GradVp = "GradientVp.sgy";
    const char * const InvertedVp = "InvertedVp.sgy";

    fopen(TrueSg, "wb");
    fopen(GradVp, "wb");
    fopen(InvertedVp, "wb");

    write_sgs_t_Data(TrueSg, (usht)Pa->Nt, ip->St[0].rn, (usht)(Pa->dt * 1000000), sgs_t, pt, *Pa, 1);
    WriteData(GradVp, Pa->Nz, Pa->Nx, Pa->dz * 1000, ip->GradVp, pt, *Pa, 0, WRITE_INTER);
    WriteData(InvertedVp, nnz, nnx, Pa->dz * 1000, ip->CurrVp, pt, *Pa, 0, WRITE_ALL);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}


