/******************************************
 * author:dwx
 * 2015.10.15
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

    ip->St = new Shot[ip->ShotN];
    memset((void *)ip->St, 0, ip->ShotN * sizeof(Shot));
    for(uint i = 0; i < ip->ShotN; ++i)
    {
        ip->St[i].rn = 510;
    }

//cout << 24324 << endl;
    uint cpu_x = 2, cpu_z = 2;
    Partition pt(Pa, ip, nnx, nnz, cpu_x, cpu_z, temph_U, temph_VW, 8, rank, p_size);//borderlength = 4
    //AFDPU2D &Pa, IP &ip, int totallength_x, int totallength_z, int sumblock_x, int sumblock_z, H_U h_U, H_VW h_VW, int rank, int size

    uint length_x = pt.getblockLength_x();
    uint length_z = pt.getblockLength_z();

    uint RL_num = pt.getRL_num();

    uint indexmin_x = pt.getindexmin_x();
    uint indexmax_x = pt.getindexmax_x();
    uint indexmin_z = pt.getindexmin_z();
    uint indexmax_z = pt.getindexmax_z();

//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "rank = " << rank << " " << pt.isfirstblock_x() << pt.isfirstblock_z() << pt.islastblock_x() << pt.islastblock_z() << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "x = " << indexmin_x << " " << indexmax_x << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "z = " << indexmin_z << " " << indexmax_z << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
    try
    {
        ip->TrueVp = new float[length_x * length_z];
        memset((void *)ip->TrueVp, 0, sizeof(float) * length_x * length_z);
        ip->CurrVp = new float[length_x * length_z];
        memset((void *)ip->CurrVp, 0, sizeof(float) * length_x * length_z);

    //    uint inside_length_x = pt.getinteriorLength_x();
    //    uint inside_length_z = pt.getinteriorLength_z();

        ip->GradVp = new float[Pa->Nz * Pa->Nx];
        memset((void *)ip->GradVp, 0, sizeof(float) * Pa->Nz * Pa->Nx);
        ip->ObjIter = new float[ip->IterN];
        memset((void *)ip->ObjIter, 0, sizeof(float) * ip->IterN);
    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: " << e.what() << endl;
    }


    // 炮信息




    char TrueVp[] = "TrueVp-510-134-marmousi.sgy";
    char InitVp[] = "InitVp-510-134-marmousi.sgy";
    // 读取速度信息
    ReadData(TrueVp, ip->TrueVp, pt, 0);
    ReadData(InitVp, ip->CurrVp, pt, 0);

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


    try
    {

    }
    catch(const std::bad_alloc& e)
    {
        cout << "Allocation failed: " << e.what() << endl;
    }

    // 给全局变量开辟空间
    MallocVariables(*Pa, ip, plan, pt);////


//cout << pt.getrank() << endl;
    // 求取NPML的参数
    GenerateNPML(*Pa, plan, pt);

    // 给有效边界的坐标赋值
    SetCoord(Pa, plan->h_Coord, pt);// h_Coord 有效边界存储策略需要存储的点的坐标
//    cout << "********************************************************************************************" << endl;
//    cout << "Doing Time domain Full Waveform Inversion ..." << endl;
//    cout << "Parameters of Inversion are as follows:" << endl;
//    cout << "\tNx = " << Pa->Nx << endl;
//    cout << "\tNz = " << Pa->Nz << endl;
//    cout << "\tdx = " << Pa->dx << "m" << endl;
//    cout << "\tdz = " << Pa->dz << "m" << endl;
//    cout << "\tNt = " << Pa->Nt << endl;
//    cout << "\tdt = " << Pa->dt << "s" << endl;
//    cout << "\tNpml = " << Pa->PMLx << endl;
//    cout << "\tf0 = " << Pa->f0 << endl;
//    cout << "\tNshot = " << ip->ShotN << endl;
//    cout << "\tIteration number = " << ip->IterN << endl;
//    cout << "********************************************************************************************" << endl;
    clock_t begin, duration;

    // 求取观测波场
//    cout << "\tCalculating the observed data..." << endl;

    begin = clock();



    CalTrueWF(*Pa, ip, plan, sgs_t, pt);//求取观测波场 plan  sgs_t
    duration = clock() - begin;

//    cout << "\tCalculating the observed data used:\t" << duration / CLOCKS_PER_SEC << "s" << endl;

    // 迭代过程中为了求取步长使用的试探步长
    float e = 24.0f;


    for (uint It = 0; It < ip->IterN; It++)
    {
        // 求取梯度
//        cout << "\n\tDoing the " << It << "th iteration" << endl;
        cout << "\tCalculating the Gradient..." << endl;

        begin = clock();
        CalGrad(*Pa, ip, plan, sgs_t, sgs_c, sgs_r, It, pt);
//cout << "********************************************************************************************" << endl;
        // 梯度后处理
        PostProcessGrad(*Pa, ip->GradVp, plan->h_Vp, pt);

        // 求取步长
        CalStepLength(*Pa, ip, plan, sgs_t, sgs_c, e, pt);
//cout << "aaaaaaaaaaaaaaa" << endl;
        duration = clock() - begin;

//        cout << "\tObjective function value:\t" << ip->ObjIter[It] << endl;
//        cout << "\tStep length:\t" << ip->Alpha << endl;
//        cout << "\tThe " << It << "th iteration used " << duration / CLOCKS_PER_SEC << "s" << endl;

        // 下一步迭代预处理
        PreProcess(*Pa, ip, plan, pt);
    }
//cout << rank << endl;
    cout << "\tWriting data to .sgy1" << endl;



//    char TrueSg[255];
//    char GradVp[255];
//    char InvertedVp[255];
//    sprintf(TrueSg, "TrueSG.sgy");
//    sprintf(GradVp, "GradientVp.sgy");
//    sprintf(InvertedVp, "InvertedVp.sgy");

    const char * const TrueSg = "TrueSG.sgy";
    const char * const GradVp = "GradientVp.sgy";
    const char * const InvertedVp = "InvertedVp.sgy";

    fopen(TrueSg, "wb");
    fopen(GradVp, "wb");
    fopen(InvertedVp, "wb");
    //WriteData(TrueSg, (usht)Pa->Nt, (usht)ip->St[0].rn, (usht)(Pa->dt * 1000000), sgs_t, 1);


//cout << "rank1" << endl;

    if(RL_num)
    {
//        MPI_Barrier(MPI_COMM_WORLD);
//        cout << "rank = " << rank << endl;
//        cout << " RL_num= " << RL_num << endl;
//        MPI_Barrier(MPI_COMM_WORLD);
        cout << rank << endl;
        write_sgs_t_Data(TrueSg, (usht)Pa->Nt, ip->St[0].rn, (usht)(Pa->dt * 1000000), sgs_t, pt, *Pa, 1);
        cout << rank << endl;//0 mei
    }

//cout << "rank3oo" << endl;

    WriteData(GradVp, Pa->Nz, Pa->Nx, Pa->dz * 1000, ip->GradVp, pt, *Pa, 0, WRITE_INTER);
cout << "rank4" << endl;
    WriteData(InvertedVp, nnz, nnx, Pa->dz * 1000, ip->CurrVp, pt, *Pa, 0, WRITE_ALL);
cout << "rank5" << endl;


    MPI_Finalize();
    return 0;
}


