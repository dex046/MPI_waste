/**************************************************************************
* 2维时间域全波形反演（Time domain Full Waveform Inversion）

* 使用非分裂完全匹配曾（NPML）技术处理吸收边界

* 使用2阶位移运动方程

* 正演使用空间8阶时间2阶精度的交错网格有限差分技术

* 作者：高照奇
**************************************************************************/

#include	"TDFWI.h"
#include        "RWsgy.h"
#include        "RWsgy.cpp"
#include        "TDFWI.cpp"

using namespace std;

int main()
{
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

	// 反演参数 
	IP *ip;
	ip = new IP[1];
	memset((void *)ip, 0, sizeof(IP));
    ip->ShotN = 2;
	ip->IterN = 1;
	ip->Alpha = 0.0f;
	ip->TrueVp = new float[nnz * nnx];	
	memset((void *)ip->TrueVp, 0, sizeof(float) * nnz * nnx);
	ip->CurrVp = new float[nnz * nnx];
	memset((void *)ip->CurrVp, 0, sizeof(float) * nnz * nnx);
	ip->GradVp = new float[Pa->Nz * Pa->Nx];
	memset((void *)ip->GradVp, 0, sizeof(float) * Pa->Nz * Pa->Nx);
	ip->ObjIter = new float[ip->IterN];
	memset((void *)ip->ObjIter, 0, sizeof(float) * ip->IterN);
	ip->St = new Shot[ip->ShotN];
	memset((void *)ip->St, 0, ip->ShotN * sizeof(Shot));

	// 炮信息
	for (uint is = 0; is < ip->ShotN; is++)
	{
		ip->St[is].rn = 510;
        ip->St[is].s.Sx = Pa->PMLx + is * 10;
		ip->St[is].s.Sz = Pa->PMLz + 2;
		ip->St[is].re = new RL[ip->St[is].rn];
		memset((void *)ip->St[is].re, 0, sizeof(RL) * ip->St[is].rn);

		// 检波器位置
		for (uint m = 0; m < ip->St[is].rn; m++)
		{
			ip->St[is].re[m].Rx = Pa->PMLx + 1 * m;
			ip->St[is].re[m].Rz = Pa->PMLz + 2;
		}
	}

	// 读取速度信息 
	char TrueVp[] = "TrueVp-510-134-marmousi.sgy";
	char InitVp[] = "InitVp-510-134-marmousi.sgy";

	ReadData(TrueVp, ip->TrueVp, 0);
	ReadData(InitVp, ip->CurrVp, 0);
//    for(int i = 0; i < nnx * nnz; ++i)
//        cout << ip->TrueVp[i];
    // 反演中使用到的全局变量
	CPUVs *plan;
	plan = new CPUVs[1];
	memset((void *)plan, 0, sizeof(CPUVs));

	float *sgs_t, *sgs_c, *sgs_r;
    sgs_t = new float[ip->ShotN * Pa->Nt * ip->St[0].rn];// 反演中的炮数 * 正演的时间步数 * 检波器个数
	memset((void *)sgs_t, 0, 
		sizeof(float) * ip->ShotN * Pa->Nt * ip->St[0].rn);
	sgs_c = new float[ip->ShotN * Pa->Nt * ip->St[0].rn];
	memset((void *)sgs_c, 0, 
		sizeof(float) * ip->ShotN * Pa->Nt * ip->St[0].rn);
	sgs_r = new float[ip->ShotN * Pa->Nt * ip->St[0].rn];
	memset((void *)sgs_r, 0, 
		sizeof(float) * ip->ShotN * Pa->Nt * ip->St[0].rn);

	// 给全局变量开辟空间
	MallocVariables(*Pa, ip, plan);

	// 求取NPML的参数 
	GenerateNPML(*Pa, plan);

	// 给有效边界的坐标赋值
    SetCoord(Pa, plan->h_Coord);// h_Coord 有效边界存储策略需要存储的点的坐标

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
	
	clock_t begin, duration;

	// 求取观测波场 
	cout << "\tCalculating the observed data..." << endl;
	
	begin = clock();
    CalTrueWF(*Pa, ip, plan, sgs_t);//求取观测波场
	duration = clock() - begin;
	
	cout << "\tCalculating the observed data used:\t" << duration / CLOCKS_PER_SEC << "s" << endl;
	
	// 迭代过程中为了求取步长使用的试探步长
    float e = 24.0f;

    for (uint It = 0; It < ip->IterN; It++)
    {
        // 求取梯度
        cout << "\n\tDoing the " << It << "th iteration" << endl;
        cout << "\tCalculating the Gradient..." << endl;
        begin = clock();
        CalGrad(*Pa, ip, plan, sgs_t, sgs_c, sgs_r, It);
		
        // 梯度后处理
        PostProcessGrad(*Pa, ip->GradVp, plan->h_Vp);

        // 求取步长
        CalStepLength(*Pa, ip, plan, sgs_t, sgs_c, e);
        duration = clock() - begin;

        cout << "\tObjective function value:\t" << ip->ObjIter[It] << endl;
        cout << "\tStep length:\t" << ip->Alpha << endl;
        cout << "\tThe " << It << "th iteration used " << duration / CLOCKS_PER_SEC << "s" << endl;

        // 下一步迭代预处理
        PreProcess(*Pa, ip, plan);
    }

//    ofstream fout0("RGradVp0.txt");
//    for(int i = 0; i < 67; ++i)
//    {
//        for(int j = 0; j < 255; ++j)
//        {
//            fout0 << *(ip->GradVp + i * Pa->Nx + j) << " ";
//        }
//        fout0 << endl;
//    }

//    ofstream fout1("RGradVp1.txt");
//    for(int i = 0; i < 67; ++i)
//    {
//        for(int j = 255; j < 510; ++j)
//        {
//            fout1 << *(ip->GradVp + i * Pa->Nx + j) << " ";
//        }
//        fout1 << endl;
//    }

//    ofstream fout2("RGradVp2.txt");
//    for(int i = 67; i < 134; ++i)
//    {
//        for(int j = 0; j < 255; ++j)
//        {
//            fout2 << *(ip->GradVp + i * Pa->Nx + j) << " ";
//        }
//        fout2 << endl;
//    }

//    ofstream fout3("RGradVp3.txt");
//    for(int i = 67; i < 134; ++i)
//    {
//        for(int j = 255; j < 510; ++j)
//        {
//            fout3 << *(ip->GradVp + i * Pa->Nx + j) << " ";
//        }
//        fout3 << endl;
//    }

    ofstream fout4("RCurrVp0.txt", ios_base::out | ios_base::trunc);
    for(int i = 0; i < 117; ++i)
    {
        for(int j = 0; j < 305; ++j)
        {
            fout4 << *(ip->CurrVp + i * nnx + j) << " ";
        }
        fout4 << endl;
    }

    ofstream fout5("RCurrVp1.txt", ios_base::out | ios_base::trunc);
    for(int i = 0; i < 117; ++i)
    {
        for(int j = 305; j < 610; ++j)
        {
            fout5 << *(ip->CurrVp + i * nnx + j) << " ";
        }
        fout5 << endl;
    }

    ofstream fout6("RCurrVp2.txt", ios_base::out | ios_base::trunc);
    for(int i = 117; i < 234; ++i)
    {
        for(int j = 0; j < 305; ++j)
        {
            fout6 << *(ip->CurrVp + i * nnx + j) << " ";
        }
        fout6 << endl;
    }

    ofstream fout7("RCurrVp3.txt", ios_base::out | ios_base::trunc);
    for(int i = 117; i < 234; ++i)
    {
        for(int j = 305; j < 610; ++j)
        {
            fout7 << *(ip->CurrVp + i * nnx + j) << " ";
        }
        fout7 << endl;
    }

    cout << "\tWriting data to .sgy" << endl;

    char TrueSg[255];
    char GradVp[255];
    char InvertedVp[255];
    sprintf(TrueSg, "TrueSG.sgy");
    sprintf(GradVp, "GradientVp.sgy");
    sprintf(InvertedVp, "InvertedVp.sgy");

    WriteData(TrueSg, (usht)Pa->Nt, (usht)ip->St[0].rn, (usht)(Pa->dt * 1000000), sgs_t, 1);
    WriteData(GradVp, Pa->Nz, Pa->Nx, Pa->dz * 1000, ip->GradVp, 0);
    WriteData(InvertedVp, nnz, nnx, Pa->dz * 1000, ip->CurrVp, 0);

	return 0;
}
