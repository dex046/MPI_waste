/**************************************************************************
* 2维时间域全波形反演（Time domain Full Waveform Inversion）

* 使用非分裂完全匹配曾（NPML）技术处理吸收边界

* 使用2阶位移运动方程

* 正演使用空间8阶时间2阶精度的交错网格有限差分技术

* 作者：高照奇
**************************************************************************/

#ifndef		_TDFWI_
#define		_TDFWI_

#include	<iostream>
#include	<stdio.h>
#include	<string.h>
#include	<math.h>
#include	<stdlib.h>
#include	<time.h>
#include	"RWsgy.h"

#define		uint			unsigned int
#define		usht			unsigned short
#define		PI				3.141592653589793f
#define		PowerGrad		2.0f

// 空间8阶交错网格有限差分系数
#define		C1_4			1.196289062506247f
#define		C2_4			-0.079752604166839f
#define		C3_4			0.009570312500021f
#define		C4_4			-0.000697544642859f

// 震源位置
typedef struct{
	uint Sz;
	uint Sx;
} SL;

// 检波器位置
typedef struct{
	uint Rz;
	uint Rx;
} RL;

// 每一炮相关信息
struct Shot{
	SL s;				// 震源位置
	RL *re;				// 检波器位置数组
	float *wf_t;		// 真实模型对应的正演数据(观测数据)
	float *wf_c;		// 当前模型对应的正演数据(合成数据)
	float *wf_r;		// 残差波场
	uint rn;			// 检波器个数
};

// 有限差分正演参数
struct AFDPU2D{
	uint Nx;			// 横向的网格数
	uint Nz;			// 纵向网格数
	uint Nt;			// 正演的时间步数
	float dx;			// 横向网格间距
	float dz;			// 纵向网格间距
	float dt;			// 采样间隔
	float f0;			// 正演子波的主频
	uint PMLx;			// 横向NPML内的网格数
	uint PMLz;			// 纵向NPML内的网格数
};

// 反演中使用到的全局变量
struct CPUVs{
	
	// 与NPML有关的变量
	float *h_dx;			// NPML横向参数
	float *h_dz;			// NPML纵向参数
	float *h_Bx;			// NPML横向参数
	float *h_Bz;			// NPML纵向参数
	float *h_PHIx_V_x;		// 波场变量V横向的卷积系数 
	float *h_PHIz_W_z;		// 波场变量W纵向的卷积系数 
	float *h_PHIx_U_x;		// 波场变量U横向的卷积系数
	float *h_PHIz_U_z;		// 波场变量U纵向的卷积系数
	float *h_PHIx_V_x_r;	// 波场变量V横向的卷积系数(反传)
	float *h_PHIz_W_z_r;	// 波场变量W纵向的卷积系数(反传)
	float *h_PHIx_U_x_r;	// 波场变量U横向的卷积系数(反传)
	float *h_PHIz_U_z_r;	// 波场变量U纵向的卷积系数(反传)

	// 波场相关变量
	float *h_U_next;		// (n+1)时刻的U波场快照
	float *h_U_now;			// (n)时刻的U波场快照
	float *h_U_past;		// (n-1)时刻的U波场快照
	float *h_V;				// (n)时刻的V波场快照
	float *h_W;				// (n)时刻的W波场快照
	float *h_U_next_r;		// (n+1)时刻的U波场快照(反传)
	float *h_U_now_r;		// (n)时刻的U波场快照(反传)
	float *h_U_past_r;		// (n-1)时刻的U波场快照(反传)
	float *h_V_r;			// (n)时刻的V波场快照(反传)
	float *h_W_r;			// (n)时刻的W波场快照(反传)

	float *h_Vp;			// 正演中使用的速度
	RL *h_re;				// 检波器位置 
	
	uint *h_Coord;			// 有效边界存储策略需要存储的点的坐标
	float *h_RU;			// 存储的波场

	float *h_TrueWF;		// 观测数据炮集
	float *h_CurrWF;		// 正演数据炮集
	float *h_ResWF;			// 残差数据炮集

	float *h_Grad;			// 迭代过程中的梯度
	float *h_U_Der;			// 正传波场对时间的2阶导

	float *h_Obj;			// 目标函数

	// 计算步长的相关变量
	float *h_TrailWF;		// 试探模型对应的波场
	float *h_SumResTrial;	// 求取步长的临时变量
	float *h_SumResCurr;	// 求取步长的临时变量
};

// 反演参数 
struct IP{
	uint IterN;				// 反演中迭代次数
	uint ShotN;				// 反演中的炮数
	Shot *St;				// 炮信息
	float *TrueVp;			// 真实速度模型
	float *CurrVp;			// 当前速度模型
	float *GradVp;			// 梯度
	float Alpha;			// 迭代步长
	float *ObjIter;			// 迭代目标函数
};

/*------------------------ 反演程序中的函数 --------------------------*/
// Ricker子波
float Ricker(float f0,
			float t);

// 从Sgy文件中读取数据
void ReadData(char FileName[],
			float *Data,
			usht flag);

// 将数据写入Sgy文件
void WriteData(char FileName[],
			usht SampleNum,
			usht TraceNum,
			usht SampleInt,
			float *data,
			usht flag);

// 在内存上为变量申请空间
void MallocVariables(AFDPU2D Pa,
					IP *ip,
					CPUVs *plan);

// 生成NPML系数
void GenerateNPML(AFDPU2D Pa,
				CPUVs *plan);

// 找出矩阵中的最大值
float MatMax(float *Mat,
			uint Len);

// 正演加震源
void AddSource(AFDPU2D Pa,
				float *h_U,
				SL s,
				float Wave,
				float *h_Vp);

// 一步更新波场U的卷积项
void StepPHIU(AFDPU2D Pa,
			float *h_U,
			float *h_PHIx_U_x,
			float *h_PHIz_U_z,
			float *h_Bx,
			float *h_Bz);

// 一步更新波场V和W的卷积项
void StepPHIVW(AFDPU2D Pa,
			float *h_V,
			float *h_W,
			float *h_PHIx_V_x,
			float *h_PHIz_W_z,
			float *h_Bx,
			float *h_Bz);

// 一步更新波场U
void StepU(AFDPU2D Pa,
			float *h_U_next,
			float *h_U_now,
			float *h_U_past,
			float *h_V,
			float *h_W,
			float *h_PHIx_V_x,
			float *h_PHIz_W_z,
			float *h_Bx,
			float *h_Bz,
			float *h_Vp);

// 一步更新波场V和W
void StepVW(AFDPU2D Pa,
			float *h_U,
			float *h_V,
			float *h_W,
			float *h_PHIx_U_x,
			float *h_PHIz_U_z,
			float *h_Bx,
			float *h_Bz);

// 一步记录检波器位置的波场值 
void StepShotGather(AFDPU2D Pa,
					float *h_U,
					float *h_SG,
					RL *h_re,
					uint nstep,
					uint rn);

// 一步记录有效边界存储策略中需要记录的波场值
void StepRecordU(AFDPU2D Pa,
				uint *h_Coord,
				uint RecNum,
				float *h_U,
				float *h_RU);

// 一步替换有效边界内的波场 
void StepReplaceU(AFDPU2D Pa,
				uint *h_Coord,
				uint RecNum,
				float *h_U,
				float *h_RU);

// 一步逆时间更新波场V和W
void StepRtVW(AFDPU2D Pa,
			float *h_U,
			float *h_V,
			float *h_W,
			float *h_PHIx_U_x,
			float *h_PHIz_U_z,
			float *h_Bx,
			float *h_Bz);

// 一步逆时间更新波场U
void StepRtU(AFDPU2D Pa,
			float *h_U_next,
			float *h_U_now,
			float *h_U_past,
			float *h_V,
			float *h_W,
			float *h_PHIx_V_x,
			float *h_PHIz_W_z,
			float *h_Bx,
			float *h_Bz,
			float *h_Vp);

// 一步计算正传波场对时间的2阶导
void StepCal2Der(AFDPU2D Pa,
				float *h_U_past,
				float *h_U_now,
				float *h_U_next,
				float *h_U_Der);

// 一步求取梯度
void StepCalGrad(AFDPU2D Pa,
				float *h_Grad,
				float *h_U_Der,
				float *h_U_r);

// 求取残差波场
void StepResidual(AFDPU2D Pa,
				float *h_TrueWF,
				float *h_CurrWF,
				float *h_ResWF,
				uint Rn);

// 残差反传加震源
void AddResidual(AFDPU2D Pa,
				float *h_ResWf,
				float *h_U_r,
				RL *h_re,
				uint nstep,
				uint Rn,
				float *h_Vp);

// 梯度后处理 
void PostProcessGrad(AFDPU2D Pa,
					float *h_Grad,
					float *h_Vp);

// 更新内部网格内的速度
void UpdateVp(AFDPU2D Pa,
			float *h_Vp,
			float *h_Grad,
			float e);

// 更新PML网格内的速度
void UpdateVpPML(AFDPU2D Pa,
			float *h_Vp,
			float *h_Grad,
			float e);

// 给有效边界存储策略需要存储的波场坐标赋值
void SetCoord(AFDPU2D *Pa,
			uint *h_Coord);

// 矩阵加减法
void MatAdd(float *Mat,
			float *Mat1,
			float *Mat2,
			uint row,
			uint col,
			float coeff1,
			float coeff2,
			uint flag);

// 矩阵点对点乘除法
void MatMul(float *Mat,
			float *Mat1,
			float *Mat2,
			uint row,
			uint col,
			float coeff1,
			float coeff2,
			uint flag);

// 求取观测波场 
void CalTrueWF(AFDPU2D Pa,
			IP *ip,
			CPUVs *plan,
			float *sgs_t);

// 求取梯度
void CalGrad(AFDPU2D Pa,
			IP *ip,
			CPUVs *plan,
			float *sgs_t,
			float *sgs_c,
			float *sgs_r,
			uint It);

// 求取步长
void CalStepLength(AFDPU2D Pa, 
				IP *ip,
				CPUVs *plan,
				float *sgs_t,
				float *sgs_c,
				float e);

// 下次迭代前的预处理 
void PreProcess(AFDPU2D Pa,
				IP *ip,
				CPUVs *plan);

#endif
