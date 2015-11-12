/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#ifndef TESTTDFWI_H
#define TESTTDFWI_H

#include	<iostream>
#include	<stdio.h>
#include	<string.h>
#include	<math.h>
#include	<stdlib.h>
#include	<time.h>
#include	"RWsgy.h"
#include    "Partition.h"
#include    "DataTran.h"





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



/*------------------------ 反演程序中的函数 --------------------------*/
// Ricker子波
float Ricker(float f0,
            float t);

// 从Sgy文件中读取数据
void ReadData(char FileName[],
            float *Data, const Partition &pt,
            usht flag);

void write_sgs_t_Data(const char * const FileName, usht SampleNum, usht TraceNum, usht SampleInt, float *data, const Partition &pt, const AFDPU2D &Pa, usht flag);


// 将数据写入Sgy文件
void WriteData(const char * const FileName,
            usht SampleNum,
            usht TraceNum,
            usht SampleInt,
            float *data, const Partition &pt, const AFDPU2D Pa,
            usht flag, unsigned short tag);

// 在内存上为变量申请空间
void MallocVariables(AFDPU2D Pa,
                    IP *ip,
                    CPUVs *plan, const Partition &pt);

// 生成NPML系数
void GenerateNPML(AFDPU2D Pa,
                CPUVs *plan, const Partition &pt);

// 找出矩阵中的最大值
float MatMax(float *Mat,
            uint Len);

// 正演加震源
void AddSource(AFDPU2D Pa,
                float *h_U,
                float Wave,
                float *h_Vp, const Partition &pt);

// 一步更新波场U的卷积项
void StepPHIU(AFDPU2D Pa,
            float *h_U,
            float *h_PHIx_U_x,
            float *h_PHIz_U_z,
            float *h_Bx,
            float *h_Bz, const Partition &pt);

// 一步更新波场V和W的卷积项
void StepPHIVW(AFDPU2D Pa,
            float *h_V,
            float *h_W,
            float *h_PHIx_V_x,
            float *h_PHIz_W_z,
            float *h_Bx,
            float *h_Bz, const Partition &pt);

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
            float *h_Vp, const Partition &pt);

// 一步更新波场V和W
void StepVW(AFDPU2D Pa,
            float *h_U,
            float *h_V,
            float *h_W,
            float *h_PHIx_U_x,
            float *h_PHIz_U_z,
            float *h_Bx,
            float *h_Bz, const Partition &pt);

// 一步记录检波器位置的波场值
void StepShotGather(AFDPU2D Pa,
                    float *h_U,
                    float *h_SG,
                    uint nstep,
                    const Partition &pt);

// 一步记录有效边界存储策略中需要记录的波场值
void StepRecordU(AFDPU2D Pa,
                uint *h_Coord,
                uint RecNum,
                float *h_U,
                float *h_RU, const Partition &pt);

// 一步替换有效边界内的波场
void StepReplaceU(AFDPU2D Pa,
                uint *h_Coord,
                uint RecNum,
                float *h_U,
                float *h_RU, const Partition &pt);

// 一步逆时间更新波场V和W
void StepRtVW(AFDPU2D Pa,
            float *h_U,
            float *h_V,
            float *h_W,
            float *h_PHIx_U_x,
            float *h_PHIz_U_z, const Partition &pt);

// 一步逆时间更新波场U
void StepRtU(AFDPU2D Pa,
            float *h_U_next,
            float *h_U_now,
            float *h_U_past,
            float *h_V,
            float *h_W,
            float *h_PHIx_V_x,
            float *h_PHIz_W_z,
            float *h_Vp,
             const Partition &pt);

// 一步计算正传波场对时间的2阶导
//extern class Partition &pt;
void StepCal2Der(AFDPU2D Pa,
                float *h_U_past,
                float *h_U_now,
                float *h_U_next,
                float *h_U_Der, const Partition &pt);

// 一步求取梯度
void StepCalGrad(AFDPU2D Pa,
                float *h_Grad,
                float *h_U_Der,
                float *h_U_r, const Partition &pt);

// 求取残差波场
void StepResidual(AFDPU2D Pa,
                float *h_TrueWF,
                float *h_CurrWF,
                float *h_ResWF,
                int RL_num);

// 残差反传加震源
void AddResidual(AFDPU2D Pa,
                float *h_ResWf,
                float *h_U_r,
                uint nstep,
                float *h_Vp, const Partition &pt);

// 梯度后处理
void PostProcessGrad(AFDPU2D Pa,
                    float *h_Grad,
                    float *h_Vp, const Partition &pt);

// 更新内部网格内的速度
void UpdateVp(AFDPU2D Pa,
            float *h_Vp,
            float *h_Grad,
            float e, const Partition &pt);

// 更新PML网格内的速度
void UpdateVpPML(AFDPU2D Pa,
            float *h_Vp,
            float *h_Grad,
            float e, const Partition &pt);

// 给有效边界存储策略需要存储的波场坐标赋值
void SetCoord(AFDPU2D *Pa,
            uint *h_Coord, const Partition &pt);

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
            float *sgs_t, const Partition &pt);

// 求取梯度
void CalGrad(AFDPU2D Pa,
            IP *ip,
            CPUVs *plan,
            float *sgs_t,
            float *sgs_c,
            float *sgs_r,
            uint It,
             const Partition &pt);

// 求取步长
void CalStepLength(AFDPU2D Pa,
                IP *ip,
                CPUVs *plan,
                float *sgs_t,
                float *sgs_c,
                float e, const Partition &pt);

// 下次迭代前的预处理
void PreProcess(AFDPU2D Pa,
                IP *ip,
                CPUVs *plan, const Partition &pt);


#endif // TESTTDFWI_H
