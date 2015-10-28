#ifndef PARTITION_H
#define PARTITION_H
//#include "testTDFWI.h"
#include <vector>
using namespace std;
#define		uint			unsigned int
#define		usht			unsigned short
#define		PI				3.141592653589793f
#define		PowerGrad		2.0f

// 空间8阶交错网格有限差分系数
#define		C1_4			1.196289062506247f
#define		C2_4			-0.079752604166839f
#define		C3_4			0.009570312500021f
#define		C4_4			-0.000697544642859f
// h_Vp_border_flag
#define     NONE_BORDER     0
#define     TOP_BORDER      1
#define     LEFT_BORDER     2
#define     BOTTOM_BORDER   3
#define     RIGHT_BORDER    4
//
class Inside;
class h_Coord;
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
//class APDPU2D;
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

class Partition{
private:
    int rank, size;
    int blockPosition_x;
    int blockPosition_z;

    int totallength_x;
    int totallength_z;

    int sumBlock_x;
    int sumBlock_z;

    int indexmin_x;
    int indexmax_x;
    int indexmin_z;
    int indexmax_z;

    int interiormin_x;
    int interiormax_x;
    int interiormin_z;
    int interiormax_z;

    int interiorLength_x;
    int interiorLength_z;

    h_Coord *h_coord;
    int h_Coord_num;
    int border_h_Coord;

    Inside *inside;
    int inside_num;

    int blockLength_x;
    int blockLength_z;

    int borderLength_x;
    int borderLength_z;

    bool iscoverbyNPML;

    int h_V_leftborder;
    int h_V_rightborder;

    int h_W_topborder;
    int h_W_bottomborder;

    int h_U_leftborder;
    int h_U_rightborder;
    int h_U_topborder;
    int h_U_bottomborder;

//    int h_Vp_topborder;
//    int h_Vp_leftborder;
//    int h_Vp_bottomborder;
//    int h_Vp_rightborder;

    int h_Vp_border[4];//top left bottom right

    vector<pair<int, int>> shot;

    vector<pair<int, int>> rl;

public:
    Partition();
    Partition(AFDPU2D &Pa, int totallength_x, int totallength_z, int borderLength_x, int borderLength_z, int sumblock_x, int sumblock_z, int rank, int size);
    ~Partition();

    int getrank() const;
    int getsize() const;

    int getborderlength_x() const;
    int getborderlength_z() const;
    int gettotallength_x() const;
    int gettotallength_z() const;
    int getsumBlock_x() const;
    int getsumBlock_z() const;
    int getblockLength_x() const;
    int getblockLength_z() const;
    int getindexmin_x() const;
    int getindexmin_z() const;
    int getindexmax_x() const;
    int getindexmax_z() const;

    int getinteriormin_x() const;
    int getinteriormax_x() const;
    int getinteriormin_z() const;
    int getinteriormax_z() const;

    int getinteriorLength_x() const;
    int getinteriorLength_z() const;

    bool isfirstblock_x() const;
    bool islastblock_x() const;
    bool isfirstblock_z() const;
    bool islastblock_z() const;

    int getinsidenum() const;
    Inside* getInside() const;
    void setInside(AFDPU2D &Pa);

    int get_h_Coord_num() const;
    h_Coord* get_h_Coord() const;
    void set_h_Coord(AFDPU2D &Pa);

    int geth_U_leftborder() const;
    int geth_U_rightborder() const;
    int geth_U_topborder() const;
    int geth_U_bottomborder() const;

    int geth_V_leftborder() const;
    int geth_V_rightborder() const;

    int geth_W_topborder() const;
    int geth_W_bottomborder() const;

    void seth_Vp_border(AFDPU2D &Pa);
    int *geth_Vp_border() const;

    int getShot_num() const;
    vector<pair<int, int>> getShot() const;

    int getRL_num() const;
    vector<pair<int, int>> getRL() const;
};

class Inside{
private:
    int indexmin_x;
    int indexmin_z;
    int indexmax_x;
    int indexmax_z;

    int length_x;
    int length_z;
public:
    int getindexmin_x() const;
    int getindexmin_z() const;
    int getindexmax_x() const;
    int getindexmax_z() const;

    int getlength_x() const;
    int getlength_z() const;

    friend void Partition::setInside(AFDPU2D &Pa);

};

class H_Coord
{
private:
    int indexmin_x;
    int indexmax_x;
    int indexmin_z;
    int indexmax_z;

    int length_x;
    int length_z;

public:
    int getindexmin_x() const;
    int getindexmax_x() const;
    int getindexmin_z() const;
    int getindexmax_z() const;

    int getlength_x() const;
    int getlength_z() const;

    friend void Partition::set_h_Coord(AFDPU2D &Pa);
};
#endif // PARTITION_H
