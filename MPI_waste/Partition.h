/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#ifndef PARTITION_H
#define PARTITION_H
//#include "testTDFWI.h"
#include <vector>
#include <limits.h>
#include <stddef.h>
using namespace std;




#define		usht			unsigned short
#define		PI				3.141592653589793f
#define		PowerGrad		2.0f

//#define		uint			unsigned int
typedef     unsigned int    uint;
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
class H_Coord;
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

class H_Border
{
public:
    uint length_x;
    uint length_z;

    uint topborder;
    uint leftborder;
    uint bottomborder;
    uint rightborder;

    H_Border();
    H_Border(uint topborder, uint leftborder, uint bottomborder, uint rightborder);
    H_Border(uint length_x, uint legnth_z, uint topborder, uint leftborder, uint bottomborder, uint rightborder);
};


class Partition{
private:
    uint rank, size;

    uint in_rank, in_size;

    //uint shot_num;

    uint blockPosition_x;
    uint blockPosition_z;

    uint in_blockPosition_x;
    uint in_blockPosition_z;

    uint totallength_x;
    uint totallength_z;

    uint sumBlock_x;
    uint sumBlock_z;

    uint in_sumBlock_x;
    uint in_sumBlock_z;

    uint indexmin_x;
    uint indexmax_x;
    uint indexmin_z;
    uint indexmax_z;

    uint blockLength_x;
    uint blockLength_z;

    uint interiormin_x;
    uint interiormax_x;
    uint interiormin_z;
    uint interiormax_z;

    uint interiorLength_x;
    uint interiorLength_z;

    H_Coord *h_coord;
    uint h_Coord_num;
    uint h_Coord_length;
    uint border_h_Coord;

    Inside *inside;
    uint inside_num;
    uint inside_length;

    bool iscoverbyNPML;

    H_Border h_U;
    H_Border h_VW;
    H_Border h_Vp;

    vector< pair<uint, uint>> shot;
    vector< pair<uint, uint>> rl;
    uint RL_beginnum;
    uint RL_endnum;

    vector<uint> trans_h_Vp;//top left bottom right

public:
    Partition();
    Partition(const AFDPU2D *Pa, const IP *ip, uint totallength_x, uint totallength_z, uint sumblock_x, uint sumblock_z, uint in_sumblock_x, uint in_sumblock_z, H_Border h_U, H_Border h_VW, uint border_h_Coord, uint rank, uint size, uint in_rank, uint in_size);
    ~Partition();


    uint getrank() const;
    uint getsize() const;

    uint get_in_rank() const;
    uint get_in_size() const;

    //uint get_shot_num() const;

    uint getblockPosition_x() const;
    uint getblockPosition_z() const;

    uint get_in_blockPosition_x() const;
    uint get_in_blockPosition_z() const;

    uint gettotallength_x() const;
    uint gettotallength_z() const;
    uint getsumBlock_x() const;
    uint getsumBlock_z() const;

    uint get_in_sumBlock_x() const;
    uint get_in_sumBlock_z() const;

    uint getblockLength_x() const;
    uint getblockLength_z() const;
    uint getindexmin_x() const;
    uint getindexmin_z() const;
    uint getindexmax_x() const;
    uint getindexmax_z() const;

    uint getinteriormin_x() const;
    uint getinteriormax_x() const;
    uint getinteriormin_z() const;
    uint getinteriormax_z() const;

    uint getinteriorLength_x() const;
    uint getinteriorLength_z() const;

    bool isfirstblock_x() const;
    bool islastblock_x() const;
    bool isfirstblock_z() const;
    bool islastblock_z() const;

    bool in_isfirstblock_x() const;
    bool in_islastblock_x() const;
    bool in_isfirstblock_z() const;
    bool in_islastblock_z() const;

    uint getinsidenum() const;
    Inside* getInside() const;
    void setInside(AFDPU2D Pa);
    uint getinside_length() const;

    uint get_h_Coord_num() const;
    H_Coord* get_h_Coord() const;
    void set_h_Coord(AFDPU2D Pa);
    uint geth_Coord_length() const;

    H_Border geth_U() const;
    H_Border geth_VW() const;
    H_Border geth_Vp() const;

    void seth_Vp_border(AFDPU2D Pa);

    uint getShot_num() const;
    vector<pair<uint, uint>> getShot() const;

    void setRL(const IP *ip, const AFDPU2D *Pa);
    uint getRL_num() const;
    vector<pair<uint, uint>> getRL() const;
    uint getRL_beginnum() const;
    uint getRL_endnum() const;

    void seth_Vp_trans(AFDPU2D Pa);
    vector<uint> gettrans_h_Vp() const;

    bool getiscoverbyNPML() const;
};

class Inside{
private:
    uint indexmin_x;
    uint indexmin_z;
    uint indexmax_x;
    uint indexmax_z;

    uint length_x;
    uint length_z;
public:
    uint getindexmin_x() const;
    uint getindexmin_z() const;
    uint getindexmax_x() const;
    uint getindexmax_z() const;

    uint getlength_x() const;
    uint getlength_z() const;

    friend void Partition::setInside(AFDPU2D Pa);

};

class H_Coord
{
private:
    uint indexmin_x;
    uint indexmax_x;
    uint indexmin_z;
    uint indexmax_z;

    uint length_x;
    uint length_z;

public:
    uint getindexmin_x() const;
    uint getindexmax_x() const;
    uint getindexmin_z() const;
    uint getindexmax_z() const;

    uint getlength_x() const;
    uint getlength_z() const;

    friend void Partition::set_h_Coord(AFDPU2D Pa);
};
#endif // PARTITION_H
