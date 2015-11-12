/******************************************
 * author:dwx
 * 2015.11.7
 ******************************************/
#include "Partition.h"
#include <iostream>
using namespace std;
Partition::Partition()
{
}

Partition::Partition(const AFDPU2D *Pa, const IP *ip, uint totallength_x, uint totallength_z, uint sumblock_x, uint sumblock_z, H_Border h_U, H_Border h_VW, uint border_h_Coord, uint rank, uint size)
{
    this->totallength_x = totallength_x;
    this->totallength_z = totallength_z;
    this->sumBlock_x = sumblock_x;
    this->sumBlock_z = sumblock_z;
    this->rank = rank;
    this->size = size;
    
    this->blockPosition_x = rank % this->getsumBlock_x();
    this->blockPosition_z = rank / this->getsumBlock_x();
    this->blockLength_x = totallength_x / this->getsumBlock_x();
    this->blockLength_z = totallength_z / this->getsumBlock_z();

    uint remainder_x = totallength_x % this->getsumBlock_x();
    uint remainder_z = totallength_z % this->getsumBlock_z();

    this->indexmin_x = blockPosition_x * blockLength_x + (blockPosition_x <= remainder_x ? blockPosition_x : remainder_x);
    this->indexmax_x = this->indexmin_x + blockLength_x - 1 + (blockPosition_x < remainder_x ? 1 : 0);
    this->indexmin_z = blockPosition_z * blockLength_z + (blockPosition_z <= remainder_z ? blockPosition_z : remainder_z);
    this->indexmax_z = this->indexmin_z + this->blockLength_z - 1 + (blockPosition_z < remainder_z ? 1 : 0);
    
    this->interiormin_x = this->indexmin_x > Pa->PMLx ? this->indexmin_x : Pa->PMLx;
    this->interiormax_x = this->indexmax_x < Pa->PMLx + Pa->Nx - 1 ? this->indexmax_x : Pa->PMLx + Pa->Nx - 1;
    this->interiormin_z = this->indexmin_z > Pa->PMLz ? this->indexmin_z : Pa->PMLz;
    this->interiormax_z = this->indexmax_z < Pa->PMLz + Pa->Nz - 1 ? this->indexmax_z : Pa->PMLz + Pa->Nz - 1;
    this->interiorLength_x = this->interiormax_x >= this->interiormin_x ? this->interiormax_x - this->interiormin_x : 0;
    this->interiorLength_z = this->interiormax_z >= this->interiormin_z ? this->interiormax_z - this->interiormin_z : 0;

    for (uint is = 0; is < ip->ShotN; is++)
    {
        int temp_x = Pa->PMLx + is * 10 + 200;
        int temp_z = Pa->PMLz + 2;
        if(temp_x >= indexmin_x && temp_x <= indexmax_x && temp_z >= indexmin_z && temp_z <= indexmax_z)
        {
            this->shot.push_back({temp_x - indexmin_x, temp_z - indexmin_z});
        }
    }

    this->border_h_Coord = border_h_Coord;

    this->h_U = h_U;
    this->h_U.length_x = this->blockLength_x + h_U.leftborder + h_U.rightborder;
    this->h_U.length_z = this->blockLength_z + h_U.topborder + h_U.bottomborder;

    this->h_VW = h_VW;
    this->h_VW.length_x = this->blockLength_x + h_VW.leftborder + h_VW.rightborder;
    this->h_VW.length_z = this->blockLength_z + h_VW.topborder + h_VW.bottomborder;

    seth_Vp_trans(*Pa);
    seth_Vp_border(*Pa);
    setInside(*Pa);
    set_h_Coord(*Pa);
    setRL(ip, Pa);
}
Partition::~Partition()
{

}

uint Partition::getrank() const
{
    return this->rank;
}
uint Partition::getsize() const
{
    return this->size;
}

uint Partition::getblockPosition_x() const
{
    return this->blockPosition_x;
}
uint Partition::getblockPosition_z() const
{
    return this->blockPosition_z;
}
uint Partition::getblockLength_x() const
{
    return this->blockLength_x;
}
uint Partition::getblockLength_z() const
{
    return this->blockLength_z;
}
uint Partition::gettotallength_x() const
{
    return this->totallength_x;
}
uint Partition::gettotallength_z() const
{
    return this->totallength_z;
}
uint Partition::getsumBlock_x() const
{
    return this->sumBlock_x;
}
uint Partition::getsumBlock_z() const
{
    return this->sumBlock_z;
}
uint Partition::getindexmin_x() const
{
    return this->indexmin_x;
}
uint Partition::getindexmax_x() const
{
    return this->indexmax_x;
}
uint Partition::getindexmin_z() const
{
    return this->indexmin_z;
}
uint Partition::getindexmax_z() const
{
    return this->indexmax_z;
}
uint Partition::getinteriormin_x() const
{
    return this->interiormin_x;
}
uint Partition::getinteriormax_x() const
{
    return this->interiormax_x;
}
uint Partition::getinteriormin_z() const
{
    return this->interiormin_z;
}
uint Partition::getinteriormax_z() const
{
    return this->interiormax_z;
}
uint Partition::getinteriorLength_x() const
{
    return this->interiorLength_x;
}
uint Partition::getinteriorLength_z() const
{
    return this->interiorLength_z;
}
bool Partition::isfirstblock_x() const
{
    return this->rank % this->getsumBlock_x() == 0;
}
bool Partition::isfirstblock_z() const
{
    return this->rank < this->getsumBlock_x();
}
bool Partition::islastblock_x() const
{
    return this->rank % this->getsumBlock_x() == this->getsumBlock_x() - 1;
}
bool Partition::islastblock_z() const
{
    return this->rank / this->getsumBlock_x() >= this->getsumBlock_z() - 1;
}
uint Partition::getinsidenum() const
{
    return this->inside_num;
}
Inside* Partition::getInside() const
{
    return this->inside;
}
uint Partition::getinside_length() const
{
    return this->inside_length;
}

uint Partition::get_h_Coord_num() const
{
    return this->h_Coord_num;
}
H_Coord* Partition::get_h_Coord() const
{
    return this->h_coord;
}
uint Partition::geth_Coord_length() const
{
    return this->h_Coord_length;
}
H_Border Partition::geth_U() const
{
    return this->h_U;
}
H_Border Partition::geth_VW() const
{
    return this->h_VW;
}
H_Border Partition::geth_Vp() const
{
    return this->h_Vp;
}

uint Partition::getShot_num() const
{
    return this->shot.size();
}
bool Partition::getiscoverbyNPML() const
{
    return this->iscoverbyNPML;
}


vector<pair<uint, uint>> Partition::getShot() const
{
    return this->shot;
}
void Partition::setRL(const IP *ip, const AFDPU2D *Pa)
{
    this->RL_beginnum = INT_MAX;
    this->RL_endnum = INT_MIN;
    for(uint is = 0; is < ip->ShotN; ++is)
    {
        for (uint m = 0; m < ip->St[is].rn; m++)
        {
            int temp_x = Pa->PMLx + 1 * m;
            int temp_z = Pa->PMLz + 2;
            if(temp_x >= indexmin_x && temp_x <= indexmax_x && temp_z >= indexmin_z && temp_z <= indexmax_z)
            {
                this->RL_beginnum = this->RL_beginnum < m ? this->RL_beginnum : m;
                this->RL_endnum = this->RL_endnum > m ? this->RL_endnum : m;
                this->rl.push_back({temp_x - indexmin_x, temp_z - indexmin_z});
            }
        }
    }
}

uint Partition::getRL_num() const
{
    return this->rl.size();
}
vector<pair<uint, uint>> Partition::getRL() const
{
    return this->rl;
}
uint Partition::getRL_beginnum() const
{
    return this->RL_beginnum;
}
uint Partition::getRL_endnum() const
{
    return this->RL_endnum;
}
void Partition::seth_Vp_trans(AFDPU2D Pa)
{
    if(this->indexmin_z <= Pa.PMLz && this->indexmax_z >= Pa.PMLz)//to top
    {
        this->trans_h_Vp.push_back(1);
    }
    else
    {
        this->trans_h_Vp.push_back(0);
    }

    if(this->indexmin_x <= Pa.PMLx && this->indexmax_x >= Pa.PMLx)//to left
    {
        this->trans_h_Vp.push_back(1);
    }
    else
    {
        this->trans_h_Vp.push_back(0);
    }

    if(this->indexmin_z < Pa.PMLz + Pa.Nz && this->indexmax_z >= Pa.PMLz + Pa.Nz - 1)//to bottom
    {
        this->trans_h_Vp.push_back(1);
    }
    else
    {
        this->trans_h_Vp.push_back(0);
    }

    if(this->indexmin_x < Pa.PMLx + Pa.Nx && this->indexmax_x >= Pa.PMLx + Pa.Nx - 1)//to right
    {
        this->trans_h_Vp.push_back(1);
    }
    else
    {
        this->trans_h_Vp.push_back(0);
    }
}
vector<uint> Partition::gettrans_h_Vp() const
{
    return this->trans_h_Vp;
}
void Partition::seth_Vp_border(AFDPU2D Pa)
{
    if(this->indexmax_x < Pa.PMLx)
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 1;
            this->h_Vp.rightborder = 1;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp.topborder = 1;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 1;
        }
        else
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 1;
        }
    }
    else if(this->indexmin_x >= Pa.PMLx + Pa.Nx)
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 1;
            this->h_Vp.bottomborder = 1;
            this->h_Vp.rightborder = 0;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp.topborder = 1;
            this->h_Vp.leftborder = 1;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 0;
        }
        else
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 1;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 0;
        }
    }
    else
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 1;
            this->h_Vp.rightborder = 0;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp.topborder = 1;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 0;
        }
        else
        {
            this->h_Vp.topborder = 0;
            this->h_Vp.leftborder = 0;
            this->h_Vp.bottomborder = 0;
            this->h_Vp.rightborder = 0;
        }
    }
    this->h_Vp.length_x = this->h_Vp.leftborder + this->blockLength_x + this->h_Vp.rightborder;
    this->h_Vp.length_z = this->h_Vp.topborder + this->blockLength_z + this->h_Vp.bottomborder;
}
void Partition::setInside(AFDPU2D Pa)
{
    if(indexmin_x > Pa.Nx + Pa.PMLx || indexmax_x < Pa.PMLx || indexmin_z > Pa.Nz + Pa.PMLz || indexmax_z < Pa.PMLz)
    {
        this->iscoverbyNPML = true;
        this->inside_num = 1;
        Inside *ins = new Inside[this->inside_num];
        this->inside = ins;
        ins->indexmin_x = this->indexmin_x;
        ins->indexmax_x = this->indexmax_x;
        ins->indexmin_z = this->indexmin_z;
        ins->indexmax_z = this->indexmax_z;
        ins->length_x = ins->indexmax_x - ins->indexmin_x;
        ins->length_z = ins->indexmax_z - ins->indexmin_z;
    }
    else
    {
        this->iscoverbyNPML = false;
        if(indexmin_x < Pa.PMLx && indexmin_z < Pa.PMLz)
        {
            if(indexmax_x < Pa.PMLx + Pa.Nx && indexmax_z < Pa.PMLz + Pa.Nz)
            {
                this->inside_num = 2;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_x = Pa.PMLx - 1;
                ins[0].indexmax_z = this->indexmax_z;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;

                ins[1].indexmin_x = Pa.PMLx;
                ins[1].indexmin_z = this->indexmin_z;
                ins[1].indexmax_x = this->indexmax_x;
                ins[1].indexmax_z = Pa.PMLz - 1;
                ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;


                this->inside = ins;
            }
            else
            {
                ///tai da meiyiyi
            }
        }
        else if(indexmin_x < Pa.PMLx && indexmax_z >= Pa.PMLz + Pa.Nz)
        {
            if(indexmax_x < Pa.PMLx + Pa.Nx && indexmin_z >= Pa.PMLz)
            {
                this->inside_num = 2;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = Pa.PMLx - 1;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_z = this->indexmax_z;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                ins[1].indexmin_x = Pa.PMLx;
                ins[1].indexmax_x = this->indexmax_x;
                ins[1].indexmin_z = Pa.Nz + Pa.PMLz;
                ins[1].indexmax_z = this->indexmax_z;
                ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                this->inside = ins;
            }
            else
            {

            }
        }
        else if(indexmax_x >= Pa.PMLx + Pa.Nx && indexmin_z < Pa.PMLz)
        {
            if(indexmin_x >= Pa.PMLx && indexmax_z < Pa.Nz + Pa.PMLz)
            {
                this->inside_num = 2;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = Pa.Nx + Pa.PMLx - 1;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_z = Pa.Nz - 1;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                ins[1].indexmin_x = Pa.Nx + Pa.PMLx;
                ins[1].indexmin_z = this->indexmin_z;
                ins[1].indexmax_x = this->indexmax_x;
                ins[1].indexmax_z = this->indexmax_z;
                ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                this->inside = ins;
            }
            else
            {

            }
        }
        else if(indexmax_x >= Pa.PMLx + Pa.Nx && indexmax_z >= Pa.PMLz + Pa.Nz)
        {
            if(indexmin_x < Pa.Nx && indexmin_z < Pa.Nz)
            {
                this->inside_num = 2;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = Pa.Nx + Pa.PMLx - 1;
                ins[0].indexmin_z = Pa.Nz + Pa.PMLz;
                ins[0].indexmax_z = this->indexmax_z;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                ins[1].indexmin_x = Pa.Nx + Pa.PMLx;
                ins[1].indexmax_x = Pa.Nx + Pa.PMLx;
                ins[1].indexmin_z = this->indexmin_z;
                ins[1].indexmax_z = this->indexmax_z;
                ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                this->inside = ins;
            }
            else
            {

            }
        }
        else if(indexmin_x < Pa.PMLx && indexmax_x >= Pa.PMLx)
        {
            if(indexmax_x < Pa.PMLz + Pa.Nz)
            {
                if(indexmin_z >= Pa.PMLz && indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->inside_num = 1;
                    Inside *ins = new Inside[this->inside_num];
                    ins[0].indexmin_x = this->indexmin_x;
                    ins[0].indexmax_x = Pa.PMLx - 1;
                    ins[0].indexmin_z = this->indexmin_z;
                    ins[0].indexmax_z = this->indexmax_z;
                    ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                    ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                    this->inside = ins;
                }
                else if(indexmin_z < Pa.PMLz && indexmax_z >= Pa.PMLz + Pa.Nz)
                {
                    this->inside_num = 3;
                    Inside *ins = new Inside[this->inside_num];
                    ins[0].indexmin_x = this->indexmin_x;
                    ins[0].indexmax_x = this->indexmax_x;
                    ins[0].indexmin_z = this->indexmin_z;
                    ins[0].indexmax_z = Pa.PMLz - 1;
                    ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                    ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                    ins[1].indexmin_x = this->indexmin_x;
                    ins[1].indexmax_x = Pa.PMLx - 1;
                    ins[1].indexmin_z = Pa.PMLz;
                    ins[1].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                    ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                    ins[2].indexmin_x = this->indexmin_x;
                    ins[2].indexmax_x = this->indexmax_x;
                    ins[2].indexmin_z = Pa.PMLz + Pa.Nz;
                    ins[2].indexmax_z = this->indexmax_z;
                    ins[2].length_x = ins[2].getindexmax_x() - ins[2].getindexmin_x() + 1;
                    ins[2].length_z = ins[2].getindexmax_z() - ins[2].getindexmin_z() + 1;

                    this->inside = ins;
                }
            }
            else
            {
                if(indexmin_z >= Pa.PMLz && indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->inside_num = 2;
                    Inside *ins = new Inside[this->inside_num];
                    ins[0].indexmin_x = this->indexmin_x;
                    ins[0].indexmax_x = Pa.PMLx - 1;
                    ins[0].indexmin_z = this->indexmin_z;
                    ins[0].indexmax_z = this->indexmax_z;
                    ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                    ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                    ins[1].indexmin_x = Pa.PMLx + Pa.Nx;
                    ins[1].indexmax_x = this->indexmax_x;
                    ins[1].indexmin_z = this->indexmin_z;
                    ins[1].indexmax_z = this->indexmax_z;
                    ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                    ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                    this->inside = ins;
                }
                else if(indexmin_z < Pa.PMLz && indexmax_z >= Pa.PMLz + Pa.Nz)
                {
                    this->inside_num = 4;
                    Inside *ins = new Inside[this->inside_num];
                    ins[0].indexmin_x = this->indexmin_x;
                    ins[0].indexmax_x = Pa.PMLx - 1;
                    ins[0].indexmin_z = this->indexmin_z;
                    ins[0].indexmax_z = this->indexmax_z;
                    ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                    ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                    ins[1].indexmin_x = Pa.PMLx;
                    ins[1].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    ins[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    ins[1].indexmax_z = this->indexmax_z;
                    ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                    ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                    ins[2].indexmin_x = Pa.PMLx + Pa.Nx;
                    ins[2].indexmax_x = this->indexmax_x;
                    ins[2].indexmin_z = this->indexmin_z;
                    ins[2].indexmax_z = this->indexmax_z;
                    ins[2].length_x = ins[2].getindexmax_x() - ins[2].getindexmin_x() + 1;
                    ins[2].length_z = ins[2].getindexmax_z() - ins[2].getindexmin_z() + 1;

                    ins[3].indexmin_x = Pa.PMLx;
                    ins[3].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    ins[3].indexmin_z = this->indexmin_z;
                    ins[3].indexmax_z = Pa.PMLz - 1;
                    ins[3].length_x = ins[3].getindexmax_x() - ins[3].getindexmin_x() + 1;
                    ins[3].length_z = ins[3].getindexmax_z() - ins[3].getindexmin_z() + 1;

                    this->inside = ins;
                }
            }
        }
        else if(indexmin_x >= Pa.PMLx && indexmin_x < Pa.PMLx + Pa.Nx && indexmax_x >= Pa.PMLx + Pa.Nx)
        {
            if(indexmin_z >= Pa.PMLz && indexmax_z < Pa.PMLz + Pa.Nz)
            {
                this->inside_num = 1;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = Pa.PMLx + Pa.Nx;
                ins[0].indexmax_x = this->indexmax_x;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_z = this->indexmax_z;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                this->inside = ins;
            }
            else if(indexmin_z < Pa.PMLz && indexmax_z >= Pa.PMLz + Pa.Nz)
            {
                this->inside_num = 3;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = this->indexmax_x;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_z = Pa.PMLz - 1;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                ins[1].indexmin_x = Pa.PMLx + Pa.Nx;
                ins[1].indexmax_x = this->indexmax_x;
                ins[1].indexmin_z = Pa.PMLz;
                ins[1].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                ins[2].indexmin_x = this->indexmin_x;
                ins[2].indexmax_x = this->indexmax_x;
                ins[2].indexmin_z = Pa.PMLz + Pa.Nz;
                ins[2].indexmax_z = this->indexmax_z;
                ins[2].length_x = ins[2].getindexmax_x() - ins[2].getindexmin_x() + 1;
                ins[2].length_z = ins[2].getindexmax_z() - ins[2].getindexmin_z() + 1;

                this->inside = ins;
            }
        }
        else if(indexmin_z < Pa.PMLz && indexmax_z >= Pa.PMLz)
        {
            if(indexmax_z >= Pa.PMLz + Pa.Nz)
            {
                if(indexmin_x >= Pa.PMLx && indexmax_x < Pa.PMLx + Pa.Nx)
                {
                    this->inside_num = 2;
                    Inside *ins = new Inside[this->inside_num];
                    ins[0].indexmin_x = this->indexmin_x;
                    ins[0].indexmax_x = this->indexmax_x;
                    ins[0].indexmin_z = this->indexmin_z;
                    ins[0].indexmax_z = Pa.PMLz - 1;
                    ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                    ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                    ins[1].indexmin_x = this->indexmin_x;
                    ins[1].indexmax_x = this->indexmax_x;
                    ins[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    ins[1].indexmax_z = this->indexmax_z;
                    ins[1].length_x = ins[1].getindexmax_x() - ins[1].getindexmin_x() + 1;
                    ins[1].length_z = ins[1].getindexmax_z() - ins[1].getindexmin_z() + 1;

                    this->inside = ins;
                }
            }
            else if(indexmax_z < Pa.PMLz + Pa.Nz)
            {
                this->inside_num = 1;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = this->indexmax_x;
                ins[0].indexmin_z = this->indexmin_z;
                ins[0].indexmax_z = Pa.PMLz - 1;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                this->inside = ins;
            }
        }
        else if(indexmin_z >= Pa.PMLz && indexmax_z >= Pa.PMLz + Pa.Nz)
        {
            if(indexmin_x >= Pa.PMLx && indexmax_x < Pa.PMLx + Pa.Nx)
            {
                this->inside_num = 1;
                Inside *ins = new Inside[this->inside_num];
                ins[0].indexmin_x = this->indexmin_x;
                ins[0].indexmax_x = this->indexmax_x;
                ins[0].indexmin_z = Pa.PMLz + Pa.Nz;
                ins[0].indexmax_z = this->indexmax_z;
                ins[0].length_x = ins[0].getindexmax_x() - ins[0].getindexmin_x() + 1;
                ins[0].length_z = ins[0].getindexmax_z() - ins[0].getindexmin_z() + 1;


                this->inside = ins;
            }
        }
        else
        {
            this->inside_num = 0;
        }
    }

    this->inside_length = 0;
    for(uint i = 0; i < this->inside_num; ++i)
    {
        this->inside_length += this->inside[i].getlength_x() * this->inside[i].getlength_z();
    }
}

void Partition::set_h_Coord(AFDPU2D Pa)
{
    if(this->indexmax_x < Pa.PMLx && (this->indexmax_z < Pa.PMLz || this->indexmin_z >= Pa.PMLz + Pa.Nz))
    {
        this->h_Coord_num = 0;
        this->h_coord = NULL;
    }
    else if(this->indexmin_x >= Pa.PMLx + Pa.Nx && (this->indexmax_z < Pa.PMLz || this->indexmin_z >= Pa.PMLz + Pa.Nz))
    {
        this->h_Coord_num = 0;
        this->h_coord = NULL;
    }
    else if(this->indexmax_z < Pa.PMLz - this->border_h_Coord && this->indexmax_x < Pa.PMLx - this->border_h_Coord && this->indexmin_z >= Pa.PMLz + Pa.Nz + this->border_h_Coord && this->indexmin_x >= Pa.PMLx + Pa.Nx + this->border_h_Coord)
    {
        this->h_Coord_num = 0;
        this->h_coord = NULL;
    }
    else if(this->indexmin_x >= Pa.PMLx && this->indexmax_x < Pa.PMLx + Pa.Nx && this->indexmin_z >= Pa.PMLz && this->indexmax_z < Pa.PMLz + Pa.Nz)
    {
        this->h_Coord_num = 0;
        this->h_coord = NULL;

    }
    else if(this->indexmin_x < Pa.PMLx && this->indexmax_x >= Pa.PMLx)
    {
        if(this->indexmax_x < Pa.PMLx + Pa.Nx)
        {
            if(this->indexmax_z < Pa.PMLz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = Pa.PMLx;
                coord[0].indexmax_x = this->indexmax_x;
                coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                coord[0].indexmax_z = this->indexmax_z;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else if(this->indexmin_z < Pa.PMLz && this->indexmax_z >= Pa.PMLz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 2;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = Pa.PMLz;
                    coord[0].indexmax_z = this->indexmax_z;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = this->indexmax_x;
                    coord[1].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[1].indexmax_z = Pa.PMLz - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 3;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = Pa.PMLz;
                    coord[0].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = this->indexmax_x;
                    coord[1].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[1].indexmax_z = Pa.PMLz - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    coord[2].indexmin_x = Pa.PMLx;
                    coord[2].indexmax_x = this->indexmax_x;
                    coord[2].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[2].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.Nz + Pa.PMLz + this->border_h_Coord - 1;
                    coord[2].length_x = coord[2].indexmax_x - coord[2].indexmin_x + 1;
                    coord[2].length_z = coord[2].indexmax_z - coord[2].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz && this->indexmin_z < Pa.PMLz + Pa.Nz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = this->indexmax_z;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 2;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = this->indexmax_x;
                    coord[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[1].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.Nz + Pa.PMLz + this->border_h_Coord - 1;;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = Pa.PMLx;
                coord[0].indexmax_x = this->indexmax_x;
                coord[0].indexmin_z = this->indexmin_z;
                coord[0].indexmax_z = Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 < this->indexmax_z ? Pa.Nz + Pa.PMLz + this->border_h_Coord - 1 : this->indexmax_z;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else
            {

            }
        }
        else if(this->indexmax_x >= Pa.PMLx + Pa.Nx)
        {
            if(this->indexmax_z < Pa.PMLz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = Pa.PMLx;
                coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                coord[0].indexmax_z = this->indexmax_z;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else if(this->indexmin_z < Pa.PMLz && this->indexmax_z >= Pa.PMLz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 3;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = Pa.PMLz;
                    coord[0].indexmax_z = this->indexmax_z;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[1].indexmin_z = Pa.PMLz - this->border_h_Coord - 1 > this->indexmin_z ? Pa.PMLz - this->border_h_Coord - 1 : this->indexmax_z;
                    coord[1].indexmax_z = Pa.PMLz;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    coord[2].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[2].indexmax_x = Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 < this->indexmax_x ? Pa.Nx + Pa.PMLx + this->border_h_Coord - 1 : this->indexmax_x;///////
                    coord[2].indexmin_z = Pa.PMLz;
                    coord[2].indexmax_z = this->indexmax_z;
                    coord[2].length_x = coord[2].indexmax_x - coord[2].indexmin_x + 1;
                    coord[2].length_z = coord[2].indexmax_z - coord[2].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 4;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = Pa.PMLz;
                    coord[0].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[1].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[1].indexmax_z = Pa.PMLz - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    coord[2].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[2].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;
                    coord[2].indexmin_z = Pa.PMLz;
                    coord[2].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[2].length_x = coord[2].indexmax_x - coord[2].indexmin_x + 1;
                    coord[2].length_z = coord[2].indexmax_z - coord[2].indexmin_z + 1;

                    coord[3].indexmin_x = Pa.PMLx;
                    coord[3].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[3].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[3].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[3].length_x = coord[3].indexmax_x - coord[3].indexmin_x + 1;
                    coord[3].length_z = coord[3].indexmax_z - coord[3].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz && this->indexmin_z < Pa.PMLz + Pa.Nz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 2;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = this->indexmax_z;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[1].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;;
                    coord[1].indexmin_z = this->indexmin_z;
                    coord[1].indexmax_z = this->indexmax_z;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 3;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x > Pa.PMLx - this->border_h_Coord ? this->indexmin_x : Pa.PMLx - this->border_h_Coord;
                    coord[0].indexmax_x = Pa.PMLx - 1;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = Pa.Nz + Pa.PMLz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx;
                    coord[1].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[1].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    coord[2].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[2].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;
                    coord[2].indexmin_z = this->indexmin_z;
                    coord[2].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[2].length_x = coord[2].indexmax_x - coord[2].indexmin_x + 1;
                    coord[2].length_z = coord[2].indexmax_z - coord[2].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
            {
                this->h_Coord_num = 1;
                H_Coord * coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = Pa.PMLx;
                coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                coord[0].indexmin_z = this->indexmin_z;
                coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else
            {
//cout << rank << endl;
            }
        }
    }
    else if(this->indexmin_x >= Pa.PMLx && this->indexmin_x < Pa.PMLx + Pa.Nx)
    {
        if(this->indexmax_x < Pa.Nx + Pa.PMLx)
        {
            if(this->indexmax_z < Pa.PMLz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = this->indexmin_x;
                coord[0].indexmax_x = this->indexmax_x;
                coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                coord[0].indexmax_z = this->indexmax_z;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else if(this->indexmin_z < Pa.PMLz && this->indexmax_z >= Pa.PMLz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = this->indexmax_x;
                    coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[0].indexmax_z = Pa.PMLz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 2;
                    H_Coord* coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = this->indexmax_x;
                    coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[0].indexmax_z = Pa.PMLz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = this->indexmin_x;
                    coord[1].indexmax_x = this->indexmax_x;
                    coord[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[1].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz)
            {
                if(this->indexmin_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = this->indexmax_x;
                    coord[0].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = this->indexmax_x;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else
            {

            }
        }
        else if(this->indexmax_x >= Pa.PMLx + Pa.Nx)
        {
            if(this->indexmax_z < Pa.PMLz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = this->indexmin_x;
                coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                coord[0].indexmax_z = this->indexmax_z;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else if(this->indexmin_z < Pa.PMLz)
            {
                if(this->indexmax_z >= Pa.PMLz && this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 2;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[0].indexmax_z = Pa.PMLz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[1].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;
                    coord[1].indexmin_z = Pa.PMLz;
                    coord[1].indexmax_z = this->indexmax_z;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 3;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[0].indexmin_z = this->indexmin_z > Pa.PMLz - this->border_h_Coord ? this->indexmin_z : Pa.PMLz - this->border_h_Coord;
                    coord[0].indexmax_z = Pa.PMLz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[1].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;
                    coord[1].indexmin_z = Pa.PMLz;
                    coord[1].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    coord[2].indexmin_x = this->indexmin_x;
                    coord[2].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[2].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[2].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[2].length_x = coord[2].indexmax_x - coord[2].indexmin_x + 1;
                    coord[2].length_z = coord[2].indexmax_z - coord[2].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz && this->indexmin_z < Pa.PMLz + Pa.Nz)
            {
                if(this->indexmax_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[0].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = this->indexmax_z;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {
                    this->h_Coord_num = 2;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = Pa.PMLx + Pa.Nx;
                    coord[0].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;;
                    coord[0].indexmin_z = this->indexmin_z;
                    coord[0].indexmax_z = Pa.PMLz + Pa.Nz - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    coord[1].indexmin_x = this->indexmin_x;
                    coord[1].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                    coord[1].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[1].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[1].length_x = coord[1].indexmax_x - coord[1].indexmin_x + 1;
                    coord[1].length_z = coord[1].indexmax_z - coord[1].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = this->indexmin_x;
                coord[0].indexmax_x = Pa.PMLx + Pa.Nx - 1;
                coord[0].indexmin_z = this->indexmin_z;
                coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                this->h_coord = coord;
            }
            else
            {

            }
        }
        else
        {

        }
    }
    else if(this->indexmin_x >= Pa.PMLx + Pa.Nx)
    {
        this->h_Coord_num = 1;
        H_Coord *coord = new H_Coord[this->h_Coord_num];
        coord[0].indexmin_x = this->indexmin_x;
        coord[0].indexmax_x = this->indexmax_x < Pa.PMLx + Pa.Nx + this->border_h_Coord - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx + this->border_h_Coord - 1;
        coord[0].indexmin_z = this->indexmin_z > Pa.PMLx ? this->indexmin_z : Pa.PMLx;
        coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz - 1;
        coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
        coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

        this->h_coord = coord;
    }
    else
    {

    }

    this->h_Coord_length = 0;
    for(uint i = 0; i < this->h_Coord_num; ++i)
    {
        this->h_Coord_length += this->h_coord[i].getlength_x() * this->h_coord[i].getlength_z();
    }
}

