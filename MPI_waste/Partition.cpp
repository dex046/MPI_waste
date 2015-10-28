#include <Partition.h>
Partition::Partition()
{
    this->borderLength_x = 0;
    this->borderLength_z = 0;
    this->blockLength_x = 0;
    this->blockLength_z = 0;
    this->indexmin_x = 0;
    this->indexmax_x = 0;
    this->indexmin_z = 0;
    this->indexmax_z = 0;
    this->totallength_x = 0;
    this->totallength_z = 0;
    this->sumBlock_x = 0;
    this->sumBlock_z = 0;
    this->blockPosition_x = 0;
    this->blockPosition_z = 0;
}

Partition::Partition(AFDPU2D &Pa, IP &ip, int totallength_x, int totallength_z, int borderLength_x, int borderLength_z, int sumblock_x, int sumblock_z, int rank, int size)
{
    this->totallength_x = totallength_x;
    this->totallength_z = totallength_z;
    this->borderLength_x = borderLength_x;
    this->borderLength_z = borderLength_z;
    this->sumBlock_x = sumblock_x;
    this->sumBlock_z = sumblock_z;
    this->rank = rank;
    this->size = size;
    
    this->blockPosition_x = rank % this->getsumBlock_x();
    this->blockPosition_z = rank / this->getsumBlock_x();
    this->blockLength_x = totallength_x / this->getsumBlock_x();
    this->blockLength_z = totallength_z / this->getsumBlock_z();
    int remainder_x = totallength_x % this->getsumBlock_x();
    int remainder_z = totallength_z % this->getsumBlock_z();
    this->indexmin_x = blockPosition_x * blockLength_x + (blockPosition_x <= remainder_x ? blockPosition_x : remainder_x);
    this->indexmax_x = this->indexmin_x + blockLength_x + (blockPosition_x < remainder_x ? 1 : 0);
    this->indexmin_z = blockPosition_z * blockLength_z + (blockPosition_z <= remainder_z ? blockPosition_z : remainder_z);
    this->indexmax_z = this->indexmin_z + this->blockLength_z + (blockPosition_z < remainder_z ? 1 : 0);
    
    this->interiormin_x = this->indexmin_x > Pa.PMLx ? this->indexmin_x : Pa.PMLx;
    this->interiormax_x = this->indexmax_x < Pa.PMLx + Pa.Nx - 1 ? this->indexmax_x : Pa.PMLx + Pa.Nx - 1;
    this->interiormin_z = this->indexmin_z > Pa.PMLz ? this->indexmin_z : Pa.PMLz;
    this->interiormax_z = this->indexmax_z < Pa.PMLz + Pa.Nz - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz - 1;
    this->interiorLength_x = this->interiormax_x >= this->interiormin_x ? this->interiormax_x - this->interiormin_x : 0;
    this->interiorLength_z = this->interiormax_z >= this->interiormin_z ? this->interiormax_z - this->interiormin_z : 0;

    for (uint is = 0; is < ip.ShotN; is++)
    {
        int temp_x = Pa->PMLx + is * 10 + 200;
        int temp_z = Pa->PMLz + 2;
        if(temp_x >= indexmin_x && temp_x <= indexmax_x && temp_z >= indexmin_z && temp_z <= indexmax_z)
        {
            this->shot.push_back({temp_x - indexmin_x, temp_z - indexmin_z});
        }
        ip->St[is].rn = 510;
        for (uint m = 0; m < ip->St[is].rn; m++)
        {
            temp_x = Pa->PMLx + 1 * m;
            temp_z = Pa->PMLz + 2;
            if(temp_x >= indexmin_x && temp_x <= indexmax_x && temp_z >= indexmin_z && temp_z <= indexmax_z)
            {
                this->rl.push_back({temp_x - indexmin_x, temp_z - indexmin_z});
            }
        }
    }

    seth_Vp_border(Pa);
    setInside(Pa);
    set_h_Coord(Pa);

//    this->insideLength_x = insidemax_x - insidemin_x;
//    this->insideLength_z = insidemax_z - insidemin_z;

    
//    this->thisBlocklength_x = islastblock_x() ? 
}
int Partition::getrank() const
{
    return this->rank;
}
int Partition::getsize() const
{
    return this->size;
}
int Partition::getborderlength_x() const
{
    return this->borderLength_x;
}
int Partition::getborderlength_z() const
{
    return this->borderLength_z;
}
int Partition::getblockLength_x() const
{
    return this->blockLength_x;
}
int Partition::getblockLength_z() const
{
    return this->blockLength_z;
}
int Partition::gettotallength_x() const
{
    return this->totallength_x;
}
int Partition::gettotallength_z() const
{
    return this->totallength_z;
}
int Partition::getsumBlock_x() const
{
    return this->sumBlock_x;
}
int Partition::getsumBlock_z() const
{
    return this->sumBlock_z;
}
int Partition::getindexmin_x() const
{
    return this->indexmin_x;
}
int Partition::getindexmax_x() const
{
    return this->indexmax_x;
}
int Partition::getindexmin_z() const
{
    return this->indexmin_z;
}
int Partition::getindexmax_z() const
{
    return this->indexmax_z;
}
int Partition::getinteriormin_x() const
{
    return this->interiormin_x;
}
int Partition::getinteriormax_x() const
{
    return this->interiormax_x;
}
int Partition::getinteriormin_z() const
{
    return this->interiormin_z;
}
int Partition::getinteriormax_z() const
{
    return this->interiormax_z;
}
int Partition::getinteriorLength_x() const
{
    return this->interiorLength_x;
}
int Partition::getinteriorLength_z() const
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
int Partition::getinsidenum() const
{
    return this->inside_num;
}
Inside* Partition::getInside() const
{
    return this->inside;
}
int Partition::geth_U_leftborder() const
{
    return this->h_U_leftborder;
}
int Partition::geth_U_rightborder() const
{
    return this->h_U_rightborder;
}
int Partition::geth_U_topborder() const
{
    return this->h_U_topborder;
}
int Partition::geth_U_bottomborder() const
{
    return this->h_U_bottomborder;
}
int Partition::geth_V_leftborder() const
{
    return this->h_V_leftborder;
}
int Partition::geth_V_rightborder() const
{
    return this->h_V_rightborder;
}
int Partition::geth_W_topborder() const
{
    return this->h_W_topborder;
}
int Partition::geth_W_bottomborder() const
{
    return this->h_W_bottomborder;
}
int Partition::getShot_num() const
{
    return this->shot.size();
}
int Partition::getRL_num() const
{
    return this->rl.size();
}
vector<pair<int, int>> Partition::getShot() const
{
    return this->shot;
}
vector<pair<int, int>> Partition::getRL() const
{
    return this->rl;
}
void Partition::seth_Vp_border(AFDPU2D &Pa)
{
    if(this->indexmax_x < Pa.PMLx)
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 1;
            this->h_Vp_border[3] = 1;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp_border[0] = 1;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 1;
        }
        else
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 1;
        }
    }
    else if(this->indexmin_x >= Pa.PMLx + Pa.Nx)
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 1;
            this->h_Vp_border[2] = 1;
            this->h_Vp_border[3] = 0;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp_border[0] = 1;
            this->h_Vp_border[1] = 1;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 0;
        }
        else
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 1;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 0;
        }
    }
    else
    {
        if(this->indexmax_z < Pa.PMLz)
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 1;
            this->h_Vp_border[3] = 0;
        }
        else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
        {
            this->h_Vp_border[0] = 1;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 0;
        }
        else
        {
            this->h_Vp_border[0] = 0;
            this->h_Vp_border[1] = 0;
            this->h_Vp_border[2] = 0;
            this->h_Vp_border[3] = 0;
        }
    }
}

int* Partition::geth_Vp_border() const
{
    return this->h_Vp_border;
}
void Partition::setInside(AFDPU2D &Pa)
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
                ins[0].indexmax_x = this->insidemax_x;
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
                ins[0].indexmin_z = this->insidemin_z;
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
//        this->insidemin_x = (indexmin_x > Pa.PMLx) ? indexmin_x : Pa.PMLx;
//        this->insidemin_z = (indexmin_z > Pa.PMLz) ? indexmin_z : Pa.PMLz;
//        this->insidemax_x = (indexmax_x < Pa.Nx + Pa.PMLx) ? indexmax_x : Pa.Nx + Pa.PMLx;
//        this->insidemax_z = (indexmax_z < Pa.Nz + Pa.PMLz) ? indexmax_z : Pa.Nz + Pa.PMLz;
    }
}

int Partition::get_h_Coord_num() const
{
    return this->h_Coord_num;
}
H_Coord* Partition::get_h_Coord() const
{
    return this->h_coord;
}
void Partition::set_h_Coord(AFDPU2D &Pa)
{
    if(this->indexmax_x < Pa.PMLx && (this->indexmax_z < Pa.PMLz || this->indexmin_z >= Pa.PMLz + Pa.Nz))
    {

    }
    else if(this->indexmin_x >= Pa.PMLx + Pa.Nx && (this->indexmax_z < Pa.PMLz || this->indexmin_z >= Pa.PMLz + Pa.Nz))
    {

    }
    else if(this->indexmax_z < Pa.PMLz - this->border_h_Coord && this->indexmax_x < Pa.PMLx - this->border_h_Coord && this->indexmin_z >= Pa.PMLz + Pa.Nz + this->border_h_Coord && this->indexmin_x >= Pa.PMLx + Pa.Nx + this->border_h_Coord)
    {

    }
    else if(this->indexmin_x >= Pa.PMLx && this->indexmax_x < Pa.PMLx + Pa.Nx && this->indexmin_z >= Pa.PMLz && this->indexmax_z < Pa.PMLz + Pa.Nz)
    {

    }
    else if(this->indexmin_x < Pa.PMLx && this->indexmax_x >= Pa.PMLx)
    {
        if(this->indexmax_x < Pa.PMLx + Pa.Nx)
        {
            if(this->indexmax_z < Pa.PMLz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = Pa.PMLz;
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
            else if(this->indexmin_z >= Pa.PMLz && this->indexmax_z > Pa.PMLz)
            {
                if(this->indexmin_z < Pa.PMLz + Pa.Nz)
                {
                    this->h_Coord_num = 1;
                    H_Coord *coord = new H_Coord[this->h_Coord_num];
                    coord[0].indexmin_x = this->indexmin_x;
                    coord[0].indexmax_x = this->indexmax_x;
                    coord[0].indexmin_z = Pa.PMLz;
                    coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
                    coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

                    this->h_coord = coord;
                }
                else
                {

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
                coord[0].indexmin_x = Pa.PMLz;
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
                    coord[2].indexmax_x = this->indexmax_x;
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
                    coord[3].indexmax_x = Pa.PMLx + Pa.Nx;
                    coord[3].indexmin_z = Pa.PMLz + Pa.Nz;
                    coord[3].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->border_h_Coord - 1;
                    coord[3].length_x = coord[3].indexmax_x - coord[3].indexmin_x + 1;
                    coord[3].length_z = coord[3].indexmax_z - coord[3].indexmin_z + 1;

                    this->h_coord = coord;
                }
            }
            else if(this->indexmin_z >= Pa.PMLz && this->indexmin_z < Pa.PMLz + Pa.Nz && this->indexmax_z > Pa.PMLz)
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
                    coord[1].indexmax_x = this->indexmax_x ;
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
                    coord[0].indexmax_z = this->indexmax_z;
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
                    coord[0].indexmax_z = this->indexmax_z;
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
                    coord[1].indexmax_z = Pa.PMLz + Pa.Nz;
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
            else if(this->indexmin_z >= Pa.PMLz + Pa.Nz)
            {
                this->h_Coord_num = 1;
                H_Coord *coord = new H_Coord[this->h_Coord_num];
                coord[0].indexmin_x = this->indexmin_x;
                coord[0].indexmax_x = Pa.PMLx + Pa.Nx;
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
        coord[0].indexmax_z = this->indexmax_z < Pa.PMLz + Pa.Nz + this->border_h_Coord - 1 ? this->indexmax_z : Pa.PMLz + Pa.Nz + this->indexmax_z;
        coord[0].length_x = coord[0].indexmax_x - coord[0].indexmin_x + 1;
        coord[0].length_z = coord[0].indexmax_z - coord[0].indexmin_z + 1;

        this->h_coord = coord;
    }
    else
    {

    }
}

