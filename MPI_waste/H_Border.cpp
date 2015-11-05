#include "Partition.h"
H_Border::H_Border()
{

}
H_Border::H_Border(uint topborder, uint leftborder, uint bottomborder, uint rightborder)
{
    this->topborder = topborder;
    this->leftborder = leftborder;
    this->bottomborder = bottomborder;
    this->rightborder = rightborder;
}
H_Border::H_Border(uint length_x, uint legnth_z, uint topborder, uint leftborder, uint bottomborder, uint rightborder)
{
    this->length_x = length_x;
    this->length_z = length_z;

    this->topborder = topborder;
    this->leftborder = leftborder;
    this->bottomborder = bottomborder;
    this->rightborder = rightborder;
}
