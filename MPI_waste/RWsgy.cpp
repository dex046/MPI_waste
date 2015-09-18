/************************************************************************/
/*
* 文件信息：本程序包含使用C语言来读写Sgy文件的函数实现

* 文件名: RWsgy.cpp

* 版本：0.1

* 作者：高照奇

* 日期：2013-12-08
*/
/************************************************************************/
#include "RWsgy.h"

using namespace std;

/*
* function name: swap
* input: tni2
* output: void
*/
// swap_short_2
void swap(short *tni2)
{
    *tni2 = (((*tni2 >> 8) & 0xff) | (*tni2 & 0xff) << 8);
}

// swap_u_short_2
void swap(unsigned short *tni2)
{
    *tni2 = (((*tni2 >> 8) & 0xff) | (*tni2 & 0xff) << 8);
}

// swap_int_4
void swap(int *tni4)
{
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
        ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

//swap_u_int_4
void swap(unsigned int *tni4)
{
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
        ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

//swap_long_4
void swap(long *tni4)
{
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
        ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

//swap_u_long_4
void swap(unsigned long *tni4)
{
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
        ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

//swap_float_4
void swap(float *tnf4)
{
    int *tni4 = (int *) tnf4;
    *tni4 = (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
        ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}

//swap_double_8
void swap(double *tndd8)
{
    char *tnd8 = (char *)tndd8;
    char tnc;

    tnc = *tnd8;
    *tnd8 = *(tnd8+7);
    *(tnd8 + 7) = tnc;

    tnc = *(tnd8 + 1);
    *(tnd8 + 1) = *(tnd8 + 6);
    *(tnd8 + 6) = tnc;

    tnc = *(tnd8 + 2);
    *(tnd8 + 2) = *(tnd8 + 5);
    *(tnd8 + 5) = tnc;

    tnc = *(tnd8 + 3);
    *(tnd8 + 3) = *(tnd8 + 4);
    *(tnd8 + 4) = tnc;
}

// IBMSwarp
float IBMF4Swap(float x)
{
    static IBMFLOAT4 f4;
    static	unsigned b;
    int bb;
    float y = 0;
    f4.a = x;
    b = f4.ch[0]&127;
    bb = b - 64;
    double yy = float(y + f4.ch[1] / 256.0 + f4.ch[2] / 65536.0 + f4.ch[3] / 16777216.0);
    y = float(yy * pow(16.0, double(bb)));
    b = (f4.ch[0] & 128);
    if(b > 0)y = -y;
    return y;
}

/* 读取Sgy文件的相关信息*/
bool InfoOfSgy(char FileName[], REEL reel, unsigned short *TraceNum, unsigned short *SampleNum,
               unsigned short *SampleInt, short *DFormat, bool *BReel, bool *BIBM)
{
    FILE *fdata;
    fdata = fopen(FileName, "rb");
    if (fdata == NULL)
    {
        printf("File Don't Exit!");
        return false;
    }
    else
    {
        // 变量定义
        unsigned char f3200[3200];
        unsigned short TempSampleInt = 0, TempSampleNum = 0, TempTraceNum = 0;
        unsigned long FileSize = 0;
        bool TempBReel, TempBIBM;
        short TempDFormat = 0;
        Head240 tH, head;
        short ISize = 4;

        // 开始读Sgy数据
        fseek(fdata, 0, 0);
        // 读前面3200字节
        fread(f3200, 3200, 1, fdata);
        // 开始读卷头400字节
        fread(&(reel).reelstruct, 400, 1, fdata);//REEL 卷头
        TempSampleInt = reel.reelstruct.hdt;//该卷采样间隔
        TempSampleNum = reel.reelstruct.hns;//该卷每道的样点数

        // 开始读道头240个字节
        fread(&(head.headstruct), 240, 1, fdata);
        fseek(fdata, 0, 0);
        // 这一个道头用于试探该Sgy有没有卷头
        fread(&(tH.headstruct), 240, 1, fdata);

        // 求文件长度
        fseek(fdata, 0, SEEK_END);
        FileSize = ftell(fdata);

        // 判断是否为工作站格式
        if (reel.reelstruct.format > 255)
        {
            TempDFormat = reel.reelstruct.format;///Struct_reelb400 reelstruct
            swap(&TempDFormat);///short format///change
        }
        else
        {
            TempDFormat = reel.reelstruct.format;
        }

        // 判断存放数据用几个字节
        if (TempDFormat == 3)	ISize = 2;
        else if (TempDFormat == 8)	ISize = 1;
        else	ISize = 4;

        if ((FileSize - 3600) % (240 + TempSampleNum * ISize) != 0 || TempSampleNum < 0 ||
            TempSampleInt < 0 || head.headstruct.ns < 0 || head.headstruct.dt < 0)
        {
            if (TempSampleNum == head.headstruct.ns && TempSampleInt == head.headstruct.dt
                && ((FileSize % (240 + tH.headstruct.ns * ISize)) != 0 || tH.headstruct.ns == 0)
                && (reel.reelstruct.format > 255))
            {
                TempBReel = true;	// 有卷头
                TempBIBM = 1;		// IBM
                TempSampleNum = head.headstruct.ns;
                swap(&TempSampleNum);
                TempSampleInt = head.headstruct.dt;
                swap(&TempSampleInt);
                TempTraceNum = (FileSize - 3600) / (240 + ISize * TempSampleNum);
            }
            else
            {
                TempSampleNum = tH.headstruct.ns;
                TempSampleInt = tH.headstruct.dt;
                TempTraceNum = FileSize / (240 + TempSampleNum * 4);
                TempBReel = false;	// 无卷头
                TempBIBM = 0;		// IEEE
                TempDFormat = 0;
            }
        }
        else
        {
            if (TempSampleNum == head.headstruct.ns && TempSampleInt == head.headstruct.dt
                && TempSampleNum > 0 && TempSampleInt > 0 && TempSampleInt < 65535)
            {
                TempDFormat = reel.reelstruct.format;
                TempBIBM = 0;		// IEEE
                TempBReel = true;	// 有卷头
                TempTraceNum = (FileSize - 3600) / (240 + TempSampleNum * ISize);
            }
            else
            {
                TempSampleNum = tH.headstruct.ns;
                TempSampleInt = tH.headstruct.dt;
                TempTraceNum = (unsigned short)(FileSize / (240 + TempSampleNum * ISize));
                TempBReel = false;	// 无卷头
                TempBIBM = 0;		// IEEE
                TempDFormat = 0;
            }
        }

        *TraceNum = TempTraceNum;//
        *SampleNum = TempSampleNum;//number of samples in this trace 该卷每道的样点数
        *SampleInt = TempSampleInt;//sample interval; in micro-seconds 该卷采样间隔
        *DFormat = TempDFormat;
        *BIBM = TempBIBM;
        *BReel = TempBReel;

        fclose(fdata);

        return true;
    }
}

/* 读取Sgy中的地震道数据 */
bool ReadSgyData(char FileName[], Trace *trace, REEL reel, unsigned short *TraceNum,
                 unsigned short *SampleNum, unsigned short *SampleInt, short *DFormat,
                 bool *BReel, bool *BIBM)
{
    float *TempData1 = NULL;
    int *TempData2 = NULL;
    short *TempData3 = NULL;
    double *TempData4 = NULL;
    unsigned char f3200[3200];

    // 根据文件中数据的存储形式来开辟空间
    if (*DFormat == 1 || *DFormat == 5)
    {
        TempData1 = new float[*SampleNum];
        memset((void *)TempData1, 0, sizeof(float) * *SampleNum);
    }
    else if (*DFormat == 2)
    {
        TempData2 = new int[*SampleNum];
        memset((void *)TempData2, 0, sizeof(int) * *SampleNum);
    }
    else if (*DFormat == 3)
    {
        TempData3 = new short[*SampleNum];
        memset((void *)TempData3, 0, sizeof(short) * *SampleNum);
    }
    else if (*DFormat == 4)
    {
        TempData4 = new double[*SampleNum];
        memset((void *)TempData4, 0, sizeof(double) * *SampleNum);
    }

    FILE *fdata;
    fdata = fopen(FileName, "rb");
    fseek(fdata, 0, 0);

    // 如果有卷头
    if (*BReel)
    {
        // 读前面3200字节
        fread(f3200, 3200, 1, fdata);
        // 开始读卷头400字节
        fread(&(reel).reelstruct, 400, 1, fdata);
    }

    for (int i = 0; i < *TraceNum; i++)
    {
        // 读道头
        fread(&trace[i].head.h2, 2, 120, fdata);//h2 h4 headstruct
        // 如果是IBM的float
        if (*BIBM)
        {
            for (int j = 0; j < 120; j++)
            {
                if (j / 2 == 7 || j / 2 == 8 || j / 2 == 17 || (j >= 44 && j < 91))
                {
                    swap((short *)&trace[i].head.h2[j]);
                }
            }
            for (int j = 0; j < 60; j++)
            {
                if (j == 7 || j == 8 || j == 17 || (j >= 22 && j < 46));
                else
                {
                    swap(&trace[i].head.h4[j]);
                }
            }
        }

        trace[i].head.h2[114] = 0;// ?
        if (trace[i].head.h2[57] != *SampleNum)
        {
            return false;
        }

        // IBM float
        if (*DFormat == 1)
        {
            fread(TempData1, 4, *SampleNum, fdata);
            for (int ii = 0; ii < *SampleNum; ii++)
            {
                if (*BIBM)
                {
                    trace[i].data[ii] = IBMF4Swap(TempData1[ii]);
                }
                else
                {
                    trace[i].data[ii] = TempData1[ii];
                }
            }
        }

        // 4字节，两互补整数
        else if (*DFormat == 2)
        {
            fread(TempData2, 4, *SampleNum, fdata);
            for (int ii = 0; ii < *SampleNum; ii++)
            {
                if (*BIBM)
                {
                    swap((int *)&TempData2[ii]);
                }
                trace[i].data[ii] = (float)TempData2[ii];
            }
        }

        // 两字节，两互补整数
        else if (*DFormat == 3)
        {
            fread(TempData3, 2, *SampleNum, fdata);
            for (int ii = 0; ii < *SampleNum; ii++)
            {
                if (*BIBM)
                {
                    swap((short *)&TempData3[ii]);
                }
                trace[i].data[ii] = TempData3[ii];
            }
        }

        // IEEE浮点
        else if (*DFormat == 5)
        {
            fread(TempData1, 4, *SampleNum, fdata);
            for (int ii = 0; ii < *SampleNum; ii++)
            {
                if (*BIBM)
                {
                    swap(&TempData1[ii]);
                }
                trace[i].data[ii] = TempData1[ii];
            }
        }

        // 无卷头
        else
        {
            fread(trace[i].data, 4, *SampleNum, fdata);
        }
    }

    fclose(fdata);

    if (*DFormat == 1 || *DFormat == 5)
    {
        delete []TempData1;
    }
    else if (*DFormat == 2)
    {
        delete []TempData2;
    }
    else if (*DFormat == 3)
    {
        delete []TempData3;
    }
    else if (*DFormat == 4)
    {
        delete []TempData4;
    }

    return true;
}

/* 将数据写到Sgy文件中:写成微机格式，IEEE的浮点类型 */
bool WriteSgy(char FileName[], unsigned char *f3200, Trace *trace, unsigned short TraceNum, unsigned short SampleNum,
              unsigned short SampleInt)
{
    FILE *fp;
    fp = fopen(FileName, "wb");

    // 写卷头前3200个字节
    fwrite(&f3200[0], 3200, 1, fp);

    // 写卷头中400个字节
    REEL reel;
    reel.reelstruct.hns = SampleNum;
    reel.reelstruct.hdt = SampleInt;
    reel.reelstruct.format = 5; // IEEE float
    reel.reelstruct.mfeet = 1;
    fwrite(&reel.reelstruct, 400, 1, fp);

    // 写每一道数据
    for (int i = 0; i < TraceNum; i++)
    {
        trace[i].head.headstruct.cdp = i;
        trace[i].head.headstruct.ns = SampleNum;
        trace[i].head.headstruct.dt = SampleInt;
        trace[i].head.headstruct.sx = 100000;
        trace[i].head.headstruct.sy = 1000000 + i * 4;
        // 写道头
        fwrite(&trace[i].head.headstruct, 240, 1, fp);
        fwrite(&trace[i].data[0], 4, SampleNum, fp);
    }
    fclose(fp);
    return true;
}
