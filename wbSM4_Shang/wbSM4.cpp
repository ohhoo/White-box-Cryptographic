//
// Created by xudong on 2020/3/16.
//

#include "wbSM4.h"
#include <fstream>
#include <string>
/**
 * 以下是对wbSM4.h中的方法的实现
 */

//实现循环移位
uint move(uint u,int n)
{
    return ((u << n) | (u >> ( UINT_LONG - n)));
}

//四个byte合并成一个uint
uint byte2uint(byte a, byte b,byte c,byte d)
{
    return (uint)a << 24 | (uint)b << 16 | (uint)c << 8 | (uint)d;
}

//用于轮密钥生成的L函数
uint L_Fun_key(uint u)
{
    return u ^ move(u,13) ^ move(u,23);
}

//对byte执行S盒查询操作，返回byte
byte S_box(byte state)
{
    int column = (uint)(state & 0xf);
    int row = (uint)(state >> 4);
    return Sbox[row * 16 + column];
}

//对uint执行S盒查询操作，返回uint
uint S_box(uint state)
{
    return byte2uint(S_box((byte)(state >> 24)), S_box((byte)((state >> 16) & 0xff)),
                     S_box((byte)((state >> 8) & 0xff)), S_box((byte)(state & 0xff)));
}

//生成轮密钥
void genRoundKey(uint* mainKey,uint* roundKey)
{
    for(int i = 0;i < 4;i++)
        mainKey[i] ^= FK[i];

    for(int i = 0;i < 32;i++)
    {
        uint state = 0;
        //后三个进行异或
        for(int j = 1;j<4;j++)
            state ^= mainKey[j];
        //异或CK常数
        state ^= CK[i];

        //S盒
        state = S_box(state);

        //L函数
        state = L_Fun_key(state);
        //与第一个进行异或
        state ^= mainKey[0];

        //移位，将新的结果加入，以进行下一轮
        for(int k = 0;k<3;k++)
            mainKey[k] = mainKey[k+1];
        mainKey[3] = state;
        roundKey[i] = state;
    }
}

//生成可逆矩阵
void creatInvMat(mat_GF2& mat,long dim)
{
    GF2 flag;
    mat_GF2 Invm;
    random(mat,dim,dim);
    inv(flag,Invm,mat);
    while(IsZero(flag))
    {
        random(mat,dim,dim);
        inv(flag,Invm,mat);
    }
}

//执行仿射变换
vec_GF2 doAffine(const affine_struct& affine,const vec_GF2& vector)
{
    vec_GF2 result;
    result = affine.matrix * vector +affine.vector;

    return result;
}

//仿射变换的逆操作
vec_GF2 doInvAffine(const affine_struct& affine,const vec_GF2& vector)
{
    vec_GF2 result;
    result = inv(affine.matrix) * (affine.vector + vector);
    return result;
}

//创建对角元素为较小可逆矩阵的对角矩阵
void creatDiagMat(mat_GF2& result,mat_GF2* a)
{
    Mat<mat_GF2> big;
    mat_GF2 middle;
    big.SetDims(4, 4);
    for (long i = 0; i < 4; i++)
    {
        for (long j = 0; j < 4; j++)
        {
            creatInvMat(middle,8);
            if (j == i)
            {
                big.put(i, j, middle);
                a[j] = middle;
            }
            else
            {
                clear(middle);
                big.put(i, j, middle);
            }
        }
    }
    random(result, 32, 32);
    long R = 0;
    for (long i = 0; i < 4; i++)
    {
        long C = 0;
        for (long j = 0; j < 4; j++)
        {
            for (long k = 0; k < 8; k++)
            {
                for (long L = 0; L < 8; L++)
                {
                    result[k + R][L + C] = (big[i][j])[k][L];
                }
            }
            C += 8;
        }
        R += 8;
    }
}

//将矩阵以若干列为单位进行切割
void matSlice(long nCol,const mat_GF2& bMat,mat_GF2* sMat)
{
    long dim = bMat.NumCols();
    if((dim % nCol) != 0)
        throw "illega!";
    for (long i = 0; i < (dim/nCol); i++)
    {
        mat_GF2 state;
        random(state, 32, 8);
        for (long r = 0; r < 32; r++)
        {
            for (long c = 0; c < 8; c++)
            {
                state[r][c] = bMat[r][c+i*8];
            }
        }
        sMat[i] = state;
    }
}

//取向量中的一部分
vec_GF2 getPartVector(const vec_GF2& v,long from,long to)
{
    vec_GF2 result;
    for(long i = from;i < to;i++)
        result.append(v[i]);
    return result;
}

//向量与uint byte数据的相互转化
uint vec2uint(const vec_GF2& vector)
{
    if (vector.length() == 32)
    {
        uint result = 0;
        for (long i = 0; i < 32; i++)
        {
            if (IsOne(vector[i]))
                result += (uint)pow(2, 31 - i);
        }
        return result;
    }
}

vec_GF2 uint2vec(const uint u)
{
    vec_GF2 state;
    for (long sti = 0; sti < 32; sti++)
    {
        GF2 flag;
        flag = (u >> (31 - sti))&(0x1);
        state.append(flag);
    }
    return state;
}

vec_GF2 byte2vec(const byte b)
{
    vec_GF2 state;
    for (long i = 0; i < 8; i++)
    {
        GF2 flag;
        flag = (b >> (7 - i))&(0x1);
        state.append(flag);
    }
    return state;
}

/*
将GF(2)8位向量转换为byte类型数据
*/
byte vec2byte(const vec_GF2& v)
{
    byte result = 0;
    if (v.length() == 8)
    {
        for (long i = 0; i < v.length(); i++)
        {
            if (IsOne(v[i]))
                result += (byte)pow(2, 7 - i);
        }
    }
    return result;
}


//向量的S盒替换函数
vec_GF2 S_box(vec_GF2 v)
{
    return byte2vec(S_box(vec2byte(v)));
}

/*
 * 根据尚培对肖雅莹白盒SMS4的改进，实现代码
 */
void writeAffineFile(affine_struct *A,ofstream & out,string affineName)
{
    uint affine_Matrix[32][32];
    uint affine_Vector[32];

    for(int i = 0;i<32;i++)
    {
        affine_Vector[i] = vec2uint(A[i].vector);
        for(int j = 0;j<32;j++)
            affine_Matrix[i][j] = vec2uint(A[i].matrix[j]);
    }

    //将向量数组写入文件
    out<<endl;
    out<<"uint "<<affineName<<"_Vector[32] = {";
    for(int i = 0;i<32;i++)
    {
        if(i != 31)
            out<<affine_Vector[i]<<",";
        else
            out<<affine_Vector[i]<<"};";
    }
    out<<endl;

    //将矩阵数组写入文件
    out<<"uint "<<affineName<<"_Matrix[32][32] = {";
    for(int i = 0;i<32;i++)
    {
        if(i!=31)
        {
            out << "{";
            for (int j = 0; j < 32; j++) {
                if (j != 31)
                    out << affine_Matrix[i][j] << ",";
                else {
                    out << affine_Matrix[i][j] << "},";
                    out<<endl;
                }
            }
        }
        else
        {
            out << "{";
            for (int j = 0; j < 32; j++)
            {
                if (j != 31)
                    out << affine_Matrix[i][j] << ",";
                else
                    out << affine_Matrix[i][j] << "}};";
            }
        }
    }
    out<<endl;
}
//将D仿射写入文件
void writeAffineDFile(affine_struct D[32][3],ofstream &outfile,string affineName)
{
    uint affine_Matrix[96][32];
    uint affine_Vector[96];

    int num = 0;
    for(int i = 0;i<32;i++)
    {
        for(int j = 0;j<3;j++)
        {
            affine_Vector[num] = vec2uint(D[i][j].vector);
            num++;
        }
    }

    int numM = 0;
    for(int i = 0;i<32;i++)
    {
        for(int j = 0;j<3;j++)
        {
            for(int k = 0;k<32;k++)
            {
                affine_Matrix[numM][k] = vec2uint(D[i][j].matrix[k]);
            }
            numM++;
        }
    }


    //将向量数组写入文件
    outfile<<endl;
    outfile<<"uint "<<affineName<<"_Vector[96] = {";
    for(int i = 0;i<96;i++)
    {
        if(i != 95)
            outfile<<affine_Vector[i]<<",";
        else
            outfile<<affine_Vector[i]<<"};";
    }
    outfile<<endl;

    //将矩阵写入数组
    outfile<<endl;
    outfile<<"uint "<<affineName<<"_Matrix[96][32] = {";
    for(int i = 0;i<96;i++)
    {
        if(i!=95)
        {
            outfile << "{";
            for (int j = 0; j < 32; j++) {
                if (j != 31)
                    outfile << affine_Matrix[i][j] << ",";
                else {
                    outfile << affine_Matrix[i][j] << "},";
                    outfile<<endl;
                }
            }
        }
        else
        {
            outfile << "{";
            for (int j = 0; j < 32; j++)
            {
                if (j != 31)
                    outfile << affine_Matrix[i][j] << ",";
                else
                    outfile << affine_Matrix[i][j] << "}};";
            }
        }
    }
    outfile<<endl;
}


void genAffineTabel(uint* mainKey)
{
    //创建文件，将仿射以及查找表以uint数组的形式存储在该文件中
    ofstream outFile;
    outFile.open("table.h",ios::app);
    outFile<<"typedef unsigned int uint;"<<endl;


    //先创建构成仿射的矩阵,向量
    mat_GF2 P[36],E[32],sE[32][4],Q[32];
    vec_GF2 p[32],e[32],q[32];

    for(int i = 0;i<36;i++)
    {
        creatInvMat(P[i],32);
        if(i>=4)
        {
            creatInvMat(Q[i-4],32);
            creatDiagMat(E[i-4],sE[i-4]);
            random(e[i-4],32);
            random(q[i-4],32);
            random(p[i-4],32);
        }
    }


    affine_struct B[32],C[32],D[32][3];

    ofstream IN_CODE,OUT_DECODE;

    IN_CODE.open("IN.h",ios::out);
    OUT_DECODE.open("OUT.h",ios::out);

    //将外部编码分别写入两个头文件中
    IN_CODE<<"typedef unsigned int uint;"<<endl<<"uint IN[4][32] = {";
    OUT_DECODE<<"typedef unsigned int uint;"<<endl<<"uint OUT[4][32] = {";
    for(int i = 0;i<4;i++)
    {
        IN_CODE<<"{";
        OUT_DECODE<<"{";
        if(i!=3)
        {
            for(int j = 0;j<32;j++)
            {
                if(j!=31)
                {
                    IN_CODE<<vec2uint(P[i][j])<<",";
                    OUT_DECODE<<vec2uint(P[32+i][j])<<",";
                } else{
                    IN_CODE<<vec2uint(P[i][j])<<"},";
                    OUT_DECODE<<vec2uint(P[32+i][j])<<"},";
                }
            }
        } else{
            for(int j = 0;j<32;j++)
            {
                if(j!=31)
                {
                    IN_CODE<<vec2uint(P[i][j])<<",";
                    OUT_DECODE<<vec2uint(P[32+i][j])<<",";
                } else{
                    IN_CODE<<vec2uint(P[i][j])<<"}}";
                    OUT_DECODE<<vec2uint(P[32+i][j])<<"}}";
                }
            }
        }
    }
    IN_CODE<<";";
    OUT_DECODE<<";";



    for(int i = 0;i < 32;i++)
    {
        //B仿射由Q仿射的逆与输出编码（仿射）构成
        B[i].matrix = P[i+4] * inv(Q[i]);
        B[i].vector = P[i+4] * inv(Q[i]) * q[i] + p[i];

        //C仿射由对上一轮输出编码的解码部分（逆矩阵）与本轮的编码部分（仿射）构成
        C[i].matrix = P[i+4] * inv(P[i]);
        C[i].vector = p[i];

        for(int j = 0;j<3;j++)
        {
            //D仿射包括对上一轮输出编码的解码（逆矩阵）与本轮的输入编码（仿射）
            D[i][j].matrix = inv(E[i]) * inv(P[i+1+j]);
            D[i][j].vector = (inv(E[i]) * e[i]);
        }
    }

    //将要保存的仿射结构存储在文件中
    writeAffineFile(B,outFile,"B");
    writeAffineFile(C,outFile,"C");
    writeAffineDFile(D,outFile,"D");


    //将M矩阵读入，该矩阵起L函数的功能
    mat_GF2 L_M;
    random(L_M,32,32);
    for(int i = 0;i<32;i++)
        L_M[i] = uint2vec(M_usefor_L[i]);




    //生成轮密钥，并将其转换为向量形式
    vec_GF2 vecKey[32];
    uint *roundKey = new uint[32];
    genRoundKey(mainKey,roundKey);
    for(int i = 0;i<32;i++)
        vecKey[i] = uint2vec(roundKey[i]);
    delete[] roundKey;

    uint TABLE[128][256];
    int tableNum = 0;
    for(int r = 0;r<32;r++)
    {
        mat_GF2 R = Q[r]*L_M;
        mat_GF2 sR[4];
        matSlice(8,R,sR);


        vec_GF2 Y;
        random(Y,32);
        clear(Y);
        for(int j = 0;j<4;j++)
        {

            byte plaintext = 0;
            for (int k = 0; k < 256; k++) {
                Y = sR[j] * S_box(sE[r][j] * byte2vec(plaintext) +
                                  getPartVector(e[r], j * 8, j * 8 + 8) + getPartVector(vecKey[r], j * 8, j * 8 + 8));
                if (j == 3)
                    Y += q[r];
                TABLE[tableNum][k] = vec2uint(Y);
                plaintext++;
            }
            tableNum++;
        }
    }

    outFile<<"uint TABLE[128][256] = {";
    for(int i = 0;i<128;i++)
    {
        if(i!=127)
        {
            outFile << "{";
            for (int j = 0; j < 256; j++) {
                if (j != 255)
                    outFile << TABLE[i][j] << ",";
                else {
                    outFile << TABLE[i][j] << "},";
                    outFile<<endl;
                }
            }
        }
        else
        {
            outFile << "{";
            for (int j = 0; j < 256; j++)
            {
                if (j != 255)
                    outFile << TABLE[i][j] << ",";
                else
                    outFile << TABLE[i][j] << "}}";
            }
        }
    }
    outFile<<";";

    outFile.close();
}

//将32个元素的数组转换为矩阵
mat_GF2 array2mat(uint* u)
{
    mat_GF2 m;
    random(m,32,32);
    for(int i = 0;i<32;i++)
        m[i] = uint2vec(u[i]);
    return m;
}
