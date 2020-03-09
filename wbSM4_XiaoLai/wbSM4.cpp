//
// Created by xudong on 2020/3/6.
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
 * 开始实现白盒算法
 * 该实现依据论文《对两个 SM4 白盒方案的分析》中
 * 对XiaoLai白盒SM4的描述而进行。
 */

//首先生成随机的可逆矩阵以及向量用以构造仿射结构，并存储相应的仿射结构以及查找表
void creatAffineTable(uint* mainKey)
{
    mat_GF2 P[36],E[32],Q[32],Ei[32][4];
    vec_GF2 p[36],e[32],q[32],p1[32];

    //生成随机可逆矩阵与随机向量
    for(int i = 0;i<32;i++)
    {
        //生成对角矩阵
        creatDiagMat(E[i],Ei[i]);
        //生成可逆矩阵
        creatInvMat(Q[i],32);

        random(e[i],32);
        random(q[i],32);
        random(p1[i],32);
    }

    for(int i = 0;i<36;i++)
    {
        creatInvMat(P[i],32);
        random(p[i],32);
    }

    //利用生成的可逆矩阵与向量复合成仿射结构
    affine_struct B[32],C[32],D[32][3];
    for(int i = 0;i<32;i++)
    {
        //仿射结构B
        B[i].matrix = P[i+4] * inv(Q[i]);
        B[i].vector = P[i+4] * inv(Q[i]) * q[i];

        //仿射结构C
        C[i].matrix = P[i+4] * inv(P[i]);
        C[i].vector = P[i+4] * inv(P[i]) * p[i] + p1[i];

        //仿射结构D
        for(int j = 0;j<3;j++)
        {
            D[i][j].matrix = E[i] * inv(P[i+j+1]);
            D[i][j].vector = E[i] * inv(P[i+j+1]) * p[i+j+1] +e[i];
        }
    }

    //将仿射写入文件
    ofstream outfile_B,outfile_C,outfile_D;
    outfile_B.open("affine_B",ios::out);
    outfile_C.open("affine_C",ios::out);
    outfile_D.open("affine_D",ios::out);

    for(int i = 0;i<32;i++)
    {
        outfile_B<<B[i].matrix<<" "<<B[i].vector;
        outfile_C<<C[i].matrix<<" "<<C[i].vector;
        for(int j = 0;j<3;j++)
            outfile_D<<D[i][j].matrix<<" "<<D[i][j].vector;
    }
    outfile_B.close();
    outfile_C.close();
    outfile_D.close();

    //将外部编码 读入文件，用矩阵数组P[]以及向量数组p的前四个组成输入外部编码，最后四个组成输出外部解码
    ofstream DECODE,INCODE;
    DECODE.open("DECODE",ios::out);
    INCODE.open("INCODE",ios::out);
    for(int i = 0;i<4;i++)
    {
        DECODE<<P[35-i]<<" "<<p[35-i]<<endl;
        INCODE<<P[i]<<" "<<p[i]<<endl;
    }
    DECODE.close();
    INCODE.close();

    //将用作L函数的矩阵读入
    mat_GF2 M;
    ifstream infile;
    infile.open("M",ios::in);
    infile >> M;
    infile.close();

    //生成轮密钥，并将其转换为向量形式
    uint * roundKey = new uint[32];
    vec_GF2 *vecKey = new vec_GF2[32];
    genRoundKey(mainKey,roundKey);
    for(int i = 0;i<32;i++)
        vecKey[i] = uint2vec(roundKey[i]);

    //生成查找表
    string NAME = "TABLE_";
    int tableNumber = 0;//生成查找表文件名

    //生成查找表
    for(long r = 0;r<32;r++)
    {
        //生成矩阵R，并将其切分
        mat_GF2 R = Q[r] * M;
        mat_GF2 sR[4];
        matSlice(8,R,sR);

        for(int i = 0;i<4;i++)
        {
            vec_GF2 useE = getPartVector(e[r],i*8,i*8+8);
            vec_GF2 useKey = getPartVector(vecKey[r],i*8,i*8+8);
            vec_GF2 cipherText[256];

            byte plaintext = 0;
            for(int k = 0;k<256;k++)
            {
                vec_GF2 state;
                state = sR[i] * S_box(inv(Ei[r][i]) * (useE + byte2vec(plaintext)) + useKey);
                if(i == 3)
                    state += q[r];
                cipherText[k] = state;
                ++plaintext;
            }

            //将得到的密文写入文件中形成查找表
            string fileName = NAME + to_string(tableNumber);
            ofstream outfile_array;
            outfile_array.open(fileName,ios::out);
            for(int j = 0;j<256;j++)
                outfile_array<<cipherText[j];
            outfile_array.close();
            ++tableNumber;
        }
    }

    delete[] roundKey;
    delete[] vecKey;
}

//对传入的明文进行加密
void wbSM4En(uint* plaintext,uint* cipher)
{
    //将仿射从文件中读入内存
    affine_struct B[32],C[32],D[32][3];
    ifstream infile_B,infile_C,infile_D;
    infile_B.open("affine_B",ios::in);
    infile_C.open("affine_C",ios::in);
    infile_D.open("affine_D",ios::in);

    for(int i = 0;i<32;i++)
    {
        infile_B>>B[i].matrix>>B[i].vector;
        infile_C>>C[i].matrix>>C[i].vector;
        for(int j = 0;j<3;j++)
            infile_D>>D[i][j].matrix>>D[i][j].vector;
    }
    infile_B.close();
    infile_C.close();
    infile_D.close();

    //将外部编码读入内存
    affine_struct IN[4],OUT[4];
    ifstream IN_CODE,OUT_CODE;
    IN_CODE.open("INCODE",ios::in);
    OUT_CODE.open("DECODE",ios::in);
    for(int i = 0;i<4;i++)
    {
        IN_CODE>>IN[i].matrix>>IN[i].vector;
        OUT_CODE>>OUT[i].matrix>>OUT[i].vector;
    }
    IN_CODE.close();
    OUT_CODE.close();

    //将查找表读入内存
    vec_GF2 table[128][256];
    string NAME = "TABLE_";
    for(int i = 0;i<128;i++)
    {
        string fileName = NAME + to_string(i);
        ifstream infile_table;
        infile_table.open(fileName,ios::in);
        for(int k = 0;k<256;k++)
            infile_table>>table[i][k];
        infile_table.close();
    }

    //将输入转换为向量，进行输入外部编码
    vec_GF2 vecPlaintext[36];
    for(int i = 0;i<4;i++)
    {
        //转换为向量
        vecPlaintext[i] = uint2vec(plaintext[i]);
        //进行编码
        vecPlaintext[i] = doAffine(IN[i],vecPlaintext[i]);
    }


    //计算密文
    int numOfTable = 0;
    for(int i = 0;i<32;i++)
    {
        vec_GF2 result;
        random(result,32);
        clear(result);
        //后三个32bit异或
        for(int j = 0;j<3;j++)
            result += doAffine(D[i][j],vecPlaintext[1+j]);

        //异或结果切割为4部分进行查表，将四个查表结果异或
        vec_GF2 z;
        z = table[numOfTable][(uint)vec2byte(getPartVector(result,0,8))] +
                table[numOfTable+1][(uint)vec2byte(getPartVector(result,8,16))] +
                table[numOfTable+2][(uint)vec2byte(getPartVector(result,16,24))] +
                table[numOfTable+3][(uint)vec2byte(getPartVector(result,24,32))];

        //将第一部分与经过处理的后三部分之和进行异或
        vec_GF2 vecCipherText;
        vecCipherText = doAffine(B[i],z) + doAffine(C[i],vecPlaintext[i]);
        vecPlaintext[i+4] = vecCipherText;

        numOfTable += 4;
    }

    //对结果进行解码
    for(int i = 0;i<4;i++)
    {
        cipher[i] = vec2uint(doInvAffine(OUT[i],vecPlaintext[32+i]));
    }
}