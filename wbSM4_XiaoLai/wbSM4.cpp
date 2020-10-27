#include "wbSM4.h"


uint shift(uint data, int n)
{
    return ((data << n) | (data >> (32 - n)));
}

uint byte2uint(byte a, byte b, byte c, byte d)
{
    return (uint)a << 24 | (uint)b << 16 | (uint)c << 8 | (uint)d;
}

uint L_Fun_key(uint u)
{
    return u ^ shift(u, 13) ^ shift(u, 23);
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
    return byte2uint(
        S_box((byte)(state >> 24)),
        S_box((byte)((state >> 16) & 0xff)),
        S_box((byte)((state >> 8) & 0xff)),
        S_box((byte)(state & 0xff))
    );
}

//生成轮密钥
void genRoundKey(uint* mainKey, uint* roundKey)
{
    for (int i = 0; i < 4; i++)
        mainKey[i] ^= FK[i];

    for (int i = 0; i < 32; i++)
    {
        uint state = 0;
        //后三个进行异或
        for (int j = 1; j < 4; j++)
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
        for (int k = 0; k < 3; k++)
            mainKey[k] = mainKey[k + 1];
        mainKey[3] = state;
        roundKey[i] = state;
    }
}

//创建可逆矩阵，mat为传入的矩阵参数，dim为矩阵维数
void creatInvMat(mat_GF2& mat, long dim)
{
    GF2 flag;
    mat_GF2 Invm;
    random(mat, dim, dim);
    inv(flag, Invm, mat);
    while (IsZero(flag))
    {
        random(mat, dim, dim);
        inv(flag, Invm, mat);
    }
}

//执行仿射变换
vec_GF2 doAffine(const affine_struct& affine, const vec_GF2& vector)
{
    vec_GF2 result;
    result = affine.matrix * vector + affine.vector;

    return result;
}

//仿射变换的逆操作
vec_GF2 doInvAffine(const affine_struct& affine, const vec_GF2& vector)
{
    vec_GF2 result;
    result = inv(affine.matrix) * (affine.vector + vector);
    return result;
}

//创建对角元素为较小可逆矩阵的对角矩阵
void creatDiagMat(mat_GF2& result, mat_GF2* a)
{
    Mat<mat_GF2> big;
    mat_GF2 middle;
    big.SetDims(4, 4);
    for (long i = 0; i < 4; i++)
    {
        for (long j = 0; j < 4; j++)
        {
            creatInvMat(middle, 8);
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
void matSlice(long nCol, const mat_GF2& bMat, mat_GF2* sMat)
{
    long dim = bMat.NumCols();
    if ((dim % nCol) != 0)
        throw "illega!";
    for (long i = 0; i < (dim / nCol); i++)
    {
        mat_GF2 state;
        random(state, 32, 8);
        for (long r = 0; r < 32; r++)
        {
            for (long c = 0; c < 8; c++)
            {
                state[r][c] = bMat[r][c + i * 8];
            }
        }
        sMat[i] = state;
    }
}

//取向量中的一部分
vec_GF2 getPartVector(const vec_GF2& v, long from, long to)
{
    vec_GF2 result;
    for (long i = from; i < to; i++)
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
        flag = (u >> (31 - sti)) & (0x1);
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
        flag = (b >> (7 - i)) & (0x1);
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

//合并两个向量
vec_GF2 comVec(vec_GF2 v1, vec_GF2 v2)
{
    vec_GF2 result = v1;
    for (int i = 0; i < v2.length(); i++)
        result.append(v2[i]);
    return result;
}

//将uint数组转换为矩阵
mat_GF2 array2mat(uint* u)
{
    mat_GF2 m;
    random(m, 32, 32);
    for (int i = 0; i < 32; i++)
        m[i] = uint2vec(u[i]);
    return m;
}
