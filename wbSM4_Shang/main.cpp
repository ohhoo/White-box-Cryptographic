#include <iostream>
#include <NTL/mat_GF2.h>
#include "wbSM4.h"
#include "table.h"
#include "IN.h"
#include "OUT.h"
using namespace NTL;
using namespace std;

void encryWBSM4(uint* ,uint*);

int main() {

    uint plaintext[4] = {0x01234567,0x89abcdef,0xfedcba98,0x76543210};
    uint mainKey[4] = {0x01234567,0x89abcdef,0xfedcba98,0x76543210};
    uint cipher[4] = {};

    encryWBSM4(plaintext,cipher);
    //genAffineTabel(mainKey);

    for(int i = 0;i<4;i++)
        printf("%08x",cipher[i]);
    return 0;
}


//加密
void encryWBSM4(uint* plain,uint* cipher)
{
    //对明文进行外部编码
    vec_GF2 plaintext[36];
    for(int i = 0;i<4;i++)
    {
        plaintext[i] = array2mat(IN[i]) * uint2vec(plain[i]);
    }

    //加密操作
    int tableNum = 0;
    for(int r = 0;r<32;r++)
    {
        affine_struct C;
        C.matrix = array2mat(C_Matrix[r]);
        C.vector = uint2vec(C_Vector[r]);

        vec_GF2 X1 = doAffine(C,plaintext[r]);

        affine_struct D[3];
        for(int j = 0;j<3;j++)
        {
            D[j].matrix = array2mat(D_Matrix[r * 3 + j]);
            D[j].vector = uint2vec(D_Vector[r * 3 + j]);
        }

        vec_GF2 X =doAffine(D[0],plaintext[r+1])+
                doAffine(D[1],plaintext[r+2])+
                doAffine(D[2],plaintext[r+3]);

        //将X分为四部分进行查表
        vec_GF2 Y;
        random(Y,32);
        clear(Y);
        for(int t = 0;t<4;t++)
        {
            Y += uint2vec( TABLE[tableNum+t][(uint)vec2byte(getPartVector(X,t*8,t*8+8))] );
        }

        affine_struct B;
        B.matrix = array2mat(B_Matrix[r]);
        B.vector = uint2vec(B_Vector[r]);

        vec_GF2 Y1 = doAffine(B,Y);

        plaintext[r+4] = Y1 + X1;

        tableNum += 4;
    }


    //对最后四个元素进行解码
    for(int i = 0;i<4;i++)
    {
        cipher[3-i] = vec2uint(inv(array2mat(OUT[i])) * plaintext[32+i]);
    }
}

