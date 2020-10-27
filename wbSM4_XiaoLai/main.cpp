#include "affineTable.h"
#include "externalEncode.h"
#include "table.h"
#include "wbSM4.h"

void encrypt(uint* plain, uint* cipher);

int main()
{
	uint plaintext[4] = { 0x01234567,0x89abcdef,0xfedcba98,0x76543210 };
	uint cipher[4] = {};
	encrypt(plaintext, cipher);
	for (int i = 0; i < 4; i++)
		printf("%08x", cipher[i]);
	cout << endl;
	return 0;
}

void encrypt(uint* plain, uint* cipher)
{
	//首先对输入的明文进行外部编码；
	vec_GF2 plaintext[36];
	for (int i = 0; i < 4; i++)
	{
		plaintext[i] = array2mat(IN[i]) * uint2vec(plain[i]) + uint2vec(IN_vec[i]);
	}

	//32轮的加密
	int numTable = 0;
	for (int r = 0; r < 32; r++)
	{
		//准备在加密过程中用到的仿射
		affine_struct D_affine[3];
		for (int j = 0; j < 3; j++)
		{
			D_affine[j].matrix = array2mat(D_matrix[r * 3 + j]);
			D_affine[j].vector = uint2vec(D_vector[r * 3 + j]);
		}

		affine_struct C_affine;
		C_affine.matrix = array2mat(C_matrix[r]);
		C_affine.vector = uint2vec(C_vector[r]);


		//开始计算
		vec_GF2 X1 = doAffine(C_affine, plaintext[r]);

		vec_GF2 X = doAffine(D_affine[0], plaintext[r + 1]) + 
			doAffine(D_affine[1], plaintext[r + 2]) + 
			doAffine(D_affine[2], plaintext[r + 3]);

		//将X分为4部分进行查表
		vec_GF2 Y;
		random(Y, 32);
		clear(Y);
		for (int t = 0; t < 4; t++)
			Y += uint2vec(TABLE[numTable + t][(uint)vec2byte(getPartVector(X, t * 8, t * 8 + 8))]);
		
		//
		affine_struct B_affine;
		B_affine.matrix = array2mat(B_matrix[r]);
		B_affine.vector = uint2vec(B_vector[r]);

		vec_GF2 Y1 = doAffine(B_affine, Y);

		plaintext[r + 4] = Y1 + X1;

		numTable += 4;
	}

	
	for (int i = 0; i < 4; i++)
		cipher[i] = vec2uint(inv(array2mat(OUT[3 - i])) * (plaintext[35 - i] + uint2vec(OUT_vec[3 - i])));
}