#include "generator.h"

void init()
{
	//首先生成可逆矩阵与向量

	//生成可逆矩阵
	for (int i = 0; i < 36; i++)
	{
		creatInvMat(P[i], 32);
		if (i < 32)
		{
			creatInvMat(Q[i], 32);
			for (int j = 0; j < 4; j++)
				creatInvMat(sE[i][j], 8);
		}
	}
	//生成对角矩阵E
	for (int i = 0; i < 32; i++)
		creatDiagMat(E[i], sE[i]);

	//生成随机向量
	for (int i = 0; i < 36; i++)
	{
		random(p[i], 32);
		p1[i] = p[i];
		if (i < 32)
		{
			random(q[i], 32);
			random(e[i], 32);

			//将e向量切割，用于sE仿射的常数部分
			se[i][0] = getPartVector(e[i], 0, 8);
			se[i][1] = getPartVector(e[i], 8, 16);
			se[i][2] = getPartVector(e[i], 16, 24);
			se[i][3] = getPartVector(e[i], 24, 32);
		}
	}

}

//将仿射结构B、C写入文件中
void writeAffineTable(ofstream& file, affine_struct* A, string name)
{
	uint affineMatrix[32][32];
	uint affineVector[32];

	//将仿射数组中的元素转换为uint数组形式
	for (int i = 0; i < 32; i++)
	{
		affineVector[i] = vec2uint(A[i].vector);
		for (int j = 0; j < 32; j++)
			affineMatrix[i][j] = vec2uint(A[i].matrix[j]);
	}

	//开始向文件写入数据
	file << endl;

	//写入向量
	file << "uint " << name << "_vector[32] = {";
	for (int i = 0; i < 32; i++)
	{
		if (i != 31)
			file << affineVector[i] << ", ";
		else
			file << affineVector[i] << "};";
	}
	file << endl;
	
	//写入矩阵
	file << "uint " << name << "_matrix[32][32] = {" << endl;
	for (int i = 0; i < 32; i++)
	{
		if (i != 31)
		{
			file << "	{";
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << affineMatrix[i][j] << ", ";
				else
				{
					file << affineMatrix[i][j] << "}," << endl;
				}
			}
		}
		else
		{
			file << "	{";
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << affineMatrix[i][j] << ", ";
				else
				{
					file << affineMatrix[i][j] << "}" << endl;
					file << "};" << endl;
				}
			}
		}
	}
	file << endl;
}

//将仿射D写入文件
void writeAffineDTable(ofstream& file, affine_struct D[32][3], string name)
{
	uint affineMartix[96][32];
	uint affineVector[96];

	//将向量转换为uint数组形式
	int num = 0;
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 3; j++)
			affineVector[num++] = vec2uint(D[i][j].vector);
	}

	//将矩阵转换为uint数组形式
	int numM = 0;
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 32; k++)
				affineMartix[numM][k] = vec2uint(D[i][j].matrix[k]);
			numM++;
		}
	}


	//将向量数组写入文件
	file << endl;
	file << "uint " << name << "_vector[96] = {";
	for (int i = 0; i < 96; i++)
	{
		if (i != 95)
			file << affineVector[i] << ", ";
		else
			file << affineVector[i] << "};" << endl;
	}

	//将矩阵写入文件
	file << endl;
	file << "uint " << name << "_matrix[96][32] = {" << endl;
	for (int i = 0; i < 96; i++)
	{
		if (i != 95)
		{
			file << "	{";
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << affineMartix[i][j] << ", ";
				else
					file << affineMartix[i][j] << "}," << endl;
			}
		}
		else
		{
			file << "    {";
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << affineMartix[i][j] << ", ";
				else
					file << affineMartix[i][j] << "}" << endl << "};";
			}
		}
	}
}

//创建仿射表，并将其写入文件中
void creatAffineTable()
{

//进行仿射表的计算
	//计算仿射B
	for (int i = 0; i < 32; i++)
	{
		B[i].matrix = P[i + 4] * inv(Q[i]);
		B[i].vector = P[i + 4] * inv(Q[i]) * q[i];
	}

	//计算仿射C
	for (int i = 0; i < 32; i++)
	{
		C[i].matrix = P[i + 4] * inv(P[i]);
		C[i].vector = P[i + 4] * inv(P[i]) * p[i] + p1[i + 4];
	}

	//计算仿射D
	for (int i = 0; i < 32; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			D[i][j].matrix = E[i] * inv(P[i + j + 1]);
			D[i][j].vector = E[i] * inv(P[i + j + 1]) * p[i + j + 1] + e[i];
		}
	}

//将仿射表保存起来
	ofstream outFile;
	outFile.open("affineTable.h",ios::app);
	outFile << "typedef unsigned int uint;" << endl;

	//将仿射B、C、D存入头文件中，以uint的形式组织矩阵与向量
	writeAffineTable(outFile, B, "B");
	writeAffineTable(outFile, C, "C");
	writeAffineDTable(outFile, D, "D");
	
	outFile.close();
}

//将外部编码写入文件中
void writeExternalEncode()
{
	ofstream file;
	//将外部编码写入到文件中
	file.open("externalEncode.h", ios::app);
	file << "typedef unsigned int uint;" << endl;

	//先写入输入编码
	file << "uint IN[4][32] = {" << endl;
	for (int i = 0; i < 4; i++)
	{
		file << "    {";
		if (i != 3)
		{
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << vec2uint(P[i][j]) << ", ";
				else
					file << vec2uint(P[i][j]) << "}," << endl;
			}
		}
		else
		{
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << vec2uint(P[i][j]) << ", ";
				else
					file << vec2uint(P[i][j]) << "}" << endl << "};";
			}
		}
	}
	file << endl;
	//写入输入编码的常量部分
	file << "uint IN_vec[4] = {";
	for (int i = 0; i < 4; i++)
	{
		if (i != 3)
			file << vec2uint(p[i]) << ", ";
		else
			file << vec2uint(p[i]) << "};" << endl;
	}

	//写入输出编码
	file << "uint OUT[4][32] = {" << endl;
	for (int i = 0; i < 4; i++)
	{
		file << "    {";
		if (i != 3)
		{
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << vec2uint(P[32 + i][j]) << ", ";
				else
					file << vec2uint(P[32 + i][j]) << "}," << endl;
			}
		}
		else
		{
			for (int j = 0; j < 32; j++)
			{
				if (j != 31)
					file << vec2uint(P[32 + i][j]) << ", ";
				else
					file << vec2uint(P[32 + i][j]) << "}" << endl << "};";
			}
		}
	}

	//写入输出编码的常量部分
	file << endl;
	file << "uint OUT_vec[4] = {";
	for (int i = 32; i < 36; i++)
	{
		if (i != 35)
			file << vec2uint(p[i]) << ", ";
		else
			file << vec2uint(p[i]) << "};" << endl;
	}

	file.close();
}

//创建查找表
void createLookUpTable(uint* mainKey)
{
	//将矩阵M读入，该矩阵起着L函数的作用
	mat_GF2 L_M;
	random(L_M, 32, 32);
	for (int i = 0; i < 32; i++)
		L_M[i] = uint2vec(M_usefor_L[i]);

	//生成轮密钥，并将其转换为向量形式参与后续运算
	vec_GF2 vecKey[32];
	uint* roundKey = new uint[32];
	genRoundKey(mainKey, roundKey);
	for (int i = 0; i < 32; i++)
		vecKey[i] = uint2vec(roundKey[i]);
	delete[] roundKey;

	uint TABLE[128][256];
	int tableNum = 0;
	for (int r = 0; r < 32; r++)
	{
		mat_GF2 R = Q[r] * L_M;
		mat_GF2 sR[4];
		matSlice(8, R, sR);

		vec_GF2 Y;
		random(Y, 32);
		clear(Y);

		for (int j = 0; j < 4; j++)
		{
			byte plaintext = 0;
			for (int k = 0; k < 256; k++)
			{
				Y = sR[j] * S_box(inv(sE[r][j]) * (byte2vec(plaintext) +
					getPartVector(e[r], j * 8, j * 8 + 8) )+
					getPartVector(vecKey[r], j * 8, j * 8 + 8));
				if (j == 3)
					Y += q[r];
				TABLE[tableNum][k] = vec2uint(Y);
				plaintext++;
			}
			tableNum++;
		}
	}




	//将生成的查找表写入文件中
	ofstream tableFile;
	tableFile.open("table.h", ios::app);
	tableFile << "typedef unsigned int uint;" << endl;


	tableFile << endl;
	tableFile << "uint TABLE[128][256] = {" << endl;
	for (int i = 0; i < 128; i++)
	{
		tableFile << "    {";
		if (i != 127)
		{
			for (int j = 0; j < 256; j++)
			{
				if (j != 255)
					tableFile << TABLE[i][j] << ", ";
				else
					tableFile << TABLE[i][j] << "}," << endl;
			}
		}
		else
		{
			for (int j = 0; j < 256; j++)
			{
				if (j != 255)
					tableFile << TABLE[i][j] << ", ";
				else
					tableFile << TABLE[i][j] << "}" << endl << "};";
			}
		}
	}
	tableFile.close();
}

void generator(uint* mainKey)
{
	//首先初始化所有的矩阵与向量
	init();

	//生成仿射表，并写入文件
	creatAffineTable();

	//将外部编码写入文件
	writeExternalEncode();

	//依据生成的矩阵、向量、仿射表以及给定的密钥生成查找表
	createLookUpTable(mainKey);
}