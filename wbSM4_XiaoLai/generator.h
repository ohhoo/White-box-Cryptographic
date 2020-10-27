#pragma once
#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "wbSM4.h"
#include <fstream>
//声明需要用到的矩阵与向量

//可逆矩阵
static mat_GF2 P[36];
static mat_GF2 Q[32];
static mat_GF2 E[32];//对角矩阵

static mat_GF2 sE[32][4];//生成E的对角上的小矩阵

//向量
static vec_GF2 p[36];
static vec_GF2 p1[36];
static vec_GF2 q[32];
static vec_GF2 e[32];
static vec_GF2 se[32][4];

//仿射表
static affine_struct B[32], C[32], D[32][3];
//初始化矩阵与向量
void init();

//生成仿射表文件
void creatAffineTable();
//以下函数负责将仿射表写入文件中
void writeAffineTable(ofstream& file, affine_struct* A,string name);
void writeAffineDTable(ofstream& file, affine_struct D[32][3], string name);

//将外部编码写入文件中
void writeExternalEncode(ofstream& file);
//生成8进32出的查找表
void creatLookUpTable(uint* mainKey);


//输入密钥，生成查找表与相关仿射表
void generator(uint* mainKey);

#endif 