本程序功能经过验证，与标准SM4加密结果相符

wbSM4.h与wbSM4.cpp是加密、生成查找表、生成仿射的过程中需要用到的函数声明与实现；

generator.cpp与gemerator.h是生成查找表与仿射的函数声明与实现；

creatTable.cpp文件执行生成查找表；

main.cpp文件执行加密；

本程序需要用到NTL库环境