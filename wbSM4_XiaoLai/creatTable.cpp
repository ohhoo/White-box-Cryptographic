#include <NTL/vec_GF2.h>
#include <NTL/mat_GF2.h>
#include <iostream>
#include <fstream>
#include "generator.h"

using namespace NTL;

int main()
{
    uint mainKey[4] = { 0x01234567,0x89abcdef,0xfedcba98,0x76543210 };
    generator(mainKey);
    return 0;
}
