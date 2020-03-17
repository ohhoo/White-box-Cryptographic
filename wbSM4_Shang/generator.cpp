#include <iostream>
#include "wbSM4.h"

int main()
{
    uint mainKey[4] = {0x01234567,0x89abcdef,0xfedcba98,0x76543210};
    genAffineTabel(mainKey);
    return 0;
}
