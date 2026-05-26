/*
For build:
mkdir -p build && g++ -o build/exp1 exp1.cpp ../*.cpp

For run:
./build/exp1

*/

// ##################################################################
#include <iostream>
#include "../geo_mag_declination.h"

int main(void)
{
    float dec = get_mag_declination_degrees_iran(35, 51);

    std::cout << "declination = " << dec << std::endl;
    return 0;
}