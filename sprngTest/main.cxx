/*    sprngTest.cxx
      Test the Ziggurat method for sprng and Normal distribution
      by Jerome Kaempf - EPFL - LESO-PB
      jerome.kaempf@epfl.ch
*/

using namespace std;

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <vector>
#include <sstream>

#include "../geometry.h"

int main(int argc, char *argv[])

{
    zigset(26041978);

    ostringstream oss1, oss2;

    for (int i=0; i<1e5; i++) {
        oss1 << randomUniform(0.0, 1.0) << "\t" << randomUniform(0.0, 1.0) << endl;
    }
    writeResultsOver("randomXY.dat", oss1);

    for (int i=0; i<1e5; i++) {
        oss2 << normallyDistributedSPRNG_Ziggurat() << endl;
    }
    writeResultsOver("Normal.dat", oss2);

    return 0;

}
