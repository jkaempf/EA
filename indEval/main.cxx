/*    genbuil.cxx
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
#include <map>

#include "../geometry.h"
#include "../problem.h"

int main(int argc, char* argv[])

{
    if (argc < 4 && argc > 7) {
        cerr << "Usage: genbuil <filenameIN> <idnumber> <filenameOUT> <problem#> <filenameSKY> <filenameSITE>" << endl;
        return 1;
    }

    vector<double> alleles;
    vector<double> fitness;

    Problem *leProb;

    // problem selection
    if (argc == 4) {

        if (string(argv[1]) == "EnergyPlus") leProb = new EnergyPlus(string(argv[2]));
        else if (string(argv[1]) == "CitySim") leProb = new CitySim(string(argv[2]));
        else {
            cout << "Missing correct arguments" << endl;
            return 1;
        }
    }
    else if (argc == 5) {
        if      (atoi(argv[4]) == 1) leProb = new Problem1();
        else if (atoi(argv[4]) == 4) leProb = new NYCity();
        else if (atoi(argv[4]) == 6) leProb = new Problem6();
        else if (atoi(argv[4]) == 7) leProb = new RandomTerracesFlatRoof();
        else if (atoi(argv[4]) == 8) leProb = new RandomSlabsSlopedRoof();
        else if (atoi(argv[4]) == 9) leProb = new RandomTerraceCourtsEW();
        else if (atoi(argv[4]) == 10) leProb = new Pavilions();
        else if (atoi(argv[4]) == 13) leProb = new RandomPavilionsFlatRoof();
        else if (atoi(argv[4]) == 14) leProb = new ProjetHuang();
        else if (atoi(argv[4]) == 15) leProb = new FourierBasis();
        else if (atoi(argv[4]) == 16) leProb = new ProjetLorenzetti();
	//else if (atoi(argv[4]) == 17) leProb = new FourierBasisMin();
	//else if (atoi(argv[4]) == 18) leProb = new FourierMinPlan();
    }
    else if (argc == 6) {
        if (argv[4] == "0") leProb = new Problem0(argv[5]);
        else if (string(argv[4]) == "17") leProb = new ProjetWagner(argv[5]);
        else if (string(argv[4]) == "EnergyPlus") leProb = new EnergyPlus(string(argv[5]));
	else if (string(argv[4]) == "CitySim") leProb = new CitySim(string(argv[5]));
    }
    else if (argc == 7) {
        if (atoi(argv[4]) == 11) leProb = new PepX(argv[5], argv[6]);
        if (atoi(argv[4]) == 12) leProb = new PepY(argv[5], argv[6]);
    }
    else return 1;

    try {
        readResultsInLine(argv[3], alleles);
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    try {
        cout << "Problem size: " << leProb->getSize() << "\t" << "Constraints: " << leProb->getnConstraints() << "\n" << "Individual not in field: " << leProb->notInField(alleles) << endl;
        fitness = leProb->evaluate(0, alleles);
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    cout << "Fitness: ";
	for (unsigned int i=0;i<fitness.size();++i)
		cout << fitness[i] << "\t";
	cout << endl;

    ostringstream ss;
    ss << setprecision(12);

	for (unsigned int i=0;i<fitness.size();++i)
		ss << fitness[i] << endl;

    writeResultsOver(string(argv[3])+".out", ss);

    return 0;

}
