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
#include "../heuristic.h"

using namespace std;

unsigned int Individual::compteur = 0;

int main(int argc, char* argv[])
{

    if (argc < 4 && argc > 7) {
        cerr << "Usage: paramEval <filenameIN> <idnumber> <filenameOUT> <problem#> <filenameSKY> <filenameSITE>" << endl;
        return 1;
    }

    vector<double> alleles;

    Problem *leProb;

    // problem selection
    if (argc == 5) {
        if      (atoi(argv[4]) == 1) leProb = new Problem1();
        else if (atoi(argv[4]) == 4) leProb = new NYCity();
        else if (atoi(argv[4]) == 6) leProb = new Problem6();
        else if (string(argv[4]) == "RandomTerracesFlatRoof") leProb = new RandomTerracesFlatRoof();
        else if (string(argv[4]) == "RandomSlabsSlopedRoof")  leProb = new RandomSlabsSlopedRoof();
        else if (string(argv[4]) == "RandomTerraceCourts") leProb = new RandomTerraceCourtsEW();
        else if (string(argv[4]) == "Pavilions") leProb = new Pavilions();
        else if (atoi(argv[4]) == 13) leProb = new RandomPavilionsFlatRoof();
        else if (atoi(argv[4]) == 14) leProb = new ProjetHuang();
        else if (atoi(argv[4]) == 15) leProb = new FourierBasis();
        else if (atoi(argv[4]) == 16) leProb = new ProjetLorenzetti();
    }
    else if (argc == 6) {
        if (argv[4] == "0") leProb = new Problem0(argv[5]);
        else if (string(argv[4]) == "17") leProb = new ProjetWagner(argv[5]);
        else if (string(argv[4]) == "EnergyPlus") leProb = new EnergyPlus(string(argv[5]));
    }
    else if (argc == 7) {
        if (atoi(argv[4]) == 11) leProb = new PepX(argv[5], argv[6]);
        if (atoi(argv[4]) == 12) leProb = new PepY(argv[5], argv[6]);
    }
    else return 1;

    try {
        readResultsInLine(argv[1], alleles);
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    // demander à l'utilisation ce qu'il veut faire

    string inputString;
    cout << "Allele to vary: " << endl;
    cin >> inputString;
    unsigned int alleleN = to<unsigned int>(inputString);

    cout << "Number of divisions: " << endl;
    cin >> inputString;
    unsigned int alleleSubDiv = to<unsigned int>(inputString);

    // tous les inputs nécessaires sont là

    // création de la famille de points - ou d'individus

    Heuristic<Individual> population(leProb); // vide, création des individus

    for (unsigned int i=0; i<alleleSubDiv+2; i++) {

        vector<double> allelesI = alleles;
        // modification de l'alleleN
        allelesI[alleleN] = leProb->getMinVector()[alleleN] + double(i)*(leProb->getMaxVector()[alleleN]-leProb->getMinVector()[alleleN])/double(alleleSubDiv+1);
        population.pushIndividual(Individual(leProb,allelesI));

    }

    population.pushIndividual(Individual(leProb,alleles));

    try {
        population.evaluate();
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    ostringstream output;
    output << setprecision(12);

    for (unsigned int i=0; i<population.getPopSize(); i++) {

        output << population.getIndividual(i).getAllele(alleleN) << "\t";
        for (unsigned int j=0; j<population.getIndividual(i).getFitnessVector().size(); j++) {
            output << population.getIndividual(i).getFitnessVector()[j] << "\t";
        }
        output << endl;

    }

    writeResultsOver(argv[3], output);

    return 0;
}
