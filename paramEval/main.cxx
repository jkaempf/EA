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
    if (argc == 1) {
        cout << "Missing arguments: Problem type and Dimension" << endl;
        return 1;
    }
    else if (argc == 2) {

        if      (atoi(argv[1]) == 4) leProb = new NYCity();
        else if (string(argv[1]) == "RandomTerracesFlatRoof") leProb = new RandomTerracesFlatRoof();
        else if (string(argv[1]) == "RandomSlablsSlopedRoof") leProb = new RandomSlabsSlopedRoof();
        else if (string(argv[1]) == "RandomTerraceCourts")    leProb = new RandomTerraceCourtsEW();
        else if (string(argv[1]) == "Pavilions")              leProb = new Pavilions();
        //else if (atoi(argv[1]) == 11) leProb = new PepX();
        //else if (atoi(argv[1]) == 12) leProb = new PepY();
        else if (atoi(argv[1]) == 13) leProb = new RandomPavilionsFlatRoof();
        else if (atoi(argv[1]) == 14) leProb = new ProjetHuang();
        else if (atoi(argv[1]) == 15) leProb = new FourierBasis();
        else if (atoi(argv[1]) == 16) leProb = new ProjetLorenzetti();
        else {
            cout << "Missing correct arguments" << endl;
            return 1;
        }
    }
    else if (argc == 4) {

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
        readResultsInLine(argv[3], alleles);
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    // demander à l'utilisation ce qu'il veut faire
    unsigned int alleleN = 0;

//    string inputString;
//    cout << "Allele to vary: " << endl;
//    cin >> inputString;
//    unsigned int alleleN = to<unsigned int>(inputString);
//
//    cout << "Number of divisions: " << endl;
//    cin >> inputString;
//    unsigned int alleleSubDiv = to<unsigned int>(inputString);

    // tous les inputs nécessaires sont là

    // création de la famille de points - ou d'individus

    Heuristic<Individual> population(leProb); // vide, création des individus

    vector<int> combination(leProb->getSize(), 0.);
    for (size_t i=0;i<combination.size();++i) {
        combination[i] = (leProb->getMaxVector()[i]-leProb->getMinVector()[i])/leProb->getStepVector()[i];
        cout << i << " " << combination[i] << endl ;
    }

    int totalCombinations=combination[0];
    for (unsigned int i=1; i<combination.size(); i++) {
            totalCombinations *= combination[i];
    }

    for(unsigned int i=0; i<totalCombinations; ++i){
            vector<int> myvec(leProb->getSize(), 0);
            for(unsigned int j=0; j<i; j++){
                for(unsigned int k=0; k<combination.size(); k++){

                    if(myvec[k] +1 < combination[k]){
                        myvec[k] += 1;
                        break;
                    }
                    else{
                        myvec[k] = 0;
                    }
                }
            }
            vector<double> finalVec(myvec.size(),0.);
            for(size_t k = 0; k<myvec.size();++k){
                finalVec[k] = myvec[k] * leProb->getStepVector()[k] + leProb->getMinVector()[k] + TINY;
                //finalVec[k] = min(leProb->getMaxVector()[k],finalVec[k]);
                cout << finalVec[k] << " ";
            }
            cout << endl;
            population.pushIndividual(Individual(leProb,finalVec));
    }

//    for (unsigned int i=0; i<alleleSubDiv+2; i++) {
//
//        vector<double> allelesI = alleles;
//        for(unsigned int j=0; j<allelesI.size();++j){
//            cout <<allelesI.at(j)<<endl;
//        }
//        // modification de l'alleleN
//        allelesI[alleleN] = leProb->getMinVector()[alleleN] + double(i)*(leProb->getMaxVector()[alleleN]-leProb->getMinVector()[alleleN])/double(alleleSubDiv+1);
//        population.pushIndividual(Individual(leProb,allelesI));
//
//    }


//    population.pushIndividual(Individual(leProb,alleles));

    try {
        population.evaluate();
    }
    catch (string msg) {

        cout << msg << endl;
        return 1;

    }

    // write the results in the Stats.txt file
    ostringstream output;
    output << setprecision(12);
    // write the header
    output << "id,";
    for (unsigned int j=0; j<leProb->getSize(); ++j)
        output << "allele[" << j << "],";
    for (unsigned int j=0; j<population.getIndividual(0).getFitnessVector().size(); j++) {
        output << "fitness[" << j << "],";
    }
    output << endl;
    // write the values
    for (unsigned int i=0; i<population.getPopSize(); i++) {
        output << population.getIndividual(i).getIdNumber() << ",";
        for (unsigned int j=0; j<leProb->getSize(); ++j)
            output << population.getIndividual(i).getAllele(j) << ",";
        for (unsigned int j=0; j<population.getIndividual(i).getFitnessVector().size(); j++) {
            output << population.getIndividual(i).getFitnessVector()[j] << ",";
        }
        output << endl;
    }
    writeResultsOver("Stats.csv", output);

    return 0;
}
