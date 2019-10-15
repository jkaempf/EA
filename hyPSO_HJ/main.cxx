/*    Genbuil Version 2.0
      Test functions & Other problems
      Evolution Algorithms
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
#include <numeric>
#include <map>
#include <algorithm>
#include <typeinfo>

#include "../geometry.h"
#include "../problem.h"
#include "../heuristic.h"

unsigned int Individual::compteur = 0;

int main(int argc, char *argv[])

{
    // timer start
    time_t start, end;

    time(&start); // start timer

    // intiate random number generator
    zigset(3021980);

    // PSO parameters
    unsigned int np = 100;
    double c1 = 2.05;
    double c2 = 2.05;
    double lambdaPSO = 0.2;
    double kappa = 1.0;

    // HJ parameters
    unsigned int r = 2;
    unsigned int s0 = 0;
    unsigned int t = 1;
    unsigned int m = 4;

    // write parameters
    ostringstream param;
    param << "PSO parameters: np=" << np << ", c1=" << c1 << ", c2=" << c2 << ", lambda=" << lambdaPSO << ", kappa=" << kappa << endl;
    param << "HJ parameters: r=" << r << ", s0=" << s0 << ", t=" << t << ", m=" << m << endl;
    writeResultsOver("param.dat", param);

for (int loopi=0; loopi < 1; loopi++) {

    // stores the performance in file
    ostringstream statRun, fileStatRunPSO, fileStatRunHJ;
    fileStatRunPSO << "StatRunPSO" << loopi << ".dat";
    fileStatRunHJ << "StatRunHJ" << loopi << ".dat";

    // Limiting parameters
    unsigned int nEvalLimit = 3000;

    // choose the problem

    Problem *leProb;

    if (argc == 1) {
        cout << "Missing arguments: Problem type and Dimension" << endl;
    }
    else if (argc == 2) {

        if      (atoi(argv[1]) == 4) leProb = new NYCity();
        else if (atoi(argv[1]) == 7) leProb = new RandomTerracesFlatRoof();
        else if (atoi(argv[1]) == 8) leProb = new RandomSlabsSlopedRoof();
        else if (atoi(argv[1]) == 9) leProb = new RandomTerraceCourtsEW();
        else if (atoi(argv[1]) == 10) leProb = new Pavilions();
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
    else if (argc == 3) {

        if (string(argv[1]) == "EnergyPlus") leProb = new EnergyPlus(string(argv[2]));
        else {
            cout << "Missing correct arguments" << endl;
            return 1;
        }
    }
    else if (argc == 5) {

        if      (string(argv[1]) == "Ackley") leProb = new Ackley(to<unsigned int>(string(argv[2])), to<double>(string(argv[3])));
        else if (string(argv[1]) == "Rastrigin") leProb = new Rastrigin(to<unsigned int>(string(argv[2])), to<double>(string(argv[3])));
        else if (string(argv[1]) == "Rosenbrock") leProb = new Rosenbrock(to<unsigned int>(string(argv[2])), to<double>(string(argv[3])));
        else if (string(argv[1]) == "Quadratic") leProb = new Quadratic(to<unsigned int>(string(argv[2])), to<double>(string(argv[3])));
        else {
            cout << "Missing correct arguments" << endl;
            return 1;
        }
        nEvalLimit = to<unsigned int>(string(argv[4]));
    }
    else {
        cout << "Missing correct arguments" << endl;
        return 1;
    }

    cout << "\nWelcome to LESOptera Version 2.0\n\n" << "Selected problem: " << typeid(*leProb).name() << endl;

    // store the results
    double bestFitness = -numeric_limits<double>::max();
    double bestFitnessPSO = -numeric_limits<double>::max();
    vector<double> bestAlleles;

    //double terminationFitness = -1e-6; //-1e-6;

//while ( bestFitness < terminationFitness ) {
//while ( leProb->getnEval() < nEvalLimit && bestFitness < -0.1 ) {

    try{

        while ( leProb->getnEval() < nEvalLimit ) {

            PSO populationPSO(leProb, np, c1, c2, lambdaPSO, kappa);

            while ( leProb->getnEval() < nEvalLimit ) {

                    populationPSO.step();

                    // write out the best fitness of this run
                    statRun.str("");
                    statRun << leProb->getnEval() << "\t" << populationPSO.getBestIndividualFitness() << endl;
                    writeResults(fileStatRunPSO.str(), statRun);

                    // update bestFitness
                    bestFitness = populationPSO.getBestIndividualFitness();
                    bestAlleles = populationPSO.getBestIndividual().getAlleles();

            }

            // save best fitness of PSO
            bestFitnessPSO = populationPSO.getBestIndividualFitness();

            Hooke_Jeeves populationHJ(leProb, (Individual&) (populationPSO.getBestIndividual()), r, s0, t, m);

            while ( populationHJ.mExceeded() == false ) {


                    populationHJ.step(); // inside this loop!

                    // write out the best fitness of this run
                    statRun.str("");
                    statRun << leProb->getnEval() << "\t" << populationHJ.getBestIndividualFitness() << endl;
                    writeResults(fileStatRunHJ.str(), statRun);


                    // update bestFitness
                    bestFitness = populationHJ.getBestIndividualFitness();
                    bestAlleles = populationHJ.getBestIndividual().getAlleles();

            }

        }
    }
    catch (string msg) {

        cout << "Exception caught:" << endl;
        cout << msg << endl;
        return 1;

    }


    time(&end);
    double simulTime = difftime(end,start);

    cout << "\nTotal number of evaluations: " << leProb->getnEval();
    cout << "\nBest fitness: " << bestFitness;
    cout << "\nBest alleles:\n[";
    for (unsigned int i=0; i<bestAlleles.size()-1;i++) cout << bestAlleles[i] << ", ";
    cout << bestAlleles[bestAlleles.size()-1] << "]";
    cout << "\nSimulation time: " << simulTime;
    cout << "\nEvaluation time: " << leProb->getTimeEval();
    cout << "\nAcceleration factor: " << (leProb->getTimeEval()/simulTime) << endl;

    ostringstream ss;
    ss << leProb->getnEval() << "\t" << bestFitnessPSO << "\t" << bestFitness << "\t" << norm(bestAlleles) << "\t";
    for (unsigned int i=0; i<bestAlleles.size(); i++) ss << bestAlleles[i] << "\t";
    ss << endl;

    if ( loopi == 0 ) writeResultsOver("Stats.dat", ss);
    else writeResults("Stats.dat", ss);

    delete leProb;

}

    return 0;

}
