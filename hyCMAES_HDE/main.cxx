/*
      Evolutionary Algorithms
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
#include <typeinfo> // to use typeid

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
    zigset(30062009);

    // DE parameters
    const double Cr = 0.1;
    const double  F = 0.3;

for (int loopi=0; loopi < 1; loopi++) {

    // stores the performance in file
    ostringstream statRun, fileStatRunCMA_ES, fileStatRunDE;
    fileStatRunCMA_ES << "StatRunCMA_ES" << loopi << ".dat";
    fileStatRunDE << "StatRunDE" << loopi << ".dat";

    // CSA_ES parameters
    double sigmaF = 0.2;
    vector< vector<double> > C;

    // Limiting parameters
    unsigned int nEvalLimit = 3000;
    const int nt = 10;

    // choose the problem

    Problem *leProb;

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
    else if (argc == 3) {

        if (string(argv[1]) == "EnergyPlus") leProb = new EnergyPlus(string(argv[2]));
        else if (string(argv[1]) == "CitySim") leProb = new CitySim(string(argv[2]));
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
		nEvalLimit = to<unsigned int>(string(argv[3]));
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

    int lambda = int(4.0 + floor(3.0*log(double(leProb->getSize()))));
    int mu = lambda/2;
    int NP = 30;

    // store the results
    vector< Individual >  popDE, popCSA_ES;
    double bestFitness  = -numeric_limits<double>::max();
    unsigned int bestId = 0;
    vector<double> bestAlleles;

    //double terminationFitness = -1e-6; //-1e-6;

//while ( bestFitness < terminationFitness ) {
//while ( leProb->getnEval() < nEvalLimit && bestFitness < -0.1 ) {

try{
        while ( leProb->getnEval() < nEvalLimit ) {

            CMA_ES *pPopulationCSA_ES = new CMA_ES (leProb, lambda, sigmaF, C, popCSA_ES);

            // output to the screen

            for (int i=0; i < nt; i++) {

                pPopulationCSA_ES->step();

                // write out the best fitness of this run
                statRun.str("");
                statRun << leProb->getnEval() << "\t";
				for (size_t nFitness=0;nFitness<pPopulationCSA_ES->getCurrentBestIndividual().getFitnessSize();++nFitness)
					statRun << pPopulationCSA_ES->getCurrentBestIndividual().getFitness(nFitness) << "\t";
				for (size_t nAlleles=0;nAlleles<pPopulationCSA_ES->getCurrentBestIndividual().getAlleles().size();++nAlleles)
					statRun << pPopulationCSA_ES->getCurrentBestIndividual().getAllele(nAlleles) << "\t";
				statRun << endl;
                writeResults(fileStatRunCMA_ES.str(), statRun);

                // update bestFitness
                if ( pPopulationCSA_ES->getBestIndividualFitness() > bestFitness ) {
                    bestFitness = pPopulationCSA_ES->getBestIndividualFitness();
                    bestAlleles = pPopulationCSA_ES->getBestIndividual().getAlleles();
                    bestId      = pPopulationCSA_ES->getBestIndividual().getIdNumber();
                }

                if ( leProb->getnEval() >= nEvalLimit ) break;
                //if (pPopulationCSA_ES->getBestIndividualFitness() > terminationFitness) break;

            }

            // clearance of popCSA_ES
            popCSA_ES.clear();

            if ( leProb->getnEval() >= nEvalLimit ) {
            //if (pPopulationCSA_ES->getBestIndividualFitness() > terminationFitness) {

                delete pPopulationCSA_ES;
                break;
            }

            pPopulationCSA_ES->getBestInPopulation(popDE);
        //    pPopulationCSA_ES->getPopulation(popCSA_ES);

            // adaptation of sigmaF and C
            sigmaF = pPopulationCSA_ES->getSigmaF();
            C      = pPopulationCSA_ES->getC();

            // start of popEA - DE
            DE *pPopulationDE = new DE (leProb, NP, popDE);

            for (int i=0; i < nt; i++) {

                pPopulationDE->step(Cr, F, 1, true);

                // write out the best fitness of this run
                statRun.str("");
                statRun << leProb->getnEval() << "\t";
				for (size_t nFitness=0;nFitness<pPopulationDE->getCurrentBestIndividual().getFitnessSize();++nFitness)
					statRun << pPopulationDE->getCurrentBestIndividual().getFitness(nFitness) << "\t";
				for (size_t nAlleles=0;nAlleles<pPopulationDE->getCurrentBestIndividual().getAlleles().size();++nAlleles)
					statRun << pPopulationDE->getCurrentBestIndividual().getAllele(nAlleles) << "\t";
				statRun << endl;
                writeResults(fileStatRunDE.str(), statRun);

                if ( leProb->getnEval() >= nEvalLimit ) break;
                //if (pPopulationDE->getBestIndividualFitness() > terminationFitness) break;

            }

            // has the DE migrated ? reinitialize sigmaF
            // if ( pPopulationDE->getMigrated() == true ) sigmaF = 0.2;

            popDE.clear(); //new extended population

        //    pPopulationDE->getPopulation(popDE);

            pPopulationDE->getPopulation( NP-nt, popDE );
            pPopulationDE->getPopulation( mu, popCSA_ES  );

            // output the sizes screen
            //cout << "Size of popCSA_ES: " << popCSA_ES.size() << "\tSize of popEA: " << popDE.size() << "\n" << endl;

            // update bestFitness
            if ( pPopulationDE->getBestIndividualFitness() > bestFitness ) {
                bestFitness = pPopulationDE->getBestIndividualFitness();
                bestAlleles = pPopulationDE->getBestIndividual().getAlleles();
                bestId      = pPopulationDE->getBestIndividual().getIdNumber();
            }

            // deletion of population
            delete pPopulationDE;
            delete pPopulationCSA_ES;

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
    cout << "\nBest ID: " << bestId;
    cout << "\nBest alleles:\n[";
    for (unsigned int i=0; i<bestAlleles.size()-1;i++) cout << bestAlleles[i] << ", ";
    cout << bestAlleles[bestAlleles.size()-1] << "]";
    cout << "\nSimulation time: " << simulTime;
    cout << "\nEvaluation time: " << leProb->getTimeEval();
    cout << "\nAcceleration factor: " << (leProb->getTimeEval()/simulTime) << endl;

    ostringstream ss;
    ss << leProb->getnEval() << "\t" << bestFitness << "\t" << norm(bestAlleles) << "\t";
    for (unsigned int i=0; i<bestAlleles.size(); i++) ss << bestAlleles[i] << "\t";
    ss << endl;

    if ( loopi == 0 ) writeResultsOver("Stats.dat", ss);
    else writeResults("Stats.dat", ss);

    delete leProb;

}

    return 0;

}
