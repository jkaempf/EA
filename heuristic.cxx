/*Algorithm definitions in all generality
  Jerome Kaempf
  LESO-PB / EPFL
  jerome.kaempf@epfl.ch

  Base class

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
#include <algorithm>
#include <map>

#include "geometry.h"
#include "heuristic.h"


// *****************************************************************
// Class Individual
// *****************************************************************

Individual::Individual(Problem *problem, vector<double> alleles)
{

    stringstream ss;

    this->problem = problem;
    simulated=false;              // not yet fitness evaluated
    fitness.assign(1,-numeric_limits<double>::max());

    this->alleles=alleles;        // initialisation of the vector alleles

    compteur++;                   // incrémentation du nombre d'individus crées
    number = compteur;
    ss << "Individual " << compteur;

    // mise en oeuvre dans le rapport
    ss << ": objet créé" << endl;
//    writeResults("report.dat", ss);

}

Individual::Individual(Problem *problem, vector<double> alleles, double fitness)
{

    stringstream ss;

    this->problem = problem;
    compteur++;                   // incrémentation du nombre d'individus crées
    number = compteur;

    simulated=true;             // fitness already evaluated
    this->fitness.assign(1,fitness);
    this->alleles=alleles;        // initialisation of the vector alleles

    // mise en oeuvre dans le rapport
    ss << "Individual " << number << ": objet créé" << endl;
//    writeResults("report.dat", ss);

}

Individual::Individual(Problem *problem, vector<double> alleles, int number, double fitness)
{

    stringstream ss;

    this->problem = problem;
    this->number = number;
    compteur = number;

    simulated=true; // fitness already evaluated
    this->fitness.assign(1,fitness);
    this->alleles=alleles;        // initialisation of the vector alleles

    // mise en oeuvre dans le rapport
    ss << "Individual " << number << ": objet créé" << endl;
//    writeResults("report.dat", ss);

}

ostream& operator<<(ostream &os, Individual obj) {

    os << "[";
    for (unsigned int i=0; i < obj.getAllele().size()-1; i++) { os << obj.getAllele(i) << ", "; }
    os << obj.getAllele(obj.getAllele().size()-1) << "]";

    return os;

}

// *****************************************************************
// Comparison of Individuals
// *****************************************************************

bool operator<(Individual &a, Individual &b) { // Define less than relative to T objects

	// test if the first individual (Pareto) dominates the other solution
	bool dominatesAll = true, dominatesOne = false;
	for (size_t i=0;i<a.getFitnessSize();++i) {
		dominatesAll &= ( a.getFitness(i) >= b.getFitness(i) );
		dominatesOne |= ( a.getFitness(i) > b.getFitness(i) );
	}
	if (dominatesAll && dominatesOne) return true; // a gets a better rank than b (<) if fitness of a is bigger than fitness b (>)
	else if ( a == b && a.getProblem() == b.getProblem() ) { // fitness égale et le même problème!!
        if ( a.getProblem()->notInField(a.getAllele()) && b.getProblem()->notInField(b.getAllele()) ) { // les deux sont infaisables

            //cout << "Les 2 sont mauvais!!" << endl;
            bool dominates = true;
            for (unsigned int k=0;k<a.getProblem()->getnConstraints();k++) {

                //cout << "g" << k+1 << "(premier): " << max(a.getProblem()->constraint(k,a.getAllele()),0.0) << "\tg" << k+1 << "(deuxieme): " << max(b.getProblem()->constraint(k,b.getAllele()), 0.0) << endl;
                dominates &= ( max(a.getProblem()->constraint(k,a.getAllele()),0.0) <= max(b.getProblem()->constraint(k,b.getAllele()), 0.0) );

            }
            //cout << a << "\ndomines: " << dominates << "\n" << b << endl;
            //wait(5);
            return dominates;

        }
        else { // egaux et dans la cible les deux
            return false;
        }
    }
    else return false; // c'est b le meilleur

}

bool operator==(Individual &a, Individual &b) { // defines perfect equality between individuals
	
	bool equal = true;
	for (size_t i=0;i<a.getFitnessSize();++i) {
		equal &= ( a.getFitness(i) == b.getFitness(i) );
	}
	return equal;
	
}

// *****************************************************************
// Class derived from Heuristic: DE
// *****************************************************************

DE::DE(Problem *problem, unsigned int size):Heuristic<Individual>(problem, size) {

    vector<double> alleles;

    for (unsigned int i=0;i<size;i++) {

        alleles.clear();

        for (unsigned int j=0 ; j < problem->getSize() ; j++) {

            alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

        }

	    pop.push_back( Individual(problem, alleles) );

    }

    migrated = false;

}

DE::DE(Problem *problem, unsigned int size, vector< Individual > &popIndividualsVector):Heuristic<Individual>(problem, size) {

    if ( popIndividualsVector.size() < size ) { // 1er lancement

        // mise en place des existants

        for (unsigned int i=0; i<popIndividualsVector.size(); i++) {

            pop.push_back( popIndividualsVector[i] );

        }

        // création de nouveaux pour compléter à NP

        vector<double> alleles;
        unsigned int nCompl = size - popIndividualsVector.size(); // nombre à compléter

        for (unsigned i=0; i < nCompl; i++) {

            alleles.clear();

            for (unsigned j=0 ; j < problem->getSize() ; j++) {

                alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

            }

            pop.push_back( Individual(problem, alleles) );

        }
    }
    else { // déjà lancé 1x

        for (unsigned int i=0; i<size; i++) {

            pop.push_back( popIndividualsVector[i] );

        }
    }

    migrated = false;

}

DE::DE(Problem *problem, string filename):Heuristic<Individual>(problem) {

    string tampon, file;
    unsigned int alleleSize;
    double fitness;
    vector<double> alleles;

    // chargement du fichier source
    ifstream input (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture
    if (!input.is_open())
      {
    cerr << "Error opening: " << filename << endl;
      }

    input >> tampon;

    setGeneration( atoi(tampon.c_str()) );

    cout << "generation: " << getGeneration() << endl;

    input >> tampon;

    setInitSize( atoi(tampon.c_str()) );

    cout << "popsize: " << atoi(tampon.c_str()) << endl;

    input >> tampon;

    alleleSize = atoi(tampon.c_str());

    cout << "alleleSize: " << alleleSize << endl;

    for (unsigned int i=0; i < getInitSize(); i++) {

        alleles.clear();

        input >> file;

        input >> tampon;
        fitness = atof(tampon.c_str());

        for (unsigned int j=0; j < alleleSize; j++) {

            input >> tampon;
            alleles.push_back(atof(tampon.c_str()));

        }

        pop.push_back( Individual(problem, alleles, atoi( (filename.substr(file.size()-1)).c_str() ), fitness) );

    }

    input.close();

    migrated = false;

}

void DE::crossover(double probability, double scale, int methodCrossover) {

    if (methodCrossover == 0) methodCrossover = (int) randomUniform(1.0, 3.0);

    for (unsigned i=0; i < getInitSize(); i++) {

        if (methodCrossover == 1) rand3(probability, scale, i);
        else rand2Best(probability, scale, i);

    }

    return;

}

void DE::migration() {

    // defines the position of the best element in the current pop
    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));

    double epsilon1 = 0.1;
    double alleleprime = 0.0;

    cout << "Population diversity: " << diversityPopulation() << endl;

    if ( diversityPopulation() < epsilon1 ) {

        // start of migration

        migrated = true;

        for (unsigned int i=0; i < getInitSize(); i++) {

//        cout << "Individual: " << i << "\tDiversity: " << diversityIndividual(i) << endl;

            if ( i!=posMax ) { // on ne change pas le meilleur...

                // boucle sur les alleles

                cout << "migration of individual: " << i << "\ttoo close from the best one: " << posMax << endl;

                for (unsigned int h=0; h<pop[i].getAllelesNumber(); h++) {

    //                    cout << "in allele: " << h << "\tvalue: " << pop[i].getAllele(h) << "\tbestValue: " << bestIndividual->getAllele(h) << endl;

                        // test si on met entre min ou max
                        if ( randomUniform(0.0, 1.0) < ( (pop[i].getAllele(h) - problem->getMinVector()[h])
                                                        / (problem->getMaxVector()[h] - problem->getMinVector()[h]) ) ) {

                            alleleprime = pop[i].getAllele(h) + randomUniform(0.0, 1.0)*(problem->getMinVector()[h]
                                                                                        - pop[posMax].getAllele(h));
                            pop[i].setAllele(h, alleleprime);

                        }
                        else {

                            alleleprime = pop[i].getAllele(h) + randomUniform(0.0, 1.0)*(problem->getMaxVector()[h]
                                                                                        - pop[posMax].getAllele(h));
                            pop[i].setAllele(h, alleleprime);

                        } // end of if/else

    //                    cout << "new allele: " << pop[i].getAllele(h) << endl;

                } // end of for loop on the alleles

                // tell him that things have changed...
                pop[i].setSimulated(false);

            } // end of if

        } // end of loop on the new population

    } // enf of if

    return;

}

bool DE::getMigrated() {

    return migrated;

}

void DE::rand2Best(double probability, double scale, int candidate) {

    // evalue l'ancienne génération
    evaluate();
    // regarde le meilleur de l'ancienne génération
    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));

    stringstream ss;

	int r1, r2;
	int n;
	int nDim = pop[candidate].getAllelesNumber();

	selectSamples(candidate,&r1,&r2);
	n = (int)randomUniform(0.0,(double)nDim);

    vector<double> trialSolution;
    trialSolution = pop[candidate].getAllele();

    ss << "candidat: " << candidate << "\tr1: " << r1 << "\tr2: " << r2 << endl;

	for (int i=0; i < nDim; i++)
	{

        // change the value of allele n, and the others with probability

	    if ( (randomUniform(0.0,1.0) <= probability) || i==n ) {

            trialSolution[i] = pop[posMax].getAllele(i)
							+ scale * ( pop[r2].getAllele(i)
							- pop[r1].getAllele(i) );
	    }

        // IF outside the limits, then put back randomly within the limits

        if ( (trialSolution[i] > problem->getMaxVector()[i]) || (trialSolution[i] < problem->getMinVector()[i]) ) {

            trialSolution[i] = randomUniform(problem->getMinVector()[i],problem->getMaxVector()[i]);

        }

	}

    pop.push_back(Individual(problem,trialSolution)); // creation du nouvel individu

//    writeResults("report.dat", ss);

    return;

}

void DE::rand3(double probability, double scale, int candidate) {

    stringstream ss;

	int r1, r2, r3;
	int n;
	int nDim = pop[candidate].getAllelesNumber();

	selectSamples(candidate,&r1,&r2,&r3);
	n = (int)randomUniform(0.0,(double)nDim);

    vector<double> trialSolution;
    trialSolution = pop[candidate].getAllele();

    ss << "candidat: " << candidate << "\tr1: " << r1 << "\tr2: " << r2 << "\tr3: " << r3 << endl;

	for (int i=0; i < nDim; i++)
	{

        // change the value of allele n, and the others (i) with probability

	    if ( (randomUniform(0.0,1.0) <= probability) || i==n ) {

            trialSolution[i] = pop[r3].getAllele(i)
							+ scale * ( pop[r1].getAllele(i)
							- pop[r2].getAllele(i) );
	    }

        // IF outside the limits, then put back randomly within the limits - or at the limit

        if ( (trialSolution[i] > problem->getMaxVector()[i]) || (trialSolution[i] < problem->getMinVector()[i]) ) {

            trialSolution[i] = randomUniform(problem->getMinVector()[i],problem->getMaxVector()[i]);
            //trialSolution[i] = min( trialSolution[i], leProb.getMaxVector()[i] );
            //trialSolution[i] = max( trialSolution[i], leProb.getMinVector()[i] );

        }

	}

    pop.push_back( Individual(problem,trialSolution) ); // creation du nouvel individu

//    writeResults("report.dat", ss);

    return;

}

void DE::selectSamples(int candidate,int *r1,int *r2,int *r3,int *r4,int *r5) {


	if (r1)
	{
		do
		{
			*r1 = (int)randomUniform(0.0,(double)getPopSize());
		}
		while (*r1 == candidate);
	}

	if (r2)
	{
		do
		{
			*r2 = (int)randomUniform(0.0,(double)getPopSize());
		}
		while ((*r2 == candidate) || (*r2 == *r1));
	}

	if (r3)
	{
		do
		{
			*r3 = (int)randomUniform(0.0,(double)getPopSize());
		}
		while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
	}

	if (r4)
	{
		do
		{
			*r4 = (int)randomUniform(0.0,(double)getPopSize());
		}
		while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
	}

	if (r5)
	{
		do
		{
			*r5 = (int)randomUniform(0.0,(double)getPopSize());
		}
		while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
													|| (*r5 == *r2) || (*r5 == *r1));
	}

	return;
}

void DE::selection() {

    stringstream ss;

    vector<int> deleteTrial, deleteTarget;

    for (unsigned int i=0; i< getInitSize(); i++) {

        if ( pop[i].getFitness() > pop[i+getInitSize()].getFitness() ) deleteTrial.push_back(i+getInitSize());
        else deleteTarget.push_back(i);

    }

    for (int i=(deleteTrial.size()-1); i >= 0 ; i--) {

        ss << "suppression de l'individu: " << pop[deleteTrial[i]].getIdNumber() << endl;
        pop.erase(pop.begin() + deleteTrial[i]); // suppression de l'individu moins bon

    }

    for (int i=(deleteTarget.size()-1); i >= 0 ; i--) {

        ss << "suppression de l'individu: " << pop[deleteTarget[i]].getIdNumber() << endl;
        pop.erase(pop.begin() + deleteTarget[i]); // suppression de l'individu moins bon

    }

//    writeResults("report.dat", ss);

    return;

}

void DE::selectionRank() {

    stringstream ss;

    vector<int> deleteTrial, deleteTarget;

    for (unsigned int i=0; i< getInitSize(); i++) {

        if ( pop[i] < pop[i+getInitSize()] ) deleteTrial.push_back(i+getInitSize());
        else deleteTarget.push_back(i);

    }

    for (int i=(deleteTrial.size()-1); i >= 0 ; i--) {

        ss << "suppression de l'individu: " << pop[deleteTrial[i]].getIdNumber() << endl;
        pop.erase(pop.begin() + deleteTrial[i]); // suppression de l'individu moins bon

    }

    for (int i=(deleteTarget.size()-1); i >= 0 ; i--) {

        ss << "suppression de l'individu: " << pop[deleteTarget[i]].getIdNumber() << endl;
        pop.erase(pop.begin() + deleteTarget[i]); // suppression de l'individu moins bon

    }

//    writeResults("report.dat", ss);

    return;

}

void DE::step(double Cr, double F, int methodCrossover, bool migrate) {

    crossover(Cr, F, methodCrossover);
    evaluate();
    selectionRank();
    statistics();

//  **** Evaluation of the diversities ****

    save("currentpopDE_" + toString(problem->getnEval()) + ".dat");
    ostringstream ss;
    ss << "bestIndividualID: " << getCurrentBestIndividual().getIdNumber();
    ss << "\tpopulationDiversity: " << diversityPopulation() << endl;
    ss << "diversities: ";
    for (unsigned int j=0; j<getInitSize(); j++) {
        ss << diversityIndividual(j) << "\t";
    }
    ss << endl;
    writeResults("currentpopDE_" + toString(problem->getnEval()) + ".dat", ss);


//  **** Migration on the top of that ****
    if (migrate == true) {
        migration();
        evaluate();
    }

    incrementGeneration();

    return;

}


// *****************************************************************
// Class derived from Heuristic: ES
// *****************************************************************


ES::ES(Problem *problem, unsigned int mu, unsigned int lambda, double sigmaF):Heuristic<IndividualES>(problem,mu) {

    this->mu = mu;
    this->lambda = lambda;

    this->sigmaF = sigmaF;

    vector<double> alleles, sigmas;

    for (unsigned int i=0;i<mu;i++) {

        alleles.clear();
        sigmas.clear();

        for (unsigned int j=0 ; j < problem->getSize() ; j++) {

            alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max
            sigmas.push_back ( sigmaF * (problem->getMaxVector()[j]-problem->getMinVector()[j]) );

        }

	    pop.push_back( IndividualES(problem, alleles, sigmas) );

    }

}

ES::ES(Problem *problem, unsigned int mu, unsigned int lambda, double sigmaF, vector< vector<double> > &allelesVector):Heuristic<IndividualES>(problem,mu) {

    this->mu = mu;
    this->lambda = lambda;

    this->sigmaF = sigmaF;

    if (allelesVector.size() < mu) { // 1er lancer

        vector<double> alleles, sigmas;

        for (unsigned int j=0 ; j < problem->getSize() ; j++) {

            sigmas.push_back ( sigmaF * (problem->getMaxVector()[j]-problem->getMinVector()[j]) );

        }

        for (unsigned int i=0; i<allelesVector.size(); i++) {

            pop.push_back( IndividualES( problem,allelesVector[i], sigmas ) );

        }

        unsigned int nCompl = mu - allelesVector.size();

        for (unsigned int i=0;i<nCompl;i++) {

            alleles.clear();

            for (unsigned int j=0 ; j < problem->getSize() ; j++) {

                alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

            }

            pop.push_back( IndividualES(problem,alleles, sigmas) );

        }

    }
    else { // déjà lancé 1x

        vector<double> sigmas;

        for (unsigned int j=0 ; j < problem->getSize() ; j++) {

            sigmas.push_back ( sigmaF * (problem->getMaxVector()[j]-problem->getMinVector()[j]) );

        }

        for (unsigned int i=0; i<mu; i++) {

            pop.push_back( IndividualES(problem,allelesVector[i], sigmas) );

        }

    }

}

ES::ES(Problem *problem, string filename, unsigned int lambda):Heuristic<IndividualES>(problem) {

    this->lambda = lambda;

    string tampon, file;
    unsigned int alleleSize;
    double fitness;
    vector<double> alleles, sigmas;

    // chargement du fichier source
    ifstream input (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture
    if (!input.is_open())
      {
    cerr << "Error opening: " << filename << endl;
      }

    input >> tampon;

    setGeneration( atoi(tampon.c_str()) );

    cout << "generation: " << getGeneration() << endl;

    input >> tampon;

    setInitSize( atoi(tampon.c_str()) );
    mu = atoi(tampon.c_str());

    cout << "popsize: " << atoi(tampon.c_str()) << endl;

    input >> tampon;

    alleleSize = atoi(tampon.c_str());

    cout << "alleleSize: " << alleleSize << endl;

    for (unsigned int i=0; i < getInitSize(); i++) {

        alleles.clear();
        sigmas.clear();

        input >> file;

        input >> tampon;
        fitness = atof(tampon.c_str());

        for (unsigned int j=0; j < alleleSize; j++) {

            input >> tampon;
            alleles.push_back(atof(tampon.c_str()));

            input >> tampon;
            sigmas.push_back(atof(tampon.c_str()));

        }

        pop.push_back( IndividualES(problem, alleles, sigmas, atoi( (filename.substr(file.size()-1)).c_str() ), fitness) );

    }

    input.close();

}

vector<unsigned int> ES::selectParents(unsigned int rho) {

    unsigned int candidate;
    vector<unsigned int> parents;
    bool bad;

    for (unsigned int i=0; i < rho; i++) {

        do {

            bad = false;
            candidate = (unsigned int)randomUniform(0.0,(double)getInitSize());

            for (unsigned int j=0; j < parents.size(); j++) {

                if ( candidate == parents[j] ) bad |= true;

            }

        }
        while ( bad == true );

        parents.push_back(candidate);

    }

    return parents;

}

void ES::globalIntermediateRecombination(double probability, unsigned int rho) {

    vector<unsigned int> parents;

    vector<double> allelesFils;
    vector<double> sigmasFils;

    stringstream ss;

    // taille de la population initiale: popInitSize()

    while ( getPopSize() < (getInitSize()+lambda) ) {        // on cherche à faire des enfants

        parents = selectParents(rho);

        if ( randomUniform(0.0, 1.0) <= probability ) {  // probabilite de croisement

            for (unsigned int i=0; i<parents.size(); i++) {

                ss << "\tparent" << i << ": " << parents[i];

            }

            ss << endl;

            unsigned int allelesNumber = pop[parents[0]].getAllelesNumber();

            allelesFils.clear();
            sigmasFils.clear();

            for (unsigned int i=0;i<allelesNumber;i++) {

                double allele=0.0, sigma=0.0;

                for (unsigned int j=0; j < parents.size(); j++ ) {

                    allele += pop[parents[j]].getAllele(i);
                    sigma  += pop[parents[j]].getSigma(i);

                }

                allelesFils.push_back( allele / (double) parents.size() );
                sigmasFils.push_back(  sigma  / (double) parents.size() );

            }

            pop.push_back(IndividualES(problem,allelesFils, sigmasFils)); // creation du nouvel individu

        }
        else { // pas de croisement (passage du parent1 à la mutation)

            pop.push_back( IndividualES(problem, pop[parents[0]].getAllele(), pop[parents[0]].getSigma()) );

        }
    }

//    writeResults("report.dat", ss);

}

void ES::intermediateRecombination(double probability) {

    vector<unsigned int> parents;
    double weight=0.5;
    vector<double> allelesFils;
    vector<double> sigmasFils;

    stringstream ss;

    // taille de la population initiale: popSize

    while ( getPopSize() < (getInitSize()+lambda) ) {        // on cherche à faire des enfants

        parents = selectParents(2);

        if ( randomUniform(0.0, 1.0) <= probability ) {  // probabilite de croisement

            ss << "parent1: " << parents[0] << "\tparent2: " << parents[1] << endl;

            int allelesNumber = pop[parents[0]].getAllelesNumber();

            allelesFils.clear();
            sigmasFils.clear();

            for (int i=0;i<allelesNumber;i++) {

                allelesFils.push_back(weight*pop[parents[0]].getAllele(i)+(1-weight)*pop[parents[1]].getAllele(i));
                sigmasFils.push_back(weight*pop[parents[0]].getSigma(i)+(1-weight)*pop[parents[1]].getSigma(i));

            }

            pop.push_back(IndividualES(problem,allelesFils, sigmasFils)); // creation du nouvel individu

        }
        else { // pas de croisement (passage du parent1 à la mutation)

            pop.push_back( IndividualES(problem,pop[parents[0]].getAllele(), pop[parents[0]].getSigma()) );

        }
    }

//    writeResults("report.dat", ss);

}

void ES::generalizedIntermediateRecombination(double probability) {

    vector<unsigned int> parents;
    double weight;
    vector<double> allelesFils;
    vector<double> sigmasFils;

    stringstream ss;

    // taille de la population initiale: popSize

    while ( getPopSize() < (getInitSize()+lambda) ) {        // on cherche à faire des enfants

        parents = selectParents(2);

        if ( randomUniform(0.0, 1.0) <= probability ) {  // probabilite de croisement

            ss << "parent1: " << parents[0] << "\tparent2: " << parents[1] << endl;

            int allelesNumber = pop[parents[0]].getAllelesNumber();

            allelesFils.clear();
            sigmasFils.clear();

            for (int i=0;i<allelesNumber;i++) {

                weight = randomUniform(0.0, 1.0);

                allelesFils.push_back(weight*pop[parents[0]].getAllele(i)+(1-weight)*pop[parents[1]].getAllele(i));
                sigmasFils.push_back(weight*pop[parents[0]].getSigma(i)+(1-weight)*pop[parents[1]].getSigma(i));

            }

            pop.push_back(IndividualES(problem,allelesFils, sigmasFils)); // creation du nouvel individu

        }
        else { // pas de croisement (passage du parent1 à la mutation)

            pop.push_back( IndividualES(problem, pop[parents[0]].getAllele(), pop[parents[0]].getSigma()) );

        }
    }

//    writeResults("report.dat", ss);

}

void ES::uniformCrossoverRecombination(double probability) {

    int parent1;

    vector<double> allelesFils;
    vector<double> sigmasFils;

    stringstream ss;

    // taille de la population initiale: popSize

    while ( getPopSize() < (getInitSize()+lambda) ) {        // on cherche à faire des enfants

        parent1 = (int)randomUniform(0.0, (double) getInitSize());

        if ( randomUniform(0.0, 1.0) <= probability ) {  // probabilite de croisement;

            int allelesNumber = pop[parent1].getAllelesNumber();

            allelesFils.clear();
            sigmasFils.clear();

            for (int i=0;i<allelesNumber;i++) {

                parent1 = (int)randomUniform(0.0, (double) getInitSize());

                allelesFils.push_back(pop[parent1].getAllele(i));
                sigmasFils.push_back(pop[parent1].getSigma(i));

            }

            pop.push_back(IndividualES(problem, allelesFils, sigmasFils)); // creation du nouvel individu

        }
        else { // pas de croisement (passage du parent1 à la mutation)

            pop.push_back( IndividualES(problem, pop[parent1].getAllele(), pop[parent1].getSigma()) );

        }
    }

//    writeResults("report.dat", ss);

}

void ES::mutationNormallyDistributed(double K) {

    stringstream ss;

    int n;
    double z0, zi, sigmaprime, sigmaMax, alleleprime;

    for (unsigned int i=getInitSize(); i < getInitSize()+lambda ;i++) {

        // initialisation du nombre d'alleles: n

        n = pop[i].getAllelesNumber();
        ss << "mutation de l'individu: " << i << endl;

        // partie commune à tous les allèles pour la mutation

        z0 = ( K / sqrt(2.0 * (double) n) ) * normallyDistributedSPRNG_Ziggurat();

        // boucle sur les alleles pour MUTATION

        for (int j=0; j < n ;j++) {

            ss << pop[i].getAllele()[j];

            // partie changeante pour chaque allele
            zi = ( K / sqrt(2.0*sqrt((double) n)) ) * normallyDistributedSPRNG_Ziggurat();

            // calcul du sigmaprime et définition dans le système
            sigmaprime = pop[i].getSigma(j)*exp( z0 + zi );

            // test si on dépasse sigmaMax (sigma0)
            sigmaMax = sigmaF * (problem->getMaxVector()[j]-problem->getMinVector()[j]);
            if ( sigmaprime > sigmaMax ) sigmaprime = sigmaMax;

            // sauvegarder le nouveau sigma
            pop[i].setSigma(j, sigmaprime);

            // calcul de l'allele'
            alleleprime = pop[i].getAllele()[j] + sigmaprime*normallyDistributedSPRNG_Ziggurat();

            // 1er check: dans les limites de la variable
            alleleprime = min( alleleprime, problem->getMaxVector()[j] - TINY );
            alleleprime = max( alleleprime, problem->getMinVector()[j] );
            pop[i].setAllele( j, alleleprime );

            ss << "/" << alleleprime << " ";

        }

        ss << endl;

        // 2ème check: dans les contraintes

        if ( problem->notInField(pop[i].getAllele()) ) {

            pop[i].setSimulated(true);  // don't simulate, out of field

        }
        else {

            pop[i].setSimulated(false);  // not yet simulated...

        }

    }

//    writeResults("report.dat", ss);

}

void ES::selectionPlus() {

    vector<double> fitnessVector;

    for (unsigned int i=0;i<getPopSize();i++) {

       fitnessVector.push_back( pop[i].getFitness() );

    }

    // lambda: offspring

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j<lambda;j++) {

        posMin = distance(fitnessVector.begin(), min_element(fitnessVector.begin(), fitnessVector.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon
        fitnessVector.erase( min_element(fitnessVector.begin(), fitnessVector.end()) );
    }

//    writeResults("report.dat", ss);

}

void ES::selectionComma() {

    vector<double> fitnessVector;

    for (unsigned int i=mu;i<getPopSize();i++) {

       fitnessVector.push_back( pop[i].getFitness() );

    }

    // mu: population size, lambda: offspring

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j < mu;j++) {     // suppression des parents

        ss << "suppression de l'individu" << pop[mu-1-j].getIdNumber() << endl;
        pop.erase(pop.begin() + mu-1-j); // suppression de l'individu

    }

    for (unsigned int j=0; j < (lambda-mu);j++) {   // selection des mu meilleurs parmi les lambda enfants restants
                                           // suppression des lambda-mu moins bons
        posMin = distance(fitnessVector.begin(), min_element(fitnessVector.begin(), fitnessVector.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon
        fitnessVector.erase( min_element(fitnessVector.begin(), fitnessVector.end()) );
    }

//    writeResults("report.dat", ss);

}

void ES::selectionPlusRank() {

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j<lambda;j++) {

        posMin = distance(pop.begin(), max_element(pop.begin(), pop.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon

    }

//    writeResults("report.dat", ss);

    return;

}

void ES::selectionCommaRank() {

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j < mu;j++) {     // suppression des parents

        ss << "suppression de l'individu" << pop[mu-1-j].getIdNumber() << endl;
        pop.erase(pop.begin() + mu-1-j); // suppression de l'individu

    }

    for (unsigned int j=0; j < (lambda-mu);j++) {   // selection des mu meilleurs parmi les lambda enfants restants
                                           // suppression des lambda-mu moins bons
        posMin = distance(pop.begin(), max_element(pop.begin(), pop.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon

    }

//    writeResults("report.dat", ss);

    return;

}

void ES::save(string filename) {

        ostringstream ss;

        ss << getGeneration() << "\t" << getPopSize() << "\t" << pop[0].getAllelesNumber() << endl;

        for (unsigned int i=0;i<getPopSize();i++) {

            ss << "individual" << pop[i].getIdNumber() << "\t";
			for (size_t nFitness=0;nFitness<pop[i].getFitnessSize();++nFitness)
				ss << pop[i].getFitness(nFitness) << "\t";

            ss << setprecision(12) << scientific;
            for (unsigned int j=0;j<pop[i].getAllelesNumber();j++) {

                ss << pop[i].getAllele(j) << "\t" << pop[i].getSigma(j) << "\t";

            }

            ss << endl;

        }

        writeResultsOver(filename, ss);

}

void ES::report() {

    stringstream ss1, ss2;

    ss1 << "\nRapport de génération n°: " << gen << endl;

    for (unsigned int i=0;i<getPopSize();i++) {

      ss1 << "Le membre n°: " << i << " (" << "individual" << pop[i].getIdNumber() << ")\tSa fitness: " << pop[i].getFitness() << endl;

      ss1 << "\talleles: " << "(";

      for (unsigned int j=0;j<pop[i].getAllelesNumber();j++) {

        ss1 << "[" << pop[i].getAllele(j) << "], ";

      }

      ss1 << "[end])" << endl;

      ss1 << "\tsigmas: " << "(";

      for (unsigned int j=0;j<pop[i].getSigmasNumber();j++) {

        ss1 << "[" << pop[i].getSigma(j) << "], ";

      }

      ss1 << "[end])" << endl;

    }

    ss1 << "\nLe meilleur : " << posMax << " (" << "building" << pop[posMax].getIdNumber() << ")\tSa fitness: " << leMax << "\nmax/min(moy): " << setprecision(3) << leMax << "/" << leMin << "(" << laMoy << ")" << endl;

    // ecriture dans le stringstream du rapport

    ss2 << gen << "\t" << laMoy << "\t" << leMax << "\t" << leMin << "\t" << "building" << pop[posMax].getIdNumber() << endl;

//    writeResults("report.dat", ss1);
//    writeResults("population.dat", ss2);

    cout << "Génération: " << gen << "\tmax/min(moy): " << leMax << "\t" << leMin << "\t" << laMoy << endl;

  }

void ES::step(double K, int methodRecombination, int methodSelection) {

    double pr = 1.0;

    if (methodRecombination == 0) methodRecombination = (int) randomUniform(1.0, 4.0);
    if (methodSelection == 0)     methodSelection     = (int) randomUniform(1.0, 3.0);

//      **** recombination ****
    if      (methodRecombination == 1) intermediateRecombination(pr);
    else if (methodRecombination == 2) generalizedIntermediateRecombination(pr);
    else if (methodRecombination == 3) uniformCrossoverRecombination(pr);
    else globalIntermediateRecombination(pr, mu);

//      **** mutation *********
    mutationNormallyDistributed(K);

//      **** evaluation *******
    evaluate();

//      **** selection ********
    if      (methodSelection == 1) selectionCommaRank();
    else if (methodSelection == 2) selectionPlusRank();
    statistics();

//  **** Evaluation of the diversities ****

    save("currentpopES_" + toString(problem->getnEval()) + ".dat");
    ostringstream ss;
    ss << "bestIndividualID: " << getCurrentBestIndividual().getIdNumber();
    ss << "\tpopulationDiversity: " << diversityPopulation() << endl;
    ss << "diversities: ";
    for (unsigned int j=0; j<getInitSize(); j++) {
        ss << diversityIndividual(j) << "\t";
    }
    ss << endl;
    writeResults("currentpopES_" + toString(problem->getnEval()) + ".dat", ss);

//      **** end of step ******
    incrementGeneration();

    return;

}

// *****************************************************************
// Class derived from ES: CSA-ES
// *****************************************************************

CSA_ES::CSA_ES(Problem *problem, unsigned int lambda, double sigmaF):Heuristic<IndividualCSA_ES>(problem,lambda/2) {

    this->lambda = lambda;
    this->mu     = lambda/2;
    this->sigmaF = sigmaF;

    cout << "mu: " << mu << "\tlambda: " << lambda << "\n" << "\tsigmaF: " << sigmaF << endl;

    N = (double) problem->getSize();

    chiN = sqrt(N)*(1.0 - 1.0/(4.0*N) + 1.0/(21.0*pow(N, 2.0)));

    cs = 1.0/sqrt(N);
    damp = sqrt(N);

    ps.assign( (int) N, 0.0);

    // individual creation

    vector<double> alleles;

    for (unsigned int i=0;i<mu;i++) {

        alleles.clear();

        for (unsigned int j=0 ; j < problem->getSize() ; j++) {

            alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

        }

	    pop.push_back( IndividualCSA_ES(problem,alleles) );

    }

}

CSA_ES::CSA_ES(Problem *problem, unsigned int lambda, double sigmaF, vector< vector<double> > &allelesVector):Heuristic<IndividualCSA_ES>(problem,lambda/2) {

    this->lambda = lambda;
    this->mu     = lambda/2;
    this->sigmaF = sigmaF;

    cout << "mu: " << mu << "\tlambda: " << lambda << "\n" << "\tsigmaF: " << sigmaF << endl;

    N = (double) problem->getSize();
    chiN = sqrt(N)*(1.0 - 1.0/(4.0*N) + 1.0/(21.0*pow(N, 2.0)));

    cs = 1.0/sqrt(N);
    damp = sqrt(N);

    ps.assign( (int) N, 0.0);

    // individual creation

    if (allelesVector.size() < mu) {

        for (unsigned int i=0; i<allelesVector.size(); i++) {

            pop.push_back( IndividualCSA_ES(problem,allelesVector[i]) );

        }

        unsigned int nCompl = mu - allelesVector.size();

        vector<double> alleles;

        for (unsigned int i=0; i<nCompl; i++) {

            alleles.clear();

            for (unsigned int j=0 ; j < problem->getSize() ; j++) {

                alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

            }

            pop.push_back( IndividualCSA_ES(problem,alleles) );

        }
    }
    else {

        for (unsigned int i=0; i<mu; i++) {

            pop.push_back( IndividualCSA_ES(problem,allelesVector[i]) );

        }

    }

}

vector<int> CSA_ES::selectParents(int rho) {

    int candidate;
    vector<int> parents;
    bool bad;

    for (int i=0; i < rho; i++) {

        do {

            bad = false;
            candidate = (int)randomUniform(0.0,(double)getInitSize());

            for (unsigned int j=0; j < parents.size(); j++) {

                if ( candidate == parents[j] ) bad |= true;

            }

        }
        while ( bad == true );

        parents.push_back(candidate);

    }

    return parents;

}

void CSA_ES::globalIntermediateRecombination() {

    vector<int> parents;
    vector<double> allelesFils;

    stringstream ss;

    // taille de la population initiale: popInitSize()

    parents = selectParents(mu);
    for (unsigned int i=0; i<parents.size(); i++) {

        ss << "\tparent" << i << ": " << parents[i];

    }
    ss << endl;

    int allelesNumber = problem->getSize();

    for (int i=0;i<allelesNumber;i++) {

        double allele=0.0;

        for (unsigned int j=0; j < parents.size(); j++ ) {

            allele += pop[parents[j]].getAllele(i);

        }

        allelesFils.push_back( allele / (double) parents.size() );

    }

    // creation de lambda nouveaux individus

    for (unsigned int i=0; i<lambda; i++) {

        pop.push_back(IndividualCSA_ES(problem,allelesFils));

    }

//    writeResults("report.dat", ss);

}

void CSA_ES::mutationNormallyDistributed() {

    stringstream ss;

    int n;
    double alleleprime;

    for (unsigned int i=mu; i < mu+lambda ;i++) {

        // initialisation du nombre d'alleles: n
        n = pop[i].getAllelesNumber();
        ss << "mutation de l'individu: " << i << endl;

        // boucle sur les alleles pour MUTATION
        for (int j=0; j < n ;j++) {

            ss << pop[i].getAllele()[j];

            // calcul de la variation dans l'allele
            pop[i].setArz(j, normallyDistributedSPRNG_Ziggurat());

            // calcul de l'allele'
            alleleprime = pop[i].getAllele(j)
                        + sigmaF*(problem->getMaxVector()[j]-problem->getMinVector()[j])*pop[i].getArz(j);

            // 1er check: dans les limites de la variable - sauvegarde de alleleprime
            alleleprime = min( alleleprime, problem->getMaxVector()[j] - TINY );
            alleleprime = max( alleleprime, problem->getMinVector()[j] );
            pop[i].setAllele( j, alleleprime );

            ss << "/" << alleleprime << " ";

        }

        ss << endl;

        // 2ème check: dans les contraintes

        if ( problem->notInField(pop[i].getAllele()) ) {

            pop[i].setSimulated(true);  // don't simulate, out of field

        }
        else {

            pop[i].setSimulated(false);  // not yet simulated...

        }

    }

//    writeResults("report.dat", ss);

}

void CSA_ES::selectionComma() {

    vector<double> fitnessVector;

    for (unsigned int i=mu;i<getPopSize();i++) {

       fitnessVector.push_back( pop[i].getFitness() );

    }

    // mu: population size, lambda: offspring

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j < mu;j++) {     // suppression des parents

        ss << "suppression de l'individu" << pop[mu-1-j].getIdNumber() << endl;
        pop.erase(pop.begin() + mu-1-j); // suppression de l'individu

    }

    for (unsigned int j=0; j < (lambda-mu);j++) {   // selection des mu meilleurs parmi les lambda enfants restants
                                           // suppression des lambda-mu moins bons
        posMin = distance(fitnessVector.begin(), min_element(fitnessVector.begin(), fitnessVector.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon
        fitnessVector.erase( min_element(fitnessVector.begin(), fitnessVector.end()) );
    }

    // mise en place du meilleur dans bestIndividual
    posMax = distance(fitnessVector.begin(), max_element(fitnessVector.begin(), fitnessVector.end()));
    bestIndividualsVector.push_back( pop[posMax] );

//    writeResults("report.dat", ss);

    return;

}

void CSA_ES::selectionCommaRank() {

    // mu: population size, lambda: offspring

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j < mu;j++) {     // suppression des parents

        ss << "suppression de l'individu" << pop[mu-1-j].getIdNumber() << endl;
        pop.erase(pop.begin() + mu-1-j); // suppression de l'individu

    }

    for (unsigned int j=0; j < (lambda-mu);j++) {   // selection des mu meilleurs parmi les lambda enfants restants
                                           // suppression des lambda-mu moins bons
        posMin = distance(pop.begin(), max_element(pop.begin(), pop.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + posMin); // suppression de l'individu moins bon

    }

    // mise en place du meilleur dans bestIndividual et ajout dans bestIndividalMM
    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));
    bestIndividualsVector.push_back( pop[posMax] );
    // à modifier pour la finale de RANK

//    writeResults("report.dat", ss);

    return;

}

void CSA_ES::adaptSigma() {

    stringstream ss;

    int allelesNumber = pop[0].getAllelesNumber();

    for (int i=0;i<allelesNumber;i++) {

        double zmeanwi=0.0;

        for (unsigned int j=0; j < mu; j++ ) {

            zmeanwi += pop[j].getArz(i);

        }

        zmeanwi /= (double) mu;

        // adaptation de sigma - ps
        ps[i] = (1-cs)*ps[i] + sqrt(mu*cs*(2.0-cs))*zmeanwi;

    }

    // adaptation de sigmaF
    sigmaF = sigmaF*exp( (norm(ps)/chiN-1)*(1.0/damp) );

    // sortie du nouveau sigmaF
    ss << "SigmaF: " << sigmaF << endl;

//    writeResults("report.dat", ss);

}

void CSA_ES::step() {

//      **** recombination ****
    globalIntermediateRecombination();

//      **** mutation *********
    mutationNormallyDistributed();

//      **** evaluation *******
    evaluate();
    //evaluateVirtual();

//      **** selection ********
    selectionCommaRank();
    statistics();

//  **** Evaluation of the diversities ****

    save("currentpopCSA_ES_" + toString(problem->getnEval()) + ".dat");
    ostringstream ss;
    ss << "bestIndividualID: " << getCurrentBestIndividual().getIdNumber();
    ss << "\tpopulationDiversity: " << diversityPopulation() << endl;
    ss << "diversities: ";
    for (unsigned int j=0; j<getInitSize(); j++) {
        ss << diversityIndividual(j) << "\t";
    }
    ss << endl;
    writeResults("currentpopCSA_ES_" + toString(problem->getnEval()) + ".dat", ss);

//      **** parameter adaptation
    adaptSigma();

//      **** end of step ******
    incrementGeneration();

    return;

}

void CSA_ES::getBestInPopulation(vector< Individual > &popIndividualsVector) {

    for (unsigned int i=0; i<bestIndividualsVector.size(); i++) {

        popIndividualsVector.push_back( (Individual) bestIndividualsVector[i] );

    }

    return;

}

double CSA_ES::getSigmaF() {

    return sigmaF;

}

// *****************************************************************
// Class derived from ES: CMA-ES
// *****************************************************************

CMA_ES::CMA_ES(Problem *problem, int lambda, double sigmaF):Heuristic<IndividualCSA_ES>(problem,lambda/2) {

    this->lambda = lambda;
    this->mu     = lambda/2;
    this->sigmaF = sigmaF;

    cout << "mu: " << mu << "\tlambda: " << lambda << "\n" << "\tsigmaF: " << sigmaF << endl;

    N = double(problem->getSize());

    chiN = sqrt(N)*(1.0 - 1.0/(4.0*N) + 1.0/(21.0*pow(N, 2.0)));

    ps.assign( (int) N, 0.0);
    pc.assign( (int) N, 0.0);

    // Compute weights and effective mu
    double sumWeight=0.0;
    weight.assign( mu, log(double(mu+1)) );

    for(unsigned int i=0; i<weight.size(); i++) {

        weight[i] -= log(double(i+1));
        sumWeight += weight[i];

    }

    muEff = 0.0;

    for(unsigned int i=0; i<weight.size(); i++) {

        weight[i] /= sumWeight;
        muEff += (weight[i] * weight[i]);

    }

    muEff = 1.0 / muEff;

    // csigma

    cs = (muEff+2.0)/(N+muEff+3.0);
    damp = 1.0 + 2.0*max( 0.0, ( sqrt( double(muEff - 1.0)/(N+1.0) ) - 1.0 ) ) + cs;
    cc = 4.0/(N+4.0);

    muCov = muEff;
    cCov = ( (1.0/muCov) * (2.0 / pow( (N+sqrt(2.0)), 2.0 ) ) )
           + ( (1.0 - (1.0/muCov) ) * min(1.0, ( ( 2.0*muEff - 1.0 ) / ( pow((N+2.0), 2.0) + muEff ) ) ) );

    // creation du vecteur deltas

    vector<double> deltas;

    for (unsigned int i=0; i<problem->getSize(); i++) {

        deltas.push_back ( problem->getMaxVector()[i]-problem->getMinVector()[i] );

    }

    // initialisation de D, B et C

    D.assign( (int) N, 0.0 );
    B.assign( (int) N, D );
    C.assign( (int) N, D );

    // creation du vecteur D

    for (unsigned int i=0; i<problem->getSize(); i++) {

        D[i] = deltas[i];

    }

    // création de la matrice B

    for (unsigned int i=0; i<problem->getSize(); i++) {
        for (unsigned int j=0; j<problem->getSize(); j++) {

            if (i==j) B[i][j]=1.0;
            else B[i][j]=0.0;
        }
    }

    // création de la matrice C

    for (unsigned int i=0; i<problem->getSize(); i++) {
        for (unsigned int j=0; j<problem->getSize(); j++) {

            for (unsigned int k=0; k<problem->getSize(); k++) C[i][j]+=B[i][k]*D[k]*D[k]*B[j][k];

        }
    }

    // individual creation

    vector<double> alleles;

    for (unsigned int i=0; i<mu; i++) {

        alleles.clear();

        for (unsigned int j=0; j<problem->getSize(); j++) {

            alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

        }

	    pop.push_back( IndividualCSA_ES(problem,alleles) );

    }

    // évaluation de la population
    evaluate();
    // classement de la population
    sort( pop.begin(), pop.end() ); // du meilleur au moins bon
    // save intial random population
    save("currentpopCMA_ES_" + toString(problem->getnEval()) + ".dat");

}

CMA_ES::CMA_ES(Problem *problem, int lambda, double sigmaF, vector< vector<double> > &C0, vector< Individual > &popIndividualsVector):Heuristic<IndividualCSA_ES>(problem,lambda/2) {

    this->lambda = lambda;
    this->mu     = lambda/2;
    this->sigmaF = sigmaF;

    cout << "mu: " << mu << "\tlambda: " << lambda << "\n" << "\tsigmaF: " << sigmaF << endl;

    N = double(problem->getSize());

    chiN = sqrt(N)*(1.0 - 1.0/(4.0*N) + 1.0/(21.0*pow(N, 2.0)));

    ps.assign( (int) N, 0.0);
    pc.assign( (int) N, 0.0);

    // Compute weights and effective mu
    double sumWeight=0.0;
    weight.assign( mu, log(double(mu+1)) );

    for(unsigned int i=0; i<weight.size(); i++) {

        weight[i] -= log(double(i+1));
        sumWeight += weight[i];

    }

    muEff = 0.0;

    for(unsigned int i=0; i<weight.size(); i++) {

        weight[i] /= sumWeight;
        muEff += (weight[i] * weight[i]);

    }

    muEff = 1.0 / muEff;

    // csigma

    cs = (muEff+2.0)/(N+muEff+3.0);
    damp = 1.0 + 2.0*max( 0.0, ( sqrt( double(muEff - 1.0)/(N+1.0) ) - 1.0 ) ) + cs;
    cc = 4.0/(N+4.0);

    muCov = muEff;
    cCov = ( (1.0/muCov) * (2.0 / pow( (N+sqrt(2.0)), 2.0 ) ) )
           + ( (1.0 - (1.0/muCov) ) * min(1.0, ( ( 2.0*muEff - 1.0 ) / ( pow((N+2.0), 2.0) + muEff ) ) ) );

    // initialisation de D, B et C

    if ( C0.size() == 0 ) {

        //cout << "INIT!!" << endl;

        // creation du vecteur deltas

        vector<double> deltas;

        for (unsigned int i=0; i<problem->getSize(); i++) {

            deltas.push_back ( problem->getMaxVector()[i]-problem->getMinVector()[i] );

        }

        // intialisation

        D.assign( (int) N, 0.0 );
        B.assign( (int) N, D );
        C.assign( (int) N, D );

        // creation du vecteur D

        for (unsigned int i=0; i<problem->getSize(); i++) {

            D[i] = deltas[i];

        }

        // création de la matrice B

        for (unsigned int i=0; i<problem->getSize(); i++) {
            for (unsigned int j=0; j<problem->getSize(); j++) {

                if (i==j) B[i][j]=1.0;
                else B[i][j]=0.0;
            }
        }

        // création de la matrice C

        for (unsigned int i=0; i<problem->getSize(); i++) {
            for (unsigned int j=0; j<problem->getSize(); j++) {

                for (unsigned int k=0; k<problem->getSize(); k++) C[i][j]+=B[i][k]*D[k]*D[k]*B[j][k];

            }
        }
    }
    else {

        //cout << "NORMAL" << endl;

        C = C0; // récupération de la matrice de covariance

        // détermination de B et D pour la suite...
        D.assign( (int) N, 0.0 );
        B = C;
        computeEigensystem(D, B);  // Principal component analysis

        // Handle too large differences in eigenvalues D
        double MaxD = *max_element(D.begin(), D.end());
        double MinD = *min_element(D.begin(), D.end());

        if(MinD <= 0) {
            for (unsigned int i=0; i<N; i++) {

                if (D[i] < 0) D[i] = 0.0;

                C[i][i] += (MaxD/1e14);
                D[i]    += (MaxD/1e14);

            }
        }
        if(MaxD > (1e14*MinD)) {
            for (unsigned int i=0; i<N; i++) {

                C[i][i] += ((MaxD/1e14)-MinD);
                D[i]    += ((MaxD/1e14)-MinD);

            }
        }
        // Adjust D as standard deviation
        for(unsigned int i=0; i<N; i++) D[i] = sqrt(D[i]);

        // transposition de B due à la méthode computeEigensystem (qui donne les vecteurs propres en ligne!)

        for(unsigned int i=0; i<N; i++) {
            for(unsigned int j=0; j<N; j++) {

                double tmp = B[i][j];
                B[i][j] = B[j][i];
                B[j][i] = tmp;

            }
        }

    }

    // écriture des exits
//    cout << "\nC(n+1): " << endl;
//    cout << C[0][0] << "\t" << C[0][1] << endl;
//    cout << C[1][0] << "\t" << C[1][1] << endl;
//
//    cout << "\nD(n+1): " << endl;
//    cout << D[0] << endl;
//    cout << D[1] << endl;
//
//    cout << "\nB(n+1): " << endl;
//    cout << B[0][0] << "\t" << B[0][1] << endl;
//    cout << B[1][0] << "\t" << B[1][1] << endl;

    // individual creation

    if (popIndividualsVector.size() < mu) {

        vector<double> alleles;

        for (unsigned int i=0; i<popIndividualsVector.size(); i++) {

            pop.push_back( popIndividualsVector[i] );

        }

        unsigned int nCompl = mu - popIndividualsVector.size();

        for (unsigned int i=0; i<nCompl; i++) {

            alleles.clear();

            for (unsigned int j=0 ; j < problem->getSize() ; j++) {

                alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

            }

            pop.push_back( IndividualCSA_ES(problem,alleles) );

        }

        // évaluation de la population
        evaluate();
        // classement de la population
        sort( pop.begin(), pop.end() ); // du meilleur au moins bon
        // save intial random population
        save("currentpopCMA_ES_" + toString(problem->getnEval()) + ".dat");

    }
    else {

        //cout << "CONSERVATION!" << endl;

        for (unsigned int i=0; i<mu; i++) {

            pop.push_back( popIndividualsVector[i] );

        }

        // classement de la population
        sort( pop.begin(), pop.end() ); // du meilleur au moins bon
        // save intial random population
        save("currentpopCMA_ES_" + toString(problem->getnEval()) + ".dat");

    }

}

void CMA_ES::globalWeightedIntermediateRecombination() {

    // la population arrive triée

    vector<double> allelesFils;

    stringstream ss;

    // début du calcul de l'endroit test "xmean"

    for (unsigned int i=0;i<problem->getSize();i++) { // boucle sur les alleles

        double allele=0.0;

        for (unsigned int j=0; j < mu; j++ ) {

            allele += weight[j]*pop[j].getAllele(i);

        }

        allelesFils.push_back( allele );

    }

    // creation de lambda nouveaux individus

    for (unsigned int i=0; i<lambda; i++) {

        pop.push_back(IndividualCSA_ES(problem,allelesFils));

    }

//    writeResults("report.dat", ss);

}

void CMA_ES::mutationNormallyDistributed() {

    stringstream ss;

    unsigned int n = problem->getSize(); // nombre d'alleles
    double alleleprime = 0.0;

    for (unsigned int i=mu; i < mu+lambda; i++) {

        // try 10 times to be inside the constraints
        unsigned int round = 10;

        do {
            // boucle sur les alleles pour DEFINITION de Zi = N(0,1)

            for (unsigned int j=0; j < n; j++) {

                pop[i].setArz(j, normallyDistributedSPRNG_Ziggurat());

            }

            // boucle sur les alleles pour MUTATION proprement dite

            for (unsigned int j=0; j < n; j++) {

                ss << pop[i].getAllele()[j];

                // calcul de BDZ

                double BDZj = 0.0;

                for (unsigned int k=0; k<n; k++) {

                    BDZj += B[j][k]*D[k]*pop[i].getArz(k);

                }

                // calcul de l'allele'
                alleleprime = pop[i].getAllele(j) + sigmaF*BDZj;

                // 1er check: dans les limites de la variable - sauvegarde de alleleprime
                alleleprime = min( alleleprime, problem->getMaxVector()[j] - TINY );
                alleleprime = max( alleleprime, problem->getMinVector()[j] );
                pop[i].setAllele( j, alleleprime );

                ss << "/" << alleleprime << " ";

            }

            ss << endl;
            round--;

        }
        while ( problem->notInField( pop[i].getAllele() ) && round > 0);

    }

//    writeResults("report.dat", ss);

}

void CMA_ES::selectionCommaRank() {

    int posMin;

    stringstream ss;

    for (unsigned int j=0; j < mu;j++) {     // suppression des parents

        ss << "suppression de l'individu" << pop[mu-1-j].getIdNumber() << endl;
        pop.erase(pop.begin() + (mu-1-j)); // suppression de l'individu

    }

    for (unsigned int j=0; j < (lambda-mu);j++) {   // selection des mu meilleurs parmi les lambda enfants restants
                                           // suppression des lambda-mu moins bons
        posMin = distance(pop.begin(), max_element(pop.begin(), pop.end()));

        ss << "suppression de l'individu" << pop[posMin].getIdNumber() << endl;
        pop.erase(pop.begin() + (posMin)); // suppression de l'individu moins bon

    }

    // mise en place du meilleur dans bestIndividual
    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));
    bestIndividualsVector.push_back( pop[posMax] );
    // à modifier pour la finale de RANK
//    writeResults("report.dat", ss);

    return;

}

void CMA_ES::adaptation() {

    stringstream ss;

    // sort population
    sort(pop.begin(), pop.end());

    // début de l'adaptation
    unsigned int n = int(N);

    vector<double> zmeanwi;
    zmeanwi.assign(n, 0.0);

    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j < mu; j++ ) {

            zmeanwi[i] += weight[j]*pop[j].getArz(i); // les poids sont normalisés

        }
    }

    vector<double> BDZm;
    BDZm.assign(n, 0.0);

    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++ ) {

            BDZm[i] += B[i][j]*D[j]*zmeanwi[j];

        }
    }

    // adaptation de la matrice de variance/covariance

    double HLeft = norm(ps) / sqrt(1.0-pow(1.0-cs, 2.0*double(getGeneration()+1)));
    double HRight = (1.5+(1.0/(N-0.5))) * chiN;

    double Hs;
    if ( HLeft < HRight ) Hs = 1.0;
    else Hs = 0.0;

    // pc
    for (unsigned int i=0; i<n; i++) {
        pc[i] = (1-cc)*pc[i] + Hs*sqrt(muEff*cc*(2.0-cc))*BDZm[i];
    }

    // écriture des exits

//    cout << "C(n): " << endl;
//    cout << C[0][0] << "\t" << C[0][1] << endl;
//    cout << C[1][0] << "\t" << C[1][1] << endl;
//
//    cout << "\nD(n): " << endl;
//    cout << D[0] << endl;
//    cout << D[1] << endl;
//
//    cout << "\nB(n): " << endl;
//    cout << B[0][0] << "\t" << B[0][1] << endl;
//    cout << B[1][0] << "\t" << B[1][1] << endl;

    // écriture de C au temps +1
    for (unsigned int i=0; i<n; i++) {
        for (unsigned int j=0; j<n; j++) {

            double sumBDZiBDZj = 0.0;

            for (unsigned int k=0; k<mu; k++) { // individu k

                // calcul de BDZ[i] et BDZ[j] pour l'individu k
                double BDZi=0.0, BDZj=0.0;

                for (unsigned int l=0; l<n; l++) {

                    BDZi += B[i][l]*D[l]*pop[k].getArz(l);
                    BDZj += B[j][l]*D[l]*pop[k].getArz(l);

                }

                sumBDZiBDZj += weight[k]*BDZi*BDZj;

            }

            C[i][j] = (1.0-cCov+(1.0-Hs)*cCov*cc*(2.0-cc)/muCov)*C[i][j]
                    + cCov*(1.0/muCov)*pc[i]*pc[j]
                    + cCov*(1.0 - 1.0/muCov)*sumBDZiBDZj;
        }
    }

    // ps
    for (unsigned int i=0; i<n; i++) {

        double BZmi = 0.0;

        for (unsigned int j=0; j<n; j++ ) {

            BZmi += B[i][j]*zmeanwi[j];

        }

        // adaptation de sigma - ps
        ps[i] = (1-cs)*ps[i] + sqrt(muEff*cs*(2.0-cs))*BZmi;

    }

    // adaptation de sigmaF
    sigmaF *= exp( (norm(ps)/chiN-1.0)*(cs/damp) );

    // Update B and D from C
    for(unsigned int i=1; i<n; ++i) {
        for(unsigned int j=0; j<i; ++j) C[i][j] = C[j][i];  // Enforce symetry
    }

    B = C;
    computeEigensystem(D, B);  // Principal component analysis

    // Handle too large differences in eigenvalues D
    double MaxD = *max_element(D.begin(), D.end());
    double MinD = *min_element(D.begin(), D.end());

    if(MinD <= 0) {
        for (unsigned int i=0; i<n; i++) {

            if (D[i] < 0) D[i] = 0.0;

            C[i][i] += (MaxD/1e14);
            D[i]    += (MaxD/1e14);

        }
    }
    if(MaxD > (1e14*MinD)) {
        for (unsigned int i=0; i<n; i++) {

            C[i][i] += ((MaxD/1e14)-MinD);
            D[i]    += ((MaxD/1e14)-MinD);

        }
    }
    // Adjust D as standard deviation
    for(unsigned int i=0; i<n; i++) D[i] = sqrt(D[i]);

    // transposition de B due à la méthode computeEigensystem (qui donne les vecteurs propres en ligne!)

    for(unsigned int i=0; i<n; i++) {
        for(unsigned int j=0; j<n; j++) {

            double tmp = B[i][j];
            B[i][j] = B[j][i];
            B[j][i] = tmp;

        }
    }

    // écriture des exits
//    cout << "\nC(n+1): " << endl;
//    cout << C[0][0] << "\t" << C[0][1] << endl;
//    cout << C[1][0] << "\t" << C[1][1] << endl;
//
//    cout << "\nD(n+1): " << endl;
//    cout << D[0] << endl;
//    cout << D[1] << endl;
//    cout << "MaxD: " << MaxD << "\tMinD: " << MinD << endl;
//
//
//    cout << "\nB(n+1): " << endl;
//    cout << B[0][0] << "\t" << B[0][1] << endl;
//    cout << B[1][0] << "\t" << B[1][1] << endl;
//
//    string essai;
//    cin >> essai;

    // Adjustments due to diversity of population
//    if (diversityPopulation() < 0.1) sigmaF *= exp(0.2+cs/damp);

    // sortie du nouveau sigmaF
    ss << "SigmaF: " << sigmaF << endl;

//    writeResults("report.dat", ss);

}

void CMA_ES::step() {

//      **** recombination ****
    globalWeightedIntermediateRecombination();

//      **** mutation *********
    mutationNormallyDistributed();

//      **** evaluation *******
    evaluate();

//      **** selection ********
    selectionCommaRank();
    statistics();

//  **** Evaluation of the diversities ****

    save("currentpopCMA_ES_" + toString(problem->getnEval()) + ".dat");
    ostringstream ss;
    ss << "bestIndividualID: " << getCurrentBestIndividual().getIdNumber();
    ss << "\tpopulationDiversity: " << diversityPopulation() << endl;
    ss << "diversities: ";
    for (unsigned int j=0; j<getInitSize(); j++) {
        ss << diversityIndividual(j) << "\t";
    }
    ss << endl;
    writeResults("currentpopCMA_ES_" + toString(problem->getnEval()) + ".dat", ss);

//      **** parameter adaptation
    adaptation();

//      **** end of step ******
    incrementGeneration();

    return;

}

void CMA_ES::getBestInPopulation(vector< Individual > &popIndividualsVector) {

    for (unsigned int i=0; i<bestIndividualsVector.size(); i++) {

        popIndividualsVector.push_back( (Individual) bestIndividualsVector[i] );

    }

    return;

}

double CMA_ES::getSigmaF() {

    return sigmaF;

}

vector< vector<double> > CMA_ES::getC() {

    return C;

}

// *****************************************************************
// Class derived from Heuristic: PSO
// *****************************************************************

PSO::PSO(Problem *problem, unsigned int size, double c1, double c2, double lambda, double kappa):Heuristic<IndividualPSO>(problem), c1(c1), c2(c2), lambda(lambda), kappa(kappa) {

    // la taille initiale devrait etre un carré
    size = (unsigned int) pow(ceil(sqrt(double(size))), 2.0);

    // sauvegarde de la taille initiale
    initSize = size;

    // initialisation of the population
    vector<double> alleles;

    for (unsigned int i=0;i<size;i++) {

        alleles.clear();

        for (unsigned int j=0; j < problem->getSize(); j++) {

            alleles.push_back( randomUniform(problem->getMinVector()[j],problem->getMaxVector()[j]) ); // entre 0 et max

        }

	    pop.push_back( IndividualPSO(problem, alleles) );

    }

    evaluate(); // evaluation of the individuals

    bestEverPop = pop;

}

void PSO::step() {

    // initialisation done with creation: randomly chosen in the domain, and speeds = 0, those are already evaluated
    updateLocalBestParticles(); // saves the local best particles in bestPopEver
    update();
    evaluateGrid();
    selection();

    statistics();
    report();

    incrementGeneration();
    return;

}

void PSO::updateLocalBestParticles() {

    for (unsigned int i=0; i<bestEverPop.size(); i++) { // comparison between the actual and old population

        if ( pop[i] < bestEverPop[i] ) bestEverPop[i]=pop[i];

    }

}

IndividualPSO* PSO::getGlobalBestParticle(unsigned int i) {

    IndividualPSO *pBestIndividual = &bestEverPop[i];

    unsigned int sqrtNP = (unsigned int)(sqrt(double(getInitSize())));

    unsigned int l = i / sqrtNP;
    unsigned int c = i % sqrtNP;

    // cas c+1
    if ( bestEverPop[ l*sqrtNP + (c+1)%sqrtNP ] < *pBestIndividual ) pBestIndividual = &bestEverPop[ l*sqrtNP + (c+1)%sqrtNP ];
    // cas c-1
    if ( bestEverPop[ l*sqrtNP + (c-1+sqrtNP)%sqrtNP ] < *pBestIndividual ) pBestIndividual = &bestEverPop[ l*sqrtNP + (c-1+sqrtNP)%sqrtNP ];
    // cas l+1
    if ( bestEverPop[ ((l+1)%sqrtNP)*sqrtNP + c ] < *pBestIndividual ) pBestIndividual = &bestEverPop[ ((l+1)%sqrtNP)*sqrtNP + c ];
    // cas l-1
    if ( bestEverPop[ ((l-1+sqrtNP)%sqrtNP)*sqrtNP + c ] < *pBestIndividual ) pBestIndividual = &bestEverPop[ ((l-1+sqrtNP)%sqrtNP)*sqrtNP + c ];

    // write out a test to see if correct!
//    cout << "index: " << i << "\t" << "l: " << l << "\tc: " << c << endl;
//    for (unsigned int ll=0; ll<sqrtNP; ll++) {
//        for (unsigned int cc=0; cc<sqrtNP; cc++) {
//            cout << bestEverPop[ ll*sqrtNP + cc%sqrtNP ].getFitness() << "\t";
//        }
//        cout << "\n";
//    }
//    cout << endl;
//    cout << "best in neighbourhood: " << (*pBestIndividual).getFitness() << endl;

    // von Neumann strategy
    return pBestIndividual;

}

void PSO::update() {

    double xsi = 0.0;
    double phi = c1 + c2;
	if (phi > 4) xsi = 2. * kappa / abs( 2. - phi - sqrt( phi * ( phi - 4. ) ) );
	else xsi = kappa;

    // update the particle location
	double rho1 = randomUniform(0.0,1.0);
	double rho2 = randomUniform(0.0,1.0);
	for (unsigned int iP = 0; iP < getInitSize(); iP++) {
	    vector<double> position,velocity;
	    position.assign(pop[iP].getAllelesNumber(), 0.0);
	    velocity.assign(pop[iP].getAllelesNumber(), 0.0);
	    // velocity of continuous parameters
	    for (unsigned int j=0; j< pop[iP].getAllelesNumber(); j++) {

            velocity[j] = pop[iP].getV(j) + c1*rho1*( bestEverPop[iP].getAllele(j) - pop[iP].getAllele(j) );
            velocity[j]+= c2*rho2*( (*getGlobalBestParticle(iP)).getAllele(j) - pop[iP].getAllele(j) );
            velocity[j]*= xsi;

            // clamp velocity
            double maxVelocity = lambda * (problem->getMaxVector()[j]-problem->getMinVector()[j]);
            // lambda \in R+, so always postive, if negative, no clamping used!
            if (maxVelocity > 0) {
                if ( velocity[j] > 0 ) velocity[j] = min(velocity[j], maxVelocity);
                else                   velocity[j] = max(velocity[j], -maxVelocity);
            }

            // update particle position, by creating a new particle
            position[j] = pop[iP].getAllele(j) + velocity[j];

            // verify position in the domain and put it back!
            double positionPre = position[j];
            do {
                positionPre = position[j];
                if (position[j] < problem->getMinVector()[j]) position[j] = 2*problem->getMinVector()[j] - position[j];
                else if (position[j] > problem->getMaxVector()[j]) position[j] = 2*problem->getMaxVector()[j] - position[j];
            }
            while (positionPre != position[j]);

	    }
	    pop.push_back( IndividualPSO(problem, position, velocity) );

	}

    return;

}

void PSO::selection() {

    pop.erase(pop.begin(), pop.begin()+getInitSize());

    return;

}

// Hooke-Jeeves

Hooke_Jeeves::Hooke_Jeeves(Problem *problem, Individual& xInit, unsigned int r, unsigned int s0, unsigned int t, unsigned int m):Heuristic<Individual>(problem), r(r),t(t),s(s0),m(m) {

    // so far only the init individual is in the population
    pop.push_back(xInit);
    pop.push_back(xInit); // the new one pop[1] is equal to the old one pop[0]

    if (problem->getStepVector().size() == 0) throw ("No step size(s), cannot launch Hooke-Jeeves");

    // modify the alleles of the first two indidividuals to put them back to the grid
    for (unsigned int i=0; i<pop.size(); i++) {

        vector<double> alleles(pop[i].getAlleles());
        // check is stepSize present and put them back to the grid
        for (unsigned int j=0; j<alleles.size(); j++) alleles[j] = problem->getStepVector()[j]*round(alleles[j]/problem->getStepVector()[j]);
        pop[i].setAlleles(alleles);

    }

    // prepare the deltas
    delta.assign( xInit.getAllelesSize(), 0.0 );
    for (unsigned int i=0; i<xInit.getAllelesSize(); i++) delta[i] = problem->getStepVector()[i]/pow(double(r),s); //pow(double,int)

    // some output for checking

//    cout << setprecision(12);
//    cout << "Alleles: " << endl;
//    for (unsigned int i=0; i<xInit.getAllelesSize()-1; i++) cout << pop[0].getAllele(i) << "\t";
//    cout << pop[0].getAllele(xInit.getAllelesSize()-1) << endl;
//    cout << "Step sizes: " << endl;
//    for (unsigned int i=0; i<xInit.getAllelesSize()-1; i++) cout << problem->getStepVector()[i] << "\t";
//    cout << problem->getStepVector()[xInit.getAllelesSize()-1] << endl;
//    cout << "Deltas: " << endl;
//    for (unsigned int i=0; i<xInit.getAllelesSize()-1; i++) cout << delta[i] << "\t";
//    cout << delta[xInit.getAllelesSize()-1] << endl;

}

void Hooke_Jeeves::globalSearch() {

    if (pop[1] < pop[0]) { // the new one is better than the old one -> make a move in the good direction!

        vector<double> alleles;
        alleles.assign( delta.size() , 0.0 );

        for (unsigned int i = 0; i < delta.size(); i++) {

           /* firstly, arrange the sign of delta[] */
           if (pop[1].getAllele(i) <= pop[0].getAllele(i))  delta[i] = - abs(delta[i]);
           else                                             delta[i] =   abs(delta[i]);
           /* now, move further in this direction */
           alleles[i] = pop[1].getAllele(i) + (pop[1].getAllele(i)-pop[0].getAllele(i));

        }

        pop.push_back( Individual(problem,alleles) );

    //    cout << "population: " << endl;
    //    for (unsigned int i=0; i<pop.size(); i++) {
    //        cout << pop[i] << "\n" << "Fitness: " << pop[i].getFitness() << endl;
    //    }
    //    cout << "\nBase individual: " << endl;
    //    cout << pop[2] << "\nFitness: " << pop[2].getFitness() << endl;

        pop.erase(pop.begin()); // delete the former best one
    }

    return;

}

void Hooke_Jeeves::localSearch() {

unsigned int iLow = pop.size()-1;

    for (unsigned int i = 0; i < delta.size(); i++) {

        // take the new one as a base
        vector<double> baseAlleles(pop[iLow].getAlleles());

        // vector for new individual
        vector<double> alleles;

        alleles=baseAlleles;
        alleles[i] += delta[i];
        pop.push_back(Individual(problem,alleles));

        evaluate();
        if (pop[pop.size()-1] < pop[iLow]) iLow = pop.size()-1;
        else { // try the other direction
            alleles=baseAlleles;
            alleles[i] -= delta[i];
            pop.push_back(Individual(problem,alleles));
            evaluate();
            if (pop[pop.size()-1] < pop[iLow]) { iLow = pop.size()-1; delta[i] = -delta[i]; }
        }

    }

}

void Hooke_Jeeves::selection() {

//    cout << "population: " << endl;
//    for (unsigned int i=0; i<pop.size(); i++) {
//        cout << pop[i] << "\n" << "Fitness: " << pop[i].getFitness() << endl;
//    }

    vector<Individual>::iterator it = min_element(pop.begin()+1, pop.end()); // find the best one along the new ones (index > 0)

//    cout << "première commande: " << endl;
    pop.erase(pop.begin()+1, it);

//    cout << "population: " << endl;
//    for (unsigned int i=0; i<pop.size(); i++) {
//        cout << pop[i] << "\n" << "Fitness: " << pop[i].getFitness() << endl;
//    }

//    cout << "deuxième commande: " << endl;
    pop.erase(pop.begin()+2, pop.end());

//    cout << "population: " << endl;
//    for (unsigned int i=0; i<pop.size(); i++) {
//        cout << pop[i] << "\n" << "Fitness: " << pop[i].getFitness() << endl;
//    }

    return; // two individuals still in the population

}

void Hooke_Jeeves::adaptation() {

    if ( !(pop[1]<pop[0]) ) { // the new one is not better than the old one -> adapt the delta

        if ( pop[1].getAlleles() == pop[0].getAlleles() ) { // copies conformes, adaptation du delta
            if (m<=0) m--; // end of the game
            else {
                cout << "No increase, reducing step size: " << endl;
                s += t;
                for (unsigned int i=0; i<delta.size(); i++) {
                    cout << delta[i] << " -> ";
                    delta[i] = sign(delta[i])*abs(problem->getStepVector()[i]/pow(double(r),s)); // pow(double,int)
                    cout << delta[i] << endl;
                }
                m--;
            }
        }

        // delete the useless new one
        pop.erase(pop.begin()+1);

        // make a copy of the old one as new one
        pop.push_back(pop[0]);

    }

}

void Hooke_Jeeves::step() {

    globalSearch();
    localSearch();
    report();
    selection();
    adaptation();

    statistics();
    report();

    incrementGeneration();

    return;

}
