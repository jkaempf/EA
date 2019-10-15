#ifndef _HEURISTIC_H
#define _HEURISTIC_H

#include "evalpthread.h"

// Individual Class definitions

#include "problem.h"
#include <numeric>

class Individual {

    private:

      Problem *problem;

    protected:

      // compteur d'individus
      static unsigned int compteur; // static car commun a toutes les instances de la classe

      // valeur associée à l'individu
      vector<double> fitness;

      // parametres
      vector<double> alleles;

      // status
      bool simulated;

      // identification
      unsigned int number;

    public:

      Individual(Problem *problem, vector<double> alleles);
      Individual(Problem *problem, vector<double> alleles, double fitness);
      Individual(Problem *problem, vector<double> alleles, int number, double fitness);

      friend bool operator<(Individual &a, Individual &b);
	  friend bool operator==(Individual &a, Individual &b);
      friend ostream& operator <<(ostream &os, Individual obj);

      Problem* getProblem() { return problem; }

      unsigned int getIdNumber() { return number; };
      bool getSimulated() { return simulated; };
      void setSimulated(bool value) { simulated=value; return; };
      double getFitness() { return fitness[0]; };
	  double getFitness(size_t index) { return fitness[index]; }
      vector<double> getFitnessVector() { return fitness; }
	  size_t getFitnessSize() { return fitness.size(); }
      void setFitness(double fitness) { this->fitness.assign(1, fitness); simulated = true; };
      void setFitness(vector<double> fitness) { this->fitness = fitness; simulated = true; };
      double getAllele(unsigned int i) { return alleles[i]; };
      void setAllele(unsigned int i, double value) { alleles[i]=value; };
      void setAlleles(vector<double> alleles) { this->alleles = alleles; };
      vector<double> getAllele() { return alleles; }; // TODO: add an "s" to Allele everywhere
      vector<double> getAlleles() { return alleles; };
      unsigned int getAllelesNumber() { return alleles.size(); }; // TODO: remove this in the code
      unsigned int getAllelesSize() { return alleles.size(); };

      static void resetCounter() { compteur = 0; return; };

};

class IndividualES : public Individual {

private:

  // parametres
  vector<double> sigmas;

public:

  IndividualES(Problem *problem, vector<double> alleles, vector<double> sigmas):Individual(problem, alleles), sigmas(sigmas) {}
  IndividualES(Problem *problem, vector<double> alleles, vector<double> sigmas, double fitness):Individual(problem, alleles, fitness), sigmas(sigmas) {}
  IndividualES(Problem *problem, vector<double> alleles, vector<double> sigmas, int number, double fitness):Individual(problem, alleles, number, fitness), sigmas(sigmas) {}

  double getSigma(unsigned int i) { return sigmas[i]; };
  void setSigma(unsigned int i, double value) { sigmas[i]=value; return; }

  unsigned int getSigmasNumber() { return sigmas.size(); };

  vector<double> getSigma() { return sigmas; };

};

class IndividualCSA_ES : public Individual {

private:

  // parametres
  vector<double> arz;

public:

  IndividualCSA_ES(Individual individual):Individual(individual) { arz.assign( alleles.size(), 0.0 ); }
  IndividualCSA_ES(Problem *problem, vector<double> alleles):Individual(problem, alleles) { arz.assign( alleles.size(), 0.0 ); }
  IndividualCSA_ES(Problem *problem, vector<double> alleles, double fitness):Individual(problem, alleles, fitness) { arz.assign( alleles.size(), 0.0); }

  double getArz(unsigned int i) { return arz[i]; };
  void setArz(unsigned int i, double value) { arz[i]=value; return; };
  unsigned int getArzNumber() { return arz.size(); };
  vector<double> getArz() { return arz; };

};

class IndividualPSO : public Individual {

private:

    // vecteur vitesses
    vector<double> v;

public:

    IndividualPSO(Problem *problem, vector<double> alleles):Individual(problem, alleles) { v.assign( alleles.size(), 0.0); }
    IndividualPSO(Problem *problem, vector<double> alleles, vector<double> velocity):Individual(problem, alleles), v(velocity) {}

    double getV(unsigned int i) { return v[i]; };
    void setV(unsigned int i, double value) { v[i]=value; return; };

};

// *************************************************************************
// Heuristic class definition - common object for all algorithms
// *************************************************************************

template <class T> class Heuristic {

    private:

        // keep the evaluated individuals in a bank for savings
        map<vector<double>, vector<double> > allelesFitnessMap;

    protected:

        Problem* problem;

        unsigned int gen;
        unsigned int initSize;

        vector<T> pop;
        T* bestIndividual;

        // statistics
        double leMax, leMin, laMoy, laSomme;
        unsigned int posMax;

    public:

        Heuristic(Problem *problem) {

          this->problem = problem;

          leMax=numeric_limits<double>::max();
          leMin=-numeric_limits<double>::max();
          laMoy=0.;
          laSomme=0.;
          posMax=0;
          gen=0;

          bestIndividual = NULL;

        }

        Heuristic(Problem *problem, int size) {

          this->problem = problem;

          leMax=numeric_limits<double>::max();
          leMin=-numeric_limits<double>::max();
          laMoy=0.;
          laSomme=0.;
          posMax=0;
          gen=0;

          bestIndividual = NULL;

          initSize = size;

        }

        virtual ~Heuristic() { delete bestIndividual; };

        void pushIndividual(T x) { pop.push_back(x); };
        T& getIndividual(unsigned int i) { return pop[i]; };
        void eraseIndividual(unsigned int i) { pop.erase(pop.begin() + i); };

        void updateBestIndividual();
        void statistics();
        virtual void report();
        void listIndividuals(string filename);
        virtual void save(string filename);

        double getMaximum() { return leMax; };
        double getMinimum() { return leMin; };
        double getMoyenne() { return laMoy; };
        unsigned int getPosMax() { return posMax; };

        void setInitSize(unsigned int size) { initSize = size; };
        unsigned int getInitSize() { return initSize; } ;

        unsigned int getPopSize() { return pop.size(); };

        vector<double> getAlleles(unsigned int id) { return pop[id].getAllele(); };
        unsigned int getIdNumber(unsigned int id) { return pop[id].getIdNumber(); };

        void incrementGeneration() { gen++; };
        int getGeneration() { return gen; };
        void setGeneration(int generationNumber) { gen = generationNumber; };

        double getBestIndividualFitness() { return bestIndividual->getFitness(); };
		double getBestIndividualFitness(size_t i) { return bestIndividual->getFitness(i); };
        T& getBestIndividual() { return *bestIndividual; };
        T& getCurrentBestIndividual() { return pop[posMax]; };
        bool compare(int n1, int n2) { return ( pop[n1].getAllele() == pop[n2].getAllele() ); /* notifie qu'il y a un doublon */ };

        void evaluateVirtual();
        void evaluate();
        void evaluateGrid();

        bool flatFitness();

        double diversityIndividual(unsigned int i);
        double diversityPopulation();

        void getPopulation(vector< Individual > &popIndividualsVector);
        void getPopulation(int size, vector< Individual > &popIndividualsVector);
        virtual void getBestInPopulation(vector<Individual> &popIndividualsVector);

};

template <class T> void Heuristic<T>::updateBestIndividual() {

// *****************************************************************
// formerly Heuristic::setBestIndividual()
// *****************************************************************

    // defines the position of the best element in the current pop
    // the comparison operator has been redefined to allow an order in the population:
    // template <class T> bool operator<(T a, T b) in this file

    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));

    // A new best one?

    if (bestIndividual == NULL) {

        bestIndividual = new T ( pop[posMax] );

    }
    else {

        if ( pop[posMax] < *bestIndividual ) {

            *bestIndividual = pop[posMax];

        }
    }

    return;

}

template <class T> void Heuristic<T>::statistics() {

    vector<double> fitnessVector;

    for (unsigned int i=0;i<pop.size();i++) {

        fitnessVector.push_back( pop[i].getFitness() );

    }

    laSomme = std::accumulate( fitnessVector.begin(), fitnessVector.end(), 0.0 ); // 0.0 valeur initiale
    leMax = *max_element(fitnessVector.begin(), fitnessVector.end());
//    posMax = distance(fitnessVector.begin(), max_element(fitnessVector.begin(), fitnessVector.end()));

    leMin = *min_element(fitnessVector.begin(), fitnessVector.end());
    laMoy = laSomme / (double) pop.size();

    updateBestIndividual();

    cout << "Génération: " << gen << "\tmax/min(moy): " << leMax << "\t" << leMin << "\t" << laMoy << endl;
    cout << "bestIndividual: id" << bestIndividual->getIdNumber() << "\tfitness: " << bestIndividual->getFitness() << endl;

    return;

}

template <class T> void Heuristic<T>::report() {

    ostringstream ss1, ss2;

    ss1 << "\nRapport de génération n°: " << gen << endl;

    for (unsigned int i=0;i<pop.size();i++) {

      ss1 << "Le membre n°: " << i << " (" << "building" << pop[i].getIdNumber() << ")\tSa fitness: " << pop[i].getFitness() << endl;

      ss1 << "\talleles: " << "(";

      for (unsigned int j=0;j<pop[i].getAllelesNumber();j++) {

        ss1 << "[" << pop[i].getAllele(j) << "], ";

      }

      ss1 << "[end])" << endl;

    }

    ss1 << "\nLe meilleur : " << posMax << " (" << "individual" << pop[posMax].getIdNumber() << ")\tSa fitness: " << leMax << "\nmax/min(moy): " << setprecision(3) << leMax << "/" << leMin << "(" << laMoy << ")" << endl;

    // ecriture dans le stringstream du rapport

    ss2 << gen << "\t" << laMoy << "\t" << leMax << "\t" << leMin << "\t" << "individual" << pop[posMax].getIdNumber() << "alleles: ";
    for (unsigned int i=0;i<pop[posMax].getAllelesSize();i++) { ss2 << pop[posMax].getAllele(i) << "\t"; }
    ss2 << endl;

    writeResults("report.dat", ss1);
    writeResults("population.dat", ss2);

    return;

}

template <class T> void Heuristic<T>::listIndividuals(string filename) {

    ostringstream ss;

    for (unsigned int i=0;i<pop.size();i++) {

        ss << "individual" << pop[i].getIdNumber() << endl;

    }

    writeResults(filename, ss);

}

template <class T> void Heuristic<T>::save(string filename) {

        ostringstream ss;

        ss << gen << "\t" << pop.size() << "\t" << pop[0].getAllelesNumber() << endl;

        for (unsigned int i=0;i<pop.size();i++) {

            ss << "individual" << pop[i].getIdNumber() << "\t";
			for (size_t nFitness=0; nFitness< pop[i].getFitnessSize(); ++nFitness)
				ss << pop[i].getFitness(nFitness) << "\t";

            ss << setprecision(12) << scientific;
            for (unsigned int j=0;j<pop[i].getAllelesNumber();j++) {

                ss << pop[i].getAllele(j) << "\t";

            }

            ss << endl;

        }

        writeResultsOver(filename, ss);

}

template <class T> void Heuristic<T>::evaluateVirtual() {

      for (int i=0;i<pop.size();i++) {

            if (!pop[i].getSimulated()) {

                pop[i].setFitness(100.e6*(i+1));

            }
      }

      return;

}

//template <class T> void Heuristic<T>::evaluate() {
//
//    for (unsigned int i=0;i<pop.size();i++) {
//
//        if ( !pop[i].getSimulated() && !problem->notInField(pop[i].getAllele()) ) {
//
//            pop[i].setFitness(problem->evaluate(pop[i].getIdNumber(), pop[i].getAllele()));
//
//        }
//    }
//
//    return;
//
//}

template <class T> void Heuristic<T>::evaluate() {

    // have a map -> allelesFitnessMap
    // see if the alleles exists in the map (find, if not found return end())

    // if they exists      -> don't put in evaluation class
    //                     -> extract result from map
    // if they don't exist -> put in evaluation class
    // launch the evaluation
    // get the results

    // create the evaluation object
    Evaluation eval(problem);

    for (unsigned int i=0;i<pop.size();i++) {

        if ( !pop[i].getSimulated() && !problem->notInField(pop[i].getAllele()) ) {

            // if equals to end, then it does not exist in the set
            if ( allelesFitnessMap.find(pop[i].getAlleles()) == allelesFitnessMap.end() ) {

//                cout << "Alleles does not exist in the set" << endl;
//
//                // show the alleles
//                cout << "[";
//                for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//                cout << alleles[alleles.size()-1] << "]" << endl;

                // add it to the evaluation -> in put to the grid form
                // insertion dans la fitnessid map de Evaluation
                eval.insert(pop[i].getIdNumber(), pop[i].getAlleles());
            }
            else {

//                cout << "Alleles exist in the set" << endl;
//
//                // show the alleles
//                cout << "[";
//                for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//                cout << alleles[alleles.size()-1] << "]" << endl;
//                // show the fitness
//                cout << "Fitness: " << allelesFitnessMap[alleles][0] << endl;

                // set the corresponding value from the map
                pop[i].setFitness(allelesFitnessMap[pop[i].getAlleles()]);

            }
        }
    }

    // début de l'évaluation
    eval.start();

    for (unsigned int i=0;i<pop.size();i++) {

        if ( !pop[i].getSimulated() && !problem->notInField(pop[i].getAllele()) ) {

            // if here then is not in the set so add it
            allelesFitnessMap.insert( pair<vector<double>, vector<double> >(pop[i].getAlleles(), eval.get(pop[i].getIdNumber())) );

            // set the fitness in the individual
            pop[i].setFitness(eval.get(pop[i].getIdNumber()));

//            cout << "Inserted in the set" << endl;
//            // show the alleles
//            cout << "[";
//            for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//            cout << alleles[alleles.size()-1] << "]" << endl;
//            // show the fitness
//            cout << "Fitness: " << eval.get(pop[i].getIdNumber())[0] << endl;

            //cout << "Received: " << pop[i].getIdNumber() << "\tfitness: " << allelesFitnessMap[pop[i].getAllele()][0] << endl;

        }
    }

    return;

}

template <class T> void Heuristic<T>::evaluateGrid() {

    // have a map -> allelesFitnessMap
    // see if the alleles exists in the map (find, if not found return end())

    // if they exists      -> don't put in evaluation class
    //                     -> extract result from map
    // if they don't exist -> put in evaluation class
    // launch the evaluation
    // get the results

    // create the evaluation object
    Evaluation eval(problem);

    for (unsigned int i=0;i<pop.size();i++) {

        if ( !pop[i].getSimulated() && !problem->notInField(pop[i].getAllele()) ) {

            // get the alleles of the current one
            vector<double> alleles = pop[i].getAlleles();

            // check is stepSize present and put them back to the grid
            for (unsigned int j=0; j<alleles.size(); j++) alleles[j] = problem->getStepVector()[j] * round( alleles[j]/problem->getStepVector()[j] );

            // if equals to end, then it does not exist in the set
            if ( allelesFitnessMap.find(alleles) == allelesFitnessMap.end() ) {

//                cout << "Alleles does not exist in the set" << endl;
//                // show the alleles
//                cout << "[";
//                for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//                cout << alleles[alleles.size()-1] << "]" << endl;

                // add it to the evaluation -> in put to the grid form
                // insertion dans la fitnessid map de Evaluation
                eval.insert(pop[i].getIdNumber(), alleles);
            }
            else {

//                cout << "Alleles exist in the set" << endl;
//                // show the alleles
//                cout << "[";
//                for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//                cout << alleles[alleles.size()-1] << "]" << endl;
//                // show the fitness
//                cout << "Fitness: " << allelesFitnessMap[alleles][0] << endl;

                // set the corresponding value from the map
                pop[i].setFitness(allelesFitnessMap[alleles]);

                // increment the evaluation
                problem->incrementnEval();

            }
        }
    }

    // début de l'évaluation
    eval.start();

    for (unsigned int i=0;i<pop.size();i++) {

        if ( !pop[i].getSimulated() && !problem->notInField(pop[i].getAllele()) ) {

            // get the alleles of the current one
            vector<double> alleles = pop[i].getAlleles();

            // check is stepSize present and put them back to the grid
            for (unsigned int j=0; j<alleles.size(); j++) alleles[j] = problem->getStepVector()[j] * round( alleles[j]/problem->getStepVector()[j] );

            // if here then is not in the set so add it
            allelesFitnessMap.insert( pair<vector<double>, vector<double> >(alleles, eval.get(pop[i].getIdNumber())) );

            // set the fitness in the individual
            pop[i].setFitness(eval.get(pop[i].getIdNumber()));

//            cout << "Inserted in the set" << endl;
//            // show the alleles
//            cout << "[";
//            for (unsigned int j=0; j<alleles.size()-1; j++) cout << alleles[j] << " ";
//            cout << alleles[alleles.size()-1] << "]" << endl;
//            // show the fitness
//            cout << "Fitness: " << eval.get(pop[i].getIdNumber())[0] << endl;

            //cout << "Received: " << pop[i].getIdNumber() << "\tfitness: " << allelesFitnessMap[pop[i].getAllele()][0] << endl;

        }
    }

    return;

}

template <class T> bool Heuristic<T>::flatFitness() {

    if ( abs(leMax - leMin) < TINY ) return true;
    else return false;

}

template <class T> void Heuristic<T>::getPopulation(vector< Individual > &popIndividualsVector) {

    for (int i=0;i<pop.size();i++) {

        popIndividualsVector.push_back( (Individual) pop[i] );

    }

    return;

}

template <class T> void Heuristic<T>::getPopulation(int size, vector< Individual > &popIndividualsVector) {

    sort(pop.begin(), pop.end()); // classe la population du rank le moins élevé au plus élevé

    for (int i=0; i < size; i++) { // prend les meilleurs

        popIndividualsVector.push_back( pop[i] );

    }

    return;

}

template <class T> void Heuristic<T>::getBestInPopulation(vector<Individual> &popIndividualsVector) {

    popIndividualsVector.push_back( *bestIndividual );

    return;

}

// *****************************************************************
// Diversity of the population defined in Chang07
// *****************************************************************

template <class T> double Heuristic<T>::diversityIndividual(unsigned int i) {

    // defines the position of the best element in the current pop
    posMax = distance(pop.begin(), min_element(pop.begin(), pop.end()));

    double epsilon2 = 0.1; // 10% or less of allele relative difference insignificant
    //double precision = 1e-3;//5e-7; // difference less than absolute value
    //precision defined by the step in the StepVector

    // calcul de rho pour l'individu i

    double khi = 0.0;
    vector<double> absDiff;
    vector<double> relDiff;

    for (unsigned int j=0;j<pop[i].getAllelesNumber();j++) {

        absDiff.push_back( abs( pop[i].getAllele(j) - pop[posMax].getAllele(j) ) );
        if ( abs(pop[posMax].getAllele(j)) > (problem->getStepVector()[j]/32.) )
        relDiff.push_back( abs( (pop[i].getAllele(j) - pop[posMax].getAllele(j)) / pop[posMax].getAllele(j) ) );
        else
        relDiff.push_back( 1.0 );

        if ( relDiff[j] > epsilon2 && absDiff[j] > (problem->getStepVector()[j]/32.) ) { khi += 1.0; }

    }

    khi /= double(pop[i].getAllelesNumber());

    //cout << "Best:[0]:" << pop[posMax].getAllele(0) << "- :[1]: " << pop[posMax].getAllele(1) << "- absDiffx: " << absDiff[0] << "- absDiffy: " << absDiff[1] << "- relDiffx: " << relDiff[0] << "- relDiffy: " << relDiff[1] << endl;
    //cout << "Individual " << i << ":" << pop[i].getIdNumber() << ":[0]: " << pop[i].getAllele(0) << "- :[1]: " << pop[i].getAllele(1) << " - diversity: " << khi << endl;
    return khi;

}

template <class T> double Heuristic<T>::diversityPopulation() {

        // calcul de rho pour tous les individus

        double rho = 0.0;

        for (unsigned int i=0;i<pop.size();i++) {

                rho += diversityIndividual(i);

        }

        rho /= double(pop.size() - 1) ;

        return rho;

}

// *********************************************************************************
// Class DE definition
// *********************************************************************************


class DE: public Heuristic<Individual> {

    private:

        bool migrated;

    protected:

        void selectSamples(int candidate,int *r1,int *r2=0,int *r3=0,int *r4=0,int *r5=0);

    public:

        DE(Problem *problem, unsigned int size);
        DE(Problem *problem, unsigned int size, vector< Individual > &popIndividualsVector);
        DE(Problem *problem, string filename);

        void crossover(double probability, double scale, int methodCrossover);
        void migration();
        bool getMigrated();

        void rand2Best(double probability, double scale, int candidate);
        void rand3(double probability, double scale, int candidate);

        void selection();
        void selectionRank();

        void step(double Cr, double F, int methodCrossover, bool migrate);

};

// *********************************************************************************
// Class ES definition
// *********************************************************************************

class ES: public Heuristic<IndividualES> {

    private:

        unsigned int mu, lambda;
        double sigmaF;
        vector<unsigned int> selectParents(unsigned int rho);

    public:

        ES(Problem *problem, unsigned int mu, unsigned int lambda, double sigma);
        ES(Problem *problem, unsigned int mu, unsigned int lambda, double sigma, vector< vector<double> > &allelesVector);
        ES(Problem *problem, string filename, unsigned int lambda);

        void globalIntermediateRecombination(double probability, unsigned int rho);
        void intermediateRecombination(double probability);
        void generalizedIntermediateRecombination(double probability);
        void uniformCrossoverRecombination(double probability);

        void mutationNormallyDistributed(double K);

        void selectionPlus ();
        void selectionComma();
        void selectionPlusRank ();
        void selectionCommaRank();

        void save(string filename);
        void report();

        void step(double K, int methodRecombination, int methodSelection);

};

// *********************************************************************************
// Class CSA-ES definition
// *********************************************************************************

class CSA_ES: public Heuristic<IndividualCSA_ES> {

    private:

        unsigned int mu, lambda;
        double N, chiN, cs, damp, sigmaF;
        vector<double> ps;

        vector< IndividualCSA_ES > bestIndividualsVector;

    protected:

        vector<int> selectParents(int rho);

    public:

        CSA_ES(Problem *problem, unsigned int lambda, double sigmaF);
        CSA_ES(Problem *problem, unsigned int lambda, double sigmaF, vector< vector<double> > &allelesVector);

        void globalIntermediateRecombination();

        void mutationNormallyDistributed();

        void selectionComma();
        void selectionCommaRank();

        void adaptSigma();

        void step();

        void getBestInPopulation( vector< Individual > &popIndividualsVector );

        double getSigmaF();

};

// *********************************************************************************
// Class CMA-ES definition
// *********************************************************************************

class CMA_ES: public Heuristic<IndividualCSA_ES> {

    private:

        unsigned int mu, lambda;
        double N, chiN, cs, damp, cc, muCov, cCov, sigmaF;

        vector<double> weight;
        double muEff;

        vector<double> ps, pc;

        vector<double> D;
        vector< vector<double> > B;
        vector< vector<double> > C;

        vector< IndividualCSA_ES > bestIndividualsVector;

    public:

        CMA_ES(Problem *problem, int lambda, double sigmaF);
        CMA_ES(Problem *problem, int lambda, double sigmaF, vector< vector<double> > &C0, vector< Individual > &popIndividualsVector);

        void globalWeightedIntermediateRecombination();
        void mutationNormallyDistributed();
        void selectionCommaRank();
        void adaptation();
        void step();

        void getBestInPopulation( vector< Individual > &popIndividualsVector );
        double getSigmaF();
        vector< vector<double> > getC();

};


class PSO: public Heuristic<IndividualPSO> {

    private:

        vector<IndividualPSO> bestEverPop;

        double c1, c2, lambda, kappa; // cognitiveAccel, socialAccel, constrictionGain


    public:

        PSO(Problem *problem, unsigned int size, double c1, double c2, double lambda, double kappa);

        void updateLocalBestParticles();
        IndividualPSO* getGlobalBestParticle(unsigned int i);
        void update();
        void selection();

        void step();

};


// NOTE: Hooke-Jeeves is not a heuristic, but I put it here for simplicity reasons
class Hooke_Jeeves: public Heuristic<Individual> {

    private:

        const unsigned int r, t;
        int s, m;
        vector<double> delta;

    public:

        Hooke_Jeeves(Problem *problem, Individual& xInit, unsigned int r, unsigned int s0, unsigned int t, unsigned int m);

        void globalSearch();
        void localSearch();
        void selection();
        void adaptation();

        void step();
        bool mExceeded() { return (m<0); }

};

#endif /* _HEURISTIC_H */
