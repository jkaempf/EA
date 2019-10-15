#ifndef _PROBLEM_H
#define _PROBLEM_H

using namespace std;

#include <cmath>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include <ctime>

#define round(x) (x<0?ceil((x)-0.5):floor((x)+0.5))
#define sign(x) x/abs(x)

const long double PI = 4.0*atan(1.0);

class Problem {

protected:

    unsigned int nEval;
    double timeEval; // time taken in the evaluations -> total cpu(s) time needed

    // constraints
    unsigned int nConstraints;

    // step size in case of real problems for rounding
    double stepSize;
    vector<double> stepVector;

    vector<double> maxVector, minVector;
    double minZ, maxZ;

public:

    // constructeur et destructeur
    Problem() { nEval = 0; timeEval = 0.0; stepSize = 0.0; nConstraints = 0; minZ = 0.0; maxZ = 0.0; };
    Problem(double stepSize): stepSize(abs(stepSize)) { nEval = 0; timeEval = 0.0; nConstraints = 0; minZ = 0.0; maxZ = 0.0; };
    Problem(vector<double> stepVector): stepVector(stepVector) { nEval = 0; timeEval = 0.0; nConstraints = 0; minZ = 0.0; maxZ = 0.0; };
    virtual ~Problem() {};

    virtual vector<double> evaluate(unsigned int id, vector<double> alleles) { vector<double> fitness; fitness.push_back(-numeric_limits<double>::max()); return fitness; }
    virtual vector<double> evaluate(vector<double> alleles) { vector<double> fitness = evaluate(0,alleles); nEval--; return fitness; }

    virtual double constraint(unsigned int k, vector<double> alleles) { return 0.0; }
    virtual bool notInField(vector<double> alleles);

    vector<double> getMaxVector() { return maxVector; };
    vector<double> getMinVector() { return minVector; };
    double getMinZ() { return minZ; };
    double getMaxZ() { return maxZ; };


    unsigned int getSize() { if (maxVector.size() == minVector.size()) return maxVector.size(); else throw(string("Error in min/max vector sizes.")); }

    unsigned int getnEval() { return nEval; };
    void resetnEval() { nEval = 0; };
    void incrementnEval() { nEval++; }
    double getTimeEval() { return timeEval; };
    unsigned int getnConstraints() { return nConstraints; };

    vector<double> getStepVector() { if (stepVector.size() == 0 && minVector.size() > 0) stepVector.assign(minVector.size(), stepSize); return stepVector; };

};

class Radiance: public Problem {

protected:

    // for Radiance like problems
    double maxDetectorArea;
    double w,l;
    unsigned int wmax, lmax;
    double hmax, hrmax;
    double ecartw, ecartl;
    string skyFile;
    string groundFile;
    string environmentFile;
    bool ground, environment, deleteFiles;

    double xinit, yinit, hbase;

public:

    Radiance():Problem() { maxDetectorArea = 4.0; ground=false; environment=false; deleteFiles=true; };
    virtual double volume(vector<double> alleles) { return 0.0; }
    virtual void runRadiance(string file);

};

class Problem0: public Radiance {

public:

    Problem0(string skyFilename);

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    bool notInField(vector<double> alleles);

};

class Problem1: public Radiance {

public:

    Problem1();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    bool notInField(vector<double> alleles, double *error);

};

class Problem2: public Radiance {

public:

    Problem2();

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Problem3: public Radiance {

protected:

    double alpha;

public:

    Problem3();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    bool notInField(vector<double> alleles, double *error=0);

};

class NYCity: public Radiance {

private:

    int floormax;
    double spacing;
    double volumeMax;

public:

    NYCity();
    vector<double> evaluate(unsigned int id, vector<double> alleles);
    double volume(vector<double> alleles);
    double constraint(unsigned int k, vector<double> alleles);

};

class Problem5: public Problem3 {

public:

    Problem5();

    vector<double> evaluate(int id, vector<double> alleles);

};

class Problem6: public Problem3 {

public:

    Problem6();

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class RandomTerracesFlatRoof: public Radiance {

    private:

        double volumeMax;

    public:

    RandomTerracesFlatRoof();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    double constraint(unsigned int k, vector<double> alleles);
    double volume(vector<double> alleles);

};

class RandomSlabsSlopedRoof: public Radiance {

    private:

        double volumeMax;

    public:

    RandomSlabsSlopedRoof();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    double constraint(unsigned int k, vector<double> alleles);
    double volume(vector<double> alleles);
    static double d(double dpara, double r1, double h2, double height);
    static double h(double dpara, double r1, double h2, double distance);
    double groundSurface(vector<double> alleles);

};

class RandomTerraceCourtsEW: public Radiance {

    private:

        vector<double> vw, vl;
        double volumeMax;

    public:

    RandomTerraceCourtsEW();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    double constraint(unsigned int k, vector<double> alleles);
    double volume(vector<double> alleles);
    static double d(double dpara, double r1, double h2, double height);
    static double h(double dpara, double r1, double h2, double distance);
    double groundSurface(vector<double> alleles);

};

class Pavilions: public Radiance {

    private:

        double volumeMax;

    public:

    Pavilions();

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class PepX: public Radiance {

    private:

    vector<double> vw, vl;

    public:

    PepX(string skyFile, string siteFile);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class PepY: public Radiance {

    private:

    vector<double> vw, vl;

    public:

    PepY(string skyFile, string siteFile);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class RandomPavilionsFlatRoof: public Radiance {

    public:

    RandomPavilionsFlatRoof();

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class ProjetHuang: public Radiance {

    public:

    ProjetHuang();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    double constraint(unsigned int k, vector<double> alleles);

};

class ProjetWagner: public Radiance {

    private:

        string varyingFile;

    public:

    ProjetWagner(string varyingFile);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class FourierBasis: public Radiance {

    private:

    unsigned int Nx,Ny;
    double h,Tx,Ty,minHorizCut,maxHorizCut,amplitude;
    double volumeMax;
    double maxGridSpacing;
    bool surroundings, deleteFiles;

    public:

    FourierBasis();

    double getA(int k, int l, vector<double> alleles);
    double getB(int k, int l, vector<double> alleles);

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    void runRadiance(string file);
    double constraint(unsigned int k, vector<double> alleles);
    double fourierSeries(double x, double y, vector<double> alleles);
    vector<double> derivativeFourierSeries(double x, double y, vector<double> alleles);
    double volumeFourierSeries(vector<double> alleles);
    double volume(vector<double> alleles);

};

class ProjetLorenzetti: public Problem {

    public:

    ProjetLorenzetti();

    vector<double> evaluate(unsigned int id, vector<double> alleles);
    void readACTrialResults(string filename, double &sumLogRmsError);

};

class EnergyPlus: public Problem {

    private:

    vector<string> labelVector, constraintVector, equalityVector, outputLabelVector;
	map<unsigned int, vector<string> > stringMap;
    vector<double> outputCoeffVector;
    string templateFile,weatherFile;

    public:

    EnergyPlus(string filename);
    double constraint(unsigned int k, vector<double> alleles);
    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class CitySim: public Problem {

    private:

    vector<string> labelVector, constraintVector, equalityVector, outputLabelVector;
	map<unsigned int, vector<string> > stringMap;
	vector<double> outputCoeffVector;
    string templateFile;

    public:

    CitySim(string filename);
    double constraint(unsigned int k, vector<double> alleles);
    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Ackley: public Problem {

public:

    Ackley();
    Ackley(unsigned int size);
    Ackley(unsigned int size, double stepSize);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Rastrigin: public Problem {

public:

    Rastrigin();
    Rastrigin(unsigned int size);
    Rastrigin(unsigned int size, double stepSize);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class RastriginConstrained: public Rastrigin {

  public:

    RastriginConstrained(unsigned int size);
    RastriginConstrained(unsigned int size, double stepSize);

    double constraint(unsigned int k, vector<double> alleles);
    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Rosenbrock: public Problem {

public:

    Rosenbrock();
    Rosenbrock(unsigned int size);
    Rosenbrock(unsigned int size, double stepSize);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Quadratic: public Problem {

public:

    Quadratic();
    Quadratic(unsigned int size);
    Quadratic(unsigned int size, double stepSize);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

class Michalewicz: public Problem {

  public:

    Michalewicz();
    Michalewicz(double stepSize);

    double constraint(unsigned int k, vector<double> alleles);

    vector<double> evaluate(unsigned int id, vector<double> alleles);

};

#endif /* _PROBLEM_H */
