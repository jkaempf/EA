#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <vector>
#include <sstream>
#include <stdint.h>

using namespace std;

const double TINY = 1e-5;

template <typename T> inline T to(const string &s) {

    T value;
    std::stringstream ss(s);
    ss >> value;

    return value;

}

template <typename T> inline string toString(const T &s) {

    std::stringstream ss;
    ss << s;

    return ss.str();

}

class Point {

private:

  vector<double> coordonnees;

public:

  Point(double x, double y);
  Point(double x, double y, double z);

  double getx();
  double gety();
  double getz();
  void setx(double x);
  void sety(double y);
  void setz(double z);
  int  getSize();
  void show();
  double operator[] (unsigned int n);
  void operator= (Point *lePoint);
  Point operator- (Point lePoint);

};

void grid(vector<Point> rectangle, vector<Point> &grille, double &nx, double &ny, double &nz, double spacing);
void gridTriangle(vector<Point> triangle, vector<Point> &grille, double &nx, double &ny, double &nz, double spacing);
void gridTriangle2(vector<Point> triangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz);
void gridRectangle(vector<Point> rectangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz);
void gridRectangle2(vector<Point> rectangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz);
void gridRectangle3(vector<Point> rectangle, double spacing, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz);
void removeVertically(Point base, vector<double> dx, vector<double> dy, vector<Point> &grille, vector<Point> &normales, vector<double> &aires);
void rotation2D(vector<Point> &vect, double phi);
void translation2D(vector<Point> &vect, double deltax, double deltay);
Point centroid2D(vector<Point> &forme);
Point centroid3D(vector<Point> &forme);
Point middlePoint(Point &point1, Point &point2);
Point baryCenter(Point &point1, Point &point2);
Point weightedCenter(Point point1, Point point2, double weight);
double norme(double vx, double vy, double vz);
double norm(vector<double> v);
void normalise(double &nx, double &ny, double &nz);
void crossProduct(double u1, double u2, double u3, double v1, double v2, double v3, double &n1, double &n2, double &n3);
void readResults(string filename, vector<double> &results);
void readResults(string filename, vector<string> &results);
void readResultsInLine(string filename, vector<double> &results);
void writeResults(string filename, ostringstream &oss);
void writeResultsOver(string filename, ostringstream &oss);
double normallyDistributedSPRNG_Sum();
double normallyDistributedSPRNG_BoxMuller();
void zigset(uint32_t jsrseed);
double normallyDistributedSPRNG_Ziggurat();
double randomUniform(double minValue,double maxValue);

void computeEigensystem(vector<double> &outValues, vector< vector<double> > &outVectors);
void tql2(vector<double> &d, vector<double> &e, vector< vector<double> > &V);
void tred2(vector<double> &d, vector<double> &e, vector< vector<double> > &V);

double equationParser(string& eq, vector<double> &alleles);
