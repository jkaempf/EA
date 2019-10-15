/* problem.cxx
   Class containing the problems solved :-)
   Jerome Kaempf
   LESO-PB / EPFL
   Idiap Research Institute
*/

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <map>
#include <numeric>

#include "geometry.h"
#include "problem.h"

// Run Radiance method

void Radiance::runRadiance(string file) {

    // paramètres de Radiance pour la simulation
    // AMBIENT: -ab 2 -ad 1024 -as 512 -aa 0
    // DIRECT : -dj 1 -ds 0.02 -dt 0 -dc 0.6 -dr 6 -dp 0

    // début de la simulation
    string leString;

    // preparation and start of OCONV
    leString = "oconv ";
    leString += "-w ";
    leString += file;
    leString += ".rad ";
    leString += skyFile;
    leString += " ";
    if (ground==true) {
        leString += groundFile;
        leString += " ";
    }
    if (environment==true) {
        leString += environmentFile;
        leString += " ";
    }
    leString += "> ";
    leString += file;
    leString += ".oct";

    // issue system command
    int errCode = system(leString.c_str());
    if ( errCode != 0 ) throw (string("Error in oconv (.oct)"));

    // preparation and start of RTRACE -I -ab 2 -ad 1000 -as 250 -dc .9 -dj .66

    leString = "rtrace ";
    leString += "-h -w -I -ab 2 -ad 1024 -as 512 -aa 0 ";
    //leString += "-h -w -I -ab 2 -ad 2048 -as 512 -aa 0.15 -ar 256 ";
    //leString += "-h -w -I -ab 2 -ad 256 -as 128 -aa 0.005 -ar 256 ";
    leString += file;
    leString += ".oct < ";
    leString += file;
    leString += ".inp > "; // white light from sky (1,1,1), multiplied by the coloured reflectance of the surfaces
    leString += file;
    leString += ".out";

    // issue system command
    errCode = system(leString.c_str());
    if ( errCode != 0 ) throw (string("Error in rtrace (.out)"));

}

// Problem::notInField returns true if outside the domain and the constraints

bool Problem::notInField(vector<double> alleles) {

    bool bad=false;

    // not in the boundaries
    for (unsigned int i=0; i < alleles.size(); i++) {

        bad |= ( (alleles[i] > maxVector[i]) || (alleles[i] < minVector[i]) );

        if ( (alleles[i] > maxVector[i]) || (alleles[i] < minVector[i]) ) cout << "allele[" << i << "] out of field" << endl;

    }

    // or not in the constraints
    for (unsigned int i=0; i < nConstraints; i++) {

        bad |= !(constraint(i,alleles) <= 0.0);

        if ( !(constraint(i,alleles) <= 0.0) ) cout << "constraint[" << i << "] not satisfied" << endl;

    }

    return bad;

}

// *************************************************************
// Zero problem: derived class Problem0 - single building
// that can rotate
// *************************************************************

Problem0::Problem0(string skyFilename):Radiance() {

    skyFile = skyFilename;
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;
    xinit = 11618.0; // 11637 centre of church, external 11813
    yinit = 68513.0; // 68569 center of church, external 68923
    hbase = 254.;

    maxDetectorArea = 1.0;

    maxVector.push_back( 180.0 ); // alphamax
    minVector.push_back(  0.0  ); // alphamin

    maxVector.push_back( 60.0 ); // wmax
    minVector.push_back(  0.0 );

    maxVector.push_back( 110.0 ); // lmax
    minVector.push_back(  0.0 );

    maxVector.push_back( 20.0  ); //hmax
    minVector.push_back(  0.0 );

    deleteFiles = false;

}

vector<double> Problem0::evaluate(unsigned int id, vector<double> alleles) {

    if (alleles.size() == 1) {

        //alleles[0]-=12.8;

        alleles.push_back(60.0);
        alleles.push_back(14.0);
        alleles.push_back(20.0);

        deleteFiles = false;

    }

    time_t start, end;
    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    // writing the alleles in a file

    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // variables needed in the method

    vector<Point> Rectangle;

 // void prepareRadiance() {

    Rectangle.clear();

    Rectangle.push_back(Point(alleles[1]/2., alleles[2]/2.));
    Rectangle.push_back(Point(-alleles[1]/2., alleles[2]/2.));
    Rectangle.push_back(Point(-alleles[1]/2., -alleles[2]/2.));
    Rectangle.push_back(Point(alleles[1]/2., -alleles[2]/2.));

    rotation2D(Rectangle, alleles[0] + 12.8); // 12.8° is the natural orientation of Matthaus

    translation2D(Rectangle, xinit, yinit);

//    writeToRad2();
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";

    // ouverture des fichiers

    fstream output1 (fileRad.c_str(), ios::binary | ios::out | ios::app );
    fstream output2 (fileInp.c_str(), ios::binary | ios::out | ios::app );
    fstream output3 (fileDet.c_str(), ios::binary | ios::out | ios::app );

    // test d'ouverture

    if (!output1.is_open()) throw (string("Error opening: " + fileRad));
    if (!output2.is_open()) throw (string("Error opening: " + fileInp));
    if (!output3.is_open()) throw (string("Error opening: " + fileDet));

    // précision de la sortie des points

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // creation du fichier .RAD
    // materiaux

    output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

    // batiments

    for (int i=0;i<4;i++) {

      output1 << "facade polygon wall" << i+1 << "\n0\n0\n12\n"
         << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
         << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
         << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[3]+hbase << "\n"
         << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[3]+hbase << "\n" << endl;

    }

    output1 << "facade polygon roof\n0\n0\n12\n";
    for (int i=0;i<4;i++) {
      output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[3]+hbase << "\n";
    }

    output1.close();

    // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;
        vector<double> detectorArea;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // points sur les 4 FACADES

        for (unsigned int id=0; id<4; id++) {

            facade.clear();
            grid1.clear();

            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+alleles[3] ) );
            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+alleles[3] ) );

            gridRectangle3(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // points sur le toit

        roof.clear();
        grid1.clear();

        for (unsigned int i=0; i<4;i++) {
            roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[3]+hbase) );
        }

        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

    // creation du fichier .inp
    // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        double gap = 0.05;

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

        output2.close();
        output3.close();

    // début de la simulation

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(string(file + ".out"), radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    cout << "totalSurface: " << totalSurface << endl;
    cout << "p1: " << 0.8*totalRadiation << "\tp2: " << Ubat*24*365*totalSurface*5.0 << endl;


    vector<double> fitness;
    fitness.assign(1,0.8*totalRadiation);

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    return fitness;

}

bool Problem0::notInField(vector<double> alleles) {

/*

    double alphac, lalpha, walpha;

    // vérification de l'angle

    if ( al == 0 ) return ( valeur > maxVector[0] );

    // vérification des paramètres

    if ( al == 2 ) {

        alphac = atan(maxVector[1]/maxVector[2])/PI*180.;  // alpha critical in degrees

        if (maxVector[0] <= alphac) {
            lalpha = maxVector[2] / cos(maxVector[0]*PI/180.);
        }
        else {
            lalpha = maxVector[1] / sin(maxVector[0]*PI/180.);
        }

        return (valeur > lalpha);

        walpha = min( (maxVector[1] - valeur*sin(maxVector[0]*PI/180.))/cos(maxVector[0]*PI/180.), (maxVector[2] - valeur*cos(maxVector[0]*PI/180.))/sin(maxVector[0]*PI/180.) );

        return (alleles[1] > walpha);

    }

    if ( al == 3 ) {

        return (valeur > ((maxVector[1]*maxVector[2]*maxVector[3]/10.0) / (alleles[1]*alleles[2])));

    }

*/

    return false;

}


// *************************************************************
// First problem: derived class Problem1 - many buildings
// that can have the height to vary
// *************************************************************


Problem1::Problem1():Radiance() {

    w = 10.0;
    l = 10.0;
    wmax=1;
    lmax=2;
    hmax=14.0;
    hrmax = 0.0;
    ecartw = (0.0/w);
    ecartl = (0.0/l);
    skyFile = "BaselSkyHH.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;
    xinit = 11637.0; // 11637 centre of church, external 11813
    yinit = 68569.0; // 68569 center of church, external 68923
    hbase = 254.;

    maxDetectorArea = 1.0;

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

    }

    deleteFiles = false;

}

vector<double> Problem1::evaluate(unsigned int id, vector<double> alleles) {

    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    // écriture du fichier des allèles

    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // variables utilisées

    vector<Point> Rectangle;
    vector<double> detectorArea;

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 12.8); // 12.8° is the natural orientation of Matthaus

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

//          void Problem::write3DShape(int hi, int hj) {

            string file3D = file + ".shape";

//          void Problem::writeRadianceInput(int hi, int hj) {

            string fileRad = file + ".rad";
            string fileInp = file + ".inp";
            string fileDet = file + ".det";

            // ouverture des fichiers

            fstream output1 (fileRad.c_str(), ios::binary | ios::out | ios::app );
            fstream output2 (fileInp.c_str(), ios::binary | ios::out | ios::app );
            fstream output3 (fileDet.c_str(), ios::binary | ios::out | ios::app );

            output1 << setprecision(12);
            output2 << setprecision(12);
            output3 << setprecision(12);

            // test d'ouverture

            if (!output1.is_open()) throw (string("Error opening: " + fileRad));
            if (!output2.is_open()) throw (string("Error opening: " + fileInp));
            if (!output3.is_open()) throw (string("Error opening: " + fileDet));

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            if ( alleles[hi*lmax+hj] > 0. ) {   // si hauteur dépasse zéro

                // POLYGON of the FACADES

                for (unsigned int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n" << endl;

                }

                // POLYGON of the FLAT ROOF

                output1 << "facade polygon wall5-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n";

                for (unsigned int i=0;i<4;i++) {
                    output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n";
                }

                output1 << endl;

            }

        output1.close();

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        if ( alleles[hi*lmax+hj] > 0.) { // si plus grand que zéro

            // points sur les 4 FACADES

            for (unsigned int id=0; id<4; id++) {

                double hrel = 0.;

                // si les batiments sont collés, seulement points visibles
                if ( ecartw == 0.0 ) {
                    if ( id == 1 ) {
                        if ( hi < (wmax-1) ) hrel = alleles[(hi+1)*lmax+hj];
                    }
                    if ( id == 3 ) {
                        //if ( hi > 0 ) hrel = alleles[(hi-1)*lmax + hj];
                        if ( hi > 0 ) hrel = alleles[hi*lmax + hj];
                    }
                }
                if ( ecartl == 0.0 ) {
                    if ( id == 0 ) {
                        //if ( hj > 0 ) hrel = alleles[hi*lmax*2+(hj-1)*2];
                        if ( hj > 0 ) hrel = alleles[hi*lmax+hj];
                    }
                    if ( id == 2 ) {
                        if ( hj < (lmax-1) ) hrel = alleles[hi*lmax + (hj+1)];
                    }
                }

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+alleles[hi*lmax+hj] ) );
                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+alleles[hi*lmax+hj] ) );

                gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

       // FLAT ROOF, no roof elevation

            roof.clear();
            grid1.clear();

            for (unsigned int i=0; i<4;i++) {
                roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax+hj]+hbase) );
            }

            gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // creation du fichier .inp
        // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        double gap = 0.05;

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

        output2.close();
        output3.close();

      } // boucle hj
    }   // boucle hi

    // début de la simulation

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(string(file + ".out"), radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    cout << "totalSurface: " << totalSurface << endl;
    cout << "p1: " << 0.8*totalRadiation << "\tp2: " << Ubat*24*365*totalSurface*5.0 << endl;


    vector<double> fitness;
    fitness.assign(1,0.8*totalRadiation);

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    return fitness;

}

bool Problem1::notInField(vector<double> alleles, double *error)
{

    *error = 0.0;

    double mean=0.0;

    // compute mean value of height

    for (unsigned int i=0; i < alleles.size(); i++) {

        mean += alleles[i];
    }

    mean /= alleles.size();

    // compute difference to reference mean value

    *error = abs(mean - 9.0);

    if ( *error > TINY ) {
        return true;
    }
    else return false;

}

// *************************************************************
// Second problem: derived class Problem2 from Problem
// *************************************************************


Problem2::Problem2():Radiance() {

    w = 10.0;
    l = 10.0;
    wmax=4;
    lmax=4;
    hmax=18.0;
    hrmax = 3.0;
    ecartw = (0.0/w);
    ecartl = (0.0/l);
    skyFile = "BaselSky.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;
    xinit = 11637.0; // 11637 centre of church, external 11813
    yinit = 68569.0; // 68569 center of church, external 68923
    hbase = 254.;

    maxDetectorArea = 4.0;

    maxVector.clear();

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

            maxVector.push_back( 1.0  );
            minVector.push_back(  0.0 );

    }

}

vector<double> Problem2::evaluate(unsigned int id, vector<double> alleles) {

    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 12.8); // 12.8° is the natural orientation of Matthaus

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

//          void Problem::write3DShape(int hi, int hj) {

            string file3D = file + ".shape";

            // ouverture des fichiers

            fstream output (file3D.c_str(), ios::binary | ios::out | ios::app );

            // test d'ouverture

            if (!output.is_open())
            {
            cerr << "Error opening: " << file << endl;
            }

                // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << endl;

            }

            //    cout << "maximum: " << wmax*lmax << endl;
            //    cout << "acces du vecteur: " << hi*lmax+hj << "\tsa valeur: " << alleles[hi*lmax+hj] << endl;

            // calcul de orient, hflat et hprime

            int orient=0;
            double hprime = 0.;
            double hflat = hbase+alleles[hi*lmax*2+hj*2];

            hprime = 4*alleles[hi*lmax*2+hj*2+1] - TINY;
            orient = (int) trunc(hprime);

            hprime -= 1.0*orient;
            hprime *= hrmax;

            if (hrmax > 0.) {

                // SLOPED ROOF PART

                output << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << endl;

                output << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << endl;

                output << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << endl;

                output << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime
                        << endl;

            }
            else {
                for (int i=0;i<4;i++) {
                    output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n";
                }
            }

            output.close();

//          void Problem::writeRadianceInput(int hi, int hj) {

            string fileRad = file + ".rad";
            string fileInp = file + ".inp";
            string fileDet = file + ".det";

            // ouverture des fichiers

            fstream output1 (fileRad.c_str(), ios::binary | ios::out | ios::app );
            fstream output2 (fileInp.c_str(), ios::binary | ios::out | ios::app );
            fstream output3 (fileDet.c_str(), ios::binary | ios::out | ios::app );

            // test d'ouverture

            if (!output1.is_open()) throw (string("Error opening: " + fileRad));
            if (!output2.is_open()) throw (string("Error opening: " + fileInp));
            if (!output3.is_open()) throw (string("Error opening: " + fileDet));

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            if ( alleles[hi*lmax*2+hj*2] > 0. ) {   // si hauteur dépasse zéro

                // POLYGON of the FACADES

                for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n" << endl;

                }
            }

            // POLYGON(S) of the ROOF

            if (hrmax > 0.) {

                // SLOPED ROOF PART

                output1 << "facade polygon roof1-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof2-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                        << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof3-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                        << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof4-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << endl;

            }
            else {

                // POLYGON of the FLAT ROOF

                output1 << "facade polygon wall5-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n";

                for (int i=0;i<4;i++) {
                    output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n";
                }

                output1 << endl;

            }

        output1.close();

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;
        vector<double> detectorArea;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        if ( alleles[hi*lmax*2+hj*2] > 0.) { // si plus grand que zéro

            // points sur les 4 FACADES

            for (int id=0; id<4; id++) {

                double hrel = 0.;

                // si les batiments sont collés, seulement points visibles
                if ( ecartw == 0.0 ) {
                    if ( id == 1 ) {
                        if ( hi < (wmax-1) ) hrel = alleles[(hi+1)*lmax*2+hj*2];
                    }
                    if ( id == 3 ) {
                        //if ( hi > 0 ) hrel = alleles[(hi-1)*lmax*2 + hj*2];
                        if ( hi > 0 ) hrel = alleles[hi*lmax*2 + hj*2];
                    }
                }
                if ( ecartl == 0.0 ) {
                    if ( id == 0 ) {
                        //if ( hj > 0 ) hrel = alleles[hi*lmax*2+(hj-1)*2];
                        if ( hj > 0 ) hrel = alleles[hi*lmax*2+hj*2];
                    }
                    if ( id == 2 ) {
                        if ( hj < (lmax-1) ) hrel = alleles[hi*lmax*2 + (hj+1)*2];
                    }
                }

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+alleles[hi*lmax*2+hj*2] ) );
                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+alleles[hi*lmax*2+hj*2] ) );

                gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }


        }

        if ( hrmax > 0. ) {

            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[(orient+1)%4][0], Rectangle[(orient+1)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+2)%4][0], Rectangle[(orient+2)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );

            gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );

            gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+1)%4][0], Rectangle[(orient+1)%4][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );

            gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[(orient+2)%4][0], Rectangle[(orient+2)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );

            gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }
        else { // FLAT ROOF, no roof elevation

            roof.clear();
            grid1.clear();

            for (int i=0; i<4;i++) {
                roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*2+hj*2]+hbase) );
            }

            gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // creation du fichier .inp
        // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*0.1 << "\t" << detectorPoint[i][1] + normalVector[i][1]*0.1 << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*0.1 << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

        output2.close();
        output3.close();


      } // boucle hj
    }   // boucle hi

    // début de la simulation

    string fileOut = file + ".out";
    string fileDet = file + ".det";

    string leString;

    // preparation and launch of Radiance

    runRadiance(file);

    // creation des vecteurs résultats

    vector<double> radianceResults, detectorAreas;

    // lecture des résultats

    radianceResults.clear(); // de RADIANCE
    readResults(fileOut, radianceResults);

    detectorAreas.clear(); // de la préparation des surfaces
    readResults(fileDet, detectorAreas);

    // détermination de la fitness de l'individu par sommation pondérée

    vector<double> fitness;
    fitness.assign(1,0.);

    for (unsigned int i=0;i<radianceResults.size();i++) {

        fitness[0]+=radianceResults[i]*detectorAreas[i];

    }

    // effacement des fichier .oct

    string fileOct = file + ".oct";

    remove(fileOct.c_str());

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    return fitness;

}

// *************************************************************
// Third problem: derived class Problem3
// *************************************************************


Problem3::Problem3():Radiance() {

    w = 10.0;
    l = 10.0;
    wmax=4;
    lmax=4;
    hmax=18.0;
    hrmax = 0.0;
    ecartw = (0.0/w);
    ecartl = (0.0/l);
    skyFile = "BaselSky.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;
    xinit = 11637.0; // 11637 centre of church, external 11813
    yinit = 68569.0; // 68569 center of church, external 68923
    hbase = 254.;

    alpha = 0.;

    maxDetectorArea = 4.0;

    for (unsigned int i=0; i < (wmax*lmax*5) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

    }
}

vector<double> Problem3::evaluate(unsigned int id, vector<double> alleles) {

    double penalty;

    // out of bounds, return very bad fitness value
    if (notInField(alleles, &penalty)) { vector<double> fitness; fitness.assign(1,-exp(penalty)); return fitness; }

    // start of evaluation
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;

    vector<double> detectorArea;

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 12.8 + alpha); // 12.8° is the natural orientation of Matthaus

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

//          void Problem::write3DShape(int hi, int hj) {

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << endl;

            }

            // calcul du centre du bâtiment

            Point centre = centroid2D(Rectangle);

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                    << centre[0] << "\t" << centre[1] << "\t" << alleles[hi*lmax*5+hj*5+4]+hbase << "\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << endl;

            }

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n" << endl;

            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {
                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                    << centre[0] << "\t" << centre[1] << "\t" << alleles[hi*lmax*5+hj*5+4]+hbase << "\n" << endl;
            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décomposition des facades

            for (int i=0;i<4;i++) {

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase ) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*5+hj*5+i]+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }



            // décompostion des toits

            for (int i=0;i<4;i++) {

                roof.clear();
                grid1.clear();

                roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*5+hj*5+i]+hbase) );
                roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase) );
                roof.push_back( Point(centre[0], centre[1], alleles[hi*lmax*5+hj*5+4]+hbase) );

                gridTriangle2(roof, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }



            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*0.1 << "\t" << detectorPoint[i][1] + normalVector[i][1]*0.1 << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*0.1 << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(file3D,  output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    vector<double> radianceResults;

    // lecture des résultats
    radianceResults.clear(); // de RADIANCE
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    vector<double> fitness;
    fitness.assign(1,0.);

    for (unsigned int i=0;i<radianceResults.size();i++) {

        fitness[0]+=radianceResults[i]*detectorArea[i];

    }

    // effacement des fichier octree Radiance

    string fileOct = file + ".oct";
    remove(fileOct.c_str());

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size() ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

bool Problem3::notInField(vector<double> alleles, double *error)
{

    *error = 0.0;

    bool bad=false;

    for (unsigned int i=0; i < alleles.size(); i++) {

        if ( ((i+1)%5) == 0 ) { // hauteur médiane qui est la plus grande

//            cout << "allele: " << num << endl;

            if ( (*max_element( alleles.begin() + i - 4, alleles.begin() + i ) > alleles[i]) ) {
//                cout << "hors champ - max: " << *max_element( alleles.begin() + num - 4, alleles.begin() + num ) << "\tvaleur: " << valeur << endl;
//                cout << "hmax: " << hmax << endl;
                *error += *max_element( alleles.begin() + i - 4, alleles.begin() + i ) - alleles[i];
                bad |= true;
            }
            else {
//                cout << "dans le champ!" << *max_element( alleles.begin() + num - 4, alleles.begin() + num ) << "\tvaleur: " << valeur << endl;
//                cout << "hmax: " << hmax << endl;
                bad |= false;
            }
        }
    }

    return bad;

}

// *************************************************************
// Fourth problem: derived class Problem4
// *************************************************************

NYCity::NYCity():Radiance() {

    // on a des tours de 369m, ce qui fait 123 étages max
    // nous allons parler en étages, ce sera plus simple et plus réaliste

    w = 60.0;
    l = 40.0;
    // les batiments font 40m*60m

    wmax=5;
    lmax=5;
    // ils sont sur un grille de 2*2

    floormax = 123;
    // 123 étages max

    ecartw = (20.0/w);
    ecartl = (20.0/l);
    // écart entre les batiments 20m et 20m, ce sont les rues

    skyFile = "BaselSky.rad";
    // ciel annuel, on peut mettre du PV avec des stores PV en facade

    groundFile = "FlatSite.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = false;

    xinit = 0.5*(double(wmax)*w+double(wmax-1)*(ecartw*w));
    yinit = 0.5*(double(lmax)*l+double(lmax-1)*(ecartl*l));
    hbase = 0.0;

    //maxDetectorArea = 10.0;
    // every 10m^2 a detector
    spacing = 3;
    // the distance in x,y on the surface is 1m at most

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( double(floormax) );
            minVector.push_back(  0.0 );

    }

    // constraints in volume... 80% of max volume
    nConstraints = 2;
    volumeMax = 0.5*wmax*lmax*w*l*floormax*3.0;

    //deleteFiles
    deleteFiles = false;

    cout << "\nProblem: NYCity" << endl;

}

vector<double> NYCity::evaluate(unsigned int id, vector<double> alleles) {

    time_t start, end;
    time(&start); // start timer

    // fichiers pour enregister
    ostringstream ossFile;
    ossFile << "building" << id;
    string file = ossFile.str();

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    // writing the alleles in a file
    ofstream outputAlleles(string(file + ".all").c_str(), ios::trunc);
    outputAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        outputAlleles << alleles[i] << "\n";
    }
    outputAlleles.flush();
    outputAlleles.close();
    //writeResultsOver(string(file + ".all"),ossAlleles);

    // choix de la précision (12 digits en tout - avant et après la virgule)
    ofstream output1(fileRad.c_str(), ios::trunc);
    ofstream output2(fileInp.c_str(), ios::trunc);
    ofstream output3(fileDet.c_str(), ios::trunc);
    //ostringstream output1, output2, output3;

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // rectangle de base pour chaque batiment
    vector<Point> Rectangle;
    vector<double> detectorArea;

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            if ( round(alleles[hi*lmax+hj]) > 0 ) {   // si nombre d'étages dépasse zéro

                // POLYGON of the FACADES

                for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << 3.0*round(alleles[hi*lmax+hj])+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << 3.0*round(alleles[hi*lmax+hj])+hbase << "\n" << endl;

                }

                // POLYGON of the FLAT ROOF

                output1 << "facade polygon wall5-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n";

                for (int i=0;i<4;i++) {
                    output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << 3.0*round(alleles[hi*lmax+hj])+hbase << "\n";
                }

                output1 << endl;

            }

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        if ( round(alleles[hi*lmax+hj]) > 0.0 ) { // si nombre d'étages > 0

            // points sur les 4 FACADES

            for (unsigned int id=0; id<4; id++) {

                double hrel = 0.;

                // si les batiments sont collés, seulement points visibles
                if ( ecartw == 0.0 ) {
                    if ( id == 1 ) {
                        if ( hi < (wmax-1) ) hrel = 3.0*round(alleles[(hi+1)*lmax+hj]);
                    }
                    if ( id == 3 ) {
                        //if ( hi > 0 ) hrel = alleles[(hi-1)*lmax + hj];
                        if ( hi > 0 ) hrel = 3.0*round(alleles[hi*lmax + hj]);
                    }
                }
                if ( ecartl == 0.0 ) {
                    if ( id == 0 ) {
                        //if ( hj > 0 ) hrel = alleles[hi*lmax*2+(hj-1)*2];
                        if ( hj > 0 ) hrel = 3.0*round(alleles[hi*lmax+hj]);
                    }
                    if ( id == 2 ) {
                        if ( hj < (lmax-1) ) hrel = 3.0*round(alleles[hi*lmax + (hj+1)]);
                    }
                }

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hrel ) );
                facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+3.0*round(alleles[hi*lmax+hj]) ) );
                facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+3.0*round(alleles[hi*lmax+hj]) ) );

                gridRectangle3(facade, spacing, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

       // FLAT ROOF, no roof elevation

            roof.clear();
            grid1.clear();

            for (int i=0; i<4;i++) {
                roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], 3.0*round(alleles[hi*lmax+hj])+hbase) );
            }

            gridRectangle3(roof, spacing, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // creation du fichier .inp
        // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        double gap = 0.05;

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t"
                    << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*gap << "\t"
                    << normalVector[i][0] << "\t"
                    << normalVector[i][1] << "\t"
                    << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

      } // boucle hj
    }   // boucle hi

    // flush des stringstreams
    output1.flush();
    output2.flush();
    output3.flush();

    // écriture dans les fichiers
    //writeResultsOver(fileRad, output1);
    //writeResultsOver(fileInp, output2);
    //writeResultsOver(fileDet, output3); // sortie pour vérifications

    output1.close();
    output2.close();
    output3.close();

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    vector<double> fitness;
    fitness.assign(1,totalRadiation);

    // writing the fitness in a file
    ofstream outputFitness(string(file + ".fit").c_str(), ios::trunc);
    outputFitness << setprecision(12);
    outputFitness << fitness[0] << "\n";
    outputFitness.flush();
    outputFitness.close();
    //writeResultsOver(string(file + ".fit"),ossFitness);

    // effacement des fichiers
    if (deleteFiles == true) {
        remove(string(file + ".oct").c_str());
        remove(fileRad.c_str());
        remove(fileInp.c_str());
        remove(fileDet.c_str());
        remove(fileOut.c_str());
    }

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

double NYCity::volume(vector<double> alleles) {

    return w*l*3.0*std::accumulate(alleles.begin(), alleles.end(), 0.0);

}

double NYCity::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        return volume(alleles) - 1.1*volumeMax;

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        return -volume(alleles) + 0.9*volumeMax;

    }
    else return 0.0;

}

// *************************************************************
// Fifth problem: derived class Problem5
// *************************************************************

Problem5::Problem5():Problem3() {

    w = 10.0;
    l = 10.0;
    wmax=1;
    lmax=1;
    hmax=18.0;
    hrmax = 0.0;
    ecartw = (0.0/w);
    ecartl = (0.0/l);
    skyFile = "BaselSky.rad";
    groundFile = "";
    environmentFile = "environment.rad";
    environment = true;
    xinit = 5.0; // 11637 centre of church, external 11813
    yinit = 5.0; // 68569 center of church, external 68923
    hbase = 0.;

    alpha = -12.8;

    maxDetectorArea = 1.0;

    maxVector.clear();

    for (unsigned int i=0; i < (wmax*lmax*5) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

    }
}

vector<double> Problem5::evaluate(int id, vector<double> alleles) {

    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;

    vector<double> detectorArea;

    vector<int> pointsFacades;

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 12.8 + alpha); // 12.8° is the natural orientation of Matthaus

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

//          void Problem::write3DShape(int hi, int hj) {

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << endl;

            }

            // calcul du centre du bâtiment

            Point centre = centroid2D(Rectangle);

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                    << centre[0] << "\t" << centre[1] << "\t" << alleles[hi*lmax*5+hj*5+4]+hbase << "\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << endl;

            }

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n" << endl;

            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {
                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*5+hj*5+i]+hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase << "\n"
                    << centre[0] << "\t" << centre[1] << "\t" << alleles[hi*lmax*5+hj*5+4]+hbase << "\n" << endl;
            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décomposition des facades

            for (int i=0;i<4;i++) {

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                pointsFacades.insert(pointsFacades.end(), grid1.size(), 1);
                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase ) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*5+hj*5+i]+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                pointsFacades.insert(pointsFacades.end(), grid1.size(), 1);

            }



            // décompostion des toits

            for (int i=0;i<4;i++) {

                roof.clear();
                grid1.clear();

                roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*5+hj*5+i]+hbase) );
                roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*5+hj*5+(i+1)%4]+hbase) );
                roof.push_back( Point(centre[0], centre[1], alleles[hi*lmax*5+hj*5+4]+hbase) );

                gridTriangle2(roof, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                pointsFacades.insert(pointsFacades.end(), grid1.size(), 0);

            }



            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*0.1 << "\t" << detectorPoint[i][1] + normalVector[i][1]*0.1 << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*0.1 << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    vector<double> radianceResults;

    // lecture des résultats
    radianceResults.clear(); // de RADIANCE
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    vector<double> fitness;
    double fitness1=0., fitness2=0.;

    for (unsigned int i=0;i<radianceResults.size();i++) {

        // fitness1: facades
        if ( pointsFacades[i] == 1 ) fitness1+=radianceResults[i]*detectorArea[i];
        // fitness2: toits
        else fitness2+=radianceResults[i]*detectorArea[i];

    }

    fitness.push_back(fitness1);
    fitness.push_back(fitness2);

    // effacement des fichier octree Radiance

    string fileOct = file + ".oct";
    remove(fileOct.c_str());

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size() ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// Sixth problem: derived class Problem
// *************************************************************

Problem6::Problem6():Problem3() {

    w = 10.0;
    l = 10.0;
    wmax=1;
    lmax=1;
    hmax  = 18.0;
    hrmax = 24.0;
    ecartw = (0.0/w);
    ecartl = (0.0/l);
    skyFile = "BaselSky.rad";
    groundFile = "";
    environmentFile = "environment.rad";
    environment = true;
    xinit = 5.0; // 11637 centre of church, external 11813
    yinit = 5.0; // 68569 center of church, external 68923
    hbase = 0.;

    alpha = -12.8;

    maxDetectorArea = 1.0;

    maxVector.clear();
    minVector.clear();

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );  //h1
            minVector.push_back(  0.0 );

            maxVector.push_back( hrmax ); //h2
            minVector.push_back(  hmax );

            maxVector.push_back( 0.0 ); //h3
            minVector.push_back( 1.0 );

            maxVector.push_back( 0.0 ); //h4
            minVector.push_back( 1.0 );

            maxVector.push_back( 0.0 ); //h5
            minVector.push_back( 1.0 );

            maxVector.push_back( 0.0 ); //h6
            minVector.push_back( 1.0 );

            maxVector.push_back( 0.0 ); //x7
            minVector.push_back( w );

            maxVector.push_back( 0.0 ); //x8
            minVector.push_back( l );

    }
}

vector<double> Problem6::evaluate(unsigned int id, vector<double> alleles) {

/*    double penalty;

    // out of bounds, return very bad fitness value
    if (notInField(alleles, &penalty)) return -exp(penalty);
*/

    // start of evaluation
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;
    vector<Point> pN;

    vector<double> detectorArea;

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

 // void prepareRadiance() {

    double h1,h2,x7,x8;
    vector<double> hN;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            h1 = alleles[hi*lmax*8+hj*8];
            h2 = alleles[hi*lmax*8+hj*8+1];
            hN.clear();
            hN.push_back( (1.0-alleles[hi*lmax*8+hj*8+2])*h1 + alleles[hi*lmax*8+hj*8+2]*h2 );
            hN.push_back( (1.0-alleles[hi*lmax*8+hj*8+3])*h1 + alleles[hi*lmax*8+hj*8+3]*h2 );
            hN.push_back( (1.0-alleles[hi*lmax*8+hj*8+4])*h1 + alleles[hi*lmax*8+hj*8+4]*h2 );
            hN.push_back( (1.0-alleles[hi*lmax*8+hj*8+5])*h1 + alleles[hi*lmax*8+hj*8+5]*h2 );
            x7 = alleles[hi*lmax*8+hj*8+6];
            x8 = alleles[hi*lmax*8+hj*8+7];

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 12.8 + alpha); // 12.8° is the natural orientation of Matthaus

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

            // calcul des points centraux...

            pN.clear();

            pN.push_back( weightedCenter(Rectangle[0], Rectangle[1], x7) );
            pN.push_back( weightedCenter(Rectangle[1], Rectangle[2], x8) );
            pN.push_back( weightedCenter(Rectangle[3], Rectangle[2], x7) );
            pN.push_back( weightedCenter(Rectangle[0], Rectangle[3], x8) );

            Point centre = weightedCenter(pN[0], pN[2], x8);

//          void Problem::write3DShape(int hi, int hj) {

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

            }

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << "\n"
                       << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;


                output << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << endl;

            }

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;
            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n9\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << endl;

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n9\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << endl;

            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décomposition des facades

            for (int i=0;i<4;i++) {

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase ) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase ) );
                facade.push_back( Point( pN[i][0], pN[i][1], hN[i]+hbase ) );

                gridTriangle2(facade, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );


            }

            // décompostion des toits

            for (int i=0;i<4;i++) {

                roof.clear();
                grid1.clear();

                roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                roof.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                roof.push_back( Point(centre[0], centre[1], h2+hbase) );

                gridTriangle2(roof, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                roof.clear();
                grid1.clear();

                roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                roof.push_back( Point(centre[0], centre[1], h2+hbase) );
                roof.push_back( Point(pN[i][0], pN[i][1], hN[i] + hbase) );

                gridTriangle2(roof, 2.0, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*0.1 << "\t" << detectorPoint[i][1] + normalVector[i][1]*0.1 << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*0.1 << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    vector<double> radianceResults;

    // lecture des résultats
    radianceResults.clear(); // de RADIANCE
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    vector<double> fitness;
    fitness.assign(1,0.);

    for (unsigned int i=0;i<radianceResults.size();i++) {

        fitness[0]+=radianceResults[i]*detectorArea[i];

    }

    // effacement des fichier octree Radiance

    string fileOct = file + ".oct";
    remove(fileOct.c_str());

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size() ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// CISBAT07 Generic Urban Forms Study: Terraces
// *************************************************************

RandomTerracesFlatRoof::RandomTerracesFlatRoof():Radiance() {

    w = 13.04;
    l = 10.0;

    wmax=5;
    lmax=5;

    hmax=14.0;

    ecartw = (0.0/w);
    ecartl = (14.59/l);

    skyFile = "BaselSkyHH.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;

    xinit = 11642.529;
    yinit = 68570.9014;
    hbase = 254.71;

    maxDetectorArea = 1.0;

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

    }

    // constraints in volume... 60% of max volume
    nConstraints = 2;
    volumeMax = 0.6*wmax*lmax*w*l*hmax;

    // delete intermediate files ?
    deleteFiles = false;

    cout << "\nProblem: RandomTerracesFlatRoof" << endl;

}

vector<double> RandomTerracesFlatRoof::evaluate(unsigned int id, vector<double> alleles) {

//    // chop precision in alleles for compatibility between MOO and our optimisation software
//    for (unsigned int i=0; i<alleles.size(); i++) {
//
//        if (alleles[i] != 0.0) {
//            double power10 = floor(log10(alleles[i]));
//            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
//        }
//        //cout << setprecision(12) << alleles[i] << endl;
//
//    }

    // timer start
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // création de la surface de base

    vector<Point> Rectangle;
    vector<double> detectorArea; // vecteur contenant les aires des surfaces de détection

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

        Rectangle.clear();

        Rectangle.push_back(Point(w/2., l/2.));
        Rectangle.push_back(Point(-w/2., l/2.));
        Rectangle.push_back(Point(-w/2., -l/2.));
        Rectangle.push_back(Point(w/2., -l/2.));

        rotation2D(Rectangle, 10.0); // 10.0 taken from Marylène's email

        translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

//          void Problem::writeRadianceInput(int hi, int hj) {

        // creation du fichier .RAD
        // materiaux

        if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

        // batiment (hi, hj)

        if ( alleles[hi*lmax+hj] > 0. ) {   // si hauteur dépasse zéro, alors les facades ont un sens

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n" << endl;

            }

        }

        // POLYGON of the FLAT ROOF

        output1 << "facade polygon wall5-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n";

        for (int i=0;i<4;i++) {
            output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax+hj]+hbase << "\n";
        }

        output1 << endl;

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // points sur les 4 FACADES

        for (unsigned int id=0; id<4; id++) {

            double hrel = 0.;

            // si les batiments sont collés, seulement points visibles
            if ( ecartw == 0.0 ) {
                if ( id == 1 ) {
                    if ( hi < (wmax-1) ) hrel = alleles[(hi+1)*lmax+hj];
                }
                if ( id == 3 ) {
                    //if ( hi > 0 ) hrel = alleles[(hi-1)*lmax + hj];
                    if ( hi > 0 ) hrel = alleles[hi*lmax + hj];
                }
            }
            if ( ecartl == 0.0 ) {
                if ( id == 0 ) {
                    //if ( hj > 0 ) hrel = alleles[hi*lmax*2+(hj-1)*2];
                    if ( hj > 0 ) hrel = alleles[hi*lmax+hj];
                }
                if ( id == 2 ) {
                    if ( hj < (lmax-1) ) hrel = alleles[hi*lmax + (hj+1)];
                }
            }

            facade.clear();
            grid1.clear();

            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hrel ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hrel ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+alleles[hi*lmax+hj] ) );
            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+alleles[hi*lmax+hj] ) );

            gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

       // FLAT ROOF, no roof elevation

        roof.clear();
        grid1.clear();

        for (int i=0; i<4;i++) {
            roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax+hj]+hbase) );
        }

        gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        // creation du fichier .inp
        // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        double gap = 0.05; // gap entre la paroi et le point de mesure (pour des raisons de précision)

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

      } // boucle hj
    }   // boucle hi

    // écriture dans les fichiers

    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    // début de la simulation
    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    //cout << "totalSurface in contact with exterior: " << totalSurface << endl;
    //cout << "totalIrradiation: " << 0.8*totalRadiation << "\tthermalLosses: " << Ubat*24*180*totalSurface*15.0 << endl;

    vector<double> fitness;
    fitness.push_back(0.8*totalRadiation - Ubat*24*180*totalSurface*15.0);
    fitness.push_back(volume(alleles));
    fitness.push_back(0.8*totalRadiation);
    fitness.push_back(Ubat*24*180*totalSurface*15.0);
    fitness.push_back(totalSurface);

    // calcul de la surface au sol
    double groundSurface = 0.0;
    for (vector<double>::iterator it=alleles.begin();it!=alleles.end(); it++)
        if ((*it) > 0.0) groundSurface += w*l;

    fitness.push_back(groundSurface);

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << "Number of measuring points: " << detectorArea.size() << "\tAverage points per second: " << detectorArea.size()/evalTime << endl;
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

double RandomTerracesFlatRoof::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        return volume(alleles) - 1.1*volumeMax;

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        return -volume(alleles) + 0.9*volumeMax;

    }
    else return 0.0;

}

double RandomTerracesFlatRoof::volume(vector<double> alleles) {

    return w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0);

}

// *************************************************************
// CISBAT07 Generic Urban Forms Study: Slabs
// *************************************************************

RandomSlabsSlopedRoof::RandomSlabsSlopedRoof():Radiance() {

    w = 10.0;
    l = 12.04;

    wmax=3;
    lmax=9;

    hmax=14.0;
    hrmax=4.0; // difference between 18m and 14m

    ecartw = (17.61/w);
    ecartl = (0.0/l);

    skyFile = "BaselSkyHH.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;

    xinit = 11642.529;
    yinit = 68570.9014;
    hbase = 254.71;

    maxDetectorArea = 1.0;

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

            maxVector.push_back( 1.0 );
            minVector.push_back( 0.0 );

    }

    // constraints in volume... 60% of max volume
    nConstraints = 2;
    volumeMax = 0.6*wmax*lmax*w*l*(hmax+hrmax/2.0); // OK

    // delete intermediate files ?
    deleteFiles = false;

    cout << "\nProblem: RandomSlabsSlopedRoof" << endl;

}

vector<double> RandomSlabsSlopedRoof::evaluate(unsigned int id, vector<double> alleles) {

//    // chop precision in alleles for compatibility between MOO and our optimisation software
//    for (unsigned int i=0; i<alleles.size(); i++) {
//
//        if (alleles[i] != 0.0) {
//            double power10 = floor(log10(alleles[i]));
//            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
//        }
//        //cout << setprecision(12) << alleles[i] << endl;
//
//    }

    // timer start
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // création de la surface de base

    vector<Point> Rectangle;
    vector<double> detectorArea; // vecteur contenant les aires des surfaces de détection

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 10.0); // 10° from Marylène information

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

            // calcul de orient, hflat et hprime

            unsigned int orient=0;
            double hprime = 0.;
            double hflat = hbase+alleles[hi*lmax*2+hj*2];

            hprime = max(4*alleles[hi*lmax*2+hj*2+1] - TINY, 0.0);
            orient = (unsigned int) trunc(hprime);

            hprime -= 1.0*orient;
            hprime *= hrmax;

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            if ( alleles[hi*lmax*2+hj*2] > 0. ) {   // si hauteur dépasse zéro

                // POLYGON of the FACADES

                for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n" << endl;

                }
            }

            // POLYGON(S) of the ROOF

            if (hrmax > 0.) {

                // SLOPED ROOF PART

                output1 << "facade polygon roof1-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof2-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                        << Rectangle[(orient+1)%4][0] << "\t" << Rectangle[(orient+1)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof3-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                        << Rectangle[(orient+2)%4][0] << "\t" << Rectangle[(orient+2)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n" << endl;

                output1 << "facade polygon roof4-" << orient << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat << "\n"
                        << Rectangle[orient][0] << "\t" << Rectangle[orient][1] << "\t" << hflat+hprime << "\n"
                        << Rectangle[(orient+3)%4][0] << "\t" << Rectangle[(orient+3)%4][1] << "\t" << hflat+hprime << "\n"
                        << endl;

            }
            else {

                // POLYGON of the FLAT ROOF

                output1 << "facade polygon wall5-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n";

                for (int i=0;i<4;i++) {
                    output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n";
                }

                output1 << endl;

            }

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        vector<Point> detectorPoint;
        vector<Point> normalVector;

        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // si les batiments sont collés, seulement points visibles
        if ( ecartw == 0.0 ) { // prochain dans la direction hi (wmax)

            double hcurrent;
            int orientnw;
            double hflatnw, hprimenw;

            if ( hi < (wmax-1) ) {
                // calcul de orient, hflat et hprime
                hflatnw = hbase+alleles[(hi+1)*lmax*2+hj*2];

                hprimenw = 4*alleles[(hi+1)*lmax*2+hj*2+1] - TINY;
                orientnw = (int) trunc(hprimenw);

                hprimenw -= 1.0*orientnw;
                hprimenw *= hrmax;
            }
            else {

                hflatnw = hbase;
                hprimenw = 0.0;
                orientnw = 0;

            }

            if ( orient == 0 || orient == 2) {

                if (orient == 0) hcurrent = hflat;
                else hcurrent = hflat + hprime;


                if (orientnw == 0 || orientnw == 2) {

                    double hrel;
                    if (orientnw == 0) hrel = hflatnw + hprimenw;
                    else hrel = hflatnw;

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hcurrent ) );
                    facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel ) );
                    facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel ) );
                    facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );

                    gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                else if (orientnw == 1 || orientnw == 3) {

                    double hrel1, hrel2;
                    if (orientnw == 1) {
                        hrel1 = hflatnw + hprimenw;
                        hrel2 = hflatnw;
                    }
                    else {
                        hrel1 = hflatnw;
                        hrel2 = hflatnw + hprimenw;
                    }

                    // compute the weight in between 1 and 2

                    double lweight = (hcurrent-hflatnw)/(hprimenw+TINY);

                    if (lweight < 1 && lweight > 0) { // entre les deux points

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( weightedCenter( Point( Rectangle[1][0], Rectangle[1][1], hcurrent), Point( Rectangle[2][0], Rectangle[2][1], hcurrent), lweight ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );
                        facade.push_back( weightedCenter( Point( Rectangle[1][0], Rectangle[1][1], hcurrent), Point( Rectangle[2][0], Rectangle[2][1], hcurrent), lweight ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else if (lweight >= 1) {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                }
            }
            else {

                    double hrel1nw, hrel2nw;
                    if (orientnw == 0) {
                        hrel1nw = hflatnw+hprimenw;
                        hrel2nw = hflatnw+hprimenw;
                    }
                    else if (orientnw == 1) {
                        hrel1nw = hflatnw + hprimenw;
                        hrel2nw = hflatnw;
                    }
                    else if (orientnw == 2) {
                        hrel1nw = hflatnw;
                        hrel2nw = hflatnw;
                    }
                    else if (orientnw == 3) {
                        hrel1nw = hflatnw;
                        hrel2nw = hflatnw + hprimenw;
                    }

                    double hrel1, hrel2;
                    if (orient == 1) {
                        hrel1 = hflat + hprime;
                        hrel2 = hflat;
                    }
                    else if (orient == 3) {
                        hrel1 = hflat;
                        hrel2 = hflat + hprime;
                    }

                    // compute the weight in between 1 and 2

                    double lweight = (hrel1-hrel1nw)/((hrel1-hrel1nw)+(hrel2nw-hrel2)+TINY);

                    //cout << "lweight: " << lweight << endl;

                    if (lweight < 1 && lweight > 0) { // entre les deux points

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[1][0], Rectangle[1][1], hrel1), Point( Rectangle[2][0], Rectangle[2][1], hrel2), lweight ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1nw ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[1][0], Rectangle[1][1], hrel1), Point( Rectangle[2][0], Rectangle[2][1], hrel2), lweight ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2nw ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else if (lweight >= 1) {

                        cout << "ici!" << endl;

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1nw ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );


                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1nw ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2nw ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1nw ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2nw ) );


                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[1][0], Rectangle[1][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2nw ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
            }
        }

        if ( ecartl == 0.0 ) { // prochain dans la direction hi (wmax)

            orient = (orient+3)%4; // transformation de l'orientation par rotation de 90°

            double hcurrent;
            int orientnl;
            double hflatnl, hprimenl;

            if ( hj < (lmax-1) ) {
                // calcul de orient, hflat et hprime
                hflatnl = hbase+alleles[hi*lmax*2+(hj+1)*2];

                hprimenl = 4*alleles[hi*lmax*2+(hj+1)*2+1] - TINY;
                orientnl = (int) trunc(hprimenl);

                hprimenl -= 1.0*orientnl;
                hprimenl *= hrmax;
            }
            else {
                hflatnl = hbase;
                hprimenl = 0.0;
                orientnl = 0;
            }

            orientnl = (orientnl+3)%4; // transformation de l'orientation par rotation de 90°

            if ( orient == 0 || orient == 2) {

                if (orient == 0) hcurrent = hflat;
                else hcurrent = hflat + hprime;


                if (orientnl == 0 || orientnl == 2) {

                    double hrel;
                    if (orientnl == 0) hrel = hflatnl + hprimenl;
                    else hrel = hflatnl;

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );
                    facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel ) );
                    facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel ) );
                    facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hcurrent ) );

                    gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                else if (orientnl == 1 || orientnl == 3) {

                    double hrel1, hrel2;
                    if (orientnl == 1) {
                        hrel1 = hflatnl + hprimenl;
                        hrel2 = hflatnl;
                    }
                    else {
                        hrel1 = hflatnl;
                        hrel2 = hflatnl + hprimenl;
                    }

                    // compute the weight in between 1 and 2

                    double lweight = (hcurrent-hflatnl)/(hprimenl+TINY);

                    if (lweight < 1 && lweight > 0) { // entre les deux points

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( weightedCenter( Point( Rectangle[2][0], Rectangle[2][1], hcurrent), Point( Rectangle[3][0], Rectangle[3][1], hcurrent), lweight ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hcurrent ) );
                        facade.push_back( weightedCenter( Point( Rectangle[2][0], Rectangle[2][1], hcurrent), Point( Rectangle[3][0], Rectangle[3][1], hcurrent), lweight ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else if (lweight >= 1) {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hcurrent ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hcurrent ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                }
            }
            else {

                    double hrel1nl, hrel2nl;
                    if (orientnl == 0) {
                        hrel1nl = hflatnl+hprimenl;
                        hrel2nl = hflatnl+hprimenl;
                    }
                    else if (orientnl == 1) {
                        hrel1nl = hflatnl + hprimenl;
                        hrel2nl = hflatnl;
                    }
                    else if (orientnl == 2) {
                        hrel1nl = hflatnl;
                        hrel2nl = hflatnl;
                    }
                    else if (orientnl == 3) {
                        hrel1nl = hflatnl;
                        hrel2nl = hflatnl + hprimenl;
                    }

                    double hrel1, hrel2;
                    if (orient == 1) {
                        hrel1 = hflat + hprime;
                        hrel2 = hflat;
                    }
                    else if (orient == 3) {
                        hrel1 = hflat;
                        hrel2 = hflat + hprime;
                    }

                    // compute the weight in between 1 and 2

                    double lweight = (hrel1-hrel1nl)/((hrel1-hrel1nl)+(hrel2nl-hrel2)+TINY);

                    //cout << "lweight: " << lweight << endl;

                    if (lweight < 1 && lweight > 0) { // entre les deux points

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[2][0], Rectangle[2][1], hrel1), Point( Rectangle[3][0], Rectangle[3][1], hrel2), lweight ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1nl ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[2][0], Rectangle[2][1], hrel1), Point( Rectangle[3][0], Rectangle[3][1], hrel2), lweight ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else if (lweight >= 1) {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );


                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2nl ) );


                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[2][0], Rectangle[2][1], hrel1 ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2nl ) );
                        facade.push_back( Point( Rectangle[3][0], Rectangle[3][1], hrel2 ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
            }

            orient = (orient+1)%4; // rétablissement de l'orientation par rotation de -90°

        }

        // points sur les 4 FACADES

        for (unsigned int id=0; id<4; id++) {

            double hrel = 0.;

            // si les batiments sont collés, seulement points visibles
            if ( ecartw == 0.0 ) {
                if ( id == 1 ) hrel = alleles[hi*lmax*2 + hj*2];
                if ( id == 3 && hi > 0 )        hrel = alleles[hi*lmax*2 + hj*2];
            }
            if ( ecartl == 0.0 ) {
                if ( id == 0 && hj > 0 )        hrel = alleles[hi*lmax*2 + hj*2];
                if ( id == 2 ) hrel = alleles[hi*lmax*2 + hj*2];
            }

            facade.clear();
            grid1.clear();

            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hrel ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hrel ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+alleles[hi*lmax*2+hj*2] ) );
            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+alleles[hi*lmax*2+hj*2] ) );

            gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // création des toits

        roof.clear();
        grid1.clear();

        roof.push_back( Point(Rectangle[(orient+1)%4][0], Rectangle[(orient+1)%4][1], hflat) );
        roof.push_back( Point(Rectangle[(orient+2)%4][0], Rectangle[(orient+2)%4][1], hflat) );
        roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );
        roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );

        gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        if ( ( (hi == 0 || ecartw > 0) && orient == 0) || ( (hi == (wmax-1) || ecartw > 0) && orient == 2) ||
              ( (hj == 0 || ecartl > 0) && orient == 1) || ( (hj == (lmax-1) || ecartl > 0) && orient == 3) ) {

            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );

            gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }
        if ( ( (hj == 0 || ecartl > 0) && orient == 0) || ( (hj == (lmax-1) || ecartl > 0) && orient == 2) ||
              ( (hi == 0 || ecartw > 0) && orient == 3) || ( (hi == (wmax-1) || ecartw > 0) && orient == 1) ) {


            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+1)%4][0], Rectangle[(orient+1)%4][1], hflat) );
            roof.push_back( Point(Rectangle[orient][0], Rectangle[orient][1], hflat+hprime) );

            gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }
        if ( ( (hj == 0 || ecartl > 0) && orient == 2) || ( (hj == (lmax-1) || ecartl > 0) && orient == 0) ||
              ( (hi == 0 || ecartw > 0) && orient == 1) || ( (hi == (wmax-1) || ecartw > 0) && orient == 3) ) {


            roof.clear();
            grid1.clear();

            roof.push_back( Point(Rectangle[(orient+2)%4][0], Rectangle[(orient+2)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat) );
            roof.push_back( Point(Rectangle[(orient+3)%4][0], Rectangle[(orient+3)%4][1], hflat+hprime) );

            gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

        // création de la grille et son output

        double gap = 0.05; // distance between surface and point

        // creation du fichier .inp
        // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

        for (unsigned int i=0;i < detectorPoint.size(); i++) {

            output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                    << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

            output3 << detectorArea[i] << endl;

        }

      } // boucle hj
    }   // boucle hi

    // écriture dans les fichiers

    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    // début de la simulation
    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    //cout << "totalSurface in contact with exterior: " << totalSurface << endl;
    //cout << "Irradiation: " << 0.8*totalRadiation << "\tThermal losses: " << Ubat*24*180*totalSurface*15.0 << endl;

    vector<double> fitness;
    fitness.push_back(0.8*totalRadiation - Ubat*24*180*totalSurface*15.0);
    fitness.push_back(volume(alleles));
    fitness.push_back(0.8*totalRadiation);
    fitness.push_back(Ubat*24*180*totalSurface*15.0);
    fitness.push_back(totalSurface);
    fitness.push_back(groundSurface(alleles));

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << "Number of measuring points: " << detectorArea.size() << "\tAverage points per second: " << detectorArea.size()/evalTime << endl;
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

double RandomSlabsSlopedRoof::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        return volume(alleles) - 1.1*volumeMax;

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        return -volume(alleles) + 0.9*volumeMax;

    }
    else return 0.0;

}

double RandomSlabsSlopedRoof::volume(vector<double> alleles) {

    // calcul du volume perdu
    vector<double> volume;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

          // pour un batiment en particulier -> h1, h2, orient
          double h1 = alleles[hi*lmax*2+hj*2];
          double h2 = 0.0;
          unsigned int orient=0;
          h2 = max(4*alleles[hi*lmax*2+hj*2+1] - TINY, 0.0);
          orient = (unsigned int) trunc(h2);
          h2 -= 1.0*orient;
          h2 *= hrmax;
          // -> r1 et h1prime
//          double h1prime = floor( h1/2.8 ) * 2.8;
//          double r1 = h1 - h1prime;
          // -> dpara et dperp
          double dpara,dperp;
          if ( orient == 0 || orient == 2 ) {
              dpara = w;
              dperp = l;
          }
          else {
              dpara = l;
              dperp = w;
          }

          double volumeTot = w*l*(h1+h2/2.0);
//          for (unsigned int i=0; i<ceil( (r1+h2)/2.8 - TINY ); i++) {
//              //cout << "id: " << i << endl;
//              double d2 = d(dpara, r1, h2, i*2.8+1.5);
//              double d1 = d(dpara, r1, h2, i*2.8);
//              double hd2 = h(dpara, r1, h2, d2);
//              double hd1 = h(dpara, r1, h2, d1);
////              cout << "d(" << i*2.8+1.5 << "): " << d2 << "\td(" << i*2.8 << "): " << d1 << endl;
////              cout << "h(d2): " << hd2 << "\th(d1): " << hd1 << endl;
////              cout << "volume: " << i << "\tdiminished by: " << (d2-d1)*((hd1+hd2)/2.0-i*2.8)*dperp << endl;
//              volumeTot -= (d2-d1)*((hd1+hd2)/2.0-i*2.8)*dperp;
//          }
          volume.push_back( volumeTot );

          // sortie des infos
//          cout << "allele1: " << alleles[hi*lmax*2+hj*2] << "\tallele2: " << alleles[hi*lmax*2+hj*2+1] << endl;
//          cout << "h1: " << h1 << "\th2: " << h2 << "\torient: " << orient << endl;
//          cout << "w: " << w << "\tl: " << l << endl;
//          cout << "volume: " << volumeTot << endl;
//          cout << "r1: " << r1 << "\tdpara: " << dpara << "\tdperp: " << dperp << endl;
//          cout << "volumeComplete: " << w*l*(h1+h2/2.0) << "\tvolumeTot: " << volumeTot << endl;

      }
    }

//    cout << "volumeTot: " << accumulate(volume.begin(), volume.end(), 0.0) << endl;

    return accumulate(volume.begin(), volume.end(), 0.0);

}

double RandomSlabsSlopedRoof::d(double dpara, double r1, double h2, double height) {

    if ( height <= r1 ) return 0.0;
    else if ( height < h2+r1 ) return dpara*(height-r1)/h2;
    else return dpara;

}

double RandomSlabsSlopedRoof::h(double dpara, double r1, double h2, double distance) {

    if ( distance <= 0.0 ) return r1;
    else if ( distance <= dpara ) return r1+h2*distance/dpara;
    else return r1+h2;

}

double RandomSlabsSlopedRoof::groundSurface(vector<double> alleles) {

    // calcul de la surface occupée par les batiments
    double groundSurface = 0.0;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

          // pour un batiment en particulier -> h1, h2, orient
          double h1 = alleles[hi*lmax*2+hj*2];
          double h2 = max(4*alleles[hi*lmax*2+hj*2+1] - TINY, 0.0);
          unsigned int orient = (unsigned int) trunc(h2);;
          h2 -= 1.0*orient;
          h2 *= hrmax;

          if ( h1 > TINY || h2 > TINY ) groundSurface += w*l;

      }
    }

    return groundSurface;

}

// *************************************************************
// CISBAT07 Generic Urban Forms Study: Terrace-courts
// *************************************************************

RandomTerraceCourtsEW::RandomTerraceCourtsEW():Radiance() {

    // création d'un vecteur intermédiare pour calculer le volume
    // la valeur 1 est donnée aux alleles présents

    vector<double> allelesPresence;
    allelesPresence.assign(32*2, 1.0);

    // insertion des points zéros dans le vecteur - pour les cours intérieures

    allelesPresence.insert( allelesPresence.begin() + 11*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 14*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 17*2, 2*2, 0.0 );

    allelesPresence.insert( allelesPresence.begin() + 21*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 24*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 27*2, 2*2, 0.0 );

    allelesPresence.insert( allelesPresence.begin() + 31*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 34*2, 2*2, 0.0 );
    allelesPresence.insert( allelesPresence.begin() + 37*2, 2*2, 0.0 );

    // dimensions de la grille

    vw.push_back(14.0);
    vw.push_back(12.4);
    vw.push_back(12.4);
    vw.push_back(12.4);
    vw.push_back(14.0);

    vl.push_back(10.0);
    vl.push_back(11.4);
    vl.push_back(11.4);
    vl.push_back(10.0);
    vl.push_back(11.4);
    vl.push_back(11.4);
    vl.push_back(10.0);
    vl.push_back(11.4);
    vl.push_back(11.4);
    vl.push_back(10.0);

    wmax=vw.size();
    lmax=vl.size();

    hmax  = 14.0;
    hrmax = 4.0;

    ecartw = (0.0);
    ecartl = (0.0);

    skyFile = "BaselSkyHH.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;

    xinit = 11642.529;
    yinit = 68570.9014;
    hbase = 254.71;

    maxDetectorArea = 1.0;

    for (int i=0; i < 32 ; i++) { // 32 buildings inside this checkerboard

            maxVector.push_back( hmax );  //h1
            minVector.push_back(  0.0 );

            maxVector.push_back( hrmax ); //h2
            minVector.push_back(  0.0 );

    }

    // constraints in volume... 60% of max volume
    nConstraints = 2;
    volumeMax = 0.0;
    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {
          volumeMax += 0.6*vw[hi]*vl[hj]*allelesPresence[hi*lmax*2+hj*2]*(hmax+hrmax/2.0);
      }
    }

    // delete intermediate files ?
    deleteFiles = false;

    cout << "\nProblem: RandomTerraceCourtsEW" << endl;

}

vector<double> RandomTerraceCourtsEW::evaluate(unsigned int id, vector<double> alleles) {

    // création du vecteur avec les hauteurs NS et EO

    double phNS[] = { 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0 };

    // création du vecteur des cours intérieures

    double court[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                       0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                       0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    // insertion des points zéros dans le vecteur - pour les cours intérieures

    alleles.insert( alleles.begin() + 11*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 14*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 17*2, 2*2, 0.0 );

    alleles.insert( alleles.begin() + 21*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 24*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 27*2, 2*2, 0.0 );

    alleles.insert( alleles.begin() + 31*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 34*2, 2*2, 0.0 );
    alleles.insert( alleles.begin() + 37*2, 2*2, 0.0 );

//    // chop precision in alleles for compatibility between MOO and our optimisation software
//    for (unsigned int i=0; i<alleles.size(); i++) {
//
//        if (alleles[i] != 0.0) {
//            double power10 = floor(log10(alleles[i]));
//            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
//        }
//        //cout << setprecision(12) << alleles[i] << endl;
//
//    }

    // montre les alleles

//    cout << "\n" << endl;
//
//    for (int j=0; j < 5; j++) {
//
//        for (int i=0; i< 20; i++) {
//
//            cout << setprecision(8) << alleles[j*20+i] << ", ";
//
//        }
//
//        cout << "\n" << endl;
//
//    }

/*    double penalty;

    // out of bounds, return very bad fitness value
    if (notInField(alleles, &penalty)) return -exp(penalty);
*/

    // start of evaluation
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;
    vector<Point> pN;

    vector<double> detectorArea;

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

 // void prepareRadiance() {

    double h1,h2,x7,x8;
    vector<double> hN;

//    ofstream dotInp(fileInp.c_str(), ios::trunc);
//    dotInp << setprecision(12);

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            h1 = alleles[hi*lmax*2+hj*2];
            h2 = h1+alleles[hi*lmax*2+hj*2+1];
            hN.clear();
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            x7 = 0.5;
            x8 = 0.5;

            Rectangle.clear();

            Rectangle.push_back(Point(vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., -vl[hj]/2.));
            Rectangle.push_back(Point(vw[hi]/2., -vl[hj]/2.));

            rotation2D(Rectangle, 10.0 ); // 10.0 donné par Marylène

            // position de la référence...
            double poshi=0.0, poshj=0.0;
            for (unsigned int i=0; i<hi; i++) poshi += vw[i];
            for (unsigned int j=0; j<hj; j++) poshj += vl[j];
            poshi /= vw[hi];
            poshj /= vl[hj];

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+poshi) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+poshj) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+poshi)	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+poshj) ) ); // posin + posbat

            // calcul des points centraux...

            pN.clear();

            pN.push_back( weightedCenter(Rectangle[0], Rectangle[1], x7) );
            pN.push_back( weightedCenter(Rectangle[1], Rectangle[2], x8) );
            pN.push_back( weightedCenter(Rectangle[3], Rectangle[2], x7) );
            pN.push_back( weightedCenter(Rectangle[0], Rectangle[3], x8) );

            Point centre = weightedCenter(pN[0], pN[2], x8);

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            if (h1 > 0.0) {
                for (int i=0;i<4;i++) {

                        output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n12\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                        output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n12\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                                << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;
                }
            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits
            //if (court[hi*lmax + hj] == 0.0) { // seulement si on est sur un bâtiment
                for (unsigned int i=0;i<4;i++) {

                    output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n9\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                            << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                            << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << endl;

                    output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n9\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                            << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << endl;

                }
            //}

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décompostion des toits
            if (court[hi*lmax + hj] == 0.0) { // seulement si on est sur un bâtiment
                for (unsigned int i=0;i<4;i++) {

                    if ( phNS[hi*lmax + hj] < 0.5 && (i == 0 || i == 2) ) {  // orientation NS

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }
                    if ( phNS[hi*lmax + hj] > 0.5 && (i == 1 || i == 3) ) {  // orientation EO

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }



                }
            }

            // décomposition des facades
            vector<double> h1n, h2n, ecartn, hNn;

            if (hi < (wmax-1)) {
                h1n.push_back(alleles[(hi+1)*lmax*2+hj*2]);
                h2n.push_back(h1n[0]+alleles[(hi+1)*lmax*2+hj*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            if (hj < (lmax-1)) {
                h1n.push_back(alleles[hi*lmax*2+(hj+1)*2]);
                h2n.push_back(h1n[1]+alleles[hi*lmax*2+(hj+1)*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            hNn.push_back( (phNS[(hi+1)*lmax + hj])*h1n[0] + (1.0-phNS[(hi+1)*lmax + hj])*h2n[0] );
            hNn.push_back( (1.0-phNS[hi*lmax + hj+1])*h1n[1] + phNS[hi*lmax + hj+1]*h2n[1] );

            ecartn.push_back(ecartw);
            ecartn.push_back(ecartl);

            for (int i=1;i<3;i++) {

                if (ecartn[i-1] < TINY) { // pas décart entre les batiments selon w

                    // calcul des cas d'intersection

                    double weight = (h1n[i-1] - h1) / ((h1n[i-1] - h1) + (h2 - hNn[i-1]));

                    if ( weight > 0 && weight < 1) { // une intersection

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }

                }
            }

            if ( hi == 0 ) {

                int i=3;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            if ( hj == 0 ) {

                int i=0;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            double gap = 0.05;  // 5cm from facade

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                       << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

                output3 << detectorArea[i] << endl;

            }

      } // boucle hj
    }   // boucle hi


    // flushing files & stringstreams
    //dotInp.flush();
    //dotInp.close();

    output1.flush();
    output2.flush();
    output3.flush();

    // writing results
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    // début de la simulation
    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    //cout << "totalSurface in contact with exterior: " << totalSurface << endl;
    //cout << "p1: " << 0.8*totalRadiation << "\tp2: " << Ubat*24*180*totalSurface*15.0 << endl;

    vector<double> fitness;
    fitness.push_back(0.8*totalRadiation - Ubat*24*180*totalSurface*15.0);
    fitness.push_back(volume(alleles));
    fitness.push_back(0.8*totalRadiation);
    fitness.push_back(Ubat*24*180*totalSurface*15.0);
    fitness.push_back(totalSurface);
    fitness.push_back(groundSurface(alleles));

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size()/3 ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << "Number of measuring points: " << detectorArea.size() << "\tAverage points per second: " << detectorArea.size()/evalTime << endl;
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

double RandomTerraceCourtsEW::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        return volume(alleles) - 1.1*volumeMax;

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        return -volume(alleles) + 0.9*volumeMax;

    }
    else return 0.0;

}

double RandomTerraceCourtsEW::volume(vector<double> alleles) {

    // création du vecteur avec les hauteurs NS et EO

//    double phNS[] = { 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
//                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
//                      0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0 };

    // insertion des points zéros dans le vecteur - pour les cours intérieures

    if (alleles.size() == 64) {
        alleles.insert( alleles.begin() + 11*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 14*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 17*2, 2*2, 0.0 );

        alleles.insert( alleles.begin() + 21*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 24*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 27*2, 2*2, 0.0 );

        alleles.insert( alleles.begin() + 31*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 34*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 37*2, 2*2, 0.0 );
    }

    // calcul du volume perdu
    vector<double> volume;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

          // pour un batiment en particulier -> h1, h2, orient
          double h1 = alleles[hi*lmax*2+hj*2];
          double h2 = alleles[hi*lmax*2+hj*2+1];
//          double dpara = (phNS[hi*lmax + hj])*vw[hi] + (1.0-phNS[hi*lmax + hj])*vl[hj];
//          double dperp = (1.0-phNS[hi*lmax + hj])*vw[hi] + phNS[hi*lmax + hj]*vl[hj];

          // -> r1 et h1prime
//          double h1prime = floor( h1/2.8 ) * 2.8;
//          double r1 = h1 - h1prime;

          double volumeTot = vw[hi]*vl[hj]*(h1+h2/2.0);
//          for (unsigned int i=0; i<ceil( (r1+h2)/2.8 - TINY ); i++) {
//              //cout << "id: " << i << endl;
//              double d2 = d(dpara, r1, h2, i*2.8+1.5);
//              double d1 = d(dpara, r1, h2, i*2.8);
//              double hd2 = h(dpara, r1, h2, d2);
//              double hd1 = h(dpara, r1, h2, d1);
////              cout << "d(" << i*2.8+1.5 << "): " << d2 << "\td(" << i*2.8 << "): " << d1 << endl;
////              cout << "h(d2): " << hd2 << "\th(d1): " << hd1 << endl;
////              cout << "volume: " << i << "\tdiminished by: " << 2.0*(d2-d1)*((hd1+hd2)/2.0-i*2.8)*dperp << endl;
//              volumeTot -= 2.0*(d2-d1)*((hd1+hd2)/2.0-i*2.8)*dperp;
//          }
          volume.push_back( volumeTot );

          // sortie des infos
//          cout << "allele1: " << alleles[hi*lmax*2+hj*2] << "\tallele2: " << alleles[hi*lmax*2+hj*2+1] << endl;
//          cout << "h1: " << h1 << "\th2: " << h2 << endl;
//          cout << "r1: " << r1 << "\tdpara: " << dpara << "\tdperp: " << dperp << endl;
//          cout << "w: " << vw[hi] << "\tl: " << vl[hj] << endl;
//          cout << "volume: " << volumeTot << endl;

      }
    }

//    cout << "Volume final: " << accumulate(volume.begin(), volume.end(), 0.0) << endl;

    return accumulate(volume.begin(), volume.end(), 0.0);

}

double RandomTerraceCourtsEW::groundSurface(vector<double> alleles) {

    // insertion des points zéros dans le vecteur - pour les cours intérieures

    if (alleles.size() == 64) {
        alleles.insert( alleles.begin() + 11*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 14*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 17*2, 2*2, 0.0 );

        alleles.insert( alleles.begin() + 21*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 24*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 27*2, 2*2, 0.0 );

        alleles.insert( alleles.begin() + 31*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 34*2, 2*2, 0.0 );
        alleles.insert( alleles.begin() + 37*2, 2*2, 0.0 );
    }

    // calcul de la surface au sol
    double groundSurface = 0.0;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

          // pour un batiment en particulier -> h1, h2, orient
          double h1 = alleles[hi*lmax*2+hj*2];
          double h2 = alleles[hi*lmax*2+hj*2+1];

          if ( h1 > TINY || h2 > TINY ) { groundSurface += vw[hi]*vl[hj]; }

      }
    }

    return groundSurface;

}

double RandomTerraceCourtsEW::d(double dpara, double r1, double h2, double height) {

    if ( height <= r1 ) return 0.0;
    else if ( height < h2+r1 ) return (dpara/2.0)*(height-r1)/h2;
    else return (dpara/2.0);

}

double RandomTerraceCourtsEW::h(double dpara, double r1, double h2, double distance) {

    if ( distance <= 0.0 ) return r1;
    else if ( distance <= dpara ) return r1+h2*distance/(dpara/2.0);
    else return r1+h2;

}

// *************************************************************
// CISBAT07 Generic Urban Forms Study: Pavilions
// *************************************************************

Pavilions::Pavilions():Radiance() {

    w = 10.0;
    l = 10.0;
    wmax=3;
    lmax=6;
    hmax=14.0;
    hrmax = 4.0;
    ecartw = (17.61/w);
    ecartl = (9.67/l);
    skyFile = "BaselSkyHH.rad";
    groundFile = "3df_matthaeus4_sit.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = true;

    xinit = 11642.529;
    yinit = 68570.9014;
    hbase = 254.71;

    maxDetectorArea = 1.0;

    for (unsigned int i=0; i < (wmax*lmax) ; i++) {

            maxVector.push_back( hmax );
            minVector.push_back(  0.0 );

            maxVector.push_back( hrmax );
            minVector.push_back(  0.0 );

    }

    // constraints in volume... 60% of max volume
    nConstraints = 2;
    volumeMax = 0.0;
    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {
          volumeMax += 0.6*w*l*(hmax+hrmax/3.0);
      }
    }

    // delete intermediate files ?
    deleteFiles = false;

    cout << "\nProblem: Pavilions" << endl;

}

vector<double> Pavilions::evaluate(unsigned int id, vector<double> alleles) {

//    // chop precision in alleles for compatibility between MOO and our optimisation software
//    for (unsigned int i=0; i<alleles.size(); i++) {
//
//        if (alleles[i] != 0.0) {
//            double power10 = floor(log10(alleles[i]));
//            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
//        }
//        //cout << setprecision(12) << alleles[i] << endl;
//
//    }

    // timer start
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    // vecteur contenant la base
    vector<Point> Rectangle;

    // vecteur des aires de détection autour des points
    vector<double> detectorArea;

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

 // void prepareRadiance() {

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            Rectangle.clear();

            Rectangle.push_back(Point(w/2., l/2.));
            Rectangle.push_back(Point(-w/2., l/2.));
            Rectangle.push_back(Point(-w/2., -l/2.));
            Rectangle.push_back(Point(w/2., -l/2.));

            rotation2D(Rectangle, 10.0 ); // 10° given by Marylène

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+hi*(1+ecartw)) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+hj*(1+ecartl)) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+hi*(1+ecartw))	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+hj*(1+ecartl))) ); // posin + posbat

            // calcul du centre du bâtiment

            Point centre = centroid2D(Rectangle);

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            if (alleles[hi*lmax*2+hj*2] > 0.0) {
                for (int i=0;i<4;i++) {

                        output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n12\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                                << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                                << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n" << endl;

                }
            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {
                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "\n0\n0\n9\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << alleles[hi*lmax*2+hj*2]+hbase << "\n"
                    << centre[0] << "\t" << centre[1] << "\t" << alleles[hi*lmax*2+hj*2+1]+alleles[hi*lmax*2+hj*2]+hbase << "\n" << endl;
            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décomposition des facades

            if (alleles[hi*lmax*2+hj*2] > 0.0) {
                for (int i=0;i<4;i++) {

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                    facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                    facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*2+hj*2]+hbase ) );
                    facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*2+hj*2]+hbase) );

                    gridRectangle3(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
            }

            // décompostion des toits

            for (int i=0;i<4;i++) {

                roof.clear();
                grid1.clear();

                roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], alleles[hi*lmax*2+hj*2]+hbase) );
                roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], alleles[hi*lmax*2+hj*2]+hbase) );
                roof.push_back( Point(centre[0], centre[1], alleles[hi*lmax*2+hj*2+1]+alleles[hi*lmax*2+hj*2]+hbase) );

                gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            double gap = 0.05; //distance between surface and point


            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    // début de la simulation
    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // évaluation de la perte d'énergie par saison

    double Ubat = 0.38;

    //cout << "totalSurface in contact with exterior: " << totalSurface << endl;
    //cout << "p1: " << 0.8*totalRadiation << "\tp2: " << Ubat*24*180*totalSurface*15.0 << endl;

    vector<double> fitness;
    fitness.push_back(0.8*totalRadiation - Ubat*24*180*totalSurface*15.0);

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size()/3 ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << "Number of measuring points: " << detectorArea.size() << "\tAverage points per second: " << detectorArea.size()/evalTime << endl;
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// Josep - Slope study: PepX
// *************************************************************

PepX::PepX(string skyFile, string siteFile):Radiance() {

    vw.push_back(12.0);
    vw.push_back(4.0);
    vw.push_back(12.0);
    vw.push_back(10.0);
    vw.push_back(6.0);
    vw.push_back(10.0);
    vw.push_back(12.0);
    vw.push_back(4.0);
    vw.push_back(12.0);
    vw.push_back(4.0);
    vw.push_back(12.0);

    vl.push_back(12.0);
    vl.push_back(4.0);
    vl.push_back(12.0);
    vl.push_back(12.0);
    vl.push_back(6.0);
    vl.push_back(6.0);
    vl.push_back(12.0);
    vl.push_back(4.0);
    vl.push_back(12.0);
    vl.push_back(4.0);
    vl.push_back(12.0);

    wmax=vw.size();
    lmax=vl.size();

    hmax  = 9.0;
    hrmax = 2.0;

    ecartw = (0.0);
    ecartl = (0.0);

    this->skyFile = skyFile;
    groundFile = siteFile;
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = false;

    xinit = 0.0;
    yinit = 0.0;
    hbase = 0.0;

    maxDetectorArea = 1.0;

    for (int i=0; i < 61 ; i++) { // 32 buildings inside this checkerboard

            maxVector.push_back( hmax );  //h1
            minVector.push_back(  0.0 );

            maxVector.push_back( hrmax ); //h2
            minVector.push_back(  0.0 );

    }

    cout << "\nProblem: PepX" << endl;

}

vector<double> PepX::evaluate(unsigned int id, vector<double> alleles) {

    // insertion des points zéros dans le vecteur - pour les cours intérieures

    alleles.insert( alleles.begin() + 7*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 9*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 11*2, 8*2, 0.0 );
    alleles.insert( alleles.begin() + 20*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 23*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 29*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 31*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 34*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 36*2, 3*2, 0.0 );
    alleles.insert( alleles.begin() + 40*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 42*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 45*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 47*2, 3*2, 0.0 );
    alleles.insert( alleles.begin() + 51*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 53*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 56*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 58*2, 3*2, 0.0 );
    alleles.insert( alleles.begin() + 62*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 64*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 67*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 73*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 75*2, 1*2, 0.0 );

    alleles.insert( alleles.begin() + 78*2, 9*2, 0.0 );

    alleles.insert( alleles.begin() + 89*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 95*2, 3*2, 0.0 );

    alleles.insert( alleles.begin() + 100*2, 9*2, 0.0 );

    alleles.insert( alleles.begin() + 111*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 117*2, 1*2, 0.0 );
    alleles.insert( alleles.begin() + 119*2, 1*2, 0.0 );

    // chop precision in alleles for compatibility between MOO and our optimisation software
    for (unsigned int i=0; i<alleles.size(); i++) {

        if (alleles[i] != 0.0) {
            double power10 = floor(log10(alleles[i]));
            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
        }
        //cout << setprecision(12) << alleles[i] << endl;

    }

    // création du vecteur avec les hauteurs NS et EO

    double phNS[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0 };


/*    double penalty;

    // out of bounds, return very bad fitness value
    if (notInField(alleles, &penalty)) return -exp(penalty);
*/

    // start of evaluation
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;
    vector<Point> pN;

    vector<double> detectorArea;

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output << setprecision(12);
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

 // void prepareRadiance() {

    double h1,h2,x7,x8;
    vector<double> hN;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            h1 = alleles[hi*lmax*2+hj*2];
            h2 = h1+alleles[hi*lmax*2+hj*2+1];
            hN.clear();
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            x7 = 0.5;
            x8 = 0.5;

            Rectangle.clear();

            Rectangle.push_back(Point(vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., -vl[hj]/2.));
            Rectangle.push_back(Point(vw[hi]/2., -vl[hj]/2.));

            rotation2D(Rectangle, 0.0 );

            // position de la référence...
            double poshi=0.0, poshj=0.0;
            for (unsigned int i=0; i<hi; i++) poshi += vw[i];
            for (unsigned int j=0; j<hj; j++) poshj += vl[j];
            poshi /= vw[hi];
            poshj /= vl[hj];

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+poshi) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+poshj) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+poshi)	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+poshj) ) ); // posin + posbat

            // calcul des points centraux...

            pN.clear();

            pN.push_back( weightedCenter(Rectangle[0], Rectangle[1], x7) );
            pN.push_back( weightedCenter(Rectangle[1], Rectangle[2], x8) );
            pN.push_back( weightedCenter(Rectangle[3], Rectangle[2], x7) );
            pN.push_back( weightedCenter(Rectangle[0], Rectangle[3], x8) );

            Point centre = weightedCenter(pN[0], pN[2], x8);

//          void Problem::write3DShape(int hi, int hj) {

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

            }

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << "\n"
                       << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;


                output << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << endl;

            }

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;
            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n9\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << endl;

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n9\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << endl;

            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décompostion des toits

            if ( h1 > 0.0 ) {

                for (int i=0;i<4;i++) {

                    if ( phNS[hi*lmax + hj] < 0.5 && (i == 0 || i == 2) ) {  // orientation NS

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }
                    if ( phNS[hi*lmax + hj] > 0.5 && (i == 1 || i == 3) ) {  // orientation EO

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }



                }

            }



            // décomposition des facades

            vector<double> h1n, h2n, ecartn, hNn;

            if (hi < (wmax-1)) {
                h1n.push_back(alleles[(hi+1)*lmax*2+hj*2]);
                h2n.push_back(h1n[0]+alleles[(hi+1)*lmax*2+hj*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            if (hj < (lmax-1)) {
                h1n.push_back(alleles[hi*lmax*2+(hj+1)*2]);
                h2n.push_back(h1n[1]+alleles[hi*lmax*2+(hj+1)*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            hNn.push_back( (phNS[(hi+1)*lmax + hj])*h1n[0] + (1.0-phNS[(hi+1)*lmax + hj])*h2n[0] );
            hNn.push_back( (1.0-phNS[hi*lmax + hj+1])*h1n[1] + phNS[hi*lmax + hj+1]*h2n[1] );

            ecartn.push_back(ecartw);
            ecartn.push_back(ecartl);

            for (int i=1;i<3;i++) {

                if (ecartn[i-1] < TINY) { // pas décart entre les batiments selon w

                    // calcul des cas d'intersection

                    double weight = (h1n[i-1] - h1) / ((h1n[i-1] - h1) + (h2 - hNn[i-1]));

                    if ( weight > 0 && weight < 1) { // une intersection

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }

                }
            }

            if ( hi == 0 ) {

                int i=3;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            if ( hj == 0 ) {

                int i=0;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            double gap = 0.05;  // 5cm from facade

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << "\t#building: " << hi+1 << "-" << hj+1 << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size();i++) {

        totalRadiation+=radianceResults[i]*detectorArea[i];
        totalSurface+=detectorArea[i];

    }

    // effacement des fichier .oct

//    string fileOct = file + ".oct";
//    remove(fileOct.c_str());
//    remove(file3D.c_str() );
//    remove(fileRad.c_str());
//    remove(fileInp.c_str());
//    remove(fileDet.c_str());
//    remove(fileOut.c_str());

    cout << "p1: " << 0.8*totalRadiation << endl;
    vector<double> fitness;
    fitness.assign(1,0.8*totalRadiation);

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size() ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// Josep - Slope study: PepY
// *************************************************************

PepY::PepY(string skyFile, string siteFile):Radiance() {

    vw.assign(8, 12.0);
    vl.assign(10, 6.0);

    wmax=vw.size();
    lmax=vl.size();

    hmax  = 9.0;
    hrmax = 2.0;

    ecartw = (0.0);
    ecartl = (0.0);

    this->skyFile = skyFile;
    groundFile = siteFile;
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = false;

    xinit = 0.0;
    yinit = 0.0;
    hbase = 0.0;

    maxDetectorArea = 1.0;

    for (unsigned int i=0; i < (wmax*lmax) ; i++) { // 32 buildings inside this checkerboard

            maxVector.push_back( hmax );  //h1
            minVector.push_back(  0.0 );

            maxVector.push_back( hrmax ); //h2
            minVector.push_back(  0.0 );

    }

    cout << "\nProblem: PepY" << endl;

}

vector<double> PepY::evaluate(unsigned int id, vector<double> alleles) {


    // chop precision in alleles for compatibility between MOO and our optimisation software
    for (unsigned int i=0; i<alleles.size(); i++) {

        if (alleles[i] != 0.0) {
            double power10 = floor(log10(alleles[i]));
            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e6 ) / 1.0e6) * pow(10.0, power10);
        }
        //cout << setprecision(12) << alleles[i] << endl;

    }

    // création du vecteur avec les hauteurs NS et EO

    vector<double> phNS;
    phNS.assign(80, 1.0);

    // start of evaluation
    time_t start, end;

    time(&start); // start timer

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    vector<Point> Rectangle;
    vector<Point> pN;

    vector<double> detectorArea;

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output << setprecision(12);
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

 // void prepareRadiance() {

    double h1,h2,x7,x8;
    vector<double> hN;

    for (unsigned int hi=0;hi<wmax;hi++) {
      for (unsigned int hj=0;hj<lmax;hj++) {

            h1 = alleles[hi*lmax*2+hj*2];
            h2 = h1+alleles[hi*lmax*2+hj*2+1];
            hN.clear();
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            hN.push_back( (1.0-phNS[hi*lmax + hj])*h1 + phNS[hi*lmax + hj]*h2 );
            hN.push_back( (phNS[hi*lmax + hj])*h1 + (1.0-phNS[hi*lmax + hj])*h2 );
            x7 = 0.5;
            x8 = 0.5;

            Rectangle.clear();

            Rectangle.push_back(Point(vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., vl[hj]/2.));
            Rectangle.push_back(Point(-vw[hi]/2., -vl[hj]/2.));
            Rectangle.push_back(Point(vw[hi]/2., -vl[hj]/2.));

            rotation2D(Rectangle, 10.0 ); // 10.0 donné par Marylène

            // position de la référence...
            double poshi=0.0, poshj=0.0;
            for (unsigned int i=0; i<hi; i++) poshi += vw[i];
            for (unsigned int j=0; j<hj; j++) poshj += vl[j];
            poshi /= vw[hi];
            poshj /= vl[hj];

            translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5+poshi) + (Rectangle[3][0]-Rectangle[0][0])*(0.5+poshj) ) , ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5+poshi)	 + (Rectangle[3][1]-Rectangle[0][1])*(0.5+poshj) ) ); // posin + posbat

            // calcul des points centraux...

            pN.clear();

            pN.push_back( weightedCenter(Rectangle[0], Rectangle[1], x7) );
            pN.push_back( weightedCenter(Rectangle[1], Rectangle[2], x8) );
            pN.push_back( weightedCenter(Rectangle[3], Rectangle[2], x7) );
            pN.push_back( weightedCenter(Rectangle[0], Rectangle[3], x8) );

            Point centre = weightedCenter(pN[0], pN[2], x8);

//          void Problem::write3DShape(int hi, int hj) {

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

            }

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << "\n"
                       << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;


                output << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                       << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                       << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                       << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << endl;

            }

//          void Problem::writeRadianceInput(int hi, int hj) {

            // creation du fichier .RAD
            // materiaux

            if ((hi == 0) && (hj == 0)) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

            // batiment (hi, hj)

            // POLYGON of the FACADES

            for (int i=0;i<4;i++) {

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;

                    output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n12\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                            << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i]+hbase << "\n"
                            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << endl;
            }

            // POLYGON(S) of the ROOF

            // dessin des triangles des toits

            for (int i=0;i<4;i++) {

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-1" << "\n0\n0\n9\n"
                        << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << h1+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2 + hbase << endl;

                output1 << "facade polygon roof" << i+1 << "-" << hi+1 << "-" << hj+1 << "-2" << "\n0\n0\n9\n"
                        << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << h1+hbase << "\n"
                        << centre[0] << "\t" << centre[1] << "\t" << h2+hbase << "\n"
                        << pN[i][0] << "\t" << pN[i][1] << "\t" << hN[i] + hbase << endl;

            }

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            // décompostion des toits

            if ( h1 > 0.0 ) {

                for (int i=0;i<4;i++) {

                    if ( phNS[hi*lmax + hj] < 0.5 && (i == 0 || i == 2) ) {  // orientation NS

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }
                    if ( phNS[hi*lmax + hj] > 0.5 && (i == 1 || i == 3) ) {  // orientation EO

                        roof.clear();
                        grid1.clear();

                        roof.push_back( Point(Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        roof.push_back( Point(Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                        roof.push_back( Point(pN[(i+1)%4][0], pN[(i+1)%4][1], hN[(i+1)%4]+hbase) );
                        roof.push_back( Point(pN[(i+3)%4][0], pN[(i+3)%4][1], hN[(i+3)%4]+hbase) );

                        gridRectangle3(roof, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );
                    }



                }

            }



            // décomposition des facades

            vector<double> h1n, h2n, ecartn, hNn;

            if (hi < (wmax-1)) {
                h1n.push_back(alleles[(hi+1)*lmax*2+hj*2]);
                h2n.push_back(h1n[0]+alleles[(hi+1)*lmax*2+hj*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            if (hj < (lmax-1)) {
                h1n.push_back(alleles[hi*lmax*2+(hj+1)*2]);
                h2n.push_back(h1n[1]+alleles[hi*lmax*2+(hj+1)*2+1]);
            }
            else {
                h1n.push_back(0.0);
                h2n.push_back(0.0);
            }

            hNn.push_back( (phNS[(hi+1)*lmax + hj])*h1n[0] + (1.0-phNS[(hi+1)*lmax + hj])*h2n[0] );
            hNn.push_back( (1.0-phNS[hi*lmax + hj+1])*h1n[1] + phNS[hi*lmax + hj+1]*h2n[1] );

            ecartn.push_back(ecartw);
            ecartn.push_back(ecartl);

            for (int i=1;i<3;i++) {

                if (ecartn[i-1] < TINY) { // pas décart entre les batiments selon w

                    // calcul des cas d'intersection

                    double weight = (h1n[i-1] - h1) / ((h1n[i-1] - h1) + (h2 - hNn[i-1]));

                    if ( weight > 0 && weight < 1) { // une intersection

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( weightedCenter( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase), Point(pN[i][0], pN[i][1], hN[i]+hbase) , weight ) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }
                    else {

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hNn[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                        facade.clear();
                        grid1.clear();

                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1+hbase) );
                        facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                        facade.push_back( Point( Rectangle[i+1][0], Rectangle[i+1][1], h1n[i-1]+hbase) );

                        gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    }

                }
            }

            if ( hi == 0 ) {

                int i=3;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            if ( hj == 0 ) {

                int i=0;

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], h1+hbase) );
                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[i][0], Rectangle[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                facade.clear();
                grid1.clear();

                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], h1+hbase) );
                facade.push_back( Point(pN[i][0], pN[i][1], hN[i]+hbase) );
                facade.push_back( Point( Rectangle[(i+1)%4][0], Rectangle[(i+1)%4][1], hbase) );

                gridTriangle2(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

            }

            double gap = 0.05;  // 5cm from facade

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << "\t#building: " << hi+1 << "-" << hj+1 << endl;

            }

      } // boucle hj
    }   // boucle hi

    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size();i++) {

        totalRadiation+=radianceResults[i]*detectorArea[i];
        totalSurface+=detectorArea[i];

    }

    // effacement des fichier .oct

//    string fileOct = file + ".oct";
//    remove(fileOct.c_str());
//    remove(file3D.c_str() );
//    remove(fileRad.c_str());
//    remove(fileInp.c_str());
//    remove(fileDet.c_str());
//    remove(fileOut.c_str());

    // évaluation de la perte d'énergie par saison

    cout << "p1: " << 0.8*totalRadiation << endl;
    vector<double> fitness;
    fitness.assign(1,0.8*totalRadiation);

    // vérification de la longueur du fichier (entrée = sortie)

    if ( detectorArea.size() != radianceResults.size() ) {

        string msg = "Error in simulating: ";
        msg += id;
        msg += " - .Det file not the same size as .Out file";

        throw (msg);

    }

    // fin de l'evaluation

    time(&end);

    cout << file << "->\ttime to simulate: " << difftime(end, start) << " s" << endl;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// Optimisation - AppliedSoftComputing
// *************************************************************

RandomPavilionsFlatRoof::RandomPavilionsFlatRoof():Radiance() {

    w = 10.0;
    l = 10.0;

    int nBuildings = 11;

    wmax = 40; //m
    lmax = 90; //m
    hmax = 20.0; //m

    skyFile = "BaselSkyHH.rad";
    groundFile = "FlatSite.rad";
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = false;

    xinit = 0.0;
    yinit = 0.0;
    hbase = 0.0;

    maxDetectorArea = 1.0;

    for (int i=0; i < nBuildings ; i++) {

            // limite en position x
            maxVector.push_back( wmax );
            minVector.push_back(  0.0 );

            // limite en position y
            maxVector.push_back( lmax );
            minVector.push_back(  0.0 );

    }

    cout << "\nProblem: RandomPavilionsFlatRoof" << endl;

}

vector<double> RandomPavilionsFlatRoof::evaluate(unsigned int id, vector<double> alleles) {

    // chop precision in alleles for compatibility between MOO and our optimisation software
    for (unsigned int i=0; i<alleles.size(); i++) {

        if (alleles[i] != 0.0) {
            double power10 = floor(log10(alleles[i]));
            alleles[i] = (floor( (alleles[i]/pow(10.0, power10)) * 1.0e3 ) / 1.0e3) * pow(10.0, power10);
        }
        //cout << setprecision(12) << alleles[i] << endl;

    }

    // timer start
    time_t start, end;

    time(&start); // start timer

    ostringstream ss;
    ss << "building" << id;

    string file = ss.str();

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    ostringstream output, output1, output2, output3;

    // choix de la précision (12 digits en tout - avant et après la virgule)

    output << setprecision(12);
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // open the output files
    //ofstream dotInp(fileInp.c_str(), ios::trunc);
    //dotInp << setprecision(12);

    // création de la surface de base

    vector<Point> Rectangle;
    vector<double> detectorArea; // vecteur contenant les aires des surfaces de détection
    vector<Point> detectorPoint;
    vector<Point> normalVector;

    // number of buildings

    int nBuildings = alleles.size()/2;

 // void prepareRadiance() {

    for (int hi=0; hi<nBuildings; hi++) {

        Rectangle.clear();

        Rectangle.push_back(Point(w/2., l/2.));
        Rectangle.push_back(Point(-w/2., l/2.));
        Rectangle.push_back(Point(-w/2., -l/2.));
        Rectangle.push_back(Point(w/2., -l/2.));

        rotation2D(Rectangle, 0.0); // no rotation

        translation2D(Rectangle, ( xinit + (Rectangle[1][0]-Rectangle[0][0])*(0.5) + (Rectangle[3][0]-Rectangle[0][0])*(0.5) + alleles[2*hi]),
                                 ( yinit + (Rectangle[1][1]-Rectangle[0][1])*(0.5) + (Rectangle[3][1]-Rectangle[0][1])*(0.5) + alleles[2*hi+1])); // posin + posbat

//          void Problem::write3DShape(int hi, int hj) {

        // POLYGON of the FACADES

        for (int i=0;i<4;i++) {

            output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
            << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hmax+hbase << "\n"
            << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hmax+hbase << endl;

        }

        // ROOF

        for (int i=0;i<4;i++) {
            output << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hmax+hbase << "\n";
        }

//          void Problem::writeRadianceInput(int hi, int hj) {

        // creation du fichier .RAD
        // materiaux

        if (hi == 0) output1 << "void plastic facade" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

        // POLYGON of the FACADES

        for (int i=0;i<4;i++) {

            output1 << "facade polygon wall" << i+1 << "-" << hi+1 << "\n0\n0\n12\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hbase << "\n"
                    << Rectangle[(i+1)%4][0] << "\t" << Rectangle[(i+1)%4][1] << "\t" << hmax+hbase << "\n"
                    << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hmax+hbase << "\n" << endl;

        }

        // POLYGON of the FLAT ROOF

        output1 << "facade polygon wall5-" << hi+1 << "\n0\n0\n12\n";

        for (int i=0;i<4;i++) {
            output1 << Rectangle[i][0] << "\t" << Rectangle[i][1] << "\t" << hmax+hbase << "\n";
        }

        output1 << endl;

//          void Problem::writeINP(int hi, int hj) {

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> facade;
        vector<Point> roof;
        vector<Point> grid1;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // création d'une place nette de points pour y mettre le batiment

        vector<double> d1, d2;

        d1.push_back(Rectangle[1][0]-Rectangle[0][0]);
        d1.push_back(Rectangle[1][1]-Rectangle[0][1]);

        d2.push_back(Rectangle[3][0]-Rectangle[0][0]);
        d2.push_back(Rectangle[3][1]-Rectangle[0][1]);

        removeVertically( Point( Rectangle[0][0], Rectangle[0][1], Rectangle[0][2] ), d1, d2, detectorPoint, normalVector, detectorArea );

        // points sur les 4 FACADES

        for (int id=0; id<4; id++) {

            facade.clear();
            grid1.clear();

            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase ) );
            facade.push_back( Point( Rectangle[(id+1)%4][0], Rectangle[(id+1)%4][1], hbase+hmax ) );
            facade.push_back( Point( Rectangle[id][0], Rectangle[id][1], hbase+hmax ) );

            gridRectangle3(facade, maxDetectorArea, grid1, pointArea, n1, n2, n3);

            detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
            normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
            detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

        }

       // FLAT ROOF

        roof.clear();
        grid1.clear();

        for (int i=0; i<4;i++) {
            roof.push_back( Point( Rectangle[i][0], Rectangle[i][1], hmax+hbase) );
        }

        gridRectangle3(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

    }   // boucle hi

    // creation du fichier .inp
    // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

    double gap = 0.05;

    for (unsigned int i=0;i < detectorPoint.size(); i++) {

        output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

        output3 << detectorArea[i] << endl;

    }

    // flushing everything
    //dotInp.flush();
    //dotInp.close();

    output.flush();
    output1.flush();
    output2.flush();
    output3.flush();

    // écriture dans les fichiers
    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResults(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size();i++) {

        totalRadiation+=radianceResults[i]*detectorArea[i];
        totalSurface+=detectorArea[i];

    }

    // effacement des fichier .oct

    string fileOct = file + ".oct";
    remove(fileOct.c_str());
    remove(file3D.c_str() );
    //remove(fileRad.c_str());
    //remove(fileInp.c_str());
    //remove(fileDet.c_str());
    //remove(fileOut.c_str());

    vector<double> fitness;
    fitness.assign(1,0.8*totalRadiation);

    // fin de l'evaluation

    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

// *************************************************************
// Optimisation - Projet Huang
// *************************************************************

ProjetHuang::ProjetHuang():Radiance() {

    skyFile = "BaselSky.rad";
    groundFile = "triangulation.rad";
    environmentFile = "lumiere.rad";
    environment = true;

    maxDetectorArea = 1.0;

    for (int i=0; i < 31 ; i++) {

            // limite en position x
            maxVector.push_back( 6.0 );
            minVector.push_back( 3.0 );
    }

    nConstraints = 32;

    cout << "\nProblem: ProjetHuang" << endl;

}

vector<double> ProjetHuang::evaluate(unsigned int id, vector<double> alleles) {

    // timer start
    time_t start, end;

    time(&start); // start timer

    // fichiers pour enregister
    ostringstream ss;
    ss << "building" << id;
    string file = ss.str();

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // choix de la précision (12 digits en tout - avant et après la virgule)
    ostringstream output, output1, output2, output3;
    output << setprecision(12);
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // open the output files
    //ofstream dotInp(fileInp.c_str(), ios::trunc);
    //dotInp << setprecision(12);

    // création de la surface de base

    vector<string>          labels;
    vector<vector<Point> >  polygons;
    vector<double> detectorArea; // vecteur contenant les aires des surfaces de détection
    vector<Point>  detectorPoint;
    vector<Point>  normalVector;

    // lecture du fichier radiance en input décomposé en triangles

    string tampon;
    char tampon1[200];

    // chargement des données a afficher
    ifstream input1 ("varying.rad", ios::binary | ios::in );

    // test d'ouverture
    if (!input1.is_open()) throw(string("Error opening: varying.rad"));

    while (!input1.eof()) {

        input1 >> tampon; // varying

        if ( tampon[0] == '#' ) input1.getline(tampon1, 200,'\n');
        else if ( tampon == "varying") {

            if (input1.eof()) break;
            input1 >> tampon; // polygon
            input1 >> tampon; // label
            labels.push_back(tampon);
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of points (3x3)
            int number = atoi(tampon.c_str())/3;

            vector<Point> polygon;
            for (int i=0; i<number; i++) {

                input1 >> tampon;
                double x = atof(tampon.c_str());
                input1 >> tampon;
                double y = atof(tampon.c_str());
                input1 >> tampon;
                int z = atoi(tampon.substr(1).c_str());

                polygon.push_back( Point(x,y,alleles[z]) );

            }
            polygons.push_back(polygon);
        }
    }

    input1.close();

 // void prepareRadiance() {

    for (unsigned int i=0; i<polygons.size(); i++) {

//          void Problem::write3DShape(int hi, int hj) {
        for (unsigned int j=0; j<polygons[i].size(); j++) output << polygons[i][j][0] << "\t" << polygons[i][j][1] << "\t" << polygons[i][j][2] << "\n";
        output << endl;

//          void Problem::writeRadianceInput(int hi, int hj) {

        // creation du fichier .RAD
        // materiaux

        //if (i == 0) output1 << "void plastic varying" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

        output1 << "varying polygon " << labels[i] << "\n0\n0\n9\n";

        for (unsigned int j=0;j<polygons[i].size();j++) output1 << polygons[i][j][0] << "\t" << polygons[i][j][1] << "\t" << polygons[i][j][2] << "\n";
        output1 << endl;

        // close the shape
        for (unsigned int j=0;j<polygons[i].size();j++) {
            output1 << "l_0 polygon " << labels[i] << "\n0\n0\n12\n";
            output1 << polygons[i][j][0] << "\t" << polygons[i][j][1] << "\t" << polygons[i][j][2] << "\n";
            output1 << polygons[i][j][0] << "\t" << polygons[i][j][1] << "\t" << "3.0" << "\n";
            output1 << polygons[i][(j+1)%(polygons[i].size())][0] << "\t" << polygons[i][(j+1)%(polygons[i].size())][1] << "\t" << "3.0" << "\n";
            output1 << polygons[i][(j+1)%(polygons[i].size())][0] << "\t" << polygons[i][(j+1)%(polygons[i].size())][1] << "\t" << polygons[i][(j+1)%(polygons[i].size())][2] << "\n";
            output1 << endl;
        }

//          void Problem::writeINP(int hi, int hj) {


        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> roof;
        vector<Point> grid1;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // calcul des points sur les surfaces
        // boucle sur les surfaces

        gridTriangle2(polygons[i], maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

    }   // boucle sur les polygons


    // creation du fichier .inp
    // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

    double gap = 0.05;

    for (unsigned int i=0;i < detectorPoint.size(); i++) {

        output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

        output3 << detectorArea[i] << endl;

    }

    // flushing everything
    //dotInp.flush();
    //dotInp.close();

    output.flush();
    output1.flush();
    output2.flush();
    output3.flush();

    // écriture dans les fichiers
    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // effacement des fichier .oct

    string fileOct = file + ".oct";
    remove(fileOct.c_str());
    //remove(file3D.c_str() );
    //remove(fileRad.c_str());
    //remove(fileInp.c_str());
    remove(fileDet.c_str());
    //remove(fileOut.c_str());

    // fin de l'evaluation

    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,totalRadiation);

    return fitness;

}

double ProjetHuang::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        //return volume(alleles) - 1.1*volumeMax;
        return alleles[0]-alleles[6];

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        //return -volume(alleles) + 0.9*volumeMax;
        return alleles[1]-alleles[6];

    }
    else if (k==2) return alleles[2]-alleles[7];
    else if (k==3) return alleles[3]-alleles[2];
    else if (k==4) return alleles[3]-alleles[4];
    else if (k==5) return alleles[4]-alleles[7];
    else if (k==6) return alleles[5]-alleles[7];
    else if (k==7) return alleles[5]-alleles[6];
    // fin du premier batiment
    else if (k==8) return alleles[8]-alleles[13];
    else if (k==9) return alleles[9]-alleles[13];
    else if (k==10) return alleles[10]-alleles[13];
    else if (k==11) return alleles[10]-alleles[14];
    else if (k==12) return alleles[11]-alleles[14];
    else if (k==13) return alleles[12]-alleles[14];
    else if (k==14) return alleles[12]-alleles[13];
    // fin du deuxième batiment
    else if (k==15) return alleles[24]-alleles[29];
    else if (k==16) return alleles[25]-alleles[29];
    else if (k==17) return alleles[25]-alleles[30];
    else if (k==18) return alleles[26]-alleles[30];
    else if (k==19) return alleles[27]-alleles[30];
    else if (k==20) return alleles[28]-alleles[30];
    else if (k==21) return alleles[28]-alleles[29];
    // fin du troisième batiment
    else if (k==22) return alleles[15]-alleles[22];
    else if (k==23) return alleles[16]-alleles[22];
    else if (k==24) return alleles[16]-alleles[23];
    else if (k==25) return alleles[17]-alleles[23];
    else if (k==26) return alleles[18]-alleles[23];
    else if (k==27) return alleles[19]-alleles[23];
    else if (k==28) return alleles[19]-alleles[22];
    else if (k==29) return alleles[20]-alleles[19];
    else if (k==30) return alleles[20]-alleles[21];
    else if (k==31) return alleles[21]-alleles[22];
    // fin du quatrième batiment
    else return 0.0;

}

// *************************************************************
// Optimisation - Projet Wagner
// *************************************************************

ProjetWagner::ProjetWagner(string varyingFile):Radiance() {

    skyFile = "BaselSkyHH.rad";
    groundFile = "sol.rad";
    environmentFile = "batimentFixes.rad";
    environment = true;
    this->varyingFile = varyingFile;

    maxDetectorArea = 3.0;

    for (int i=0; i < 11 ; i++) {

            // limite en position x
            maxVector.push_back( 10.0 );
            minVector.push_back( 0.0 );
    }

    nConstraints = 0;

    cout << "\nProblem: ProjetWagner" << endl;

}

vector<double> ProjetWagner::evaluate(unsigned int id, vector<double> alleles) {

    // chop precision in alleles for having an integer
    for (unsigned int i=0; i<alleles.size(); i++) alleles[i] = floor( alleles[i] ); // on garde que la partie entière

    // timer start
    time_t start, end;
    time(&start); // start timer

    // fichiers pour enregister
    ostringstream ss;
    ss << "building" << id;
    string file = ss.str();

    string file3D =  file + ".shape";
    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // choix de la précision (12 digits en tout - avant et après la virgule)
    ostringstream output, output1, output2, output3;
    output << setprecision(12);
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);

    // open the output files
    //ofstream dotInp(fileInp.c_str(), ios::trunc);
    //dotInp << setprecision(12);

    // création de la surface de base

    map<string, vector<Point> >  quad, quadMes, tri, triMes;

    vector<double> detectorArea; // vecteur contenant les aires des surfaces de détection
    vector<Point>  detectorPoint;
    vector<Point>  normalVector;

    // lecture du fichier radiance en input décomposé en triangles

    string tampon;
    char tampon1[200];

    // chargement des données a afficher
    ifstream input1 (varyingFile.c_str(), ios::binary | ios::in );

    // test d'ouverture
    if (!input1.is_open()) throw(string("Error opening: " + varyingFile));

    while (!input1.eof()) {

        input1 >> tampon; // varying

        if ( tampon[0] == '#' ) input1.getline(tampon1, 200,'\n');
        else if ( tampon == "varyingMeasured") {

            string label;
            vector<Point> polygon;

            if (input1.eof()) break;
            input1 >> tampon; // polygon
            input1 >> tampon; // label
            label = tampon;
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of coordinates (3x3)


            unsigned int number = atoi(tampon.c_str())/3; // number of Points
            for (unsigned int i=0; i<number; i++) {

                input1 >> tampon;
                double x = atof(tampon.c_str());
                input1 >> tampon;
                double y = atof(tampon.c_str());
                input1 >> tampon;
                int z = atoi(tampon.substr(1).c_str());

                polygon.push_back( Point(x,y,alleles[z]) );

            }
            if (number == 4) quadMes.insert( pair<string,vector<Point> > (label, polygon) );
            else if (number == 3) triMes.insert( pair<string,vector<Point> > (label, polygon) );

        }
        else if ( tampon == "varyingNonMeasured") {

            string label;
            vector<Point> polygon;

            if (input1.eof()) break;
            input1 >> tampon; // polygon
            input1 >> tampon; // label
            label = tampon;
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of coordinates (3x3)


            unsigned int number = atoi(tampon.c_str())/3; // number of Points
            for (unsigned int i=0; i<number; i++) {

                input1 >> tampon;
                double x = atof(tampon.c_str());
                input1 >> tampon;
                double y = atof(tampon.c_str());
                input1 >> tampon;
                int z = atoi(tampon.substr(1).c_str());

                polygon.push_back( Point(x,y,alleles[z]) );

            }
            if (number == 4) quad.insert( pair<string,vector<Point> > (label, polygon) );
            else if (number == 3) tri.insert( pair<string,vector<Point> > (label, polygon) );

        }
        else if ( tampon == "nonvaryingMeasured") {

            string label;
            vector<Point> polygon;

            if (input1.eof()) break;
            input1 >> tampon; // polygon
            input1 >> tampon; // label
            label = tampon;
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of points (3x3)

            unsigned int number = atoi(tampon.c_str())/3;
            for (unsigned int i=0; i<number; i++) {

                input1 >> tampon;
                double x = atof(tampon.c_str());
                input1 >> tampon;
                double y = atof(tampon.c_str());
                input1 >> tampon;
                double z = atof(tampon.c_str());

                polygon.push_back( Point(x,y,z) );

            }
            if (number == 4) quadMes.insert( pair<string,vector<Point> > (label, polygon) );
            else if (number == 3) triMes.insert( pair<string,vector<Point> > (label, polygon) );
        }
        else if ( tampon == "nonvaryingNonMeasured") {

            string label;
            vector<Point> polygon;

            if (input1.eof()) break;
            input1 >> tampon; // polygon
            input1 >> tampon; // label
            label = tampon;
            input1 >> tampon; // 0
            input1 >> tampon; // 0
            input1 >> tampon; // number of points (3x3)

            unsigned int number = atoi(tampon.c_str())/3;
            for (unsigned int i=0; i<number; i++) {

                input1 >> tampon;
                double x = atof(tampon.c_str());
                input1 >> tampon;
                double y = atof(tampon.c_str());
                input1 >> tampon;
                double z = atof(tampon.c_str());

                polygon.push_back( Point(x,y,z) );

            }
            if (number == 4) quad.insert( pair<string,vector<Point> > (label, polygon) );
            else if (number == 3) tri.insert( pair<string,vector<Point> > (label, polygon) );
        }
    }

    input1.close();

 // void prepareRadiance() {

    output1 << "void plastic measured" << "\n0\n0\n5 0.07 0.07 0.07 0.07 0\n\n" << endl;
    output1 << "void plastic nonMeasured" << "\n0\n0\n5 0.2 0.2 0.2 0 0\n\n" << endl;

    for (map<string,vector<Point> >::iterator it=triMes.begin() ; it != triMes.end(); it++ ) {

    //cout << (*it).first << " => " << (*it).second << endl;
        // écriture dans le fichier .shape
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output << (*it).first << "\n";
            output << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n\n";
        }
        output << endl;

        // creation du fichier .RAD
        output1 << "measured polygon " << (*it).first << "\n0\n0\n9\n";
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output1 << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n";
        }
        output1 << endl;

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> roof;
        vector<Point> grid1;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // calcul des points sur les surfaces
        // boucle sur les surfaces

        gridTriangle2((*it).second, maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

    }   // boucle sur les triangles variables
    for (map<string,vector<Point> >::iterator it=tri.begin() ; it != tri.end(); it++ ) {

        // creation du fichier .RAD
        output1 << "nonMeasured polygon " << (*it).first << "\n0\n0\n9\n";
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output1 << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n";
        }
        output1 << endl;

    }   // boucle sur les triangles variables
    for (map<string,vector<Point> >::iterator it=quadMes.begin() ; it != quadMes.end(); it++ ) {

    //cout << (*it).first << " => " << (*it).second << endl;
        // écriture dans le fichier .shape
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output << (*it).first << "\n";
            output << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n\n";
        }
        output << endl;

        // creation du fichier .RAD
        output1 << "measured polygon " << (*it).first << "\n0\n0\n12\n";
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output1 << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n";
        }
        output1 << endl;

        // calcul des points de surfaces (facade, roof): determination de la grille (detectorPoint)

        vector<Point> roof;
        vector<Point> grid1;
        double n1=0., n2=0., n3=0.;
        double pointArea=0.0;

        // calcul des points sur les surfaces
        // boucle sur les surfaces

        gridRectangle2((*it).second, maxDetectorArea, grid1, pointArea, n1, n2, n3);

        detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
        normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
        detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

    }   // boucle sur les triangles variables
    for (map<string,vector<Point> >::iterator it=quad.begin() ; it != quad.end(); it++ ) {

        // creation du fichier .RAD
        output1 << "nonMeasured polygon " << (*it).first << "\n0\n0\n12\n";
        for (unsigned int j=0; j<(*it).second.size(); j++) {
            output1 << (*it).second[j][0] << "\t" << (*it).second[j][1] << "\t" << (*it).second[j][2] << "\n";
        }
        output1 << endl;

    }   // boucle sur les triangles variables


    // creation du fichier .inp
    // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

    double gap = 0.05;

    for (unsigned int i=0;i < detectorPoint.size(); i++) {

        output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;

        output3 << detectorArea[i] << endl;

    }

    // flushing everything
    //dotInp.flush();
    //dotInp.close();

    output.flush();
    output1.flush();
    output2.flush();
    output3.flush();

    // écriture dans les fichiers
    writeResultsOver(file3D, output);
    writeResultsOver(fileRad, output1);
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications

    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats

    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // détermination de la fitness de l'individu par sommation pondérée

    double totalRadiation=0., totalSurface=0.;

    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // effacement des fichier .oct

    string fileOct = file + ".oct";
    remove(fileOct.c_str());
    //remove(file3D.c_str() );
    //remove(fileRad.c_str());
    //remove(fileInp.c_str());
    //remove(fileDet.c_str());
    //remove(fileOut.c_str());

    // fin de l'evaluation

    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,totalRadiation);

    return fitness;

}
//
//double ProjetHuang::constraint(unsigned int k, vector<double> alleles) {
//
//    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%
//
//        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
//        //return volume(alleles) - 1.1*volumeMax;
//        return alleles[0]-alleles[6];
//
//    }
//    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%
//
//        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
//        //return -volume(alleles) + 0.9*volumeMax;
//        return alleles[1]-alleles[6];
//
//    }
//    else if (k==2) return alleles[2]-alleles[7];
//    else if (k==3) return alleles[3]-alleles[2];
//    else if (k==4) return alleles[3]-alleles[4];
//    else if (k==5) return alleles[4]-alleles[7];
//    else if (k==6) return alleles[5]-alleles[7];
//    else if (k==7) return alleles[5]-alleles[6];
//    // fin du premier batiment
//    else if (k==8) return alleles[8]-alleles[13];
//    else if (k==9) return alleles[9]-alleles[13];
//    else if (k==10) return alleles[10]-alleles[13];
//    else if (k==11) return alleles[10]-alleles[14];
//    else if (k==12) return alleles[11]-alleles[14];
//    else if (k==13) return alleles[12]-alleles[14];
//    else if (k==14) return alleles[12]-alleles[13];
//    // fin du deuxième batiment
//    else if (k==15) return alleles[24]-alleles[29];
//    else if (k==16) return alleles[25]-alleles[29];
//    else if (k==17) return alleles[25]-alleles[30];
//    else if (k==18) return alleles[26]-alleles[30];
//    else if (k==19) return alleles[27]-alleles[30];
//    else if (k==20) return alleles[28]-alleles[30];
//    else if (k==21) return alleles[28]-alleles[29];
//    // fin du troisième batiment
//    else if (k==22) return alleles[15]-alleles[22];
//    else if (k==23) return alleles[16]-alleles[22];
//    else if (k==24) return alleles[16]-alleles[23];
//    else if (k==25) return alleles[17]-alleles[23];
//    else if (k==26) return alleles[18]-alleles[23];
//    else if (k==27) return alleles[19]-alleles[23];
//    else if (k==28) return alleles[19]-alleles[22];
//    else if (k==29) return alleles[20]-alleles[19];
//    else if (k==30) return alleles[20]-alleles[21];
//    else if (k==31) return alleles[21]-alleles[22];
//    // fin du quatrième batiment
//    else return 0.0;
//
//}

// *************************************************************
// Optimisation - Fourrier basis
// *************************************************************

FourierBasis::FourierBasis():Radiance() {

    h  = 10.0; // 10m mean height
    Tx = 20.0;
    Ty = 30.0;
    Nx = 5;
    Ny = 5;

    minHorizCut = 0.0;
    maxHorizCut = 1000.0;

    amplitude = 2.0;

    skyFile = "BaselSky.rad";
    groundFile = "FlatSite.rad";
    ground = true;
    environmentFile = "3df_matthaeus4_bld2.rad";
    environment = false;

    deleteFiles = false;
    surroundings = true;
    if (surroundings == true) ground = false;

    xinit = 0.0;
    yinit = 0.0;
    hbase = 0.0;

    maxDetectorArea = 0.01; // for the subdivision of surface in measurement points
    maxGridSpacing =  0.5; // discretisation of the continuous function by steps of

    if (surroundings == false) {
        maxVector.push_back( h ); // basic height (A00)
        minVector.push_back( 0.0 );
    }

    for (unsigned int i=0; i < (Nx*(Ny-1)/2)+(Nx-1)/2 ; i++) { // matrice A

            maxVector.push_back( amplitude );
            minVector.push_back( -amplitude );

    }

//    maxVector.push_back( 1.0 ); // element that defines the orientation
//    minVector.push_back( 0.0 );

    for (unsigned int i=0; i < (Nx*(Ny-1)/2)+(Nx-1)/2 ; i++) { // matrice B

            maxVector.push_back( amplitude );
            minVector.push_back( -amplitude );

    }

    // constraints in volume... 80% of max volume
    nConstraints = 2; //2; // can change
    if (surroundings == false) volumeMax = 0.8*h*Tx*Ty;
    else volumeMax = 0.2*h*(double(Nx)*(Tx/double(Nx-1)))*(double(Ny)*(Ty/double(Ny-1)));

    cout << "\nProblem: FourierBasis" << endl;

}

vector<double> FourierBasis::evaluate(unsigned int id, vector<double> alleles) {

    // start of evaluation
    time_t start, end;
    time(&start); // start timer

    // verification de la taille de alleles
    if (alleles.size() == (Nx*Ny-1)) alleles.insert(alleles.begin(), 0.0);

    vector<double> detectorArea;

    stringstream ss;
    ss << "building" << id;

    string file = ss.str();

    string fileRad = file + ".rad";
    string fileInp = file + ".inp";
    string fileDet = file + ".det";
    string fileOut = file + ".out";
    string fileDat = file + ".dat";

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    ostringstream output1, output2, output3, output4;

    // choix de la précision (12 digits en tout - avant et après la virgule)
    output1 << setprecision(12);
    output2 << setprecision(12);
    output3 << setprecision(12);
    output4 << setprecision(12);

    unsigned int nbDivx, nbDivy;
    double deltax, deltay;

    // use of maxGridSpacing
    if ( surroundings == true ) { //:D case with a cut, no repeating
        nbDivx = (unsigned int) ceil((double(Nx)*(Tx/double(Nx-1)))/maxGridSpacing);
        nbDivy = (unsigned int) ceil((double(Ny)*(Ty/double(Ny-1)))/maxGridSpacing);
        deltax = ((double(Nx)*(Tx/double(Nx-1)))/nbDivx);
        deltay = ((double(Ny)*(Ty/double(Ny-1)))/nbDivy);
    }
    else { //:J case with a repeat of the forms
        nbDivx = (unsigned int) ceil(Tx/maxGridSpacing);
        nbDivy = (unsigned int) ceil(Ty/maxGridSpacing);
        deltax = (Tx/nbDivx);
        deltay = (Ty/nbDivy);
    }

    // création des points de mesure
    for (unsigned int v = 0; v <= nbDivy; v++) {
        for(unsigned int u = 0; u <= nbDivx; u++) {

            // création du fichier .dat pour gensurf
            output4 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            vector<Point> facade;
            vector<Point> roof;
            vector<Point> grid1;
            vector<Point> detectorPoint;
            vector<Point> normalVector;
            double n1=0., n2=0., n3=0.;
            double pointArea=0.0;

            if ( surroundings == false ) {
                if ( (u<nbDivx) && (v == 0) ) {

                    // décomposition des facades

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );
                    facade.push_back( Point( u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                if ( (u<nbDivx) && (v == nbDivy) ) {

                    // décomposition des facades

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );
                    facade.push_back( Point( (u+1)*deltax, v*deltay, 0.0 ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                if ( (v<nbDivy) && (u == 0) ) {

                    // décomposition des facades

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, 0.0 ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                if ( (v<nbDivy) && (u == nbDivx) ) {

                    // décomposition des facades

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    facade.clear();
                    grid1.clear();

                    facade.push_back( Point( u*deltax, v*deltay, 0.0 ) );
                    facade.push_back( Point( u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );
                    facade.push_back( Point( u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );

                    gridTriangle2(facade, sqrt(maxDetectorArea), grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
            }

            if (v != nbDivy && u != nbDivx) { // pas de points dans les coins...

                // décompostion des toits

                if ( v & 1 ) {

                    roof.clear();
                    grid1.clear();

                    roof.push_back( Point(u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point(u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );

                    if (roof[0][2] > 0.0 || roof[1][2] > 0.0 || roof[2][2] > 0.0)
                    gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    roof.clear();
                    grid1.clear();

                    roof.push_back( Point(u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, (v+1)*deltay, fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) ) );

                    if (roof[0][2] > 0.0 || roof[1][2] > 0.0 || roof[2][2] > 0.0)
                    gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }
                else {

                    roof.clear();
                    grid1.clear();

                    roof.push_back( Point(u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, (v+1)*deltay, fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) ) );
                    roof.push_back( Point(u*deltax, (v+1)*deltay, fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) ) );

                    if (roof[0][2] > 0.0 || roof[1][2] > 0.0 || roof[2][2] > 0.0)
                    gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                    roof.clear();
                    grid1.clear();

                    roof.push_back( Point(u*deltax, v*deltay, fourierSeries(double(u*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, v*deltay, fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) ) );
                    roof.push_back( Point((u+1)*deltax, (v+1)*deltay, fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) ) );

                    if (roof[0][2] > 0.0 || roof[1][2] > 0.0 || roof[2][2] > 0.0)
                    gridTriangle2(roof, maxDetectorArea, grid1, pointArea, n1, n2, n3);

                    detectorArea.insert( detectorArea.end(), grid1.size(), pointArea );
                    normalVector.insert( normalVector.end(), grid1.size(), Point(n1, n2, n3) );
                    detectorPoint.insert( detectorPoint.end(), grid1.begin(), grid1.end() );

                }

            } // end if les points sont les derniers

            double gap = 0.05; //distance between surface and point

            // creation du fichier .inp
            // celui-ci contient les points de mesure ainsi que les normales au faces en ces points

            for (unsigned int i=0;i < detectorPoint.size(); i++) {

                output2 << detectorPoint[i][0] + normalVector[i][0]*gap << "\t" << detectorPoint[i][1] + normalVector[i][1]*gap << "\t"
                        << detectorPoint[i][2] + normalVector[i][2]*gap << "\t" << normalVector[i][0] << "\t" << normalVector[i][1] << "\t" << normalVector[i][2] << endl;
                output3 << detectorArea[i] << endl;

            }

/* utilisation de la derivee de la fonction pour le vector normal, implique des petites différences...
            double gap = 0.1;

            if (!(v&1) && u != nbDivx && v != nbDivy) { // sortie des points de mesures (un par triangle)
                vector<double> n1 = derivativeFourierSeries(double((u+2./3.)*deltax),double((v+1./3.)*deltay),alleles);
                output3 << gap*n1[0]+(u+2./3.)*deltax << "\t" << gap*n1[1]+(v+1./3.)*deltay << "\t" << gap*n1[2]+fourierSeries(double((u+2./3.)*deltax),double((v+1./3.)*deltay),alleles) << "\t" << n1[0] << "\t" << n1[1] << "\t" << n1[2] << endl;
                vector<double> n2 = derivativeFourierSeries(double((u+1./3.)*deltax),double((v+2./3.)*deltay),alleles);
                output3 << gap*n2[0]+(u+1./3.)*deltax << "\t" << gap*n2[1]+(v+2./3.)*deltay << "\t" << gap*n2[2]+fourierSeries(double((u+1./3.)*deltax),double((v+2./3.)*deltay),alleles) << "\t" << n2[0] << "\t" << n2[1] << "\t" << n2[2] << endl;
            }
            if ((v&1) && u != nbDivx && v != nbDivy) { // sortie des points de mesures (un par triangle)
                vector<double> n1 = derivativeFourierSeries(double((u+1./3.)*deltax),double((v+1./3.)*deltay),alleles);
                output3 << gap*n1[0]+(u+1./3.)*deltax << "\t" << gap*n1[1]+(v+1./3.)*deltay << "\t" << gap*n1[2]+fourierSeries(double((u+1./3.)*deltax),double((v+1./3.)*deltay),alleles) << "\t" << n1[0] << "\t" << n1[1] << "\t" << n1[2] << endl;
                vector<double> n2 = derivativeFourierSeries(double((u+2./3.)*deltax),double((v+2./3.)*deltay),alleles);
                output3 << gap*n2[0]+(u+2./3.)*deltax << "\t" << gap*n2[1]+(v+2./3.)*deltay << "\t" << gap*n2[2]+fourierSeries(double((u+2./3.)*deltax),double((v+2./3.)*deltay),alleles) << "\t" << n2[0] << "\t" << n2[1] << "\t" << n2[2] << endl;
            }
*/
        }
    }

    output2.flush();
    output3.flush();

    // ecriture des .inp, .det, .dat
    writeResultsOver(fileInp, output2);
    writeResultsOver(fileDet, output3); // sortie pour vérifications
    writeResultsOver(fileDat, output4);

    // écriture du fichier .rad2
    stringstream leSS;
    if ( surroundings == true ) {
        leSS << "gensurf roof roof1 " << fileDat << " " << fileDat << " " << fileDat << " " << (unsigned int) ceil((double(Nx)*(Tx/double(Nx-1)))/maxGridSpacing) << " " << (unsigned int) ceil((double(Ny)*(Ty/double(Ny-1)))/maxGridSpacing) << " > " << fileRad << "2" << endl;
    }
    else {
        leSS << "gensurf roof roof1 " << fileDat << " " << fileDat << " " << fileDat << " " << (unsigned int) ceil(Tx/maxGridSpacing) << " " << (unsigned int) ceil(Ty/maxGridSpacing) << " > " << fileRad << "2" << endl;
    }
    // issue system command
    if ( system(leSS.str().c_str()) != 0 ) throw (string("Error in gensurf (.rad2)"));

    // écriture du fichier .rad
    output1 << "void plastic roof\n" << "0\n0\n" << "5\t0.2\t0.2\t0.2\t0.0\t0.0\n" << endl;
    if ( surroundings == false ) {
        for (unsigned int v = 0; v < nbDivy; v++) {
            for(unsigned int u = 0; u < nbDivx; u++) {

                if ( (u<nbDivx) && (v == 0) ) {
                    output1 << "roof polygon facade1-" << v << "-" << u << endl;
                    output1 << "0\n0\n12" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << 0.0 << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << 0.0 << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                }
                if ( (u<nbDivx) && (v == nbDivy-1) ) {
                    output1 << "roof polygon facade2-" << v << "-" << u << endl;
                    output1 << "0\n0\n12" << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << 0.0 << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << 0.0 << endl;
                }
                if ( (v<nbDivy) && (u == 0) ) {
                    output1 << "roof polygon facade3-" << v << "-" << u << endl;
                    output1 << "0\n0\n12" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << 0.0 << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << 0.0 << endl;
                }
                if ( (v<nbDivy) && (u == nbDivx-1) ) {
                    output1 << "roof polygon facade4-" << v << "-" << u << endl;
                    output1 << "0\n0\n12" << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << 0.0 << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << 0.0 << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                }
                // triangle des toits
                if ( v & 1 ) {
                    output1 << "roof polygon roof-" << v << "-" << u << "-0" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << "\nroof polygon roof-" << v << "-" << u << "-1" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << "\n" << endl;
                }
                else {
                    output1 << "roof polygon roof-" << v << "-" << u << "-0" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << "\nroof polygon roof-" << v << "-" << u << "-1" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << "\n" << endl;
                }

            }
        }
    }
    else { // :J case, pas de base, seulements des "toits"
        for (int v = -nbDivy; v != int(2*nbDivy); v++) {
            for(int u = -nbDivx; u != int(2*nbDivx); u++) {
                // triangle des toits
                if ( v & 1 ) {
                    output1 << "roof polygon roof-" << v << "-" << u << "-0" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << "\nroof polygon roof-" << v << "-" << u << "-1" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << "\n" << endl;
                }
                else {
                    output1 << "roof polygon roof-" << v << "-" << u << "-0" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << u*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double(u*deltax),double((v+1)*deltay),alleles) << endl;
                    output1 << "\nroof polygon roof-" << v << "-" << u << "-1" << endl;
                    output1 << "0\n0\n9" << endl;
                    output1 << u*deltax << "\t" << v*deltay << "\t" << fourierSeries(double(u*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << v*deltay << "\t" << fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) << endl;
                    output1 << (u+1)*deltax << "\t" << (v+1)*deltay << "\t" << fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) << "\n" << endl;
                }

            }
        }

    }
    //output1 << "!gensurf roof roof1 " << fileDat << " " << fileDat << " " << fileDat << " " << (unsigned int) ceil((double(size-1)*(Tx/double(size)))/maxGridSpacing) << " " << (unsigned int) ceil((double(size-1)*(Ty/double(size)))/maxGridSpacing) << endl;
    output1.flush();

    // writing .rad file
    writeResultsOver(fileRad, output1);

    // merge des deux fichiers
//    leSS.str("");
//    leSS << "cat " << fileRad3 << " " << fileRad2 << " > " << fileRad << endl;
//    system(leSS.str().c_str());

    // run the simulation on the files
    runRadiance(file);

    // creation des vecteurs résultats
    // lecture des résultats
    vector<double> radianceResults;
    readResultsInLine(fileOut, radianceResults);

    // effacement des fichiers
    remove(string(file + ".oct").c_str());
    remove(string(file + ".det").c_str());
    if (deleteFiles == true) {
        remove(string(file + ".rad").c_str());
        remove(string(file + ".inp").c_str());
        remove(string(file + ".out").c_str());
        remove(string(file + ".dat").c_str());
    }

    // détermination de la fitness de l'individu par sommation pondérée
    double totalRadiation=0., totalSurface=0.;
    for (unsigned int i=0;i<radianceResults.size()/3;i++) {

        // RED channel
        totalRadiation+=radianceResults[i*3]*0.265*detectorArea[i];

        // GREEN channel
        totalRadiation+=radianceResults[i*3+1]*0.67*detectorArea[i];

        // BLUE channel
        totalRadiation+=radianceResults[i*3+2]*0.065*detectorArea[i];

        // total surface
        totalSurface+=detectorArea[i];

    }

    // sortie du volume
    if ( surroundings == false ) {
        cout << "volume (prismes): " << volume(alleles) << endl;
        cout << "volume (theorique): " << volumeFourierSeries(alleles) << endl;
    }

    // évaluation du gain d'energie
    vector<double> fitness;
    fitness.assign(1,1.0*totalRadiation);
    cout << "fitness: " << fitness[0] << endl;

    // writing the fitness in a file
    ostringstream ossFitness;
    ossFitness << setprecision(12);
    ossFitness << fitness[0] << "\n";
    ossFitness.flush();
    writeResultsOver(string(file + ".fit"),ossFitness);

    // vérification de la longueur du fichier (entrée = sortie)
    if ( detectorArea.size() != radianceResults.size()/3 ) {

        throw (string("Error in simulating: " + string(""+id) + " - .Det file not the same size as .Out file"));

    }

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    return fitness;

}

void FourierBasis::runRadiance(string file) {

    // début de la simulation
    string leString;

    // preparation and start of OCONV
    leString = "oconv ";
    leString += "-w ";
    leString += file;
    leString += ".rad ";
    leString += skyFile;
    leString += " ";
    if (ground==true) {
        leString += groundFile;
        leString += " ";
    }
    if (environment==true) {
        leString += environmentFile;
        leString += " ";
    }
    leString += "> ";
    leString += file;
    leString += ".oct";

    // issue system command
    int errCode = system(leString.c_str());
    if ( errCode != 0 ) throw (string("Error in oconv (.oct)"));

    // preparation and start of RTRACE -I -ab 2 -ad 1000 -as 250 -dc .9 -dj .66

    leString = "rtrace ";
    leString += "-h -w -I -ab 1 -ad 1024 -as 512 -aa 0.1 -ar 256 ";
    //leString += "-h -w -I -ab 2 -ad 2048 -as 512 -aa 0.15 -ar 256 ";
    //leString += "-h -w -I -ab 2 -ad 256 -as 128 -aa 0.005 -ar 256 ";
    leString += file;
    leString += ".oct < ";
    leString += file;
    leString += ".inp > "; // white light from sky (1,1,1), multiplied by the coloured reflectance of the surfaces
    leString += file;
    leString += ".out";

    // issue system command
    errCode = system(leString.c_str());
    if ( errCode != 0 ) throw (string("Error in rtrace (.out)"));

}

double FourierBasis::getA(int k, int l, vector<double> alleles) { return alleles[ abs(k + l*int(Nx)) ]; }

double FourierBasis::getB(int k, int l, vector<double> alleles) {
    if (k==0 && l==0) return 0.0;
    else if ( k < 0 && l==0) return -alleles[ abs(k + l*int(Nx)) + (Ny-1)/2*Nx + (Nx-1)/2 ];
    else return alleles[ abs(k + l*int(Nx)) + (Ny-1)/2*Nx + (Nx-1)/2 ];
}

double FourierBasis::fourierSeries(double x, double y, vector<double> alleles) {

    // verification de la taille de alleles
    if (alleles.size() == (Nx*Ny-1)) alleles.insert(alleles.begin(), 0.0);

    // Pi
    double Pi = 4.0*atan(1.0);

    // compute the fourier series at point (x,y)
    double series = 0.0;
    for (int k=-(int(Nx)-1)/2;k!=(int(Nx)-1)/2+1;k++) {
        for (int l=0;l!=(int(Ny)-1)/2+1;l++) {
            //cout << "A[" << k << "," << l << "]: " << getA(k,l,alleles) << "\t" << "B[" << k << "," << l << "]: " << getB(k,l,alleles) << endl;
            series += getA(k,l,alleles)*cos(2.0*Pi*k*(x/Tx)*((double(Nx)-1.0)/double(Nx)) +2*Pi*l*(y/Ty)*((double(Ny)-1.0)/double(Ny))) + getB(k,l,alleles)*sin(2.0*Pi*k*(x/Tx)*((double(Nx)-1.0)/double(Nx)) +2.0*Pi*l*(y/Ty)*((double(Ny)-1.0)/double(Ny)));
        }
    }

    return min(max(series, minHorizCut), maxHorizCut);

}

vector<double> FourierBasis::derivativeFourierSeries(double x, double y, vector<double> alleles) {

    // verification de la taille de alleles
    if (alleles.size() == (Nx*Ny-1)) alleles.insert(alleles.begin(), 0.0);

    // matrix sizes
    const unsigned int size = (unsigned int) sqrt(alleles.size()/2.0);

    // look at the sign of x and y using the first element of matrix B which is of no importance in Fourier Series
    unsigned int B00i = alleles.size()/2;
    double signCoeff =  alleles[B00i]; // between 0 and 1
    double signX, signY, transX, transY;
    if (signCoeff < 0.25) { signX = 1.0; signY = 1.0; transX = 0.0; transY = 0.0; }
    else if (signCoeff < 0.5) { signX = -1.0; signY = 1.0; transX = Tx/double(size); transY = 0.0; }
    else if (signCoeff < 0.75) { signX = 1.0; signY = -1.0; transX = 0.0; transY = Ty/double(size); }
    else { signX = -1.0; signY = -1.0; transX = Tx/double(size); transY = Ty/double(size); }

    // NOT CORRECTED!!!!

    // note: the values inside the alleles are the A and B matrix elements
    double A[size][size], B[size][size];

    for (unsigned int i=0; i<size; i++) {
        for (unsigned int j=0;j<size; j++) {
            A[i][j]=alleles[i*size + j];
            B[i][j]=alleles[size*size +i*size +j];
        }
    }

    vector<double> normal;
    normal.assign(3, 0.0);

    // Pi
    double Pi = 4.0*atan(1.0);

    // compute the derivative with respect to x of the fourier series at point (x,y)
    for (unsigned int n=0;n<size;n++) {
        for (unsigned int m=0;m<size;m++) {
            normal[0] += - (-A[n][m]*sin(2*Pi/Tx*n*signX*(x+transX) +2*Pi/Ty*m*signY*(y+transY))*2*Pi/Tx*n + B[n][m]*cos(2*Pi/Tx*n*signX*(x+transX) +2*Pi/Ty*m*signY*(y+transY))*2*Pi/Tx*n);
        }
    }

    // compute the derivative with respect to y of the fourier series at point (x,y)
    for (unsigned int n=0;n<size;n++) {
        for (unsigned int m=0;m<size;m++) {
            normal[1] += - (-A[n][m]*sin(2*Pi/Tx*n*signX*(x+transX) +2*Pi/Ty*m*signY*(y+transY))*2*Pi/Ty*m + B[n][m]*cos(2*Pi/Tx*n*signX*(x+transX) +2*Pi/Ty*m*signY*(y+transY))*2*Pi/Ty*m);
        }
    }

    normal[2] = 1.0;

    // normalisation
    double norm = sqrt(pow(normal[0], 2.0) + pow(normal[1], 2.0) + pow(normal[2], 2.0));
    normal[0]/=norm;
    normal[1]/=norm;
    normal[2]/=norm;

    return normal;

}

double FourierBasis::volumeFourierSeries(vector<double> alleles) {

//    // matrix sizes
//    const unsigned int size = (unsigned int) sqrt(alleles.size()/2.0);
//
//    // look at the sign of x and y using the first element of matrix B which is of no importance in Fourier Series
//    unsigned int B00i = alleles.size()/2;
//    double signCoeff = alleles[B00i]; // between 0 and 1
//    double signX, signY, transX, transY;
//    if (signCoeff < 0.25) { signX = 1.0; signY = 1.0; transX = 0.0; transY = 0.0; }
//    else if (signCoeff < 0.5) { signX = -1.0; signY = 1.0; transX = Tx/double(size); transY = 0.0; }
//    else if (signCoeff < 0.75) { signX = 1.0; signY = -1.0; transX = 0.0; transY = Ty/double(size); }
//    else { signX = -1.0; signY = -1.0; transX = Tx/double(size); transY = Ty/double(size); }
//
//
//// NOT CORRECTED!!!!
//
//
//    // domain boundaries
//    double x1,x2;
//    if (signX < 0) {
//        x1 = -Tx;
//        if (surroundings == false ) x2 = -Tx/double(size);
//        else x2 = 0.0;
//    }
//    else {
//        x1 = 0.0;
//        if (surroundings == false) x2 = double(size-1)*(Tx/double(size));
//        else x2 = Tx;
//    }
//
//    double y1,y2;
//    if (signY < 0) {
//        y1 = -Ty;
//        if (surroundings == false) y2 = -Ty/double(size);
//        else y2 = -Ty;
//    }
//    else {
//        y1 = 0.0;
//        if (surroundings == false) y2 = double(size-1)*(Ty/double(size));
//        else y2 = Ty;
//    }
//
//    //cout << "signX: " << signX << "\tsignY: " << signY << endl;
//    //cout << "x1: " << x1 << "\tx2: " << x2 << "\ty1: " << y1 << "\ty2: " << y2 << endl;

    // verification de la taille de alleles
    if (alleles.size() == (Nx*Ny-1)) alleles.insert(alleles.begin(), 0.0);

    // Pi
    double Pi = 4.0*atan(1.0);

    double x1,x2,y1,y2;
    if (surroundings == false) {
        x1 = 0.0;
        x2 = Tx;
        y1 = 0.0;
        y2 = Ty;
    }
    else {
        x1 = 0.0;
        x2 = Tx * (double(Nx)/(double(Nx)-1));
        y1 = 0.0;
        y2 = Ty * (double(Ny)/(double(Ny)-1));
    }

    // compute the double integral of fourier series over the domain of the building
    double volume = 0.0;
    for (int k=-(int(Nx)-1)/2;k!=(int(Nx)-1)/2+1;k++) {
        for (int l=0;l!=(int(Ny)-1)/2+1;l++) {

            if (k==0 && l==0) volume += getA(k,l,alleles)*(x1 - x2)*(y1 - y2);
            else if (k==0)    volume += (double(Ny)*Ty*(-x1 + x2)*(getB(k,l,alleles)*cos((2.0*l*(-1.0 + double(Ny))*Pi*y1)/(double(Ny)*Ty)) -
                                        getB(k,l,alleles)*cos((2.0*l*(-1. + double(Ny))*Pi*y2)/(double(Ny)*Ty)) +
                                        getA(k,l,alleles)*(-sin((2.*l*(-1. + double(Ny))*Pi*y1)/(double(Ny)*Ty)) + sin((2.*l*(-1. + double(Ny))*Pi*y2)/(double(Ny)*Ty)))))/
                                        (2.*l*(-1. + double(Ny))*Pi);
            else if (l==0)    volume += -(double(Nx)*Tx*(y1 - y2)*(getB(k,l,alleles)*cos((2.*k*(-1. + double(Nx))*Pi*x1)/(double(Nx)*Tx)) -
                                        getB(k,l,alleles)*cos((2.*k*(-1. + double(Nx))*Pi*x2)/(double(Nx)*Tx)) +
                                        getA(k,l,alleles)*(-sin((2.*k*(-1. + double(Nx))*Pi*x1)/(double(Nx)*Tx)) + sin((2.*k*(-1. + double(Nx))*Pi*x2)/(double(Nx)*Tx)))))/
                                        (2.*k*(-1. + double(Nx))*Pi);
            else              volume += (double(Nx)*double(Ny)*Tx*Ty*(-(getA(k,l,alleles)*cos(2.0*Pi*((k*(-1.0 + double(Nx))*x1)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y1)/(double(Ny)*Ty)))) +
                                        getA(k,l,alleles)*cos(2.0*Pi*((k*(-1.0 + double(Nx))*x2)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y1)/(double(Ny)*Ty))) +
                                        getA(k,l,alleles)*cos(2.0*Pi*((k*(-1.0 + double(Nx))*x1)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y2)/(double(Ny)*Ty))) -
                                        getA(k,l,alleles)*cos(2.0*Pi*((k*(-1.0 + double(Nx))*x2)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y2)/(double(Ny)*Ty))) -
                                        getB(k,l,alleles)*sin(2.0*Pi*((k*(-1.0 + double(Nx))*x1)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y1)/(double(Ny)*Ty))) +
                                        getB(k,l,alleles)*sin(2.0*Pi*((k*(-1.0 + double(Nx))*x2)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y1)/(double(Ny)*Ty))) +
                                        getB(k,l,alleles)*sin(2.0*Pi*((k*(-1.0 + double(Nx))*x1)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y2)/(double(Ny)*Ty))) -
                                        getB(k,l,alleles)*sin(2.0*Pi*((k*(-1.0 + double(Nx))*x2)/(double(Nx)*Tx) + (l*(-1.0 + double(Ny))*y2)/(double(Ny)*Ty)))))/
                                        (4.*k*l*(-1.0 + double(Nx))*(-1.0 + double(Ny))*pow(Pi,2.0));
        }
    }

    return volume;

}

double FourierBasis::volume(vector<double> alleles) {

    // verification de la taille de alleles
    if (alleles.size() == (Nx*Ny-1)) alleles.insert(alleles.begin(), 0.0);

    unsigned int nbDivx, nbDivy;
    double deltax, deltay;

    // use of maxGridSpacing
    if ( surroundings == true ) { //:D case with a cut, no repeating
        nbDivx = (unsigned int) ceil((double(Nx)*(Tx/double(Nx-1)))/maxGridSpacing);
        nbDivy = (unsigned int) ceil((double(Ny)*(Ty/double(Ny-1)))/maxGridSpacing);
        deltax = ((double(Nx)*(Tx/double(Nx-1)))/nbDivx);
        deltay = ((double(Ny)*(Ty/double(Ny-1)))/nbDivy);
    }
    else { //:J case with a repeat of the forms
        nbDivx = (unsigned int) ceil(Tx/maxGridSpacing);
        nbDivy = (unsigned int) ceil(Ty/maxGridSpacing);
        deltax = (Tx/nbDivx);
        deltay = (Ty/nbDivy);
    }

    // compute the volume over the domain of the building
    double volume = 0.0;
    for (unsigned int v = 0; v < nbDivy; v++) {
        for(unsigned int u = 0; u < nbDivx; u++) {

            double n1, n2, n3;
            double u1,u2,u3,v1,v2,v3;
            double d;

            // décompostion des toits en prismes
            if ( v & 1 ) {

                // calcul de la normale au plan
                u1 = deltax;
                u2 = 0.0;
                u3 = (fourierSeries(double((u+1)*deltax),double(v*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                v1 = 0.0;
                v2 = deltay;
                v3 = (fourierSeries(double(u*deltax),double((v+1)*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                crossProduct(u1,u2,u3,v1,v2,v3,n1,n2,n3);

                // calcul de d du plan
                d = -n1*u*deltax -n2*v*deltay -n3*fourierSeries(double(u*deltax),double(v*deltay),alleles);
                volume += (deltax*deltay/2.)*(-d - n1*(2*u*deltax+(u+1)*deltax)/3.0 - n2*(2*v*deltay+(v+1)*deltay)/3.0)/n3;

                // calcul de la normale au plan
                u1 = deltax;
                u2 = -deltay;
                u3 = (fourierSeries(double((u+1)*deltax),double(v*deltay),alleles) - fourierSeries(double(u*deltax),double((v+1)*deltay),alleles));
                v1 = deltax;
                v2 = 0.0;
                v3 = (fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles) - fourierSeries(double(u*deltax),double((v+1)*deltay),alleles));
                crossProduct(u1,u2,u3,v1,v2,v3,n1,n2,n3);

                // calcul de d du plan
                d = -n1*u*deltax -n2*v*deltay -n3*fourierSeries(double(u*deltax),double(v*deltay),alleles);
                volume += (deltax*deltay/2.)*(-d - n1*(2*u*deltax+(u+1)*deltax)/3.0 - n2*(2*v*deltay+(v+1)*deltay)/3.0)/n3;

            }
            else {

                // calcul de la normale au plan
                u1 = deltax;
                u2 = deltay;
                u3 = (fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                v1 = 0.0;
                v2 = deltay;
                v3 = (fourierSeries(double(u*deltax),double((v+1)*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                crossProduct(u1,u2,u3,v1,v2,v3,n1,n2,n3);

                // calcul de d du plan
                d = -n1*u*deltax -n2*v*deltay -n3*fourierSeries(double(u*deltax),double(v*deltay),alleles);
                volume += (deltax*deltay/2.)*(-d - n1*(2*u*deltax+(u+1)*deltax)/3.0 - n2*(2*v*deltay+(v+1)*deltay)/3.0)/n3;

                // calcul de la normale au plan
                u1 = deltax;
                u2 = 0.0;
                u3 = (fourierSeries(double((u+1)*deltax),double(v*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                v1 = deltax;
                v2 = deltay;
                v3 = (fourierSeries(double((u+1)*deltax),double((v+1)*deltay),alleles)-fourierSeries(double(u*deltax),double(v*deltay),alleles));
                crossProduct(u1,u2,u3,v1,v2,v3,n1,n2,n3);

                // calcul de d du plan
                d = -n1*u*deltax -n2*v*deltay -n3*fourierSeries(double(u*deltax),double(v*deltay),alleles);
                volume += (deltax*deltay/2.)*(-d - n1*(2*u*deltax+(u+1)*deltax)/3.0 - n2*(2*v*deltay+(v+1)*deltay)/3.0)/n3;

            }

        }
    }

    return volume;

}

double FourierBasis::constraint(unsigned int k, vector<double> alleles) {

    if ( k == 0 ) { // première contrainte: vcal_alleles <= volumeMax*110%

        //cout << "Contrainte1: " << w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " - " << 1.1*volumeMax << endl;
        return volume(alleles) - 1.1*volumeMax;

    }
    else if ( k == 1 ) { // deuxième contrainte: vcal_alleles >= volumeMax*90%

        //cout << "Contrainte2: " << -w*l*std::accumulate(alleles.begin(), alleles.end(), 0.0) << " + " << 0.9*volumeMax << endl;
        return -volume(alleles) + 0.9*volumeMax;

    }
    else return 0.0;

}

// *************************************************************
// Test problem: derived class Problem - Projet Lorenzetti
// *************************************************************

ProjetLorenzetti::ProjetLorenzetti():Problem() {


    minVector.push_back( 0.15 );
    maxVector.push_back( 0.6 );

    minVector.push_back( 10000.0 );
    maxVector.push_back( 18000.0 );

    minVector.push_back( 0.00001 );
    maxVector.push_back( 1.0 );

    minVector.push_back( 0.00001 );
    maxVector.push_back( 1.0 );

    nConstraints = 0;

    cout << "\nProblem: ProjetLorenzetti" << endl;

}

vector<double> ProjetLorenzetti::evaluate(unsigned int id, vector<double> alleles) {

    // start of evaluation
    time_t start, end;
    time(&start); // start timer

    stringstream ss;
    ss << "ACTrial" << id;

    string file = ss.str();
    string fileOut = file + ".out";

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    string leString;
    ostringstream alleles0, alleles1, alleles2, alleles3;

    alleles0 << alleles[0];
    alleles1 << alleles[1];
    alleles2 << alleles[2];
    alleles3 << alleles[3];

    // Assemble inputs to batch file.
    leString = "ac-trials.bat";
    leString += " -w";
    leString += alleles0.str();
    leString += " -e";
    leString += alleles1.str();
    leString += " -x";
    leString += alleles2.str();
    leString += " -i";
    leString += alleles3.str();
    // Redirect batch file output.
    leString += " > ";
    leString += fileOut;


    // issue system command

    int errCode = system(leString.c_str());
    if ( errCode != 0 ) throw (string("Error in launching system command."));


    // Scrape output.
    double sumLogRmsError = 0.0;
    readACTrialResults(fileOut, sumLogRmsError);

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,-sumLogRmsError);

    return fitness;

}

void ProjetLorenzetti::readACTrialResults(string filename, double &sumLogRmsError)
{

    string tampon;
    ifstream input1 (filename.c_str(), ios::binary | ios::in );


    // test d'ouverture
    if (!input1.is_open()) throw (string("Error opening: " + filename));

    // Read file
    do {

        input1 >> tampon;
        if( tampon == "" ) break;
        if( tampon == "Sum" ) {
            input1 >> tampon;
            if( tampon == "RMS" ) {
                input1 >> tampon;
                if( tampon == "log" ) {
                    input1 >> tampon;
                    if( tampon == "error:" ) {
                        input1 >> tampon;
                        sscanf(tampon.c_str(), "%le", &sumLogRmsError);
                    }
                }
            }
        }

      } while( !input1.eof() );

    input1.close();

    // Check result.
    if( sumLogRmsError <= 0.0 ) throw (string("Negative or missing RMS log error."));

    return;

}

// *************************************************************
// Test problem: derived class Problem - Projet EnergyPlus
// *************************************************************

EnergyPlus::EnergyPlus(string filename):Problem() {

    // welcome
    cout << "\nProblem: EnergyPlus, configuration file: " << filename << endl;

    // tampon for saving the lines
    string tampon;
    // read configuration filename
    ifstream input1 (filename.c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input1.is_open()) throw (string("Error opening: " + filename));

    // Read file
    do {

        input1 >> tampon;
        if( tampon == "" ) break;
        else if (tampon == "Template:") {
            input1 >> tampon;
            templateFile = tampon;
        }
        else if (tampon == "Weather:") {
            input1 >> tampon;
            weatherFile = tampon;
        }
        else if( tampon == "Parameters:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nParams = to<unsigned int>(tampon);
//            cout << "Parameters number: " << nParams << endl;

            for (unsigned int i=0; i<nParams; i++) { // read the parameters
                input1 >> tampon;
                labelVector.push_back(tampon); // label for the parameters

                input1 >> tampon;
                minVector.push_back( to<double>(tampon) );

                input1 >> tampon;
                maxVector.push_back( to<double>(tampon) );

                input1 >> tampon;
                stepVector.push_back( to<double>(tampon) );

                input1 >> tampon; // comments about the parameters

            }
        }
        else if ( tampon == "Constraints:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            nConstraints = to<unsigned int>(tampon);
//            cout << "Constraints number: " << nConstraints << endl;

            input1 >> tampon; // comment on the formation of constraints

            for (unsigned int i=0; i<nConstraints; i++) { // read the parameters
                input1 >> tampon;
                constraintVector.push_back(tampon); // label for the parameters

            }

        }
        else if ( tampon == "Equalities:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nEqualities = to<unsigned int>(tampon);
            cout << "Equalities number: " << nEqualities << endl;

            for (unsigned int i=0; i<nEqualities; i++) { // read the parameters
                input1 >> tampon;
                equalityVector.push_back(tampon); // label for the parameters
                input1 >> tampon; // comments on the parameter
            }

        }
        else if ( tampon == "Strings:" ) {

            input1 >> tampon; // get the number of parameters

            // convert it to int

            unsigned int nStrings = to<unsigned int>(tampon);

            cout << "Strings number: " << nStrings << endl;

            for (unsigned int i=0; i<nStrings; i++) { // read the parameters

                //input1 >> tampon;

                do {

                    getline(input1,tampon,'\n');

                }

                while (tampon.empty());

                //cout << "String: " << to<unsigned int>(tampon.substr(1,tampon.find("=")-1)) << "\tValue: " << tampon.substr(tampon.find("=")+1) << endl;

                vector<string> stringi;

                size_t pos = tampon.find("=");

                do {

                    //cout << tampon.substr(pos+1,tampon.find(",",pos+1)-pos-1) << "\t";

                    stringi.push_back(tampon.substr(pos+1,tampon.find(",",pos+1)-pos-1));

                    pos = tampon.find(",",pos+1);

                }

                while (pos!=string::npos);

                //cout << endl;



                stringMap.insert( pair<unsigned int,vector<string> >(to<unsigned int>(tampon.substr(1,tampon.find("=")-1)), stringi) );



//                cout << "Accessing elements..." << endl;

//                for (unsigned int j=0; j < minVector.size(); j++) {

//                    if ( stringMap.find(j) != stringMap.end() ) { cout << j << "\t";

//                        for (unsigned int k=0; k<stringMap.find(j)->second.size(); k++) cout << stringMap.find(j)->second[k] << "\t";

//                        cout << endl;

//                    }

//                }

            }
        }
        else if ( tampon == "Outputs:") {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nOutput = to<unsigned int>(tampon);

            for (unsigned int i=0; i<nOutput; i++) { // read the parameters
                input1 >> tampon;
                outputLabelVector.push_back(tampon); // label for the parameters
                input1 >> tampon;
                outputCoeffVector.push_back( to<double>(tampon) );

            }

        }

      } while( !input1.eof() );

    input1.close();

}

double EnergyPlus::constraint(unsigned int k, vector<double> alleles) {

    // display for debugging reasons
    //cout << "Constraint " << k << ": " << constraintVector[k] << endl;
    //cout << "Alleles: [";
    //for (unsigned int i=0; i<alleles.size()-1; i++) cout << alleles[i] << "\t";
    //cout << alleles[alleles.size()-1] << "]" << endl;
    //cout << "Computed constraint: " << equationParser(constraintVector[k],alleles) << endl;

    // we are looking at constraint k, g_k(x) <= 0
    // interpret constraint k
    return equationParser(constraintVector[k],alleles);

}

vector<double> EnergyPlus::evaluate(unsigned int id, vector<double> alleles) {

    // start of evaluation
    time_t start, end;
    time(&start); // start timer

    stringstream ss;
    ss << "Individual" << id;
    string file = ss.str();

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // extension of the alleles due to the equalities
    vector<double> newAlleles;
    for (unsigned int i=0; i<equalityVector.size(); i++) {

        string secondMember = equalityVector[i].substr(equalityVector[i].find("=")+1); // from posEqual to the end of string
        newAlleles.push_back(equationParser(secondMember, alleles));

    }

    // merge the alleles
    if (newAlleles.size() > 0 ) alleles.insert(alleles.end(), newAlleles.begin(), newAlleles.end());

    // open template
    // tampon for saving the lines
    string tampon;
    // read configuration filename
    ifstream input1 (templateFile.c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input1.is_open()) throw (string("Error opening: " + templateFile));
    // create a ostringstream for the output
    ostringstream oss;

    do {
        getline(input1,tampon,'\n');
        //if (tampon.empty()) break;

        string::size_type first = tampon.find_first_of("%", 0);
        string::size_type second = tampon.find_first_of("%",first+1);
        if ( first != string::npos && second != string::npos ) {
            //cout << "\nTampon: " << tampon << endl;
            unsigned int number = to<unsigned int>(tampon.substr(first+1,second-first-1));
            //cout << "substr: " << tampon.substr(first+1,second-first-1) << "\tnumber: " << number << endl;
            //cout << "toString(alleles): " << toString(alleles[number]) << "x" << endl;
            if ( stringMap.find(number) != stringMap.end() )
                 tampon.replace(first,second-first+1,stringMap.find(number)->second[static_cast<unsigned int>(trunc(alleles[number]-TINY))],0,stringMap.find(number)->second[static_cast<unsigned int>(trunc(alleles[number]-TINY))].length());
            else tampon.replace(first,second-first+1,toString(alleles[number]),0,toString(alleles[number]).length());     			
            //cout << "\nTampon: " << tampon << endl;
            oss << tampon << endl;
        }
        else { oss << tampon << endl; }
    }
    while (!input1.eof());
    input1.close();

    // write output in directory under in.idf
    writeResultsOver(string(file + ".idf"),oss);

    // system the simulation
    int err = system(string("RunEPlus2 " + file + " " + weatherFile + " > " + file + ".log").c_str());
    if ( err != 0 ) throw ("Cannot start EnergyPlus.");

	// cleans the directory
	err = system(string("CleanEPlus2 " + file).c_str());
    if ( err != 0 ) throw ("Cannot clean directory.");
	
    // read the output file
    // read configuration filename
    ifstream input2(string(file + ".eso").c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input2.is_open()) throw (string("Error opening: " + file + ".eso"));

    do {
        // read dictionnary
        getline(input2,tampon,'\n');

    } while ( tampon.substr(0,22) != "End of Data Dictionary" );

    vector<double> fitness;
    fitness.assign(2,0.0); // two objectives

    do {

        getline(input2,tampon,'\n');
        if (tampon == "") break;

        string::size_type first = tampon.find_first_of(",", 0);
        string::size_type second = tampon.find_first_of(",",first+1);

        // line label
        for (unsigned int i=0;i<outputLabelVector.size();i++) {
            if ( tampon.substr(0,first) == outputLabelVector[i] ) { // electric consumption
                fitness[i] -= outputCoeffVector[i]*to<double>(tampon.substr(first+1,second-first-1));
                //cout << "\nline id: " << tampon.substr(0,first) << "\tconsumption: " << tampon.substr(first+1,second-first-1);
            }
        }

        // // ajout de particuliers
        // if ( tampon.substr(0,first) == "1073" ) { // cooling electric consumption
            // fitness[1] += to<double>(tampon.substr(first+1,second-first-1));
            // //cout << "\nline id: " << tampon.substr(0,first) << "\tconsumption: " << tampon.substr(first+1,second-first-1);
        // }
        // if ( tampon.substr(0,first) == "1310" ) { // fan electric consumption
            // fitness[2] += to<double>(tampon.substr(first+1,second-first-1));
            // //cout << "\nline id: " << tampon.substr(0,first) << "\tconsumption: " << tampon.substr(first+1,second-first-1);
        // }

    }
    while (!input2.eof());
    input2.close();

    // // read the ssz file (system sizing)
    // ifstream input3(string(file + ".ssz").c_str(), ios::binary | ios::in );
    // // test d'ouverture
    // if (!input3.is_open()) throw (string("Error opening: " + file + ".ssz"));

    // do {

        // getline(input3,tampon,'\n');
        // if (tampon == "") break;

        // string::size_type first = tampon.find_first_of(",", 0);
        // string::size_type second = tampon.find_first_of(",",first+1);
        // string::size_type third = tampon.find_first_of(",",second+1);

        // // line label
        // if ( tampon.substr(0,first) == "NonCoinc Peak" ) { // system design mass flow rate
            // fitness[3] += to<double>(tampon.substr(second+1,third-second-1));
            // //cout << "\nline id: " << tampon.substr(0,first) << "\tconsumption: " << tampon.substr(first+1,second-first-1);
        // }

    // }
    // while (!input3.eof());
    // input3.close();

    // destroys the files
    //remove(string(file + ".idf").c_str());
    //remove(string(file + ".eso").c_str());

    // writing the fitness in a file
    ostringstream ossFitness;
    ossFitness << setprecision(12);
    for (vector<double>::iterator it=fitness.begin(); it!=fitness.end(); it++) {
        ossFitness << (*it) << endl;
    }
    writeResultsOver(string(file + ".fit"),ossFitness);

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    // returns the final value
    return fitness;

}

// *************************************************************
// Test problem: derived class Problem - Projet CitySim
// *************************************************************

CitySim::CitySim(string filename):Problem() {

    // welcome
    cout << "\nProblem: CitySim, configuration file: " << filename << endl;

    // tampon for saving the lines
    string tampon;
    // read configuration filename
    ifstream input1 (filename.c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input1.is_open()) throw (string("Error opening: " + filename));

    // Read configuration file
    do {

        input1 >> tampon;
        if( tampon == "" ) break;
        else if (tampon == "Template:") {
            input1 >> tampon;
            templateFile = tampon;
        }
        else if( tampon == "Parameters:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nParams = to<unsigned int>(tampon);
//            cout << "Parameters number: " << nParams << endl;

            for (unsigned int i=0; i<nParams; i++) { // read the parameters
                input1 >> tampon;
                labelVector.push_back(tampon); // label for the parameters

                input1 >> tampon;
                minVector.push_back( to<double>(tampon) );

                input1 >> tampon;
                maxVector.push_back( to<double>(tampon) );

                input1 >> tampon;
                stepVector.push_back( to<double>(tampon) );

                input1 >> tampon; // comments about the parameters

            }
        }
        else if ( tampon == "Constraints:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            nConstraints = to<unsigned int>(tampon);
//            cout << "Constraints number: " << nConstraints << endl;

            input1 >> tampon; // comment on the formation of constraints

            for (unsigned int i=0; i<nConstraints; i++) { // read the parameters
                input1 >> tampon;
                constraintVector.push_back(tampon); // label for the parameters
                input1 >> tampon; // comments on the parameter
            }

        }
        else if ( tampon == "Equalities:" ) {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nEqualities = to<unsigned int>(tampon);
            cout << "Equalities number: " << nEqualities << endl;

            for (unsigned int i=0; i<nEqualities; i++) { // read the parameters
                input1 >> tampon;
                equalityVector.push_back(tampon); // label for the parameters
                input1 >> tampon; // comments on the parameter
            }

        }
        else if ( tampon == "Strings:" ) {

            input1 >> tampon; // get the number of parameters

            // convert it to int

            unsigned int nStrings = to<unsigned int>(tampon);

            cout << "Strings number: " << nStrings << endl;

            for (unsigned int i=0; i<nStrings; i++) { // read the parameters

                //input1 >> tampon;

                do {

                    getline(input1,tampon,'\n');

                }

                while (tampon.empty());

                //cout << "String: " << to<unsigned int>(tampon.substr(1,tampon.find("=")-1)) << "\tValue: " << tampon.substr(tampon.find("=")+1) << endl;

                vector<string> stringi;

                size_t pos = tampon.find("=");

                do {

                    //cout << tampon.substr(pos+1,tampon.find(",",pos+1)-pos-1) << "\t";

                    stringi.push_back(tampon.substr(pos+1,tampon.find(",",pos+1)-pos-1));

                    pos = tampon.find(",",pos+1);

                }

                while (pos!=string::npos);

                //cout << endl;



                stringMap.insert( pair<unsigned int,vector<string> >(to<unsigned int>(tampon.substr(1,tampon.find("=")-1)), stringi) );



//                cout << "Accessing elements..." << endl;

//                for (unsigned int j=0; j < minVector.size(); j++) {

//                    if ( stringMap.find(j) != stringMap.end() ) { cout << j << "\t";

//                        for (unsigned int k=0; k<stringMap.find(j)->second.size(); k++) cout << stringMap.find(j)->second[k] << "\t";

//                        cout << endl;

//                    }

//                }

            }



        }
        else if ( tampon == "Outputs:") {
            input1 >> tampon; // get the number of parameters
            // convert it to int
            unsigned int nOutput = to<unsigned int>(tampon);

            for (unsigned int i=0; i<nOutput; i++) { // read the parameters
                input1 >> tampon;
                outputLabelVector.push_back(tampon); // label for the parameters
                input1 >> tampon;
                outputCoeffVector.push_back( to<double>(tampon) );

            }

        }

      } while( !input1.eof() );

    input1.close();

}

double CitySim::constraint(unsigned int k, vector<double> alleles) {

    // display for debugging reasons
    //cout << "Constraint " << k << ": " << constraintVector[k] << endl;
    //cout << "Alleles: [";
    //for (unsigned int i=0; i<alleles.size()-1; i++) cout << alleles[i] << "\t";
    //cout << alleles[alleles.size()-1] << "]" << endl;
    //cout << "Computed constraint: " << equationParser(constraintVector[k],alleles) << endl;

    // we are looking at constraint k, g_k(x) <= 0
    // interpret constraint k
    return equationParser(constraintVector[k],alleles);

}

vector<double> CitySim::evaluate(unsigned int id, vector<double> alleles) {

    // start of evaluation
    time_t start, end;
    time(&start); // start timer

    stringstream ss;
    ss << "Individual" << id;
    string file = ss.str();

    // writing the alleles in a file
    ostringstream ossAlleles;
    ossAlleles << setprecision(12);
    for (unsigned int i=0; i<alleles.size(); i++) {
        ossAlleles << alleles[i] << "\n";
    }
    ossAlleles.flush();
    writeResultsOver(string(file + ".all"),ossAlleles);

    // extension of the alleles due to the equalities
    vector<double> newAlleles;
    for (unsigned int i=0; i<equalityVector.size(); i++) {

        string secondMember = equalityVector[i].substr(equalityVector[i].find("=")+1); // from posEqual to the end of string
        newAlleles.push_back(equationParser(secondMember, alleles));

    }

    // merge the alleles
    if (newAlleles.size() > 0 ) alleles.insert(alleles.end(), newAlleles.begin(), newAlleles.end());

    // open template
    // tampon for saving the lines
    string tampon;
    // read configuration filename
    ifstream input1 (templateFile.c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input1.is_open()) throw (string("Error opening: " + templateFile));
    // create a ostringstream for the output
    ostringstream oss;

    do {
        getline(input1,tampon,'\n');
        //if (tampon.empty()) break;

        string::size_type first = tampon.find_first_of("%", 0);
        string::size_type second = tampon.find_first_of("%",first+1);
        if ( first != string::npos && second != string::npos ) {
            //cout << "\nTampon: " << tampon << endl;
            unsigned int number = to<unsigned int>(tampon.substr(first+1,second-first-1));
            //cout << "substr: " << tampon.substr(first+1,second-first-1) << "\tnumber: " << number << endl;
            //cout << "toString(alleles): " << toString(alleles[number]) << "x" << endl;
            if ( stringMap.find(number) != stringMap.end() )
                 tampon.replace(first,second-first+1,stringMap.find(number)->second[static_cast<unsigned int>(trunc(alleles[number]-TINY))],0,stringMap.find(number)->second[static_cast<unsigned int>(trunc(alleles[number]-TINY))].length());
            else tampon.replace(first,second-first+1,toString(alleles[number]),0,toString(alleles[number]).length());     
			//cout << "\nTampon: " << tampon << endl;
            oss << tampon << endl;
        }
        else { oss << tampon << endl; }
    }
    while (!input1.eof());
    input1.close();

    // write output in directory under in.idf
    writeResultsOver(string(file + ".xml"),oss);

    // system the simulation
    int err = system(string("CitySim " + file + ".xml > " + file + ".log").c_str());
    if ( err != 0 ) throw ("Cannot start CitySim.");

    // read the output file
    // read configuration filename
    ifstream input2(string(file + "_YearlyResults.out").c_str(), ios::binary | ios::in );
    // test d'ouverture
    if (!input2.is_open()) throw (string("Error opening: " + file + "_YearlyResults.out"));

    vector<double> fitness;
    fitness.assign(outputLabelVector.size(),0.0);  // lecture de 1'objectif

    do {

        getline(input2,tampon,'\n');
        if (tampon == "") break;

        string::size_type first = tampon.find_first_of(" ", 0);
        string::size_type second = tampon.find_first_of(":",first+1);
        string::size_type third = tampon.find_first_of("\n",second+1);

        // line label
        for (unsigned int i=0;i<outputLabelVector.size();i++) {
            if ( tampon.substr(0,first) == outputLabelVector[i] ) {
                fitness[i] -= outputCoeffVector[i]*to<double>(tampon.substr(second+1,third-second-1));
                cout << "\nline id: " << tampon.substr(0,first) << "\tconsumption: " << tampon.substr(second+1,third-second-1);
            }
        }
    }
    while (!input2.eof());
    input2.close();

    // destroys the files
    //remove(string(file + ".idf").c_str());
    //remove(string(file + ".eso").c_str());

    // writing the fitness in a file
    ostringstream ossFitness;
    ossFitness << setprecision(12);
    for (vector<double>::iterator it=fitness.begin(); it!=fitness.end(); it++) {
        ossFitness << (*it) << endl;
    }
    writeResultsOver(string(file + ".fit"),ossFitness);

    // fin de l'evaluation
    time(&end);

    double evalTime = difftime(end, start);
    cout << file << "->\ttime to simulate: " << evalTime << " s" << endl;
    timeEval += evalTime;

    // incrémente le compteur de simulations
    nEval++;

    // returns the final value
    return fitness;

}

// *************************************************************
// Test problem: derived class Problem - Ackley Function
// *************************************************************

Ackley::Ackley():Problem() {}

Ackley::Ackley(unsigned int size):Problem() {

    minZ = -20.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back(  32.768 ); // changed for CMA-ES normally 32.768/-32.768
        minVector.push_back( -32.768 );

    }

}

Ackley::Ackley(unsigned int size, double stepSize):Problem(stepSize) {

    minZ = -20.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back(  32.768 ); // changed for CMA-ES normally 32.768/-32.768
        minVector.push_back( -32.768 );

    }

}

vector<double> Ackley::evaluate(unsigned int id, vector<double> alleles) {

    double squares = 0.0, cosines = 0.0;

    unsigned int n = alleles.size();

    for (unsigned int i=0; i<n; i++) {

        squares += pow( alleles[i], 2.0 );
        cosines += cos( 2.0 * PI * alleles[i] );

    }

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,-(-20.0*exp( -0.2 * sqrt( (1.0/ (double) n) * squares ) ) - exp( (1.0/ (double) n) * cosines ) + 20.0 + exp(1.0)));

    return fitness;

}

// *************************************************************
// Test problem: derived class Problem - Rastrigin Function
// *************************************************************

Rastrigin::Rastrigin():Problem() {}

Rastrigin::Rastrigin(unsigned int size):Problem() {

    minZ = -80.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 5.12 );
        minVector.push_back(-5.12 );

    }

}

Rastrigin::Rastrigin(unsigned int size, double stepSize):Problem(stepSize) {

    minZ = -80.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 5.12 );
        minVector.push_back(-5.12 );

    }

}

vector<double> Rastrigin::evaluate(unsigned int id, vector<double> alleles) {

    double A = 10.0;

    int n = alleles.size();

    double sum = A * (double) n;

    for (int i=0; i<n; i++) {

        sum += pow( alleles[i], 2.0 ) - A * cos( 2.0 * PI * alleles[i] );

    }

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,-sum);

    return fitness;
}

// ***********************************************************************
// Test problem: derived class Problem - Rastrigin Constrained Function
// ***********************************************************************

RastriginConstrained::RastriginConstrained(unsigned int size):Rastrigin() {

    minZ = -80.0;
    maxZ = 0.0;

    nConstraints = 1;

    maxVector.push_back(5.12);
    minVector.push_back(0.0);

    for (unsigned int i=1; i < size ; i++) {

        maxVector.push_back( 5.12);
        minVector.push_back(-5.12);

    }

}

RastriginConstrained::RastriginConstrained(unsigned int size, double stepSize):Rastrigin(stepSize) {

    minZ = -80.0;

    maxZ = 0.0;

    nConstraints = 1;

    maxVector.push_back(5.12);
    minVector.push_back(0.0);

    for (unsigned int i=1; i < size ; i++) {

        maxVector.push_back( 5.12);

        minVector.push_back(-5.12);
		
    }

}

double RastriginConstrained::constraint(unsigned int k, vector<double> alleles) {

    double sum = 0.;

    if (k == 0) { // first constraint

        sum = -4.; // constrained in a sphere of radius 2

        sum += pow(alleles[0]-2+0.02,2);

        for (unsigned int i=1; i<alleles.size(); i++) sum += pow(alleles[i],2);

    }

    return sum;

}

vector<double> RastriginConstrained::evaluate(unsigned int id, vector<double> alleles) {

    double A = 10.0;

    int n = alleles.size();

    double sum = A * (double) n;

    for (int i=0; i<n; i++) sum += pow( alleles[i], 2.0 ) - A * cos( 2.0 * PI * alleles[i] );

    // incrémente le compteur de simulations

    nEval++;

    vector<double> fitness;

    fitness.assign(1,-sum);

    return fitness;

}
// *************************************************************
// Test problem: derived class Problem - Rosenbrock Function
// *************************************************************

Rosenbrock::Rosenbrock():Problem() {}

Rosenbrock::Rosenbrock(unsigned int size):Problem() {

    minZ = -3900.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 2.048 );
        minVector.push_back(-2.048 );

    }

}

Rosenbrock::Rosenbrock(unsigned int size, double stepSize):Problem(stepSize) {

    minZ = -3900.0;
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 2.048 );
        minVector.push_back(-2.048 );

    }

}

vector<double> Rosenbrock::evaluate(unsigned int id, vector<double> alleles) {

    int n = alleles.size();

    double sum = 0.0;

    for (int i=0; i<(n-1); i++) {

        sum += 100.0*pow( ( alleles[i]*alleles[i] - alleles[i+1] ) , 2.0) + pow( (1.0-alleles[i]) , 2.0);

    }

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,-sum);

    return fitness;
}

// *************************************************************
// Test problem: derived class Problem - Quadratic Function
// *************************************************************

Quadratic::Quadratic():Problem() {}

Quadratic::Quadratic(unsigned int size):Problem() {

    minZ = -1.0*double(size);
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 1.0 );
        minVector.push_back(-1.0 );

    }

}

Quadratic::Quadratic(unsigned int size, double stepSize):Problem(stepSize) {

    minZ = -1.0*double(size);
    maxZ = 0.0;

    for (unsigned int i=0; i < size ; i++) {

        maxVector.push_back( 1.0 );
        minVector.push_back(-1.0 );

    }

}

vector<double> Quadratic::evaluate(unsigned int id, vector<double> alleles) {

    unsigned int n = alleles.size();

    double sum = 0.0;

    for (unsigned int i=0; i<n; i++) {

        sum += alleles[i]*alleles[i];

    }

    // incrémente le compteur de simulations
    nEval++;

    vector<double> fitness;
    fitness.assign(1,-sum);

    return fitness;
}

// *************************************************************

// Test problem: derived class Problem - Michalewicz Function

// *************************************************************



Michalewicz::Michalewicz():Problem() {



    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.01);

    stepVector.push_back(0.1);

    stepVector.push_back(0.1);

    stepVector.push_back(0.1);

    stepVector.push_back(0.01);



    nConstraints = 9;



    // box constraints definition



    for (unsigned int i=0; i < 9 ; i++) {

        maxVector.push_back( 1.0 );

        minVector.push_back( 0.0 );

    }



    for (unsigned int i=9; i < 12 ; i++) {

        maxVector.push_back( 100.0 );

        minVector.push_back( 0.0 );

    }

    for (unsigned int i=12; i < 13 ; i++) {

        maxVector.push_back( 1.0 );

        minVector.push_back( 0.0 );

    }





}



Michalewicz::Michalewicz(double stepSize):Problem(stepSize) {



    nConstraints = 9;



    // box constraints definition



    for (unsigned int i=0; i < 9 ; i++) {

        maxVector.push_back( 1.0 );

        minVector.push_back( 0.0 );

    }



    for (unsigned int i=9; i < 12 ; i++) {

        maxVector.push_back( 100.0 );

        minVector.push_back( 0.0 );

    }

    for (unsigned int i=12; i < 13 ; i++) {

        maxVector.push_back( 1.0 );

        minVector.push_back( 0.0 );

    }



}



vector<double> Michalewicz::evaluate(unsigned int id, vector<double> alleles) {



    double sum = 5.*alleles[0] + 5.*alleles[1] + 5.*alleles[2] + 5.*alleles[3];



    for (unsigned int i=0; i<4; i++) sum -= 5.*alleles[i]*alleles[i];



    for (unsigned int i=4; i<13; i++) sum -= alleles[i];



    sum += 15.0; // to get a minimum at zero



    // incrémente le compteur de simulations

    nEval++;

    cout << "Evaluation: " << nEval << endl;



    vector<double> fitness;

    fitness.assign(1,-sum);



    return fitness;

}



double Michalewicz::constraint(unsigned int k, vector<double> alleles) {



    double sum = 0.;



    if (k == 0)      sum = 2.*alleles[0] + 2.*alleles[1] + alleles[9] + alleles[10] - 10.;

    else if (k == 1) sum = 2.*alleles[0] + 2.*alleles[2] + alleles[9] + alleles[11] -10.;

    else if (k == 2) sum = 2.*alleles[1] + 2.*alleles[2] + alleles[10] + alleles[11] -10.;

    else if (k == 3) sum = -8.*alleles[0] + alleles[9];

    else if (k == 4) sum = -8.*alleles[1] + alleles[10];

    else if (k == 5) sum = -8.*alleles[2] + alleles[11];

    else if (k == 6) sum = -2.*alleles[3] - alleles[4] + alleles[9];

    else if (k == 7) sum = -2.*alleles[5] - alleles[6] + alleles[10];

    else if (k == 8) sum = -2.*alleles[7] - alleles[8] + alleles[11];



    return sum;



}