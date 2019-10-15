/*Geometry.cxx
   Jerome Kaempf
   LESO-PB
   2006
*/

#include "geometry.h"

//calculates absolute value
uint32_t Abs(int32_t a){
 if(a < 0)
  return -1*a;
 else
  return a;
}

Point::Point(double x, double y) {

    coordonnees.push_back(x);
    coordonnees.push_back(y);

}

Point::Point(double x, double y, double z) {

    coordonnees.push_back(x);
    coordonnees.push_back(y);
    coordonnees.push_back(z);

}

double Point::getx() {
  return coordonnees[0];
}

double Point::gety() {
    return coordonnees[1];
}

double Point::getz() {
    return coordonnees[2];
}

void Point::setx(double x) {

    coordonnees[0]=x;
}

void Point::sety(double y) {

    coordonnees[1]=y;
}

void Point::setz(double z) {

    coordonnees[2]=z;
}

int Point::getSize() {

    return coordonnees.size();

}

void Point::show() {

    if ( coordonnees.size() == 2 ) {
        cout << "[" << coordonnees[0] << "," << coordonnees[1] << "]";
    }
    if ( coordonnees.size() == 3 ) {
        cout << "[" << coordonnees[0] << "," << coordonnees[1] << "," << coordonnees[2] << "]";
    }

}

double Point::operator[] (unsigned int n) {

    return coordonnees[n];

}

void Point::operator= (Point *lePoint) {

    coordonnees[0]=lePoint->getx();
    coordonnees[1]=lePoint->gety();

}

Point Point::operator- (Point lePoint) {

    coordonnees[0]-=lePoint.getx();
    coordonnees[1]-=lePoint.gety();

    return *this;

}

void grid(vector<Point> rectangle, vector<Point> &grille, double &nx, double &ny, double &nz, double spacing) {

    // obsolète!!!

    double l1, l2, dl1, dl2;

    vector<double> d1;
    vector<double> d2;
    vector<double> n;

    // calcul de l1 et l2

    l1 = sqrt( pow( (rectangle[1][0]-rectangle[0][0]), 2.)
                + pow( (rectangle[1][1]-rectangle[0][1]), 2.)
                + pow( (rectangle[1][2]-rectangle[0][2]), 2.) );

    l2 = sqrt( pow( (rectangle[3][0]-rectangle[0][0]), 2.)
                + pow( (rectangle[3][1]-rectangle[0][1]), 2.)
                + pow( (rectangle[3][2]-rectangle[0][2]), 2.) );

    // calcul de dl1 et dl2

    dl1 = (l1 - spacing*int(l1/spacing+TINY) )/2.;
    dl2 = (l2 - spacing*int(l2/spacing+TINY) )/2.;

    //cout << "dl1: " << dl1 << "\tdl2: " << dl2 << endl;

    // calcul de d1 et d2

    d1.push_back( (rectangle[1][0]-rectangle[0][0]) / l1 );
    d1.push_back( (rectangle[1][1]-rectangle[0][1]) / l1 );
    d1.push_back( (rectangle[1][2]-rectangle[0][2]) / l1 );

    d2.push_back( (rectangle[3][0]-rectangle[0][0]) / l2 );
    d2.push_back( (rectangle[3][1]-rectangle[0][1]) / l2 );
    d2.push_back( (rectangle[3][2]-rectangle[0][2]) / l2 );

    // calcul de n

    nx = (d1[1]*d2[2] - d1[2]*d2[1]);
    ny = (d1[2]*d2[0] - d1[0]*d2[2]);
    nz = (d1[0]*d2[1] - d1[1]*d2[0]);

    // calcul des points de la grille

    for (int i=0; i < int(l1/spacing+TINY); i++) {
        for (int j=0; j < int(l2/spacing+TINY); j++) {

                grille.push_back( Point( rectangle[0][0] + d1[0]*(dl1+spacing*(0.5+i)) + d2[0]*(dl2+spacing*(0.5+j)),
                                         rectangle[0][1] + d1[1]*(dl1+spacing*(0.5+i)) + d2[1]*(dl2+spacing*(0.5+j)),
                                         rectangle[0][2] + d1[2]*(dl1+spacing*(0.5+i)) + d2[2]*(dl2+spacing*(0.5+j)) ) );

        }
    }

}

void gridTriangle(vector<Point> triangle, vector<Point> &grille, double &nx, double &ny, double &nz, double spacing) {

    // décomposition du triangle selon deux directions privilégiées
    // regarde si les points sont dedans

    double l1, l2, l3, H, dl1, dl2;
    int cas1;

    vector<double> d1;
    vector<double> d2;

    // calcul de l1, l2 et l3

    l1 = sqrt( pow( (triangle[1][0]-triangle[0][0]), 2.)
                + pow( (triangle[1][1]-triangle[0][1]), 2.)
                + pow( (triangle[1][2]-triangle[0][2]), 2.) );

    l2 = sqrt( pow( (triangle[2][0]-triangle[1][0]), 2.)
                + pow( (triangle[2][1]-triangle[1][1]), 2.)
                + pow( (triangle[2][2]-triangle[1][2]), 2.) );

    l3 = sqrt( pow( (triangle[0][0]-triangle[2][0]), 2.)
                + pow( (triangle[0][1]-triangle[2][1]), 2.)
                + pow( (triangle[0][2]-triangle[2][2]), 2.) );

    //cout << "l1: " << l1 << "\tl2: " << l2 << "\tl3: " << l3 << endl;

    // distinction des cas

        if ( (l2 > l1) && (l2 > l3) ) {  // cas 1

                // calcul de d1 et d2

                d1.push_back( (triangle[1][0]-triangle[0][0]) / l1 );
                d1.push_back( (triangle[1][1]-triangle[0][1]) / l1 );
                d1.push_back( (triangle[1][2]-triangle[0][2]) / l1 );

                d2.push_back( (triangle[2][0]-triangle[0][0]) / l3 );
                d2.push_back( (triangle[2][1]-triangle[0][1]) / l3 );
                d2.push_back( (triangle[2][2]-triangle[0][2]) / l3 );

                H = l3; // hauteur l3

                cas1 = 1;

                //cout << "cas1" << endl;

                // calcul de n

                nx = (d1[1]*d2[2] - d1[2]*d2[1]);
                ny = (d1[2]*d2[0] - d1[0]*d2[2]);
                nz = (d1[0]*d2[1] - d1[1]*d2[0]);

        }
        else if ( (l3 > l1) && (l3 > l2) ) {  // cas 2

                // calcul de d1 et d2

                d1.push_back( (triangle[0][0]-triangle[1][0]) / l1 );
                d1.push_back( (triangle[0][1]-triangle[1][1]) / l1 );
                d1.push_back( (triangle[0][2]-triangle[1][2]) / l1 );

                d2.push_back( (triangle[2][0]-triangle[1][0]) / l2 );
                d2.push_back( (triangle[2][1]-triangle[1][1]) / l2 );
                d2.push_back( (triangle[2][2]-triangle[1][2]) / l2 );

                H = l2; // hauteur l2

                cas1 = 0;

                //cout << "cas2" << endl;

                // calcul de n -> negatif du au sens de d1

                nx = -(d1[1]*d2[2] - d1[2]*d2[1]);
                ny = -(d1[2]*d2[0] - d1[0]*d2[2]);
                nz = -(d1[0]*d2[1] - d1[1]*d2[0]);

        }
        else return;

    // calcul de dl1 et dl2

    dl1 = (l1 - spacing*int(l1/spacing+TINY) )/2.;
    dl2 = (H - spacing*int(H/spacing+TINY) )/2.;

    // calcul des points de la grille

    for (int i=0; i < int(l1/spacing+TINY); i++) {

        l2 = H*(1.0 - (dl1+spacing*(0.5+i))/l1);

        for (int j=0; j < int( l2/spacing + TINY); j++) {

                grille.push_back( Point( cas1*triangle[0][0] + (1-cas1)*triangle[1][0] + d1[0]*(dl1+spacing*(0.5+i)) + d2[0]*(dl2+spacing*(0.5+j)),
                                         cas1*triangle[0][1] + (1-cas1)*triangle[1][1] + d1[1]*(dl1+spacing*(0.5+i)) + d2[1]*(dl2+spacing*(0.5+j)),
                                         cas1*triangle[0][2] + (1-cas1)*triangle[1][2] + d1[2]*(dl1+spacing*(0.5+i)) + d2[2]*(dl2+spacing*(0.5+j)) ) );

        }
    }

}

void gridTriangle2(vector<Point> triangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz) {

    // décomposition du triangle en 4 sous-triangles et ainsi de suite

    double l1, l2, l3, aire;

    vector<double> d1;
    vector<double> d2;

    vector<Point> subdiv, subdiv2, triangle2;

    // calcul de l1, l2 et l3

    l1 = sqrt( pow( (triangle[1][0]-triangle[0][0]), 2.)
                + pow( (triangle[1][1]-triangle[0][1]), 2.)
                + pow( (triangle[1][2]-triangle[0][2]), 2.) );

    l2 = sqrt( pow( (triangle[2][0]-triangle[1][0]), 2.)
                + pow( (triangle[2][1]-triangle[1][1]), 2.)
                + pow( (triangle[2][2]-triangle[1][2]), 2.) );

    l3 = sqrt( pow( (triangle[0][0]-triangle[2][0]), 2.)
                + pow( (triangle[0][1]-triangle[2][1]), 2.)
                + pow( (triangle[0][2]-triangle[2][2]), 2.) );

//    cout << "l1: " << l1 << "\tl2: " << l2 << "\tl3: " << l3 << endl;

    // calcul de l'aire

    aire = 0.25*sqrt( (l1 + l2 + l3)*(-l1 + l2 + l3)*(l1 - l2 + l3)*(l1 + l2 - l3) );

 //   cout << "aire: " << aire << endl;

    if (aire < TINY) return;

     // calcul de d1 et d2

    d1.push_back( (triangle[1][0]-triangle[0][0]) / l1 );
    d1.push_back( (triangle[1][1]-triangle[0][1]) / l1 );
    d1.push_back( (triangle[1][2]-triangle[0][2]) / l1 );

    d2.push_back( (triangle[2][0]-triangle[0][0]) / l3 );
    d2.push_back( (triangle[2][1]-triangle[0][1]) / l3 );
    d2.push_back( (triangle[2][2]-triangle[0][2]) / l3 );

    // calcul de n

    nx = (d1[1]*d2[2] - d1[2]*d2[1]);
    ny = (d1[2]*d2[0] - d1[0]*d2[2]);
    nz = (d1[0]*d2[1] - d1[1]*d2[0]);

    // définition de l'aire des Points primaire (modifiée lors des subdivisions)

    airePoint = aire;

    // distinction des cas

    if ( aire <= aireDetecteur ) {

            grille.clear();
            grille.push_back( centroid3D(triangle) );

    }
    else {

        unsigned int niveau;

        // préparation de subdiv: les 3 points du triangle

        subdiv.push_back(triangle[0]);
        subdiv.push_back(triangle[1]);
        subdiv.push_back(triangle[2]);

        // calcul du niveau de division à obtenir

        niveau = (unsigned int) int( ceil( log( aire / aireDetecteur )/log(4.0) ) );

//        cout << "niveau: " << niveau << endl;

        for (unsigned int i=0; i < niveau; i++) {

            grille.clear();
            subdiv2.clear();

            airePoint /= 4.; // subdivision de la surface par point

            for (unsigned int j=0; j < (subdiv.size()/3); j++) {

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3],   subdiv[j*3+1]) );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                triangle2.push_back( middlePoint(subdiv[j*3+2], subdiv[j*3]) );
                grille.push_back( centroid3D(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( subdiv[j*3] );
                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+1]) );
                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+2]) );
                grille.push_back( centroid3D(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+1]) );
                triangle2.push_back( subdiv[j*3+1] );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                grille.push_back( centroid3D(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

                triangle2.clear();

                triangle2.push_back( middlePoint(subdiv[j*3], subdiv[j*3+2]) );
                triangle2.push_back( middlePoint(subdiv[j*3+1], subdiv[j*3+2]) );
                triangle2.push_back( subdiv[j*3+2] );
                grille.push_back( centroid3D(triangle2) );

                subdiv2.insert( subdiv2.end(), triangle2.begin(), triangle2.end() );

            }

            subdiv = subdiv2;

        }

    }

//    cout << "airePoint: " << airePoint << endl;

    return;

}

void gridRectangle(vector<Point> rectangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz) {

    // décomposition du rectangle en 9 sous-rectangles (réutilisation du point central)

    double l1, l2, aire;

    vector<double> d1;
    vector<double> d2;

    vector<Point> subdiv, subdiv2, rectangle2;

    // calcul de l1, l2 et l3

    l1 =    sqrt( pow( (rectangle[1][0]-rectangle[0][0]), 2.)
                + pow( (rectangle[1][1]-rectangle[0][1]), 2.)
                + pow( (rectangle[1][2]-rectangle[0][2]), 2.) );

    l2 =    sqrt( pow( (rectangle[2][0]-rectangle[1][0]), 2.)
                + pow( (rectangle[2][1]-rectangle[1][1]), 2.)
                + pow( (rectangle[2][2]-rectangle[1][2]), 2.) );

//    cout << "l1: " << l1 << "\tl2: " << l2 << endl;

    // calcul de l'aire

    aire = l1*l2;

//    cout << "aire: " << aire << endl;

    if (aire < TINY) return; // aire valant zéro: retour sans points...

     // calcul de d1 et d2

    d1.push_back( (rectangle[1][0]-rectangle[0][0]) / l1 );
    d1.push_back( (rectangle[1][1]-rectangle[0][1]) / l1 );
    d1.push_back( (rectangle[1][2]-rectangle[0][2]) / l1 );

    d2.push_back( (rectangle[2][0]-rectangle[1][0]) / l2 );
    d2.push_back( (rectangle[2][1]-rectangle[1][1]) / l2 );
    d2.push_back( (rectangle[2][2]-rectangle[1][2]) / l2 );

    // calcul de n

    nx = (d1[1]*d2[2] - d1[2]*d2[1]);
    ny = (d1[2]*d2[0] - d1[0]*d2[2]);
    nz = (d1[0]*d2[1] - d1[1]*d2[0]);

//    cout << "normale: " << nx << " " << ny << " " << nz << endl;

    // définition de l'aire des Points primaire (modifiée lors des subdivisions)

    airePoint = aire;

    grille.clear();
    grille.push_back( centroid3D(rectangle) );

    // distinction des cas

    if ( aire > aireDetecteur ) {

        int niveau;

        // préparation de subdiv: les 3 points du triangle

        subdiv.push_back(rectangle[0]);
        subdiv.push_back(rectangle[1]);
        subdiv.push_back(rectangle[2]);
        subdiv.push_back(rectangle[3]);

        // calcul du niveau de division à obtenir

        niveau = int( ceil( log( aire / aireDetecteur )/log(9.0) ) );

//        cout << "niveau: " << niveau << endl;

        for (int i=0; i < niveau; i++) {

            subdiv2.clear();

            airePoint /= 9.; // subdivision de la surface par point

            for (unsigned int j=0; j < (subdiv.size()/4); j++) {

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4],   subdiv[j*4+2]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+3]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+1]) );
                // point deja existant dans grille

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( subdiv[j*4] );
                rectangle2.push_back( baryCenter(subdiv[j*4],   subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+2]) );
                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+3]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+3]) );
                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+2]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4]) );
                rectangle2.push_back( subdiv[j*4+1] );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+2]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+3]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+3]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+1], subdiv[j*4+2]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4+1]) );
                rectangle2.push_back( subdiv[j*4+2] );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4+3]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+2], subdiv[j*4+3]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+2]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+2]) );
                rectangle2.push_back( subdiv[j*4+3] );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+3]) );
                rectangle2.push_back( baryCenter(subdiv[j*4], subdiv[j*4+2]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4+1]) );
                rectangle2.push_back( baryCenter(subdiv[j*4+3], subdiv[j*4]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

            }

            subdiv = subdiv2;

        }

    }

//    cout << "airePoint: " << airePoint << endl;

    return;

}

void gridRectangle2(vector<Point> rectangle, double aireDetecteur, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz) {

    // décompose le rectangle en 4 parties (perte du point central)

    double l1, l2, aire;

    vector<double> d1;
    vector<double> d2;

    vector<Point> subdiv, subdiv2, rectangle2;

    // calcul de l1, l2 et l3

    l1 =    sqrt( pow( (rectangle[1][0]-rectangle[0][0]), 2.)
                + pow( (rectangle[1][1]-rectangle[0][1]), 2.)
                + pow( (rectangle[1][2]-rectangle[0][2]), 2.) );

    l2 =    sqrt( pow( (rectangle[2][0]-rectangle[1][0]), 2.)
                + pow( (rectangle[2][1]-rectangle[1][1]), 2.)
                + pow( (rectangle[2][2]-rectangle[1][2]), 2.) );

//    cout << "l1: " << l1 << "\tl2: " << l2 << endl;

    // calcul de l'aire

    aire = l1*l2;

//    cout << "aire: " << aire << endl;

    if (aire < TINY) return; // aire valant zéro: retour sans points...

     // calcul de d1 et d2

    d1.push_back( (rectangle[1][0]-rectangle[0][0]) / l1 );
    d1.push_back( (rectangle[1][1]-rectangle[0][1]) / l1 );
    d1.push_back( (rectangle[1][2]-rectangle[0][2]) / l1 );

    d2.push_back( (rectangle[2][0]-rectangle[1][0]) / l2 );
    d2.push_back( (rectangle[2][1]-rectangle[1][1]) / l2 );
    d2.push_back( (rectangle[2][2]-rectangle[1][2]) / l2 );

    // calcul de n

    nx = (d1[1]*d2[2] - d1[2]*d2[1]);
    ny = (d1[2]*d2[0] - d1[0]*d2[2]);
    nz = (d1[0]*d2[1] - d1[1]*d2[0]);

//    cout << "normale: " << nx << " " << ny << " " << nz << endl;

    // définition de l'aire des Points primaire (modifiée lors des subdivisions)

    airePoint = aire;

    // distinction des cas

    if ( aire <= aireDetecteur ) {

            grille.clear();
            grille.push_back( centroid3D(rectangle) );

    }
    else {

        int niveau;

        // préparation de subdiv: les 4 points du rectangle

        subdiv.push_back(rectangle[0]);
        subdiv.push_back(rectangle[1]);
        subdiv.push_back(rectangle[2]);
        subdiv.push_back(rectangle[3]);

        // calcul du niveau de division à obtenir -> par 4

        niveau = int( ceil( log( aire / aireDetecteur )/log(4.0) ) );

//        cout << "niveau: " << niveau << endl;

        for (int i=0; i < niveau; i++) {

            grille.clear();
            subdiv2.clear();

            airePoint /= 4.; // subdivision de la surface par point

            for (unsigned int j=0; j < (subdiv.size()/4); j++) {

                rectangle2.clear();

                rectangle2.push_back( subdiv[j*4] );
                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+1]) );
                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+2]) );
                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+3]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+1]) );
                rectangle2.push_back( subdiv[j*4+1] );
                rectangle2.push_back( middlePoint(subdiv[j*4+1],   subdiv[j*4+2]) );
                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+2]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+2]) );
                rectangle2.push_back( middlePoint(subdiv[j*4+1], subdiv[j*4+2]) );
                rectangle2.push_back( subdiv[j*4+2] );
                rectangle2.push_back( middlePoint(subdiv[j*4+2], subdiv[j*4+3]) );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

                rectangle2.clear();

                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+3]) );
                rectangle2.push_back( middlePoint(subdiv[j*4], subdiv[j*4+2]) );
                rectangle2.push_back( middlePoint(subdiv[j*4+2], subdiv[j*4+3]) );
                rectangle2.push_back( subdiv[j*4+3] );
                grille.push_back( centroid3D(rectangle2) );

                subdiv2.insert( subdiv2.end(), rectangle2.begin(), rectangle2.end() );

            }

            subdiv = subdiv2;

        }

    }

//    cout << "airePoint: " << airePoint << endl;

    return;

}

void gridRectangle3(vector<Point> rectangle, double spacing, vector<Point> &grille, double &airePoint, double &nx, double &ny, double &nz) {

    // deux directions privilégiées de recherche, points dedans ou dehors (éliminés)
    // méthode tirée de PPF, pour compatibilité avec les données de Marylène

    double l1, l2, aire;
    double lx, ly, max_u, max_v;

    vector<double> d1;
    vector<double> d2;

    vector<Point> subdiv, subdiv2, rectangle2;

    // calcul de l1, l2 et l3

    l1 =    sqrt( pow( (rectangle[1][0]-rectangle[0][0]), 2.)
                + pow( (rectangle[1][1]-rectangle[0][1]), 2.)
                + pow( (rectangle[1][2]-rectangle[0][2]), 2.) );

    l2 =    sqrt( pow( (rectangle[2][0]-rectangle[1][0]), 2.)
                + pow( (rectangle[2][1]-rectangle[1][1]), 2.)
                + pow( (rectangle[2][2]-rectangle[1][2]), 2.) );

//    cout << "l1: " << l1 << "\tl2: " << l2 << endl;

    // calcul de l'aire

    aire = l1*l2;

//    cout << "aire: " << aire << endl;

    if (aire < TINY) return; // aire valant zéro: retour sans points...

     // calcul de d1 et d2

    d1.push_back( (rectangle[1][0]-rectangle[0][0]) / l1 );
    d1.push_back( (rectangle[1][1]-rectangle[0][1]) / l1 );
    d1.push_back( (rectangle[1][2]-rectangle[0][2]) / l1 );

    d2.push_back( (rectangle[2][0]-rectangle[1][0]) / l2 );
    d2.push_back( (rectangle[2][1]-rectangle[1][1]) / l2 );
    d2.push_back( (rectangle[2][2]-rectangle[1][2]) / l2 );

    // calcul de n

    nx = (d1[1]*d2[2] - d1[2]*d2[1]);
    ny = (d1[2]*d2[0] - d1[0]*d2[2]);
    nz = (d1[0]*d2[1] - d1[1]*d2[0]);

//    cout << "normale: " << nx << " " << ny << " " << nz << endl;

    max_u=ceil(l1/spacing);
    max_v=ceil(l2/spacing);

    lx = l1 / max_u;
    ly = l2 / max_v;

    // calcul des points de la grille

    for (int i=0; i < max_u; i++) {
        for (int j=0; j < max_v; j++) {

                grille.push_back( Point( rectangle[0][0] + d1[0]*(lx*(0.5+i)) + d2[0]*(ly*(0.5+j)),
                                         rectangle[0][1] + d1[1]*(lx*(0.5+i)) + d2[1]*(ly*(0.5+j)),
                                         rectangle[0][2] + d1[2]*(lx*(0.5+i)) + d2[2]*(ly*(0.5+j)) ) );

        }
    }

    // calcul du volume par point

    airePoint = aire / (max_u*max_v);

    return;

}

void removeVertically(Point base, vector<double> dx, vector<double> dy, vector<Point> &grille, vector<Point> &normales, vector<double> &aires) {

    for (int i=grille.size()-1; i>=0; i--) {

        double px = grille[i][0];
        double py = grille[i][1];

        px -= base[0];
        py -= base[1];

        // ps entre (px,py) et dx

        double psx = (px*dx[0]+py*dx[1])/sqrt(dx[0]*dx[0] + dx[1]*dx[1]);

        // ps entre (px,py) et dy

        double psy = (px*dy[0]+py*dy[1])/sqrt(dy[0]*dy[0] + dy[1]*dy[1]);

        if ( ( psx < sqrt(dx[0]*dx[0] + dx[1]*dx[1])+0.1 && psx > -0.1 ) && ( psy < sqrt(dy[0]*dy[0] + dy[1]*dy[1])+0.1 && psy > -0.1 ) ) {

            // en dehors de la zone donnée

           grille.erase(grille.begin() + i);
           normales.erase(normales.begin() + i);
           aires.erase(aires.begin() + i);

        }

    }

}

void rotation2D(vector<Point> &vect, double phi) {

  double xprime, yprime;
  double phiRad = phi*M_PI/180.;

  int vectsize = vect.size(); // initial vector size

  for (int i=0;i<vectsize;i++) {
    xprime = vect[i][0]*cos(phiRad) - vect[i][1]*sin(phiRad);
    yprime = vect[i][0]*sin(phiRad) + vect[i][1]*cos(phiRad);

    vect.push_back(Point(xprime, yprime));
  }

  vect.erase(vect.begin(), vect.begin() + vectsize);

}

void translation2D(vector<Point> &vect, double deltax, double deltay) {

  double xprime,yprime;
  int vectsize = vect.size(); // initial vector size

  for (int i=0;i<vectsize;i++) {
    xprime=vect[i][0]+deltax;
    yprime=vect[i][1]+deltay;
    vect.push_back(Point(xprime,yprime));
  }

  vect.erase(vect.begin(), vect.begin() + vectsize);

}

Point centroid2D(vector<Point> &forme) {

  double xprime=0.,yprime=0.;

  for (unsigned int i=0;i<forme.size();i++) {

    xprime += forme[i][0];
    yprime += forme[i][1];

  }

  xprime/=forme.size();
  yprime/=forme.size();

  return Point(xprime,yprime);

}

Point centroid3D(vector<Point> &forme) {

  double xprime=0.,yprime=0.,zprime=0.;

  for (unsigned int i=0;i<forme.size();i++) {

    xprime += forme[i][0];
    yprime += forme[i][1];
    zprime += forme[i][2];

  }

  xprime/=forme.size();
  yprime/=forme.size();
  zprime/=forme.size();

  return Point(xprime,yprime,zprime);

}

Point middlePoint(Point &point1, Point &point2) {

    double xprime, yprime, zprime;

    xprime = 0.5*(point1[0]+point2[0]);
    yprime = 0.5*(point1[1]+point2[1]);
    zprime = 0.5*(point1[2]+point2[2]);

    return Point(xprime, yprime,zprime);

}

Point baryCenter(Point &point1, Point &point2) {

    double xprime, yprime, zprime;

    xprime = (2./3.)*point1[0] + (1./3.)*point2[0];
    yprime = (2./3.)*point1[1] + (1./3.)*point2[1];
    zprime = (2./3.)*point1[2] + (1./3.)*point2[2];

    return Point(xprime, yprime,zprime);

}

Point weightedCenter(Point point1, Point point2, double weight) {

    double xprime, yprime, zprime;

    if ( (point1.getSize() == 2) && (point2.getSize() == 2) ) {

        xprime = (1.0 - weight)*point1[0] + weight*point2[0];
        yprime = (1.0 - weight)*point1[1] + weight*point2[1];

        return Point(xprime, yprime);
    }
    else if ( (point1.getSize() == 3) && (point2.getSize() == 3) ) {

        xprime = (1.0 - weight)*point1[0] + weight*point2[0];
        yprime = (1.0 - weight)*point1[1] + weight*point2[1];
        zprime = (1.0 - weight)*point1[2] + weight*point2[2];

        return Point(xprime, yprime,zprime);
    }
    else throw ("Error in the two points sizes\n");

}

double norme(double vx, double vy, double vz) {

  return sqrt(pow(vx, 2.)+pow(vy, 2.)+pow(vz,2.));

}

double norm(vector<double> v) {

    double accumulateur=0.0;

    for (unsigned int i=0; i < v.size(); i++) {

        accumulateur += pow(v[i], 2.0);

    }

    return sqrt(accumulateur);

}

void normalise(double &nx, double &ny, double &nz) {

    double laNorme;

    laNorme = norme(nx, ny, nz);

    nx /= laNorme;
    ny /= laNorme;
    nz /= laNorme;

    return;

}

void crossProduct(double u1, double u2, double u3, double v1, double v2, double v3, double &n1, double &n2, double &n3) {

    n1 = u2*v3-u3*v2;
    n2 = u3*v1-u1*v3;
    n3 = u1*v2-u2*v1;

    normalise(n1,n2,n3);

    return;

}

void readResults(string filename, vector<double> &results)
{
  double mesure;
  string tampon;
  char tampon1[200];

    // chargement des données a afficher
  ifstream input1 (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture

    if (!input1.is_open())
      {
    cerr << "Error opening: " << filename << endl;
      }

    results.clear();

    while (!input1.eof()) {
      input1.getline(tampon1, 200,'\n');
      if (input1.eof()) break;
      tampon=tampon1;
      if (tampon.empty()) break;
      sscanf(tampon1, "%le", &mesure);
      results.push_back(mesure);

    }

    input1.close();

}

void readResults(string filename, vector<string> &results)
{
  string tampon;
  char tampon1[200];

    // chargement des données a afficher
  ifstream input1 (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture

    if (!input1.is_open())
      {
    cerr << "Error opening: " << filename << endl;
      }

    results.clear();

    while (!input1.eof()) {
      input1.getline(tampon1, 200,'\n');
      if (input1.eof()) break;
      tampon=tampon1;
      if (tampon.empty()) break;
      results.push_back(tampon);
    }

    input1.close();

}

void readResultsInLine(string filename, vector<double> &results)
{
    double mesure;
    string tampon;

    // chargement des données a afficher
    ifstream input1 (filename.c_str(), ios::binary | ios::in );

    // test d'ouverture

    if (!input1.is_open()) throw (string("Error opening: " + filename));

    results.clear();

    do {

        input1 >> tampon;

	if (tampon == "") break;

        sscanf(tampon.c_str(), "%le", &mesure);

        results.push_back(mesure);

	tampon = "";

    } while (!input1.eof());

    input1.close();

}

void writeResults(string filename, ostringstream &oss) {

  ofstream outputInFile(filename.c_str(), ios::binary | ios::app );

    // test d'ouverture

    if (!outputInFile.is_open())
      {
            string msg = "Error opening: ";
            msg += filename;

            throw (msg);
      }

    outputInFile << oss.str();

    outputInFile.flush();
    outputInFile.close();

}

void writeResultsOver(string filename, ostringstream &oss) {

    ofstream outputInFile(filename.c_str(), ios::binary | ios::trunc );

    // test d'ouverture

    if (!outputInFile.is_open())
      {
            string msg = "Error opening: ";
            msg += filename;

            throw (msg);
      }

    string unCutString = oss.str();
    // decoupage du string -> écriture facilitée
    unsigned int maxSize = 10000;
    for (unsigned int i=0; i<(unCutString.size()/maxSize+1); i++) {
        outputInFile << unCutString.substr(i*maxSize, maxSize);// oss.str();
    }
    outputInFile.close();

    return;

}

double normallyDistributedSPRNG_Sum() {

  double utot=0.;

  for (int i=0;i<12;i++) {

    utot += randomUniform(0.0, 1.0);

  }

  return (utot-6.);

}

double normallyDistributedSPRNG_BoxMuller() {

    double u1,u2;

    u1 = randomUniform(0.0, 1.0);
    u2 = randomUniform(0.0, 1.0);

    if ( randomUniform(0.0, 1.0) > 0.5 ) {

        return (sqrt(-2.0*log(1.0-u1)) * cos( 2*M_PI*u2 ));

    }
    else {

        return (sqrt(-2.0*log(1.0-u1)) * sin( 2*M_PI*u2 ));

    }

}

// ******************************** Ziggurat Method for Randn **********************************
// changed for all platforms -> uint32_t, JK - 22/02/09
// note: on a 32-bits machine, the unsigned long int and the unsigned int are of the same size
// (maximum of the 32-bits available = 0..2^32-1)
// HOWEVER, on a 64-bits machine, the unsigned int is still in the range 0..2^32-1, but the
// unsigned long int is now between 0..2^64-1
//
// THEREFORE, we need to use the following definition for the unsigned int:
// typedef unsigned long int uint32_t
// as this method was designed using unsigned int in 32 bits
// *********************************************************************************************

static uint32_t jz,jsr=123456789;

#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3*.2328306e-9)
#define IUNI SHR3

static int32_t hz;
static uint32_t iz, kn[128];
static float wn[128],fn[128];

#define RNOR (hz=SHR3, iz=hz&127, (Abs(hz) < kn[iz])? hz*wn[iz] : nfix())

/* nfix() generates variates from the residue when rejection in RNOR occurs. */

float nfix(void)
{
const float r = 3.442620f;     /* The start of the right tail */
static float x, y;
 for(;;)
  {  x=hz*wn[iz];      /* iz==0, handles the base strip */
     if(iz==0)
       { do{ x=-log(UNI)*0.2904764; y=-log(UNI);}	/* .2904764 is 1/r */
        while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
       }
                         /* iz>0, handle the wedges of other strips */
      if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;

     /* initiate, try to exit for(;;) for loop*/
      hz=SHR3;
      iz=hz&127;
      if( Abs(hz) < kn[iz] ) return (hz*wn[iz]);
  }

}

/*--------This procedure sets the seed and creates the tables------*/

void zigset(uint32_t jsrseed)
{  const double m1 = 2147483648.0, m2 = 4294967296.;
   double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
   double de=7.697117470131487, te=de, ve=3.949659822581572e-3;

   jsr^=jsrseed;

/* Set up tables for RNOR */
   q=vn/exp(-.5*dn*dn);
   kn[0]=(uint32_t) ((dn/q)*m1);
   kn[1]=0;

   wn[0]=q/m1;
   wn[127]=dn/m1;

   fn[0]=1.;
   fn[127]=exp(-.5*dn*dn);

    for(unsigned short int i=126;i>=1;i--) // put the index i in the loop declaration, JK - 22/02/09
    {dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
     kn[i+1]=(uint32_t) ((dn/tn)*m1);
     tn=dn;
     fn[i]=exp(-.5*dn*dn);
     wn[i]=dn/m1;
    }
}

// ******************************** End of Ziggurat Method for Randn *************************** //

double normallyDistributedSPRNG_Ziggurat() {

    return RNOR;

}

double randomUniform(double minValue,double maxValue)
{
	return minValue + UNI * (maxValue - minValue);
}

// ******** Maths extensions ***********

void computeEigensystem(vector<double> &outValues, vector< vector<double> > &outVectors)
{

	// Computer eigenvectors/eigenvalues using Triagonal QL method
	vector<double> lE;
	lE.assign(outVectors.size(), 0.0); // imaginary part
	tred2(outValues, lE, outVectors);
	tql2(outValues, lE, outVectors);

	// Sort by eigenvalues.
	for(unsigned int j = 0; j < outValues.size(); ++j) {
		double lMax=outValues[j];
		unsigned int lMaxArg=j;
		for(unsigned int l = j+1; l<outValues.size(); ++l) {
			if(outValues[l] > lMax) {
				lMax=outValues[l];
				lMaxArg=l;
			}
		}
		if(lMaxArg != j) {
			for(unsigned int r = 0; r < outVectors[j].size(); ++r) {
				double lTmp = outVectors[r][j];
				outVectors[r][j] = outVectors[r][lMaxArg];
				outVectors[r][lMaxArg] = lTmp;
			}
			double lTmp = outValues[j];
			outValues[j] = outValues[lMaxArg];
			outValues[lMaxArg] = lTmp;
		}
	}
}

/*!
 *  \param d Real part of eigenvalues computed from the matrix.
 *  \param e Imaginary part of eigenvalues computed from the matrix.
 *  \param V Eigenvectors computed from the matrix.
 *
 *  This method is derived from procedure tql2 of the Java package JAMA,
 *  which is itself derived from the Algol procedures tql2, by
 *  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
 *  Fortran subroutine in EISPACK.
 */
void tql2(vector<double> &d, vector<double> &e, vector< vector<double> > &V)
{
	const unsigned int n=V.size();
	for(unsigned int i = 1; i < n; i++) e[i-1] = e[i];
	e[n-1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = std::pow(2.0,-52.0);
	for(unsigned int l = 0; l < n; l++) {
		// Find small subdiagonal element
		tst1 = max(tst1, abs(d[l]) + abs(e[l]));
		unsigned int m=l;
		while((m+1) < n) {
			if(std::abs(e[m]) <= eps*tst1) break;
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.
		if(m > l) {
			unsigned int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)
													// Compute implicit shift
				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot(p,1.0);
				if(p < 0) r = -r;
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for(unsigned int i = l+2; i < n; i++) d[i] -= h;
				f = f + h;

				// Implicit QL transformation.
				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for(unsigned int i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.
					for(unsigned int k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
					if(i == 0) break;
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.
			} while (std::abs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}
}

/*!
 *  \param d Real part of eigenvalues computed from the matrix.
 *  \param e Imaginary part of eigenvalues computed from the matrix.
 *  \param V Eigenvectors computed from the matrix.
 *
 *  This method is derived from procedure tred2 of the Java package JAMA,
 *  which is itself derived from the Algol procedures tred2, by
 *  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
 *  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
 *  Fortran subroutine in EISPACK.
 */
void tred2(vector<double> &d, vector<double> &e, vector< vector<double> > &V)
{
	const unsigned int n=V.size();
	for(unsigned int j = 0; j < n; ++j) d[j] = V[n-1][j];

	// Householder reduction to tridiagonal form.
	for(unsigned int i = n-1; i > 0; --i) {

		// Scale to avoid under/overflow.
		double scale = 0.0;
		double h = 0.0;
		for(unsigned int k = 0; k < i; ++k) scale += abs(d[k]);
		if(scale == 0.0) {
			e[i] = d[i-1];
			for(unsigned int j = 0; j < i; ++j) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} else {
			// Generate Householder vector.
			for(unsigned int k=0; k<i; ++k) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if(f > 0.0) g = -g;
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for(unsigned int j = 0; j < i; j++) e[j] = 0.0;

			// Apply similarity transformation to remaining columns.
			for(unsigned int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for(unsigned int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for(unsigned int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for(unsigned int j=0; j<i; j++) e[j] -= hh * d[j];
			for(unsigned int j=0; j<i; j++) {
				f = d[j];
				g = e[j];
				for(unsigned int k = j; k <= i-1; k++) V[k][j] -= (f * e[k] + g * d[k]);
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.
	for(unsigned int i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if(h!=0.0) {
			for(unsigned int k=0; k<=i; k++) d[k] = V[k][i+1] / h;
			for(unsigned int j=0; j<=i; j++) {
				double g = 0.0;
				for(unsigned int k=0; k<=i; k++) g += V[k][i+1] * V[k][j];
				for(unsigned int k=0; k<=i; k++) V[k][j] -= g * d[k];
			}
		}
		for(unsigned int k=0; k<=i; k++) V[k][i+1] = 0.0;
	}
	for(unsigned int j=0; j<n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
}

// Equation parser based upon one found  *  at www.wizardscripts.com/forum

double equationParser(string& eq, vector<double>& alleles)
{
   double result=0.0;
   string::size_type curr_index = 0;

   //loop while there are tokens
   do
   {
         string::size_type start = curr_index;
         string::size_type end = eq.find_first_of("+-*/^",start+1);

         if ( eq[curr_index] == '+' ) { // addition
             // of what?
             if ( eq[curr_index+1] == 'x' ) { // allele
                 result += alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
//                 cout << "Add: " << alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
             }
             else { // number
                 result += to<double>(eq.substr(start+1,end-start-1));
//                 cout << "Add: " << to<double>(eq.substr(start+1,end-start-1));
             }
         }
         else if ( eq[curr_index] == '-' ) { // substraction
             // of what?
             if ( eq[curr_index+1] == 'x' ) { // allele
                 result -= alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
//                 cout << "Substract: " << alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
             }
             else { // number
                 result -= to<double>(eq.substr(start+1,end-start-1));
//                 cout << "Substract: " << to<double>(eq.substr(start+1,end-start-1));
             }
         }
         else if ( eq[curr_index] == '*' ) { // addition
             // of what?
             if ( eq[curr_index+1] == 'x' ) { // allele
                 result *= alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
             }
             else { // number
                 result *= to<double>(eq.substr(start+1,end-start-1));
             }
         }
         else if ( eq[curr_index] == '/' ) { // addition
             // of what?
             if ( eq[curr_index+1] == 'x' ) { // allele
                 result /= alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ];
             }
             else { // number
                 result /= to<double>(eq.substr(start+1,end-start-1));
             }
         }
         else if ( eq[curr_index] == '^' ) { // addition
             // of what?
             if ( eq[curr_index+1] == 'x' ) { // allele
                 result = pow(result, alleles[ to<unsigned int>(eq.substr(start+2,end-start-2)) ]);
             }
             else { // number
                 result = pow(result, to<double>(eq.substr(start+1,end-start-1)));
             }
         }
         else throw (string("Constraint parser: unrecognized operation (not in the set)"));
//         cout << endl;

         curr_index = end;

   }
   while (curr_index != string::npos);

   return result;

}
