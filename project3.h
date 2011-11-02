/* Finding a smallest circle */
#include <cmath>
#include <float.h>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
//enable this library & namespace below for printing the output in command line
#include <iostream>
//#include "project3.h"
using namespace std;

//point class definition
class Point {
    private:
	    double x, y;
    public:
	    Point(): x(-1.0), y(-1.0) {};
	    Point(double init_x, double init_y) { x = init_x; y = init_y; };
	    double get_x() {return x;};
	    double get_y() {return y;};
            void print() { printf("x: %f y: %f",x,y); };
            bool isEqual(Point p) {
                double dx = x - p.x; double dy = y - p.y;
                if (fabs(dx*dx + dy*dy) <= DBL_EPSILON)
                    return true;
                else
                    return false;
            }
};

//circle class definition
class Circle {
    private:
        double sqRadius;     //the square radius of smallest circle
        Point center;     //the center of smallest circle
    public:    
        Circle(){};
         //compute circle for single point
        Circle(Point p) { center = p; sqRadius = 0.0; };
        //compute circle for double point
        Circle(Point p1, Point p2) {
            double xDiff =p2.get_x() - p1.get_x();
            double yDiff =p2.get_y() - p1.get_y();
            center = Point( (p1.get_x() + p2.get_x())*0.5, (p1.get_y() + p2.get_y())*0.5);
            sqRadius = (xDiff*xDiff + yDiff*yDiff) * 0.25;
        };
        //compute circle for triple point
        Circle(Point p1, Point p2, Point p3) {
            Point o21 = Point(p2.get_x() - p1.get_x(), p2.get_y() - p1.get_y());
            Point o31 = Point(p3.get_x() - p1.get_x(), p3.get_y() - p1.get_y());
            
            double l21 = (o21.get_x()*o21.get_x()+o21.get_y()*o21.get_y()) * 0.5;
            double l31 = (o31.get_x()*o31.get_x()+o31.get_y()*o31.get_y()) * 0.5;
            double detMatrix = o21.get_x()*o31.get_y() - o21.get_y()*o31.get_x();
            if (fabs(detMatrix) > DBL_EPSILON) {
                double inv = 1.0/detMatrix;
                double dx = (o31.get_y()*l21 - o21.get_y()*l31) * inv;
                double dy = (o21.get_x()*l31 - o31.get_x()*l21) * inv;
                center = Point(dx+p1.get_x(), dy+p1.get_y());
                sqRadius = dx*dx+dy*dy;
            }
        };
        double getSqRadius() { return sqRadius; };
        double getRadius() { return sqrt(sqRadius);};
        Point getCenter() { return center;};
        
        void printCircle() {
            cout << "The circle center => ";
            center.print();
            printf("  radius => %f", getRadius());
        }

        //check if p is within circle ... euclidean distance - radius <=0
        bool contains(Point p) {
            double xDiff = p.get_x() - center.get_x();
            double yDiff = p.get_y() - center.get_y();
            double distance = xDiff*xDiff + yDiff*yDiff;
            return (distance - sqRadius <= 0.0);
        };
};

class SmallestCircle {
    private:
        int boundarySize;
        Point bp0, bp1, bp2;//max we have 3 boundary points

        int nbPoints;       //number of point to be enclosed
        Point *points;      //pointer to list of points
        Circle minCircle;

        bool hasBoundaryPoint(Point p) {
            bool contain = false;
            switch(boundarySize) {
                case 1:
                    contain = p.isEqual(bp0);
                    break;
                case 2:
                    contain = (p.isEqual(bp0) && p.isEqual(bp1));
                    break;
                case 3:
                    contain = (p.isEqual(bp0) && p.isEqual(bp1) && p.isEqual(bp2));
                    break;
            }
            return contain;
        };

        void updateCircle(Point p) {
            switch(boundarySize) {
                case 1: {
                        bp1 = p;
                        minCircle = Circle(bp0, p);
                        boundarySize = 2;
                        break;
                    }
                case 2: {
                        double minSqRadius = DBL_MAX;
                        Point p0 = bp0; 
                        Point p1 = bp1;
                        int circleTest = -1;
                        Circle c0 = Circle(p0, p);
                        if (c0.contains(p1)) {
                            circleTest = 0;
                            minSqRadius = c0.getSqRadius();
                        }
                        Circle c1 = Circle(p1, p);
                        if ((c1.getSqRadius() < minSqRadius) && c1.contains(p0))
                            circleTest = 1;
                        if (circleTest>-1) {
                            if (circleTest==0) {
                                minCircle = c0;
                                bp1 = p;
                            }
                            else {
                                minCircle = c1;
                                bp0 = p;
                            }
                            boundarySize = 2;
                        }
                        else {
                            minCircle = Circle(p0, p1, p);
                            bp2 = p;
                            boundarySize = 3;
                        }
                        break;
                    }
                case 3: {
                        Point p0 = bp0;
                        Point p1 = bp1;
                        Point p2 = bp2;
                        int circleTest = -1;
                        double minSqRadius = DBL_MAX;
                        //circle of (p0,p) -> p1, p2 are inside
                        Circle c0 = Circle(p0,p);
                        if (c0.contains(p1) && c0.contains(p2)) {
                            circleTest = 0;
                            minSqRadius = c0.getSqRadius();
                        }
                        //circle of (p1,p) -> p0, p2 are inside
                        Circle c1 = Circle(p1,p);
                        if (c1.getSqRadius() < minSqRadius && c1.contains(p0) && c1.contains(p2)) {
                            circleTest = 1;
                            minSqRadius = c1.getSqRadius();
                        }
                        //circle of (p2,p) -> p0, p1 are inside
                        Circle c2 = Circle(p2,p);
                        if (c2.getSqRadius() < minSqRadius && c2.contains(p0) && c2.contains(p1)) {
                            circleTest = 2;
                            minSqRadius = c2.getSqRadius();
                        }
                        //circle of (p0,p1,p) -> p2 is inside
                        Circle c3 = Circle(p0,p1,p);
                        if (c3.getSqRadius() < minSqRadius && c3.contains(p2)) {
                            circleTest = 3;
                            minSqRadius = c3.getSqRadius();
                        }
                        //circle of (p0,p2,p) -> p1 is inside
                        Circle c4 = Circle(p0,p2,p);
                        if (c4.getSqRadius() < minSqRadius && c4.contains(p1)) {
                            circleTest = 4;
                            minSqRadius = c4.getSqRadius();
                        }
                        //circle of (p1,p2,p) -> p0 is inside
                        Circle c5 = Circle(p1,p2,p);
                        if (c5.getSqRadius() < minSqRadius && c5.contains(p0)) {
                            circleTest = 5;
                            minSqRadius = c5.getSqRadius();
                        }
                        switch(circleTest) {
                            case 0:
                                //keep p0
                                minCircle = c0;
                                bp1 = p; bp2 = Point();
                                boundarySize = 2;
                                break;
                            case 1:
                                //keep p1
                                minCircle = c1;
                                bp0 = p; bp2 = Point();
                                boundarySize = 2;
                                break;
                            case 2:
                                //keep p2
                                minCircle = c2;
                                bp0 = bp2; bp1=p; bp2 = Point();
                                boundarySize = 2;
                                break;
                            case 3:
                                //keep p0 and p1
                                minCircle = c3;
                                bp2 = p;
                                boundarySize = 3;
                                break;
                            case 4:
                                //keep p0 and p2
                                minCircle = c4;
                                bp1 = p;
                                boundarySize = 3;
                                break;
                            case 5:
                                //keep p1 and p2
                                minCircle = c5;
                                bp0 = p;
                                boundarySize = 3;
                                break;
                        }
                        break;
                    }
            }
        };

    public:
        SmallestCircle(){};
        SmallestCircle(int n, Point &p) { nbPoints = n; points = &p; 
            boundarySize = 0;
        };
        double getRadius() {
            return minCircle.getRadius();
        };
        Point getCenter() {
            return minCircle.getCenter();
        }
        //use iteration, we know the number of points to be enclosed
        void getMinCircle() {
            /* shuffling process has already done out side circumcircle */
            //shuffle the list of points
            /*
            srand(time(0));
            int *shuffledIndex = new int[nbPoints];
            for (int i=0; i<nbPoints; i++) { shuffledIndex[i] = i; }
            int tmpshuff;
            for (int i=0; i<nbPoints; i++) {
                int r = rand() % nbPoints;
                int tmpshuff = shuffledIndex[i]; shuffledIndex[i] = shuffledIndex[r]; shuffledIndex[r] = tmpshuff;
            }*/

            //initialize circle with first point
            Point p = points[0];
            boundarySize = 1; bp0 = p;
            minCircle = Circle(p);
            int iterIndex = 1;
            while (iterIndex < nbPoints) { 
                p = points[iterIndex];
                if (! hasBoundaryPoint(p)) {
                    if (! minCircle.contains(p)) {
                        updateCircle(p);               
                    }
                }
                iterIndex++;
            } 
        };
};

//the points will all lie in the unit circle (x*x + y*y <= 1)
void circumcircle(int n, double *pts, double *radius, double *center) {
    //copy input into Point class structure
    Point *point = new Point[n];
    for (int i = 0;i < n; i++) {
        point[i] = Point (pts[i*2], pts[i*2 + 1]);
    }
    //clock_t tStart = clock();

    SmallestCircle c = SmallestCircle(n, *point);
    c.getMinCircle();
    *radius = c.getRadius();
    Point p = c.getCenter();
    center[0] = p.get_x();
    center[1] = p.get_y();
    //printf("\nTime taken: %fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
};
