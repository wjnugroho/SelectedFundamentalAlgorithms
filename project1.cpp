/* compiles with command line  g++ project1.cpp -lX11 -lm -L/usr/X11R6/lib 
   Find the shortest path between two input points with multiple line obstacles
*/
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
//enable this library & namespace below for printing the output in command line
#include <iostream>
using namespace std;

#define INFINITY_VALUE 1000000000
#define NO_INTERSECTION 0
#define HAS_INTERSECTION 1
#define COLLINEAR 2

/*classes definition*/
//point class definition
class Point {
	private:
		int x, y;
	public:
		Point(): x(-1), y(-1) {};
		Point(int init_x, int init_y);
		int get_x() {return x;};
		int get_y() {return y;};
		bool isEqual(int a, int b);
};

// type definition of line object
class Line {
	private:
		Point* pstart;
		Point* pend;
	public:
		Line() {};
		Line(Point &p1, Point &p2);
		//method
		Point &getStart() { return *pstart; };
		Point &getEnd() { return *pend; };
		void draw(Display* ptr, Window win, GC gc_line);
};

//a pseudo-multidimensional array for storing the distance table
class Matrix {
	private:
		int row, col;
		float* value;
	public:
		Matrix() {};
		Matrix(int r, int c) {
			row=r; col=c;
			value = new float[row*col];
		};
		float get(int r, int c) {
			return value[r * col + c];
		}; 
		void set(int r, int c, float v) {
			value[r * col + c] = v;		
		};
};

class DistanceTableBuilder {
	private:
		Point* p_ptr; // points pointer
		Line* l_ptr; // lines pointer
		int source, target;
		int nbPoints;
		int getOrientation (Point a, Line l);
		double calcEuclideanDistance(Point a, Point b);
		bool hasObstacle(Point a, Point b);
		//storage for distanceTable: pseudo-multidimensional array
		Matrix distanceTable;
		int determineLineIntersection(Point a, Point b, Line l);
		bool equalSign(int a, int b);
		
	public:
		DistanceTableBuilder() {};
		DistanceTableBuilder(int s, int t, int n, Point *&p, Line *&l);
		void createTable();
		void print();
		Matrix getTable() { return distanceTable; };
};

//find the shortest path from start to target point  using Dijkstra's algorithm
class ShortPathFinder {
	private:
		int source, target; //index of source and target
		Point* p_ptr;
		Matrix distanceTable;
		float *distance; //path length from the source to node[i] of pathLength[i]
		int *predecessor; // shortest path from the source to node[i] of pred[i]
		int nbPoints;
		void printPath(int dest, int prev, bool redraw);
		Display* display_ptr;
		Window oWin;
		GC gc_green;

	public:
		ShortPathFinder() {};
		ShortPathFinder(int s, int t, int n, Point *&p, Matrix table);
		void find(); //dijkstra algorithm implementation
		void print();
		void draw(Display *&disp, Window win, Colormap color, bool redraw); //draw the path
};

class PWindow {
	private:
		Display *display_ptr;
		int border_width;
		int win_x, win_y;
		unsigned int display_width, display_height;
		XWMHints *wm_hints;
		XClassHint *class_hints;
		XSizeHints *size_hints;
		XTextProperty win_name, icon_name;
		
	public:
		Screen *screen_ptr;
		int screen_num;
		Window win;
		unsigned int win_width, win_height;
		Colormap color_map;
		
		//constructor
		PWindow(Display *ptr);
		
		//method
		void createWindow(int argc, char **argv);
};		

/* Below are the implementation of all the classes members*/
Point::Point(int init_x, int init_y) {
	x = init_x;
	y = init_y;
};

bool Point::isEqual(int a, int b){
	if ((x==a) && (y==b))
		return true;
	else
		return false;
};

Line::Line(Point &p1, Point &p2) {
	pstart = &p1;
	pend = &p2;
};

void Line::draw(Display * ptr, Window win, GC gcLine) {
	XDrawLine(ptr, win, gcLine, (*pstart).get_x(), (*pstart).get_y(), (*pend).get_x(), (*pend).get_y() );
};

DistanceTableBuilder::DistanceTableBuilder(int s, int t, int n, Point *&p, Line *&l) {
	//set the fields
	source = s;
	target = t;
	nbPoints = n;
	p_ptr = *&p;
	l_ptr = *&l;
	
	cout << "List of Points:\n";
	for (int i=0;i<nbPoints;i++) {
		cout << i << ": x=" << p_ptr[i].get_x() << " y=" << p_ptr[i].get_y() << "\n";
	}
	cout << "Source Point=" << s;
	cout << ", Target Point=" << t << "\n";
	distanceTable = Matrix(nbPoints, nbPoints);
};

// return euclidean distance between two points
double DistanceTableBuilder::calcEuclideanDistance(Point a, Point b) {
	float diffx = (float) (b.get_x() - a.get_x());
	float diffy = (float) (b.get_y() - a.get_y());
	float distance = sqrt( diffx * diffx + diffy * diffy );	
	return distance;
};

// calculate the triangle orientation of l.start->l.end->a
int DistanceTableBuilder::getOrientation (Point a, Line l) {
	float determ = 0; 
	Point s = l.getStart();
	Point e = l.getEnd();
	
	float x1,x2,x3, y1,y2,y3;
	x1 = s.get_x();
	y1 = s.get_y();
	x2 = e.get_x();
	y2 = e.get_y();
	x3 = a.get_x();
	y3 = a.get_y();
	
	determ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	
	if (determ > 0.0)
		return 1;	//counter clock wise
	if (determ < 0.0)
		return -1;	//clock wise
	return 0; //collinear
};

bool DistanceTableBuilder::equalSign(int a, int b) {
	if ( (a^b) >= 0 )
		return true;
	else
		return false;
};

// determine if there is any intersection between segment ab and line l
int DistanceTableBuilder::determineLineIntersection(Point a, Point b, Line l) {
	int x1,y1,	//point a 
            x2,y2,	//point b
	    x3,y3,	//line l.start
	    x4,y4;	//line l.end
	x1 = a.get_x(); y1 = a.get_y();
	x2 = b.get_x(); y2 = b.get_y();
	x3 = (l.getStart()).get_x(); y3 = (l.getStart()).get_y();
	x4 = (l.getEnd()).get_x(); y4 = (l.getEnd()).get_y();
	
	//calculate coeff line for segment ab =>  a1x + b1y + c1 = 0
	int a1 = y2 - y1;
	int b1 = x1 - x2;
        int c1 = x2 * y1 - x1 * y2;
	//calculate the orientation signs
	int o1 = a1 * x3 + b1 * y3 + c1;
	int o2 = a1 * x4 + b1 * y4 + c1;
	//compare o1 and o2 => if they have the same sign then there is no intersection
	if (o1!=0 && o2!=0 && equalSign(o1, o2))
		return NO_INTERSECTION;

	//calculate coeff line for segment l => a2x + b2y + c2 = 0
	int a2 = y4 - y3;
	int b2 = x3 - x4;
	int c2 = x4 * y3 - x3 * y4;
	//calculate the orientation sign
	int o3 = a2 * x1 + b2 * y1 + c2;
	int o4 = a2 * x2 + b2 * y2 + c2;
	//compare o3 and o4 => if they have the same sign then there is no intersection
	if (o3!=0 && o4!=0 && equalSign(o3,o4))
		return NO_INTERSECTION;

	if ( (a1 * b2 - a2 * b1) == 0)
		return COLLINEAR;

	//must have an intersection.
	return HAS_INTERSECTION;
};

// do we have any obstacle between a and b
bool DistanceTableBuilder::hasObstacle(Point a, Point b) {
	bool hasIntersection = false;
	for (int i=0;i<(nbPoints/2)-1;i++) {

		int ao = getOrientation(a,l_ptr[i]);
		int bo = getOrientation(b,l_ptr[i]);
		
		if ((ao!=0) && (bo!=0)) {
			if (ao!=bo) {	// a & b have different orientation, but do we have any line intersection between segment ab and line l? ... check further
				int flagIntersection = determineLineIntersection(a,b,l_ptr[i]);
				if ( (flagIntersection==HAS_INTERSECTION) || (flagIntersection==COLLINEAR) ) {
					hasIntersection = true;
					break;
				}
			}
		}
		else {
			if (ao==bo) { // colinear, a & b is part of the segment
				hasIntersection = true;
				break;
			}
		}
	}
	return hasIntersection;
};

//create distance table
void DistanceTableBuilder::createTable() {
	for (int i=0;i<nbPoints;i++) {
		for (int j=0;j<nbPoints;j++) {
			//ignore it, if they are the same point
			if (i==j) {
				distanceTable.set(i,j,INFINITY_VALUE);
			}
			else {
				//ignore it, if the segment ij has any of the obstacle segments
				if (hasObstacle(p_ptr[i],p_ptr[j])) {
					distanceTable.set(i,j,INFINITY_VALUE);
				}
				else 
					//set distance as euclidean distance
					distanceTable.set(i,j, calcEuclideanDistance(p_ptr[i],p_ptr[j]));
			}
		}
	}
        //cout << "end of create table";
};

void DistanceTableBuilder::print() {
	cout << "Distance Table (print only pairs which have euclidean distance)" <<"\n";
	for (int i=0; i < nbPoints; i++) {
		for (int j=0; j <nbPoints; j++) {
			float distVal = distanceTable.get(i,j);
			if (distVal!=INFINITY_VALUE)
				cout << " [" << i << "][" << j << "] = " <<distVal << "\n";
		}
	}
};

ShortPathFinder::ShortPathFinder(int s, int t, int n, Point *&p, Matrix table)  { 
	p_ptr = *&p;	
	distanceTable = table;
	source = s;
	target = t;
	nbPoints = n;
	distance = new float[nbPoints];
	predecessor = new int[nbPoints];
	
};

void ShortPathFinder::print() {
	for (int i=0; i < nbPoints; i++) {
		for (int j=0; j < nbPoints; j++) {
			cout << " [" << i << "][" << j << "] = " << distanceTable.get(i,j);
		}
		cout << "\n";
	}
};

void ShortPathFinder::printPath(int dest, int next, bool redraw) {
	if (predecessor[dest] != -1)
		printPath(predecessor[dest],dest,redraw);
        if (!redraw)
		cout << "  ===> " << dest << " , cost=" <<distance[dest] <<"\n";
	if (next!=-1) {
		Point p1 = p_ptr[next];
		Point p2 = p_ptr[dest];
		XDrawLine(display_ptr, oWin, gc_green, p1.get_x(), p1.get_y(), p2.get_x(), p2.get_y());	
	}
};

void ShortPathFinder::draw( Display *&disp, Window win, Colormap colormap, bool redraw)
{
	//set private variable
	display_ptr = *&disp;
	gc_green;
	oWin = win;

	unsigned long valuemask = 0;
	XGCValues gc_green_values;
	XColor tmp_color1, tmp_color2;

	gc_green = XCreateGC( display_ptr, win, valuemask, &gc_green_values);
  	XSetLineAttributes(display_ptr, gc_green, 6, LineSolid,CapRound, JoinRound);
  	if( XAllocNamedColor( display_ptr, colormap, "green", &tmp_color1, &tmp_color2 ) == 0 )
    		{printf("failed to get color green\n"); exit(-1);} 
  	else
    		XSetForeground( display_ptr, gc_green, tmp_color1.pixel );
	if (!redraw)
		cout << "Shortest Path with Djikstra Algorithm " <<"(from=" <<source <<" to=" <<target <<"):\n";
        printPath(target, -1, redraw);
};

//find short path using djikstra algorithm
void ShortPathFinder::find() {
	const bool visited = true;
	const bool unvisited = false;
	const int default_min_vertex = -1;

	int i; 
	int j;
	int minVertex; //store vertex whose current path length is minimal
	bool checkList[nbPoints]; 
	
	//initialize
	for (i = 0; i < nbPoints;i++) {
		distance[i] = INFINITY_VALUE;
		predecessor[i] = -1;
		checkList[i] = unvisited;
	}
	distance[source] = 0;	//distance source to source is 0
	
        for (i=0;i<nbPoints;i++) {
		minVertex = default_min_vertex;
		for (j=0;j<nbPoints;j++) {
			if (!checkList[j]) {
				if ((minVertex == default_min_vertex) || (distance[j] < distance[minVertex])) {
					minVertex = j;
					//cout << "selected min-vertex i=" <<i << " j=" << j << "\n";
				}
			}	
		}
		checkList[minVertex]=visited;
		for (j=0;j<nbPoints;j++) {
			float distVertex = distanceTable.get(minVertex,j);
			if (distVertex + distance[minVertex] < distance[j]) {
				//cout << "distance value j=" << j << " " << distVertex; 	
				distance[j] = distance[minVertex] + distVertex;
				predecessor[j] = minVertex; 
			}		
		}
	}
};

PWindow::PWindow(Display *ptr) {
	display_ptr = ptr;
};

void PWindow::createWindow(int argc, char **argv) {
	char *str_win_name = (char *) "Homework Project 1";
	char *icon_name_string = (char *) "Icon for Homework Project 1";
			
	screen_num = DefaultScreen( display_ptr );
  	screen_ptr = DefaultScreenOfDisplay( display_ptr );
  	color_map  = XDefaultColormap( display_ptr, screen_num );
  	display_width  = DisplayWidth( display_ptr, screen_num );
  	display_height = DisplayHeight( display_ptr, screen_num );

  	printf("Width %d, Height %d, Screen Number %d\n", 
           		display_width, display_height, screen_num);
  	/* creating the window */
  	border_width = 10;
  	win_x = 0; win_y = 0;
  	win_width = display_width/2;
  	win_height = (int) (win_width / 1.7); /*rectangular window*/
  
  	win= XCreateSimpleWindow( display_ptr, RootWindow( display_ptr, screen_num),
                            win_x, win_y, win_width, win_height, border_width,
                            BlackPixel(display_ptr, screen_num),
                            WhitePixel(display_ptr, screen_num) );
  			/* now try to put it on screen, this needs cooperation of window manager */
  	size_hints = XAllocSizeHints();
  	wm_hints = XAllocWMHints();
  	class_hints = XAllocClassHint();
  	if( size_hints == NULL || wm_hints == NULL || class_hints == NULL )
    	{ printf("Error allocating memory for hints. \n"); exit(-1);}

  	size_hints -> flags = PPosition | PSize | PMinSize  ;
  	size_hints -> min_width = 60;
  	size_hints -> min_height = 60;

  	XStringListToTextProperty( &str_win_name,1,&win_name);
  	XStringListToTextProperty( &icon_name_string,1,&icon_name);
  
  	wm_hints -> flags = StateHint | InputHint ;
  	wm_hints -> initial_state = NormalState;
  	wm_hints -> input = False;

  	class_hints -> res_name = (char *) "x_use_project1";
  	class_hints -> res_class = (char *) "project1";
  	
  	XSetWMProperties( display_ptr, win, &win_name, &icon_name, argv, argc,
                      size_hints, wm_hints, class_hints );

  	/* what events do we want to receive */
  	XSelectInput( display_ptr, win, 
            ExposureMask | StructureNotifyMask | ButtonPressMask );
  
  	/* finally: put window on screen */
  	XMapWindow( display_ptr, win );

  	XFlush(display_ptr);
}

/* utility function for main program */
void printEdges(Display *display_ptr, Window win, GC gc_grey, Matrix table, Point *points,int nbPoints) {
    for (int i=0; i < nbPoints; i++) {
	for (int j=0; j <nbPoints; j++) {
		if  (table.get(i,j)!=INFINITY_VALUE)
			XDrawLine(display_ptr, win, gc_grey, points[i].get_x(), points[i].get_y(), points[j].get_x(), points[j].get_y());	
			//cout << " [" << i << "][" << j << "] = " <<distVal << "\n";
		}
	}
	

};

/* window and graphic content variables */
Display *display_ptr;
char *display_name = NULL;
PWindow objWin = NULL;

XEvent report;

GC gc, gc_yellow, gc_red, gc_grey;
unsigned long valuemask = 0;
XGCValues gc_values, gc_yellow_values, gc_red_values, gc_grey_values;
XColor tmp_color1, tmp_color2;

/* core variables */
const int line_array_size = 100; //maximum there are 100 line obstacles	
//default value for capturing the points, it should reallocate as 2*current array size (NOT YET), max input 200 points + source + target = 202 points
const int point_array_size = 202;

Point *point_ptr = new Point[point_array_size]; //list of all captured points(x,y) from mouse click
Line *line_ptr = new Line[line_array_size];
int nbPoints = 0; //number of points from user click events
int source = -1;  //index of the source in list of point (point_ptr)
int target = -1;  //index of the target in list of point (point_ptr)

DistanceTableBuilder dtBuilder; /*the generator of distance table*/
ShortPathFinder pathFinder;	/*the finder of shortest path in distance table*/

/* main program */
int main(int argc, char **argv)
{
  bool endCaptureObstaclesPhase = false; //end of capturing obstacles phase
 	
  /* opening display: basic connection to X Server */
  if( (display_ptr = XOpenDisplay(display_name)) == NULL )
    { printf("Could not open display. \n"); exit(-1);}
  printf("Connected to X server  %s\n", XDisplayName(display_name) );
  objWin = PWindow(display_ptr);
  objWin.createWindow(argc, argv);

  /* create graphics context, so that we may draw in this window */
  gc = XCreateGC( display_ptr, objWin.win, valuemask, &gc_values);
  XSetForeground( display_ptr, gc, BlackPixel( display_ptr, objWin.screen_num ) );
  XSetLineAttributes( display_ptr, gc, 4, LineSolid, CapRound, JoinRound);
  /* and three other graphics contexts, to draw in yellow and red and grey*/
  
  gc_yellow = XCreateGC( display_ptr, objWin.win, valuemask, &gc_yellow_values);
  XSetLineAttributes(display_ptr, gc_yellow, 6, LineSolid,CapRound, JoinRound);
  if( XAllocNamedColor( display_ptr, objWin.color_map, "yellow", 
			&tmp_color1, &tmp_color2 ) == 0 )
    {printf("failed to get color yellow\n"); exit(-1);} 
  else
    XSetForeground( display_ptr, gc_yellow, tmp_color1.pixel );
  gc_red = XCreateGC( display_ptr, objWin.win, valuemask, &gc_red_values);
  XSetLineAttributes( display_ptr, gc_red, 6, LineSolid, CapRound, JoinRound);
  if( XAllocNamedColor( display_ptr, objWin.color_map, "red", 
			&tmp_color1, &tmp_color2 ) == 0 )
    {printf("failed to get color red\n"); exit(-1);} 
  else
    XSetForeground( display_ptr, gc_red, tmp_color1.pixel );
  gc_grey = XCreateGC( display_ptr, objWin.win, valuemask, &gc_grey_values);
  if( XAllocNamedColor( display_ptr, objWin.color_map, "light grey", 
			&tmp_color1, &tmp_color2 ) == 0 )
    {printf("failed to get color grey\n"); exit(-1);} 
  else
    XSetForeground( display_ptr, gc_grey, tmp_color1.pixel );
  

  /* and now it starts: the event loop */
  while((1))
    { XNextEvent( display_ptr, &report );
      switch( report.type ) {
	case Expose: {
		if (endCaptureObstaclesPhase) {
			if (source>-1) {
				XFillArc( display_ptr, objWin.win, gc_yellow, (point_ptr[source]).get_x()-objWin.win_height/40,  
					  (point_ptr[source]).get_y()-objWin.win_height/40, objWin.win_height/20, objWin.win_height/20, 0, 360*64);
			}
			if (target>-1) {
				XFillArc( display_ptr, objWin.win, gc_red, (point_ptr[target]).get_x()-objWin.win_height/40,  
					  (point_ptr[target]).get_y()-objWin.win_height/40, objWin.win_height/20, objWin.win_height/20, 0, 360*64);
				pathFinder.draw(display_ptr, objWin.win, objWin.color_map, true);
			}
		}
		
		if (nbPoints/2 > 0) {
			if (target > -1){
				for (int i=0; i < (nbPoints/2)-1;i++) {
					line_ptr[i].draw(display_ptr,objWin.win,gc);
				}
                        	//print all path posibilities with grey
				 printEdges(display_ptr, objWin.win, gc_grey, dtBuilder.getTable(), point_ptr, nbPoints);	
			}
			else
				for (int i=0; i < (nbPoints/2);i++)
					line_ptr[i].draw(display_ptr,objWin.win,gc);
		}
	}
            break;
        case ConfigureNotify: {
        	objWin.win_width = report.xconfigure.width;
        	objWin.win_height = report.xconfigure.height;
	}
            break;
        case ButtonPress:
        {  
            int x, y;
  	    x = report.xbutton.x;
            y = report.xbutton.y;
            if (report.xbutton.button == Button1 ) {
            	if (! endCaptureObstaclesPhase) {	
			nbPoints++;
			if (nbPoints>200) {//maximum 100 obstacles or 200 points
				cout << "Error: you can only have 100 obstacles\n";
				exit(0);
			}
	        	point_ptr[nbPoints-1] = Point(x,y);
	        	//if it's a pair then add it into list of segments
	        	if ((nbPoints & 1) == 0) { //nb of lines = nbPoints/2
	        		line_ptr[(nbPoints/2)-1] = Line(point_ptr[nbPoints-2], point_ptr[nbPoints-1]);
	        		line_ptr[(nbPoints/2)-1].draw(display_ptr,objWin.win,gc);
        		}
	       	}
	        else {
	        	if (source==-1) {
	        		nbPoints++;
	        		point_ptr[nbPoints-1] = Point(x,y);
	        		source = nbPoints - 1;
	        		XFillArc( display_ptr, objWin.win, gc_yellow, x - objWin.win_height/40, y - objWin.win_height/40,
					  objWin.win_height/20, objWin.win_height/20, 0, 360*64);
	        	}
	        	else {
	        		if (target==-1) {
	        			nbPoints++;
	        			point_ptr[nbPoints-1] =  Point(x,y);
	        			target = nbPoints - 1;
	        			//XFillArc( display_ptr, objWin.win, gc_red, x, y,objWin.win_height/20, objWin.win_height/20, 0, 360*64);
	        			XFillArc ( display_ptr, objWin.win, gc_red, x - objWin.win_height/40, y - objWin.win_height/40, 
						   objWin.win_height/20,  objWin.win_height/20, 0, 360*64 );
						
	        			//build the distance matrix after the users finish clicking event
	        			dtBuilder =  DistanceTableBuilder(source, target, nbPoints, point_ptr, line_ptr);
	        			dtBuilder.createTable();
	        			dtBuilder.print();
	        		        
					//print all path posibilities with grey
				 	printEdges(display_ptr, objWin.win, gc_grey, dtBuilder.getTable(), point_ptr, nbPoints);	
								
					//call the path finder
					pathFinder = ShortPathFinder(source, target, nbPoints, point_ptr, dtBuilder.getTable());
					pathFinder.find();
					pathFinder.draw(display_ptr, objWin.win, objWin.color_map, false);
	        		}
	        	}
	        }
            }
            if ((report.xbutton.button == Button2) || (report.xbutton.button == Button3))
            	endCaptureObstaclesPhase = true;
        }
        	break;
        default: 
          	break;
	} //end switch

    } //end while
  exit(0);
}
