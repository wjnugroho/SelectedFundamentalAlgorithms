/* compiles with command line  g++ project4.cpp -lX11 -lm -L/usr/X11R6/lib 
   Implementation of randomized server algorithm.
   
   my strategy   = a randomized server selection based on the probability of inverse distance (1)
   compared with = a work-balance (select the server which did travel least)                  (2)
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

PWindow::PWindow(Display *ptr) {
	display_ptr = ptr;
};

void PWindow::createWindow(int argc, char **argv) {
	char *str_win_name = (char *) "Homework Project 4";
	char *icon_name_string = (char *) "Icon for Homework Project 4";
			
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

  	class_hints -> res_name = (char *) "x_use_project4";
  	class_hints -> res_class = (char *) "project4";
  	
  	XSetWMProperties( display_ptr, win, &win_name, &icon_name, argv, argc,
                      size_hints, wm_hints, class_hints );

  	/* what events do we want to receive */
  	XSelectInput( display_ptr, win, 
            ExposureMask | StructureNotifyMask | ButtonPressMask );
  
  	/* finally: put window on screen */
  	XMapWindow( display_ptr, win );

  	XFlush(display_ptr);
}

/*** new classes for homework 4 ***/
class Server {
    private:
        //drawing variable
        Display *d_ptr; Window win; //unsigned int wWidth; unsigned int wHeight;
        GC gc_grey;
        GC gc_red;

        Point wbInitPosition;       //initial position, required for redrawing
        Point wbPosition;           //current position by work balance strategy
        double wbTravelDistance;    //total distance travelled by work balance strategy
        
        Point myInitPosition;       //initial position
        Point myPosition;           //current position by selected strategy
        double myTravelDistance;    //total distance travelled by selected strategy
        int nbMyTravel;             //number of total travel
        
        void drawPath(Point p1, Point p2) {
            XDrawLine(d_ptr, win, gc_grey, p1.get_x(), p1.get_y(), p2.get_x(), p2.get_y());	
        };
       
        int calcEuclideanDistance(Point a, Point b) {
	    double diffx = (double) (b.get_x() - a.get_x());
	    double diffy = (double) (b.get_y() - a.get_y());
	    double distance = sqrt( diffx * diffx + diffy * diffy );	
	    return distance;
        };

    public:
        Server() {};            //random initial position
        Server(Display *display_ptr, Window w, GC gc_grey_in, GC gc_red_in, Point p) {
            myInitPosition = wbInitPosition = p; //set init position ...
            myPosition = wbPosition = p;            
            myTravelDistance = wbTravelDistance = 0.0;
            nbMyTravel =0;
            d_ptr = display_ptr;
            win = w;
            gc_grey = gc_grey_in;
            gc_red = gc_red_in;
        };
        //print my strategy only
        void drawServer() {
            int x = myPosition.get_x();
            int y = myPosition.get_y();
            XFillArc( d_ptr, win, gc_red, x - 4, y - 4, 8, 8, 0, 360*64);
        };
        void moveByWbStrategy(Point &p) { 
            wbTravelDistance+=calcEuclideanDistance(wbPosition,p);
            wbPosition = p;         
        };

        void moveByMyStrategy(Point &p) {
            nbMyTravel+=1;
            drawPath(myPosition, p);
            myTravelDistance+=calcEuclideanDistance(myPosition,p);
            drawServer();
            myPosition = p;
            drawServer();
        };

        double getWbDistanceTravelled() { return wbTravelDistance; };
        double getWbCost(Point p) { return calcEuclideanDistance(wbPosition, p); }

        double getMyDistanceTravelled() { return myTravelDistance; };
        double getMyCost(Point p) { return calcEuclideanDistance(myPosition,p); }
        double getMyWeighted() {
            double w = 1.0;
            if (nbMyTravel==0)
                w = 1.0 / nbMyTravel;
            return w;
        }
        void reset() {
            myPosition = myInitPosition;
            wbPosition = wbInitPosition;          
            myTravelDistance = wbTravelDistance = 0.0;
            nbMyTravel =0;
        }
  };

class KServer {
    private:
        int nbServers;
        Server *servers;
        double *weighting;
        double wbTotalCost;
        double myTotalCost;
    public:
        KServer() {};
        KServer(int n) { nbServers=n; servers = new Server[n]; wbTotalCost=myTotalCost=0.0; weighting=new double[n]; };
        void initialize(Display *display_ptr, Window w, unsigned int width, unsigned int height, GC gc_p, GC gc_s) {
            //for (int i=0;i<nbServers;i++) {
            int i = 0;
            while (i < nbServers) {
                srand(time(NULL));
                int x = rand() % (width - (i*10)) + i;  //pick random position between 10 and max(width/height - 10)
                int y = rand() % (height - (i*10)) + i;
                if (x>1 && x<(width-1) && y>1 && y<(height-1)) { // 1<x<width, 1<y<height
                    servers[i] = Server(display_ptr, w, gc_p, gc_s, Point(x,y));
                    i++;
                }
            }
        };
        void draw() {
            for (int i=0;i<nbServers;i++) {
                servers[i].drawServer();
            }
        };

        //use the work-balance strategy ... always moves the server who did travel least
        void selectWbServer(Point p) {
            int cserver = 0;
            double ctravel = servers[cserver].getWbDistanceTravelled();
            double cpoint = servers[cserver].getWbCost(p);
            for (int i=0;i<nbServers;i++) {
                double itravel = servers[i].getWbDistanceTravelled();
                double ipoint  = servers[i].getWbCost(p);
                if (ctravel > itravel) {
                    cserver = i;
                    ctravel = itravel;
                    cpoint  = ipoint;
                }
            }
            wbTotalCost += cpoint;
            servers[cserver].moveByWbStrategy(p);
        };
        
        //use random selection with the probability of inverse distance        
        void selectMyServer(Point p) {
            int cserver=0;
            bool found = false;
            //assign inverse distance as probability factor
            double totalWeight=0.0;
            for (int i=0;i<nbServers;i++) {
                double cpoint = servers[i].getMyCost(p);
                if (cpoint==0.0) {
                    cserver=i;
                    found = true;
                    break;
                }
                else
                    weighting[i]=1/cpoint;
                totalWeight+=weighting[i];
                //cout << i << "-wo-" << weighting[i] << "-co-" << cpoint << "\n";
            }
            //cout << "total" << totalWeight << "\n";
            //if not found get random the server based on the probability value
            if (!found) {
                int sortedWeight[nbServers];
                int temp;
                //initialize sortedWeight & normalize weighting
                for (int i=0;i<nbServers;i++) {
                    weighting[i]=weighting[i]/totalWeight;
                    //cout << i <<"-no-" << weighting[i] << "\n";
                    sortedWeight[i]=i;
                }
                //sort weight in ascending order
                for (int i=0;i<(nbServers-1);i++) {
                    for(int j =(i+1);j<nbServers;j++)
                    {
                        if (weighting[sortedWeight[i]] > weighting[sortedWeight[j]]) // ascending order
                        {
                            temp=sortedWeight[i];          // swap
                            sortedWeight[i] = sortedWeight[j];
                            sortedWeight[j] = temp;
                        }
                    }
                }
                double rval = rand() / double(RAND_MAX);
                // cout << "random val" << rval;
                totalWeight =0.0;
                for(int i=0;i<nbServers;i++) {
                    cserver=sortedWeight[i];
                    totalWeight += weighting[cserver];
                    //cout << "totalw:" << totalWeight << "\n";
                    if (rval<=totalWeight)
                        break;
                }
                //cout << "rand:" << rval << " server:" << cserver << "\n";
            }
            myTotalCost +=servers[cserver].getMyCost(p);
            servers[cserver].moveByMyStrategy(p);
        };

        void printTotalDistanceTravel() {
            cout << "\n==>Distance travel:\n";
            cout << "work-balance strategy:\n";
            cout << "[";
            for(int i=0;i<nbServers;i++) {
                if (i==0)
                    cout << i << ":" << servers[i].getWbDistanceTravelled();
                else
                    cout << "  " << i << ":" << servers[i].getWbDistanceTravelled();
            }
            cout << "] =>Total:" << wbTotalCost;
            cout << "\n";
            cout << "my strategy: randomized server with the probability of inverse distance\n";
            cout << "[";
            for(int i=0;i<nbServers;i++) {
                if (i==0)
                    cout << i << ":" << servers[i].getMyDistanceTravelled();
                else
                    cout << "  " << i << ":" << servers[i].getMyDistanceTravelled();
            }
            cout << "] =>Total:" << myTotalCost;
            cout << "\n";

        }
        //call this before redrawing
        void reset() {
            wbTotalCost=myTotalCost=0.0;
            for (int i =0; i < nbServers; i++) {
                servers[i].reset();
            }
        }
};

/* simple linkedlist for storing list of point with insertion only */
struct Node {
    Point point;
    Node *next;
};

struct Node *head;
struct Node *tail;
void addNode(Point p) {
    Node *newNode = new Node;
    newNode->point = p;
    newNode->next = NULL;
    if (head==NULL) {
        tail = head = newNode;
    }
    else {
        Node *curr = tail;
        tail = newNode;
        curr->next = tail;
    }
    return;
};


/* window and graphic content variables */
Display *display_ptr;
char *display_name = NULL;
PWindow objWin = NULL;

int nServers;
KServer kserver;

XEvent report;

GC gc, gc_yellow, gc_red, gc_grey;
unsigned long valuemask = 0;
XGCValues gc_values, gc_yellow_values, gc_red_values, gc_grey_values;
XColor tmp_color1, tmp_color2;

/* core variables */
//const int line_array_size = 100; //maximum there are 100 line obstacles	
//default value for capturing the points, it should reallocate as 2*current array size (NOT YET), max input 200 points + source + target = 202 points
//const int point_array_size = 202;

//Point *point_ptr = new Point[point_array_size]; //list of all captured points(x,y) from mouse click
//Line *line_ptr = new Line[line_array_size];
//int nbPoints = 0; //number of points from user click events
//int source = -1;  //index of the source in list of point (point_ptr)
//int target = -1;  //index of the target in list of point (point_ptr)

//DistanceTableBuilder dtBuilder; /*the generator of distance table*/
//ShortPathFinder pathFinder;	/*the finder of shortest path in distance table*/

/* main program */
int main(int argc, char **argv)
{
  //bool endCaptureObstaclesPhase = false; //end of capturing obstacles phase
 	
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
  
  if (argc==1) {
      cout << "please add the number of servers in the parameter";
      exit(-1); 
  }
  nServers = atoi(argv[1]);

  cout << "Nb Servers:" << nServers << "\n";
  kserver = KServer(nServers);
  kserver.initialize(display_ptr, objWin.win, objWin.win_width, objWin.win_height, gc_grey, gc_red);

  /* and now it starts: the event loop */
  while((1))
    { XNextEvent( display_ptr, &report );
      switch( report.type ) {
	case Expose: {
	
                kserver.draw();
                kserver.reset();
                Node *cur = head;
                while (cur!=NULL) {
                    kserver.selectWbServer(cur->point);
                    kserver.selectMyServer(cur->point);
                    cur = cur->next;
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

            addNode(Point(x,y));  //for redrawing
            kserver.selectWbServer(Point(x,y));
            kserver.selectMyServer(Point(x,y));
            kserver.printTotalDistanceTravel();
        }
            break;
        default: 
            break;
	} //end switch

    } //end while
  exit(0);
}
