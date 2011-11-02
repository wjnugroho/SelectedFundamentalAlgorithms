#include <iostream> //provide cout
#include <cstdlib>  //provide EXIT_SUCCESS
#include <cassert>  //provide assert function
using namespace std;

//a pseudo-multidimensional array for storing coefficient matrix
class Matrix {
  private: 
    int row, col;
    double* value;
  public:
    Matrix() {};
    Matrix(int r, int c) {
      row = r; col = c;
      value = new double[row * col];
    }
    Matrix(int r, int c, double &A) {
      row = r; col = c;
      value = &A;
    }
    double get(int r, int c) {
      return value[r*col + c];
    }
    void set(int r, int c, double v) {
      value[r*col + c] = v;
    }
    int getNbRow() { return row;}
    int getNbCol() { return col;}
    double *getMatrix() {
      return value;
    }
    void print() {
      for (int i=0;i<row;i++) {
        //print header
	if (i==0) {
		for (int k =0; k < col; k++) {
			if (k==col-1)
				cout << "rhs";
			else
				cout << "r" << k << "\t";
		}
		cout << "\n";
	} 
	//print values
	for (int j=0;j<col;j++) {
		cout << get(i,j) << "\t";
	}
	cout << "\n";
      }
    }
    void freeMatrix() {
      delete[] value; //it's done, free memory pointed by value
      value = NULL;   //clear value to prevent using invalid memory reference
    }
};

class Tableau {
  private:
    Matrix matrix;
    int nbVars;
    int nbInEqualities;
    int* activeVars;
    int pivotCol;
    int pivotRow;
    int nbSteps;
  public:
    bool isOptimal;
    Tableau() {};
    Tableau(Matrix coeff, double &r, double &o) {
	double* rhs = &r;
	double* obj = &o;
	int i,j;
	isOptimal = false;
        nbSteps = 0;
	nbVars = coeff.getNbCol();
        nbInEqualities = coeff.getNbRow();

        //table column = nb-inequalities + nb-variable + 1 coeffobj + 1 result
	int col = coeff.getNbRow() + coeff.getNbCol() + 2;
	//table row    = nb-inequalities + 1 obj formula
	int row = coeff.getNbRow() + 1;
	matrix = Matrix(row, col);
	//set list of active variables
	activeVars = new int[row];	
	for (i=0;i<row;i++) {
	  activeVars[i] = nbVars+i;
          //cout << "var: " << i-coeff.getNbRow() << "=" << i;
	}	
	//cout << " end init " << "\n";
	//copy coeff into tableau
	for (i=0;i<coeff.getNbRow();i++) {
	  for (j=0;j<coeff.getNbCol();j++) {
		matrix.set(i, j, coeff.get(i, j));              
	  }
	}

	//copy obj into tableau
        for (i=0;i<coeff.getNbCol();i++)
	  matrix.set(row-1, i, obj[i]*(-1.0));

	//set slack variables
	for (i=0;i<row;i++) {
          for (j=coeff.getNbCol();j<col-1;j++) {
	    if (i==(j-coeff.getNbCol()))
	      matrix.set(i, j, 1.0);
            else
              matrix.set(i, j, 0.0);
          }
        }

        //copy right handside into tableau
        //cout << "rhs" << col << "\n";
	for (i=0;i<coeff.getNbRow();i++)
	  matrix.set(i, col-1, rhs[i]);	
    }
    //print the tableau into screen
    void print() {
	matrix.print();
    }
    int getNbSteps() { return nbSteps; }
    //set the pivotCol and pivotRow
    void selectPivot(){
	//a pivot column = the largest magnitude column value in the bottom row (excluding the right hand side column/last column)
        pivotCol=0;
        int obj = matrix.getNbRow()-1;
	int rhs = matrix.getNbCol()-1; 
	for (int i=1;i<rhs;i++) {
		if (matrix.get(obj, pivotCol) > matrix.get(obj, i))
			pivotCol = i;
	}
	nbSteps++;
	if (matrix.get(obj, pivotCol) >= 0.0) { //hooray no negative values, it's optimal solution!
		isOptimal = true;
		//cout << "[ solution is optimal ]\n";
	}
	else {
		//a pivot row = smallest ratio of rhs/matrix.get(row, pivotCol)
		double ratio1, ratio2 = 0.0;
		pivotRow = -1;
		for (int j=0;j<obj;j++) {
			if (matrix.get(j, pivotCol) > 0) // only take positive row
				if (pivotRow==-1) {
					ratio1 = ((double) matrix.get(j, rhs)/matrix.get(j, pivotCol));
					pivotRow = j;
				}
				else {
					ratio2 = ((double) matrix.get(j, rhs)/matrix.get(j, pivotCol));
					if (ratio2<ratio1) {
						ratio1 = ratio2;
						pivotRow = j;
					}
				}
		}	
		//cout << "[ pivot row:" << pivotRow << " col:" << pivotCol << " ]\n";
		//change activeVars
		activeVars[pivotRow] = pivotCol;
	}
    }

    //use the pivot to clear the pivot column
    void clearColumn() {
	double pivotDivFactor = 1.0/matrix.get(pivotRow, pivotCol);
	for (int i=0;i<matrix.getNbCol();i++) {
		matrix.set(pivotRow, i, matrix.get(pivotRow, i)*pivotDivFactor);
	}

	double r, p, v;
	for (int j=0;j<matrix.getNbRow();j++) {
		if (j!=pivotRow) {
			r = matrix.get(pivotRow, pivotCol);
			p = matrix.get(j, pivotCol);
			//if (r*p > 0.0)
			//	p = p * (-1.0);
			if (p!=0.0) {
				for (int k=0;k<matrix.getNbCol();k++) {
					v = matrix.get(j, k);
					matrix.set(j, k, v - p * matrix.get(pivotRow, k));
				}
			}	
		}
	}
    }
    
    double* getResult() {
	double* r=new double[nbVars];
	for (int i=0;i<nbVars;i++)
		r[i]=0.0;

	for (int i=0;i<matrix.getNbRow()-1;i++) {
		//cout << "active var:" << i << "value:" << activeVars[i] << "\n";
		if (activeVars[i] < nbVars) {
			int var=activeVars[i];
                        //cout << "var=" << var << "atas:" <<  matrix.get(i, matrix.getNbCol() - 1) << "bawah:" << matrix.get(i, activeVars[i]) << "\n";
			r[var] = matrix.get(i, matrix.getNbCol() - 1)/matrix.get(i, activeVars[i]) ;
			//cout << "var: " << var << " value: " << r[var] <<"\n";
		}
	}
        return r;
    }
};

// d=number of variables
// n=number of inequalities
// A=coefficient matrix
// b=right-hand side
// c=coefficient objective function
void simplex(int d, int n, double *A, double *b, double *c, double* result, int *steps) {
  //use pseudo matrix ... for accessible access
  Matrix m = Matrix(n, d, *A);

  // form a tableau corresponding to a basic feasible solution (BFS) of matrix (n+1) x (d+1) 
  Tableau tableau = Tableau(m, *b, *c);
  //cout << "==> Initial Tableau:  \n";
  while (! tableau.isOptimal) {
	//tableau.print();
	tableau.selectPivot();
	if (! tableau.isOptimal) {
		tableau.clearColumn();
	}
  }

  //get the steps
  *steps = tableau.getNbSteps();
 
  //copy the result into vector result which contains the optimum values for the d variables
  double* tmpResult = tableau.getResult(); 
  for (int i=0;i<d;i++) {
  	result[i] = tmpResult[i];
  }
  
};
