#include<stdio.h>
//#include<conio.h>
#include<curses.h>
#include<iostream>
#include<fstream>
#include<iosfwd>
#include<string>
#include<deque>
#include<sstream>
#include<time.h>
#include<stdlib.h>
#include<filesystem>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilosys.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

#define TIMLIM 3600.0
#define min(a, b)  (((a) < (b)) ? (a) : (b)) 
#define max(a, b)  (((a) > (b)) ? (a) : (b)) 
#define EPSILON 1e-09
#define KK 0
#define LL 1
void write_report_file(string execname,
                       string filename,
                       IloNum capmrktshr,
                       IloNum lstmrktshr,
                       IloNum capmrkt,
                       IloNum lstmrkt,
                       IloNum cpu,
                       IloNum relgap,
                       IloInt nnodes,
                       IloNumArray y);

typedef struct PIKL
{
  IloInt k;
  IloInt l;
  IloNum x;
} PIKL;

struct comparer{
  inline bool operator()(const PIKL& a, const PIKL& b)
  {
     return a.x > b.x;
  }
};

ILOSTLBEGIN
IloEnv env;

typedef IloArray<IloNumArray> Num2DMatrix;
typedef IloArray<Num2DMatrix>Num3DMatrix;
typedef IloArray<Num3DMatrix>Num4DMatrix;
typedef IloArray<IloBoolVarArray> BoolVar2DMatrix;

typedef IloArray<IloBoolArray> Bool2DMatrix;
typedef IloArray<Bool2DMatrix> Bool3DMatrix;

typedef IloArray<IloNumVarArray> NumVar2DMatrix;
typedef IloArray<NumVar2DMatrix> NumVar3DMatrix;
typedef IloArray<NumVar3DMatrix> NumVar4DMatrix;

typedef IloArray<IloNumArray> Num2DMatrix;

typedef IloArray<IloIntArray> Int2DMatrix;
typedef IloArray<Int2DMatrix> Int3DMatrix;
typedef IloArray<Int3DMatrix> Int4DMatrix;

typedef vector<PIKL> PIKLArray;
typedef vector<PIKLArray> PIKL2DArray;
typedef vector<PIKL2DArray> PIKL3DArray;


class CData{
  /*
   * class for instance data
   *
   */
  public:
  double alpha;
  int p;
  IloInt N,c;
  IloIntArray ycomp;
  IloIntArray yc;
  Int2DMatrix x_comp;
  Num2DMatrix uc;

  Num2DMatrix h; // flow i-j
  Num2DMatrix dist_xy;
  Num4DMatrix u;
  Num4DMatrix pi;
 
  Num4DMatrix cost_xy;
  Num4DMatrix time_xy;

  PIKL3DArray sortedPi;
  IloNum total_flow;
  string filename; 
  CData(IloNum alpha=0.2)
  {
     this->p = 3;                        // number of hubs to install
     this->alpha = alpha;                // discount factor 0 < alpha < 1
     this->h = Num2DMatrix(env);         // demand i-j
     this->dist_xy = Num2DMatrix(env);   // distance i-j   
     this->yc = IloIntArray(env);
  }
  ~CData()
  {
     IloInt N = this->N;
     for (int i = 0;i < N;++i)
     {
        this->h[i].end();
        this->dist_xy[i].end();
        this->uc[i].end();
        for (int j = 0;j < N;++j)
        {
            for (int k = 0;k < N;++k)
            {
                this->cost_xy[i][j][k].end();
                this->time_xy[i][j][k].end();
            }
            this->cost_xy[i][j].end();
            this->time_xy[i][j].end();
        }
        this->x_comp[i].end();
        this->cost_xy[i].end();
        this->time_xy[i].end();
     }
     this->h.end();
     this->ycomp.end();
     this->yc.end();
     this->dist_xy.end();
     this->cost_xy.end();
     this->x_comp.end();
     this->uc.end();
     this->time_xy.end();
  }
  /*
   * read data file using standard io 
   *
   */
  void read_data(const char* filename)
  {
      this->filename.append(filename);
      if (std::filesystem::exists(this->filename.c_str()) == false)
      {
               cerr << endl << endl;
         cerr << "Please, provide a valid data file name " << endl;
         throw(-1);
      }
      fstream datafile;
      datafile.open(filename, ios::in);
      Num2DMatrix xytemp(env);
      datafile >> this->N >> this->yc >> xytemp >> this->h;
      IloInt nrows = xytemp.getSize();
      IloInt ncols = xytemp[0].getSize();
      this->dist_xy = Num2DMatrix(env,nrows);
      for (int i = 0; i < nrows; ++i)
      {
          this->dist_xy[i] = IloNumArray(env,nrows);
          for (int j = 0; j < nrows; ++j)
          {
              if (ncols == nrows)
              {
                 this->dist_xy[i][j] = xytemp[i][j];
              }
              else
              {
                 this->dist_xy[i][j] = pow( (xytemp[i][0] - xytemp[j][0]) * (xytemp[i][0] - xytemp[j][0]) + (xytemp[i][1] - xytemp[j][1]) * (xytemp[i][1] - xytemp[j][1]), 2);
              }
          }
      }

      for (int i = 0; i < nrows; ++i)
          xytemp[i].end();
      xytemp.end();


      this->p = this->yc.getSize();  
      this->ycomp = IloIntArray(env,this->N);
      for(int i = 0; i < this->N; ++i)
         this->ycomp[i] = 0;
      for(int i = 0; i < yc.getSize(); ++i)
         this->ycomp[this->yc[i]-1] = 1;
      this->total_flow = 0.0;
      for (int i = 0; i < this->N; ++i)
      {
          for (int j = 0; j < this->N; ++j)
          {
             this->total_flow += this->h[i][j];
          }
      }
      this->set_problem_data();
      this->find_competitors_attractiveness();
      this->set_pi_values();
      this->sort_pi_values();
  }
  /*
   * create utility data for each path i-k-l-j
   *
   */
  void set_problem_data()
  {
    IloInt N = this->N;
    IloInt A = 1000;
    IloNum be = 1;
    IloNum de = 1;
    IloNum ga = .5;

    IloNum c_collect = 1;
    IloNum c_distribute = 1;
    this->cost_xy = Num4DMatrix(env, N);
    this->time_xy = Num4DMatrix(env, N);
    this->u = Num4DMatrix(env,N); 
    this->pi = Num4DMatrix(env,N); 
    // create dynamic matrices
    for (int i = 0; i < N; i++)
    {
        this->cost_xy[i] = Num3DMatrix(env, N);
        this->time_xy[i] = Num3DMatrix(env, N);
        this->u[i] = Num3DMatrix(env, N);
        this->pi[i] = Num3DMatrix(env, N);
        for (int j = 0; j < N; j++)
        {
            this->cost_xy[i][j] = Num2DMatrix(env, N);
            this->time_xy[i][j] = Num2DMatrix(env, N);
            this->u[i][j] = Num2DMatrix(env, N);
            this->pi[i][j] = Num2DMatrix(env, N);
        
            for (int k = 0; k < N; k++)
            {
                this->cost_xy[i][j][k] = IloNumArray(env, N);
                this->time_xy[i][j][k] = IloNumArray(env, N);
                this->u[i][j][k] = IloNumArray(env, N);
                this->pi[i][j][k] = IloNumArray(env, N);
            }
        }
    }
    // create cost and time matrices
    for (int i = 0; i<N; i++)
    {
        for (int j = 0; j<N; j++)
        {
            for (int k = 0; k<N; k++)
            {
                for (int l = 0; l<N; l++)
                {
                    if (i != j)
                    {
                       this->cost_xy[i][j][k][l] = c_collect*dist_xy[i][k] + this->alpha*dist_xy[k][l] + c_distribute*dist_xy[l][j];
                       this->time_xy[i][j][k][l] = dist_xy[i][k] + dist_xy[k][l] + dist_xy[l][j];
                       if (k != l)
                       {
                          this->u[i][j][k][l] = A / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
                       }
                       else 
                       {
                          this->u[i][j][k][l] = A*1.25 / (ga*pow(time_xy[i][j][k][l], be) + (1 - ga)*pow(cost_xy[i][j][k][l], de));
                       }
                    }
                    else 
                    {
                       this->cost_xy[i][j][k][l] = 0;
                       this->time_xy[i][j][k][l] = 0;
                       this->u[i][j][k][l] = EPSILON;
                    }
               }
            }
        }
    }
  }
  /*
   *
   * compute the competitors attractiveness for each i-j
   *
   */
  void find_competitors_attractiveness()
  {
     IloInt N = this->N;
     IloNum min_cost;
     this->x_comp = Int2DMatrix(env,N); 
     this->uc = Num2DMatrix(env,N);
     IloNum objValue = 0.0; 
     for (int i = 0; i < N; ++i)
     {
        this->x_comp[i] = IloIntArray(env,N); 
        this->uc[i] = IloNumArray(env,N); 
        for (int j = 0; j < N; ++j)
	{
            this->x_comp[i][j] = 0;
	    this->uc[i][j] = 0.0;
	}
     }

     for (int hubk = 0; hubk < this->p; ++hubk)
     {
         for (int hubl = 0; hubl < this->p; ++hubl)
	 {
             IloInt k = this->yc[hubk] - 1; 
             IloInt l = this->yc[hubl] - 1; 
             this->x_comp[k][l] = 1;
             for (int i = 0; i < N; ++i)
                 for (int j = 0; j < N; ++j)
		     if (i != j)
			this->uc[i][j] += this->u[i][j][k][l];
         }
     }
  }
  /*
   * set pi values
   *
   */
  void set_pi_values()
  {
     IloInt N = this->N;
     for (int i = 0; i < N; ++i)
     {
        for (int j = 0; j < N; ++j)
        {
            for (int k = 0; k < N; ++k)
            {
                for (int l = 0; l < N; ++l)
                {
                   this->pi[i][j][k][l] = this->u[i][j][k][l]/this->uc[i][j];
                }
            }
        }
     }
  }
  /*
   * sort pi values for each i-j
   *
   */
  void sort_pi_values()
  {
      IloInt N = this->N;
      IloInt N2 = N * N;

      this->sortedPi = PIKL3DArray(N);
      for (int i = 0; i < N; ++i)
      {
          this->sortedPi[i] = PIKL2DArray(N);
          for (int j = 0; j < N; ++j)
          {
              this->sortedPi[i][j] = PIKLArray(N2);
              IloInt h = 0;
              for (int k = 0; k < N; ++k)
              {
                  for (int l = 0; l < N; ++l)
                  {
                      this->sortedPi[i][j][h].k = k;
                      this->sortedPi[i][j][h].l = l;
                      this->sortedPi[i][j][h].x = this->pi[i][j][k][l];
                      ++h;
                  }
              }
              sort(this->sortedPi[i][j].begin(),this->sortedPi[i][j].begin() + N2, comparer());
          }
      }
  }
};


void SOCMIN(const CData& dat)
{
  IloModel model(env);
  IloCplex cplex(model);
  IloInt N = dat.N;

  IloBoolVarArray Y(env, N);
  IloNumArray Y_val(env, N);

  NumVar2DMatrix X(env, N); 
  Num2DMatrix X_val(env, N); 

  NumVar2DMatrix V(env, N);
  Num2DMatrix V_val(env, N);

  NumVar2DMatrix R(env, N);
  Num2DMatrix R_val(env, N);
  // setting dynamic matrices
  for (int i = 0; i < N; i++)
  {
      X[i] = IloNumVarArray(env, N, 0, 1.0, ILOFLOAT);
      X_val[i] = IloNumArray(env, N);

      V[i] = IloNumVarArray(env, N, 0, 1.0, ILOFLOAT);
      V_val[i] = IloNumArray(env, N);
      
      R[i] = IloNumVarArray(env, N, 0, IloInfinity, ILOFLOAT);
      R_val[i] = IloNumArray(env, N);
  }
  // objective function
  char varname[250];
  IloExpr objective(env);
  for (int i = 0; i<N; ++i)
  {
      sprintf(varname,"y(%d)",i);
      Y[i].setName(varname);
      for (int j = 0; j<N; ++j)
      {
          sprintf(varname,"X(%d,%d)",i,j);
          X[i][j].setName(varname);
          sprintf(varname,"V(%d,%d)",i,j);
          V[i][j].setName(varname);
          sprintf(varname,"R(%d,%d)",i,j);
          R[i][j].setName(varname);
          objective += dat.h[i][j] * V[i][j];  
      }
  }
  model.add(IloMinimize(env,objective));
  objective.end();
  // only one path for each i-j
  for (int i = 0; i<N; ++i)
  {
      for (int j = 0; j<N; ++j)
      {
          IloExpr rykla(env);
	  rykla += X[i][j] - Y[i];
          model.add(rykla <= 0);
          rykla.end();

          IloExpr ryklb(env);
	  ryklb += X[i][j] - Y[j];
          model.add(ryklb <= 0);
          ryklb.end();
      }
  }
  // install up to p hubs
  IloExpr sump(env);
  for (int k = 0; k<N; ++k)
  {
      sump += Y[k];
  }
  model.add(sump == dat.p);
  sump.end();
  
  // R_ij = sum(k,l) pi_ijkl x_kl+1 
  for (int i = 0; i<N; ++i)
  {
      for (int j = 0; j<N; ++j)
      {
          if (i != j)
          {
             IloExpr rR(env);
	     rR += R[i][j];
             for (int k = 0; k<N; ++k)
             {
                 for (int l=0; l<N; ++l)
                 {
                     rR -= dat.pi[i][j][k][l] * X[k][l];
                 }
             }
             model.add(rR == 1);
             rR.end();
          }
      }
   }
   // conic constraint 1 
   for (int i = 0; i<N; ++i)
   {
       for (int j = 0; j<N; ++j)
       {
           if (i != j)
           {
              IloExpr conic1(env);
              conic1 += V[i][j] * R[i][j];
              model.add(conic1 >= 1.0);
              conic1.end();
           }
       }
   }
   // i = j
   for (int i = 0; i<N; ++i)
   {
           IloExpr varr(env);
	   varr += R[i][i];
	   model.add(varr == 0);
	   varr.end();

           IloExpr varv(env);
	   varv += V[i][i];
	   model.add(varv == 0);
	   varv.end();
   }
   //cplex.exportModel("soc2.lp");
   cplex.setParam(IloCplex::Param::TimeLimit, TIMLIM);
   cplex.setParam(IloCplex::Param::Threads, 1);
   IloTimer crono(env);
   cplex.setOut(env.getNullStream());
   cplex.setWarning(env.getNullStream());

   crono.start();
   cplex.solve();
   crono.stop();

   if (cplex.getStatus() == IloAlgorithm::Infeasible) // if the problem is infeasible
   {
      cerr << "Problem Infeasible" << endl;
   }
   IloNum UB = cplex.getObjValue();
   for (int k = 0; k < N; ++k)
   { 
       Y_val[k] = cplex.getValue(Y[k]);
       if (Y_val[k] > 1e-5)
	  Y_val[k] = 1;
       else
	  Y_val[k] = 0;
   }
   IloNum marketshare_lost = 100.0 * UB / dat.total_flow;

   write_report_file("mpmin",
                      dat.filename,
                      100.0 - marketshare_lost,
                      marketshare_lost,
                      dat.total_flow - UB,
                      UB,
                      crono.getTime(),
                      100.0 * cplex.getMIPRelativeGap(), 
                      cplex.getNnodes(),
                      Y_val);
   cout << endl;
   cout << endl;
   cout << "Second order conic formulation: minimization (ours)" << endl;
   cout << endl;
   cout << "data file            " << dat.filename << endl;
   cout << "Captured marketshare " << fixed << setprecision(2) << 100.0 - marketshare_lost << " %" << endl;
   cout << "Lost marketshare     " << fixed << setprecision(2) << marketshare_lost << " %" << endl;
   cout << "Captured market      " << setprecision(2) << dat.total_flow - UB << endl;
   cout << "Lost market          " << setprecision(2) << UB << endl;
   cout << "runTime              " << setprecision(2) << crono.getTime() << endl;
   cout << "MIP relative gap     " << setprecision(2) << 100.0 * cplex.getMIPRelativeGap() << "%" << endl;
   cout << "Number of bb nodes   " << cplex.getNnodes() << endl;
   cout << setprecision(0) << Y_val << endl;
   
}

void write_report_file(string execname,
                       string filename,
                       IloNum capmrktshr,
                       IloNum lstmrktshr,
                       IloNum capmrkt,
                       IloNum lstmrkt,
                       IloNum cpu,
                       IloNum relgap,
                       IloInt nnodes,
                       IloNumArray y)
{
   ofstream outfile;
   outfile.open("mpmin.dat",std::ofstream::app);
   outfile << execname;
   outfile << ","; 
   outfile << filename;
   outfile << ","; 
   outfile << fixed << setprecision(2) << capmrktshr;
   outfile << ","; 
   outfile << fixed << setprecision(2) << lstmrktshr;
   outfile << ","; 
   outfile << setprecision(2) << capmrkt;
   outfile << ","; 
   outfile << setprecision(2) << lstmrkt;
   outfile << ","; 
   outfile << setprecision(2) << cpu;
   outfile << ","; 
   outfile << setprecision(2) << relgap;
   outfile << ","; 
   outfile << nnodes; 
   //outfile << ","; 
   //outfile << setprecision(0) << y;
   outfile << endl;
   outfile.close();
}
                       

int main(int argc, char **argv)
{
  try
  {
     if (argc > 2)
     {
        IloNum alpha = atof(argv[2]);
        CData dat(alpha);

        dat.read_data(argv[1]);

	      SOCMIN(dat);
        
     }
     else
     {
        cerr << "usage:   " << argv[0] << " <datafile>  <alpha>" << endl;
     }
  }
  catch (IloException& ex)
  {
        cerr << "Error: " << ex << endl;
  }
  catch (...)
  {
        cerr << "Error" << endl;
  }
  return 0;
}
