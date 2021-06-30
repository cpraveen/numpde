#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define SIGN(a) (((a) < 0) ? -1 : 1)
#define Cp      (GAMMA * gas_const / (GAMMA - 1.0))

// Coefficients for RK3 scheme
const double arks[] = {0.0, 3.0/4.0, 1.0/3.0};
const double brks[] = {1.0, 1.0/4.0, 2.0/3.0};

// These values are set below based on test case type
double GAMMA;
double gas_const;

enum class SchemeOrder { first, second };
enum class Limiter { minmod, vanleer };

using namespace std;

//------------------------------------------------------------------------------
// To store some scheme and testcase parameters
//------------------------------------------------------------------------------
struct Parameters
{
   int test_case;
   SchemeOrder order;
   Limiter limiter;
};

//------------------------------------------------------------------------------
// Compute temperature given primitive variables
//------------------------------------------------------------------------------
double temperature(const vector<double>& prim)
{
   return prim[2] / (gas_const * prim[0]);
}

//------------------------------------------------------------------------------
// Minmod of three numbers
//------------------------------------------------------------------------------
double minmod (const double& a,
               const double& b,
               const double& c)
{
   int sa = SIGN(a);
   int sb = SIGN(b);
   int sc = SIGN(c);

   double result;
   if (sa == sb && sb == sc)
   {
      result  = min( min(fabs(a), fabs(b)), fabs(c) );
      result *= SIGN(a);
   }
   else
      result = 0.0;

   return result;

}

//------------------------------------------------------------------------------
// vanleer limiter
//------------------------------------------------------------------------------
double vanleer (const double& a,
                const double& b)
{
   double du;

   if(fabs(a*b) > 0.0)
      du = (SIGN(a)+SIGN(b))*fabs(a)*fabs(b)/(fabs(a) + fabs(b));
   else
      du = 0.0;

   return du;
}
//------------------------------------------------------------------------------
//   |  ul  |  uc  |  ur  |
//                 ^
//                 |
// Reconstruct left state of right face
//------------------------------------------------------------------------------
vector<double> muscl (const Limiter         limiter,
                      const vector<double>& ul,
                      const vector<double>& uc,
                      const vector<double>& ur)
{
   unsigned int n = ul.size();
   vector<double> result (n);
   double dul, duc, dur;
   // beta = 1.0 is TVD, stronger limiter
   // beta = 2.0 should give more accurate solution
   const double beta = 2.0;

   for(unsigned int i=0; i<n; ++i)
   {
      dul = uc[i] - ul[i];
      dur = ur[i] - uc[i];
      duc = 0.5 * (ur[i] - ul[i]);
      if(limiter == Limiter::vanleer)
         result[i] = uc[i] + 0.5 * vanleer (dul, dur);
      else
         result[i] = uc[i] + 0.5 * minmod (beta*dul, duc, beta*dur);
   }

   return result;
}

//------------------------------------------------------------------------------
// Main class of the problem
//------------------------------------------------------------------------------
class FVProblem
{
   public:
      FVProblem (const Parameters param);
      void run ();

   private:
      void make_grid_and_dofs ();
      void initial_condition ();
      void compute_dt ();
      void con_to_prim ();
      void reconstruct (const unsigned int face,
                        vector<double>& left,
                        vector<double>& right) const;
      void kfvs_split_flux (const vector<double>& prim,
                            const int sign,
                            vector<double>& flux) const;
      void num_flux (const vector<double>&,
                     const vector<double>&,
                     vector<double>&) const;
      void compute_residual ();
      void update_solution (const unsigned int rk);
      void output ();

      SchemeOrder order;
      Limiter     limiter;
      int         test_case;
      double d_left, u_left, p_left;
      double d_right, u_right, p_right;
      vector<double> prim_left;
      vector<double> prim_right;
      unsigned int n_var;
      unsigned int n_cell;
      unsigned int n_face;
      double xmin;
      double xmax;
      double xmid;
      double dx;
      double dt;
      double cfl;
      double final_time;
      vector<double> xc;
      vector<double> xf;

      vector< vector<double> > primitive;     // density, velocity, pressure
      vector< vector<double> > residual;
      vector< vector<double> > conserved;     // solution at time level n+1
      vector< vector<double> > conserved_old; // solution at time level n
};

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
FVProblem::FVProblem (const Parameters param)
{
   n_var  = 3;

   test_case = param.test_case;
   order = param.order;
   limiter = param.limiter;

   if(test_case == 1)
   {
      // Sod shock tube case
      GAMMA = 1.4;
      gas_const = 1.0;
      final_time = 0.2;
      n_cell = 100;
      cfl    = 0.8;

      xmin    = 0.0;
      xmax    = 1.0;
      xmid    = 0.3;

      d_left  = 1.0;
      d_right = 0.125;

      u_left  = 0.75;
      u_right = 0.0;

      p_left  = 1.0;
      p_right = 0.1;
   }

   prim_left.resize (n_var);
   prim_right.resize (n_var);

   prim_left[0] = d_left;
   prim_left[1] = u_left;
   prim_left[2] = p_left;

   prim_right[0] = d_right;
   prim_right[1] = u_right;
   prim_right[2] = p_right;

}

//------------------------------------------------------------------------------
// Allocate memory for grid and create grid
//------------------------------------------------------------------------------
void FVProblem::make_grid_and_dofs ()
{
   cout << "Making grid and allocating memory ...\n";

   n_face = n_cell + 1;
   dx = (xmax - xmin) / n_cell;
   xc.resize (n_cell);
   xf.resize (n_face);

   // Make grid
   for(unsigned int i=0; i<n_face; ++i)
      xf[i] = xmin + i * dx;
   for(unsigned int i=0; i<n_cell; ++i)
      xc[i] = 0.5 * (xf[i] + xf[i+1]);

   primitive.resize (n_cell, vector<double>(n_var));
   residual.resize (n_cell, vector<double>(n_var));
   conserved.resize (n_cell, vector<double>(n_var));
   conserved_old.resize (n_cell, vector<double>(n_var));
}

//------------------------------------------------------------------------------
// Set initial condition
//------------------------------------------------------------------------------
void FVProblem::initial_condition ()
{
   cout << "Setting initial conditions ...\n";

   // Set initial condition
   for(unsigned int i=0; i<n_cell; ++i)
   {
      if(xf[i+1] <= xmid)
      {
         conserved[i][0] = d_left;
         conserved[i][1] = d_left * u_left;
         conserved[i][2] = p_left/(GAMMA-1.0) + 0.5 * d_left * pow(u_left,2);
      }
      else
      {
         conserved[i][0] = d_right;
         conserved[i][1] = d_right * u_right;
         conserved[i][2] = p_right/(GAMMA-1.0) + 0.5 * d_right * pow(u_right,2);
      }
   }
}

//------------------------------------------------------------------------------
// Compute time step
//------------------------------------------------------------------------------
void FVProblem::compute_dt ()
{
   dt = 1.0e20;
   for(unsigned int i=0; i<n_cell; ++i)
   {
      double speed = fabs(primitive[i][1]) +
                     sqrt(GAMMA * primitive[i][2] / primitive[i][0]);
      dt = min (dt, dx/speed);
   }

   dt *= cfl;
}

//------------------------------------------------------------------------------
// Convert conserved to primitive
//------------------------------------------------------------------------------
void FVProblem::con_to_prim ()
{
   for(unsigned int i=0; i<n_cell; ++i)
   {
      primitive[i][0] = conserved[i][0];
      primitive[i][1] = conserved[i][1]/conserved[i][0];
      primitive[i][2] = (GAMMA-1.0) * (conserved[i][2] -
                           0.5 * pow(conserved[i][1], 2.0) / conserved[i][0]);
   }
}

//------------------------------------------------------------------------------
// Reconstruct left/right state at a face
//------------------------------------------------------------------------------
void FVProblem::reconstruct (const unsigned int face,
                             vector<double>&    left,
                             vector<double>&    right) const
{
   if(order == SchemeOrder::first)
   {
      left  = primitive[face-1];
      right = primitive[face];
      return;
   }

   if(face==1)
   {
      left = primitive[0];
      right = muscl (limiter, primitive[2], primitive[1], primitive[0]);
   }
   else if(face==n_face-2)
   {
      left  = muscl (limiter, primitive[n_cell-3], primitive[n_cell-2],
                     primitive[n_cell-1]);
      right = primitive[n_cell-1];
   }
   else
   {
      left  = muscl (limiter, primitive[face-2], primitive[face-1], primitive[face]);
      right = muscl (limiter, primitive[face+1], primitive[face], primitive[face-1]);
   }
}

//------------------------------------------------------------------------------
// KFVS split fluxes: sign=+1 give positive flux and
// sign=-1 gives negative flux
//------------------------------------------------------------------------------
void FVProblem::kfvs_split_flux (const vector<double>& prim,
                                 const int             sign,
                                 vector<double>&       flux) const
{
   double beta, s, A, B, E, fact;

   beta = 0.5 * prim[0] / prim[2];
   s    = prim[1] * sqrt(beta);
   A    = 0.5 * (1.0 + sign * erf(s));
   B    = sign * 0.5 * exp(-s * s) / sqrt(beta * M_PI);
   E    = prim[2]/(GAMMA-1.0) + 0.5 * prim[0] * pow(prim[1], 2);
   fact = prim[1] * A + B;

   // inviscid flux
   flux[0] = prim[0] * fact;
   flux[1] = (prim[2] + prim[0] * pow(prim[1], 2)) * A +
             prim[0] * prim[1] * B;
   flux[2] = prim[1] * (E + prim[2]) * A +
             (E + 0.5 * prim[2]) * B;
}

//------------------------------------------------------------------------------
// Numerical flux function
//------------------------------------------------------------------------------
void FVProblem::num_flux(const vector<double>& left,
                         const vector<double>& right,
                         vector<double>&       flux) const
{
   vector<double> flux_pos(n_var);
   vector<double> flux_neg(n_var);

   kfvs_split_flux (left,  +1, flux_pos);
   kfvs_split_flux (right, -1, flux_neg);

   for(unsigned int i=0; i<n_var; ++i)
      flux[i] = flux_pos[i] + flux_neg[i];


}

//------------------------------------------------------------------------------
// Compute finite volume residual
//------------------------------------------------------------------------------
void FVProblem::compute_residual ()
{
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         residual[i][j] = 0.0;

   vector<double> flux (n_var);
   vector<double> left (n_var);
   vector<double> right(n_var);

   // Flux through left boundary
   num_flux (prim_left, primitive[0], flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[0][j] -= flux[j];

   // Flux through interior faces
   for(unsigned int i=1; i<n_face-1; ++i)
   {
      reconstruct (i, left, right);
      num_flux (left, right, flux);
      for(unsigned int j=0; j<n_var; ++j)
      {
         residual[i-1][j] += flux[j];
         residual[i][j]   -= flux[j];
      }
   }

   // Flux through right boundary
   num_flux (primitive[n_cell-1], prim_right, flux);
   for(unsigned int j=0; j<n_var; ++j)
      residual[n_cell-1][j] += flux[j];
}

//------------------------------------------------------------------------------
// Compute finite volume residual
// Apply one stage of SSPRK3
//------------------------------------------------------------------------------
void FVProblem::update_solution (const unsigned int rk)
{
   const double lam = dt / dx;
   for(unsigned int i=0; i<n_cell; ++i)
      for(unsigned int j=0; j<n_var; ++j)
         conserved[i][j] = arks[rk] * conserved_old[i][j] +
                           brks[rk] * (conserved[i][j] - lam * residual[i][j]);

}

//------------------------------------------------------------------------------
// Save solution to file
//------------------------------------------------------------------------------
void FVProblem::output ()
{
   cout << "Saving solution to sol.dat\n";

   ofstream fo("sol.dat");
   for(unsigned int i=0; i<n_cell; ++i)
      fo << xc[i] << " "
         << primitive[i][0] << " "
         << primitive[i][1] << " "
         << primitive[i][2] << "\n";
   fo.close ();
}

//------------------------------------------------------------------------------
// Start the computations
//------------------------------------------------------------------------------
void FVProblem::run ()
{
   make_grid_and_dofs ();
   initial_condition ();
   con_to_prim ();

   double time = 0.0;
   unsigned int iter = 0;
   while (time < final_time)
   {
      conserved_old = conserved;
      compute_dt ();
      if(time+dt > final_time) dt = final_time - time;
      for(unsigned int rk=0; rk<3; ++rk)
      {
         compute_residual ();
         update_solution (rk);
         con_to_prim ();
      }
      time += dt;
      ++iter;
      if(iter % 1000 == 0)
      cout << "Iter = " << iter << " Time = " << time << endl;
      if(iter % 1000 == 0) output ();
   }

   con_to_prim ();
   output ();
}

//------------------------------------------------------------------------------
void  get_command_line(int argc, char *argv[],
                       Parameters& param)
{
   if(argc == 1)
   {
      cout << "Not enough arguments\n";
      cout << "  first order: " << argv[0] << " testcase\n";
      cout << "  high order : " << argv[0] << " testcase  limiter\n";
      cout << "  Testcases  : sod\n";
      cout << "  Limiters   : minmod, vanleer\n";
      exit(0);
   }

   if(string(argv[1]) == "sod")
   {
      param.test_case = 1;
   }
   else
   {
      cout << "Unknown test case = " << argv[1] << endl;
      exit(0);
   }

   if(argc == 2)
   {
      param.order = SchemeOrder::first;
      cout << "First order scheme\n";
   }
   else if(argc == 3)
   {
      param.order = SchemeOrder::second;
      cout << "Second order scheme\n";
      string val = string(argv[2]);
      if(val == "minmod")
         param.limiter = Limiter::minmod;
      else if(val == "vanleer")
         param.limiter = Limiter::vanleer;
      else
      {
         cout << "Unknown limiter = " << val << endl;
         exit(0);
      }
   }
}

//------------------------------------------------------------------------------
// This is where it all starts
//------------------------------------------------------------------------------
int main (int argc, char *argv[])
{
   Parameters param;
   get_command_line(argc, argv, param);
   FVProblem fv_problem(param);
   fv_problem.run ();

   return 0;
}
