#include "SDKTSolver.h"
#include "SDKTSolverSWPP.h"
#include "SWFourEqn.h"
#include "RK2TimeStepper.h"

#include <cmath>
#include <fenv.h>
#include <sys/time.h>
#include <iomanip>
#include <vector>
#include <sstream>
#include <functional>


class ChannelRollWave
{
	typedef SWMuIEqn1DRhoVary<0> Eqn;
	// typedef SWMuIvEqn1DViscous VisEqn;
	typedef LimiterWENO LIMITER; 
	typedef RK2TimeStepper TIMESTEPPER;
	typedef SDKTSolverSWPP<Eqn, LIMITER> Solver;

public:
	ChannelRollWave(double h0_, double domainLength_) 
	{
		h0 = h0_;
		domainLength = domainLength_;
	}

	void Run(int n) //, double finalTheta, double change_t
	{
		Eqn eqn(9.81, 15, 0.0, 1e-4, 1e-5); //, finalTheta, change_t
		eqn.SetMuIParams(BoyerRockWater);
		eqn.EnableStoppedMaterialHandling();
		// std::string limiterName = LIMITER::Name();
		// double CFL_in = 0.03;
		// eqn.RegisterParameter("CFL", Parameter(&CFL_in, true));
		eqn.EnableInDirectoryName("theta");

		// eqn.EnableInDirectoryName("initTheta");
		// eqn.EnableInDirectoryName("finalTheta");
		// eqn.EnableInDirectoryName("change_t");

		eqn.EnableInDirectoryName("tau0");
		double lambda_in =  domainLength/h0;
		// eqn.RegisterParameter("lambda", Parameter(&lambda_in, true));
		double u0, phi0, pbterm0;
		eqn.SteadyUniformU(h0,u0,phi0,pbterm0);
		std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		
		// Don't change equation below here, change it in the solver via solver.EquationPtr()!
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]
		

		solver.SetPeriodicBoundaryConditions();
		Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		eqnpntr->solverPtr = &solver;

		TimeStepper *TSpntr = dynamic_cast<TimeStepper *>(solver.TimeStepperPtr());
		// TSpntr->SetCFL(CFL_in);
		// using namespace std::placeholders;
		// solver.LoadInitialConditions("time_d_load_h.txt",0);
		// solver.LoadInitialConditions("time_d_load_hu.txt",1);
		// solver.LoadInitialConditions("time_d_load_hphi.txt",2);
		// solver.LoadInitialConditions("time_d_load_pbh.txt",3);
		solver.SetInitialConditions([this,u0,phi0,pbterm0](double *u, double x, double y)
									{
										u[Eqn::H]=h0*(1+1e-2*(sin(2.0*M_PI*x/domainLength))); //+sin(6.0*M_PI*x/domainLength)+sin(10.0*M_PI*x/domainLength)));
										u[Eqn::HU]=h0*u0;
										u[Eqn::HPHI]=h0*phi0;
										u[Eqn::PBH]=pbterm0;
									});
		// solver.Run(10.0,10); // Integrate to t=100.0, outputting 100 times
		// Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		// eqnpntr->SwitchTheta();
		// eqn.SwitchTheta();
		solver.Run(100.0,200);

	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 
	// double lambda_ = std::stod(argv[1]);
	// double change_t_ = std::stod(argv[2]);
	// std::cout << finalTheta_ << " " << change_t_ << std::endl;

	int npts = 4000;
	{

		ChannelRollWave crw(0.02,50*0.02); //250*0.0129
		// crw.Run(npts, finalTheta_, change_t_);
		crw.Run(npts); // finalTheta_, change_t_
	}

	return 0;
}
