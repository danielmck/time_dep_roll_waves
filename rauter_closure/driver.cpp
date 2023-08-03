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
	typedef SWMuIvEqn1DRauter Eqn;
	// typedef SWMuIvEqn1DViscous VisEqn;
	typedef LimiterMinMod2 LIMITER; 
	typedef RK2TimeStepper TIMESTEPPER;
	typedef SDKTSolverSWPP<Eqn, LIMITER> Solver;

public:
	ChannelRollWave(double h0_, double domainLength_) 
	{
		h0 = h0_;
		domainLength = domainLength_;
	}

	void Run(int n)
	{
		Eqn eqn(9.81, 12, 0, 1e-4, 130,0.53,0.63);
		eqn.SetMuIvParams(BoyerRockWater);
		eqn.EnableStoppedMaterialHandling();
		std::string limiterName = LIMITER::Name();
		// eqn.RegisterParameter("limiter", Parameter(&limiterName, true));
		eqn.EnableInDirectoryName("theta");
		// eqn.EnableInDirectoryName("initTheta");
		// eqn.EnableInDirectoryName("finalTheta");
		eqn.EnableInDirectoryName("tau0");
		double u0, phi0;
		eqn.SteadyUniformU(h0,u0,phi0);
		std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]

		solver.SetPeriodicBoundaryConditions();
		
		// using namespace std::placeholders;
		// solver.LoadInitialConditions("time_d_load_h.txt",0);
		// solver.LoadInitialConditions("time_d_load_hu.txt",1);
		// solver.LoadInitialConditions("time_d_load_hphi.txt",2);
		// solver.LoadInitialConditions("time_d_load_pbh.txt",3);
		solver.SetInitialConditions([this,u0,phi0](double *u, double x, double y)
									{
										u[Eqn::H]=h0*(1+1e-2*(sin(2.0*M_PI*x/domainLength))); //+sin(6.0*M_PI*x/domainLength)+sin(10.0*M_PI*x/domainLength)));
										u[Eqn::HU]=h0*u0;
										u[Eqn::HPHI]=h0*phi0;
									});
		// solver.Run(10.0,10); // Integrate to t=100.0, outputting 100 times
		// Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		// eqnpntr->SwitchTheta();
		// eqn.SwitchTheta();
		solver.Run(50.0,100);
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 

	int npts = 3000;
	{
		ChannelRollWave crw(0.0071,0.0071*12);
		crw.Run(npts);
	}

	return 0;
}
