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
	// typedef SWMuIvEqn1DRhoVary Eqn;
	typedef SWMuIvEqnFullViscousRhoVary<0> VisEqn;
	typedef LimiterWENO LIMITER;
	typedef RK2TimeStepper TIMESTEPPER;
	typedef SDKTSolver<VisEqn, LIMITER> Solver;

public:
	ChannelRollWave(double h0_, double domainLength_) 
	{
		h0 = h0_;
		domainLength = domainLength_;
	}

	void Run(int n)
	{
		VisEqn eqn(9.81, 10, 0, 1e-4, 1e-5);
		eqn.SetMuIvParams(BoyerRockWater);
		// eqn.EnableStoppedMaterialHandling();
		eqn.EnableInDirectoryName("theta");
		eqn.EnableInDirectoryName("tau0");
		double u0, phi0, pbterm0;
		eqn.SteadyUniformU(h0,u0,phi0,pbterm0);
		std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]

		solver.SetPeriodicBoundaryConditions();
		
		using namespace std::placeholders;
		solver.SetInitialConditions([this,u0,phi0,pbterm0](double *u, double x, double y)
									{
										u[VisEqn::H]=h0*(1+1e-2*(sin(2.0*M_PI*x/domainLength))); //+sin(6.0*M_PI*x/domainLength)+sin(10.0*M_PI*x/domainLength)));
										u[VisEqn::HU]=h0*u0;
										u[VisEqn::HPHI]=h0*phi0;
										u[VisEqn::PBH]=pbterm0;
									});
		solver.Run(50.0,100); // Integrate to t=100.0, outputting 100 times
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 

	int npts = 2000;
	{
		ChannelRollWave crw(0.0188,25*0.0188);
		crw.Run(npts);
	}

	return 0;
}
