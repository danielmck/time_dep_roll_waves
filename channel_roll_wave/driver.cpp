#include "SDKTSolverSWPP.h"
#include "SWMuIvEqn.h"
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
	typedef SWMuIvEqn1D<0> Eqn;
	typedef LimiterWENO LIMITER;
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
		Eqn eqn(9.81, 12.0, 0.0);
		eqn.SetMuIvParams(BoyerRockWater);
		eqn.EnableStoppedMaterialHandling();
		eqn.EnableInDirectoryName("theta");
		// eqn.EnableInDirectoryName("initTheta");
		// eqn.EnableInDirectoryName("finalTheta");
		// eqn.EnableInDirectoryName("change_t");

		eqn.EnableInDirectoryName("tau0");
		double lambda_in =  domainLength/h0;
		eqn.RegisterParameter("lambda", Parameter(&lambda_in, true));
		double u0 = eqn.SteadyUniformU(h0);
		std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0) << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]

		Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		eqnpntr->solverPtr = &solver;
		solver.SetPeriodicBoundaryConditions();

		using namespace std::placeholders;
		// solver.LoadInitialConditions("time_d_load_h.txt",0);
		// solver.LoadInitialConditions("time_d_load_hu.txt",1);
		solver.SetInitialConditions([this,u0](double *u, double x, double y)
									{
										u[Eqn::H]=h0*(1+1e-2*sin(2.0*M_PI*x/domainLength));
										u[Eqn::HU]=h0*u0;
									});

		// solver.Run(10.0,10); // Integrate to t=100.0, outputting 100 times
		// Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		// eqnpntr->SwitchTheta();
		solver.Run(50.0,200); // Integrate to t=100.0, outputting 100 times
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 
	// double lambda_ = std::stod(argv[1]);
	int npts = 1500;
	{
		ChannelRollWave crw(0.0076,0.0076*15);
		crw.Run(npts);
	}

	return 0;
}
