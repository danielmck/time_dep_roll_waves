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

double gh0, gu0, gphi0, gpbterm0;

void InflowFunction(double *u, double *extras, double x, double y, double t)
{
	u[SWMuIvEqn1DRhoVary::H] = gh0*(1+0.001*sin(t/5));
	u[SWMuIvEqn1DRhoVary::HU] = gh0*gu0;
	u[SWMuIvEqn1DRhoVary::HPHI]=gh0*gphi0;
	u[SWMuIvEqn1DRhoVary::PBH]=gpbterm0;
}

class ChannelRollWaveInflow
{
	typedef SWMuIvEqn1DRhoVary Eqn;
	typedef LimiterWENO LIMITER;
	typedef RK2TimeStepper TIMESTEPPER;
	typedef SDKTSolverSWPP<Eqn, LIMITER> Solver;

public:
	ChannelRollWaveInflow(double h0_, double domainLength_) 
	{
		h0 = h0_;
		domainLength = domainLength_;
	}

	void Run(int n)
	{
		Eqn eqn(9.81, 10, 40, 1e-4, 1e-5);
		eqn.SetMuIvParams(BoyerRockWater);
		eqn.EnableStoppedMaterialHandling();
		eqn.EnableInDirectoryName("theta");
		eqn.EnableInDirectoryName("tau0");
		double u0, phi0, pbterm0, uin, phiin, pbtermin;
		eqn.SteadyUniformU(h0,u0,phi0,pbterm0);
		// eqn.SteadyUniformUTheta(10,h0,uin,phiin,pbtermin);
		std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]


		solver.SetBoundaryConditionType(HyperbolicSolver::East,
										{Eqn::H, Eqn::HU, Eqn::HPHI, Eqn::PBH},
										HyperbolicSolver::GradientExtrapolate0,
										HyperbolicSolver::FluxUseExtrapolated);
		
		solver.SetBoundaryConditionType(HyperbolicSolver::West,
										{Eqn::H, Eqn::HU, Eqn::HPHI, Eqn::PBH},
										HyperbolicSolver::GradientExtrapolate0,
										HyperbolicSolver::FluxSetValue);

		gh0 = h0;
		gu0 = u0;
		gphi0 = phi0;
		gpbterm0 = pbterm0;
		solver.SetBoundaryConditionFunctions(HyperbolicSolver::West, InflowFunction);
											 
		
		using namespace std::placeholders;
		solver.SetInitialConditions([this,u0,phi0,pbterm0](double *u, double x, double y)
									{
										u[Eqn::H]=h0;
										u[Eqn::HU]=h0*u0;
										u[Eqn::HPHI]=h0*phi0;
										u[Eqn::PBH]=pbterm0;
									});

		solver.Run(300.0,100); // Integrate to t=1.0, outputting 100 times
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 

	int npts = 6000;
	{
		ChannelRollWaveInflow crw(0.1,800);
		crw.Run(npts);
	}

	return 0;
}

