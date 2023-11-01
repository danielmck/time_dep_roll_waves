#include "SDKTSolver.h"
#include "SWMuIvEqn.h"
#include "RK2TimeStepper.h"

#include <cmath>
#include <fenv.h>
#include <sys/time.h>
#include <iomanip>
#include <vector>
#include <sstream>
#include <functional>


// class ChannelRollWave
// {
// 	typedef SWMuIvEqn1D Eqn;
// 	typedef LimiterWENO LIMITER;
// 	typedef RK2TimeStepper TIMESTEPPER;
// 	typedef SDKTSolver<Eqn, LIMITER> Solver;

// public:
	// ChannelRollWave(double h0_, double domainLength_) 
	// {
	// 	h0 = h0_;
	// 	domainLength = domainLength_;
	// }

	// void Run(int n)
	// {
	// 	Eqn eqn(9.81, 8);
	// 	eqn.SetMuIvParams(BoyerRockWater);
	// 	eqn.EnableInDirectoryName("theta");
	// 	double u0 = eqn.SteadyUniformU(h0);
	// 	std::cout << "u0=" << u0 << ", Fr0="<<eqn.SteadyUniformFr(h0) << std::endl;
	// 	Solver solver(n, eqn, new TIMESTEPPER());
	// 	solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]

	// 	solver.SetPeriodicBoundaryConditions();

	// 	using namespace std::placeholders;
	// 	solver.SetInitialConditions([this,u0](double *u, double x, double y)
	// 								{
	// 									u[Eqn::H]=h0*(1+1e-2*sin(2.0*M_PI*x/domainLength));
	// 									u[Eqn::HU]=h0*0.1;//u0;
	// 								});

	// 	solver.Run(100.0,100); // Integrate to t=100.0, outputting 100 times
	// }
// private:
// 	double u0, h0, domainLength;
// };


double gh0, gu0;
unsigned const int n_uns = 1;

void InflowFunction(double *u, double *extras, double x, double y, double t)
{
	u[SWMuIvEqn1D<n_uns>::H] = gh0*(1+0.001*sin(1*t));
	u[SWMuIvEqn1D<n_uns>::HU] = gh0*gu0;
}

class ChannelRollWaveInflow
{
	typedef SWMuIvEqn1D<n_uns> Eqn;
	typedef LimiterWENO LIMITER;
	typedef RK2TimeStepper TIMESTEPPER;
	typedef SDKTSolver<Eqn, LIMITER> Solver;

public:
	ChannelRollWaveInflow(double h0_, double domainLength_) 
	{
		h0 = h0_;
		domainLength = domainLength_;
	}

	void Run(int n)
	{
		Eqn eqn(9.81, 7.0, 0.0);
		eqn.SetMuIvParams(BoyerRockWater);
		eqn.EnableInDirectoryName("theta");
		double u0 = eqn.SteadyUniformU(h0);
		// std::cout << "u0=" << u0 << ", Fr0="<<eqn.SteadyUniformFr(h0) << std::endl;
		double h0frac = 0.1, uin;
		// eqn.SteadyUniformUTheta(10.0,h0,u0,phi0,pbterm0);
		uin = eqn.SteadyUniformU(h0frac*h0);
		// std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		std::cout << "uin=" << uin << ", Frin=" << eqn.SteadyUniformFr(h0frac*h0) << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		solver.SetDomain(0.0, domainLength); // Domain is x in [0, domainLength]


		solver.SetBoundaryConditionType(HyperbolicSolver::East,
										{Eqn::H, Eqn::HU},
										HyperbolicSolver::GradientExtrapolate0,
										HyperbolicSolver::FluxUseExtrapolated);
		
		solver.SetBoundaryConditionType(HyperbolicSolver::West,
										{Eqn::H, Eqn::HU},
										HyperbolicSolver::GradientExtrapolate0,
										HyperbolicSolver::FluxSetValue);

		gh0 = h0;
		// gu0 = u0;
		gu0 = eqn.SteadyUniformUTheta(10.0,gh0);
		std::cout << "gu0=" << gu0 << ", Fr0="<< eqn.SteadyUniformFrTheta(10.0,h0) << std::endl;
		solver.SetBoundaryConditionFunctions(HyperbolicSolver::West, InflowFunction);
											 
		
		using namespace std::placeholders;
		solver.SetInitialConditions([this,uin,h0frac](double *u, double x, double y)
									{
										u[Eqn::H]=h0frac*h0;
										u[Eqn::HU]=h0frac*h0*uin;
										double con_len = 30.0, ch_len = 30.0;
										if (x<con_len) {
											u[Eqn::THETA] = 10.0;
										}
										if (x<con_len+ch_len) {
											u[Eqn::THETA] = 10.0-(10.0-7.0)*(x-con_len)/(ch_len);
										}
										else {
											u[Eqn::THETA] = 7.0;
										}
									});

		solver.Run(300.0,100); // Integrate to t=1.0, outputting 100 times
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 

	int npts = 10000;
	{
		ChannelRollWaveInflow crw(0.0246,100);
		crw.Run(npts);
	}



	return 0;
}
