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
unsigned int const  n_uns = 0;

void InflowFunction(double *u, double *extras, double x, double y, double t)
{
	u[SWMuIvEqn1DRhoVary<n_uns>::H] = gh0*(1+0.001*sin(t*5.0));
	u[SWMuIvEqn1DRhoVary<n_uns>::HU] =  gh0*gu0;
	u[SWMuIvEqn1DRhoVary<n_uns>::HPHI]= gh0*gphi0;
	u[SWMuIvEqn1DRhoVary<n_uns>::PBH]= gpbterm0;
}

class ChannelRollWaveInflow
{
	typedef SWMuIvEqn1DRhoVary<n_uns> Eqn;
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
		Eqn eqn(9.81, 10.0, 0.0, 1e-4, 1e-5);
		eqn.SetMuIvParams(BoyerRockWater);
		eqn.EnableStoppedMaterialHandling();
		eqn.EnableInDirectoryName("theta");
		eqn.EnableInDirectoryName("tau0");
		double u0, phi0, pbterm0, h0frac = 1.0, uin, phiin, pbtermin;
		eqn.SteadyUniformU(h0*h0frac,u0,phi0,pbterm0);
		// eqn.SteadyUniformUTheta(5.0,h0frac*h0,uin,phiin,pbtermin);
		// double Fr_in = 3.0;
		// eqn.RegisterParameter("Fr", Parameter(&Fr_in, true));
		// std::cout << "u0=" << u0 << ", Fr0=" << eqn.SteadyUniformFr(h0,u0) << " , phi0=" << phi0 << std::endl;
		// std::cout << "uin=" << uin << ", Frin=" << eqn.SteadyUniformFr(h0frac*h0,uin) << " , phiin=" << phiin << std::endl;
		Solver solver(n, eqn, new TIMESTEPPER());
		Eqn *eqnpntr = dynamic_cast<Eqn *>(solver.EquationPtr());
		eqnpntr->solverPtr = &solver;
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
		// eqn.SteadyUniformP(gh0,Fr_in, gu0, gphi0, gpbterm0);
		eqn.SteadyUniformUTheta(10.0,gh0,gu0,gphi0,gpbterm0);
		// std::cout << "gu0=" << gu0 << ", gFr0=" << eqn.SteadyUniformFr(gh0,gu0) << " , gphi0=" << gphi0 << std::endl;
		
		solver.SetBoundaryConditionFunctions(HyperbolicSolver::West, InflowFunction);
											 
		
		using namespace std::placeholders;
		// solver.LoadInitialConditions("time_d_load_h.txt",0);
		// solver.LoadInitialConditions("time_d_load_hu.txt",1);
		// solver.LoadInitialConditions("time_d_load_hphi.txt",2);
		// solver.LoadInitialConditions("time_d_load_pbh.txt",3);
		solver.SetInitialConditions([this,h0frac,u0,phi0,pbterm0](double *u, double x, double y)
									{
										u[Eqn::H]=h0*h0frac;
										u[Eqn::HU]=h0*h0frac*u0;
										u[Eqn::HPHI]=h0*h0frac*phi0;
										u[Eqn::PBH]=pbterm0;
										// u[Eqn::THETA]=10;
										// double con_len = 30.0, ch_len = 30.0;
										// if (x<con_len) {
										// 	u[Eqn::THETA] = 10.0;
										// }
										// if (x<con_len+ch_len) {
										// 	u[Eqn::THETA] = 10.0-(10.0-7.0)*(x-con_len)/(ch_len);
										// }
										// else {
										// 	u[Eqn::THETA] = 7.0;
										// }
									});
		// solver.LoadInitialConditions("theta_profile.txt",4);
		solver.Run(100.0,100); // Integrate to t=1.0, outputting 100 times
	}
private:
	double u0, h0, domainLength;
};

int main(int argc, char *argv[])
{
	feenableexcept( FE_INVALID | FE_DIVBYZERO); 

	int npts = 5000;
	{
		ChannelRollWaveInflow crw(0.0041,0.0041*1000);
		crw.Run(npts);
	}

	return 0;
}

