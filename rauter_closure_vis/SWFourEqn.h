#ifndef SWPOULIQUEN_EQN_H
#define SWPOULIQUEN_EQN_H

/*
  1-D Shallow water equations, MuIv drag

  dh/dt + d/dx (hu) = 0

  d/dt (hu) + d/dx (hu^2 + (g cos(theta)/2) h^2) = h g sin(theta) - mu g h cos(theta) u/|u|

*/


struct MuIvParams
{
	MuIvParams()
	{
		mu1 = mu2 = I0 = phim = phi = rhog = rhof = etaf = 0.0;
	}

	MuIvParams(double mu1_, double mu2_, double I0_, double phim_, double phi_, double rhog_, double rhof_, double etaf_)
	{
		mu1 = mu1_;
		mu2 = mu2_;
		I0 = I0_;
		phi = phi_;
		phim = phim_;
		rhog = rhog_;
		rhof = rhof_;
		etaf = etaf_;
	}

	double mu1, mu2, I0, phim, phi, rhog, rhof, etaf;
};

// Boyer et al mu parameters, with grain densities corresponding to rock and water, roughly, and water viscosity
MuIvParams BoyerRockWater = MuIvParams(0.32, 0.7, 0.005, 0.585, 0.585, 2500, 1000, 1.0016e-3);

template <unsigned DIM, unsigned N_UNSOLVED=0>
class SWMuIvEqnBase : public Equation
{
public:
	/// Constructor, provide g, theta in degrees, and struct with MuIv friction law params
	SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, MuIvParams pp_, int nExtras_ = 0) 
		: SWMuIvEqnBase(g_, thetaDeg_, tau0_, nExtras_)
	{
		SetMuIvParams(pp_);
	}

	/// Provide g, theta in degrees.
	SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, double finalTheta_=-1, int nExtras_ = 0) :
		Equation(DIM+1, N_UNSOLVED, nExtras_),  zeroHeightThreshold(1e-7), huThreshold(1e-14)
	{
		SetGTheta(g_, thetaDeg_, tau0_,finalTheta_);
		stoppedMaterialHandling = false;
		this->RegisterParameter("smh", Parameter(&stoppedMaterialHandling));
	}

	void SetMuIvParams(MuIvParams pp_)
	{
		pp = pp_;
		this->RegisterParameter("mu1", Parameter(&pp.mu1));
		this->RegisterParameter("mu2", Parameter(&pp.mu2));
		this->RegisterParameter("I0", Parameter(&pp.I0));
		this->RegisterParameter("phi", Parameter(&pp.phi));
		this->RegisterParameter("phim", Parameter(&pp.phim));
		this->RegisterParameter("rhog", Parameter(&pp.rhog));
		this->RegisterParameter("rhof", Parameter(&pp.rhof));
		this->RegisterParameter("etaf", Parameter(&pp.etaf));

		P = (rhoBulk-pp.rhof)/rhoBulk;
	}

	const MuIvParams &GetMuIvParams() const
	{
		return pp;
	}
	
	double SteadyUniformU(double h)
	{
		double max=1e8, min=0;
		while (max-min > 1e-14)
		{
			((MuIv(Iv(0.5*(max+min),h))-tantheta/P+tau0/((rhoBulk-pp.rhof)*gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		return 0.5*(max+min);
	}
	double SteadyUniformFr(double h)
	{
		return SteadyUniformU(h)/sqrt(this->gcostheta*h);
	}
	double SteadyUniformFr(double h, double u)
	{
		return std::abs(u)/sqrt(this->gcostheta*h);
	}
		  
	void SetGTheta(double g_, double thetaDeg_, double tau0_, double finalThetaDeg_ = -1.0)
	{
		const double pi = 3.14159265358979323846264338327950288;
		thetaDeg = thetaDeg_;
		tau0 = tau0_;
		initThetaDeg = thetaDeg;
		finalThetaDeg = finalThetaDeg_;
		theta = thetaDeg*pi/180.0;
		
		g = g_;
		gcostheta = g_*cos(theta);
		gsintheta = g_*sin(theta);
		tantheta = tan(theta);

		this->RegisterParameter("theta", Parameter(&thetaDeg));
		this->RegisterParameter("g", Parameter(&g));
		this->RegisterParameter("tau0", Parameter(&tau0));
		if (finalThetaDeg >= 0.0)
		{
			this->RegisterParameter("initTheta", Parameter(&initThetaDeg));
			this->RegisterParameter("finalTheta", Parameter(&finalThetaDeg));
		}
	}

	void SwitchTheta()
	{
		const double pi = 3.14159265358979323846264338327950288;

		thetaDeg = finalThetaDeg;
		theta = thetaDeg*pi/180.0;
		gcostheta = g*cos(theta);
		gsintheta = g*sin(theta);
		tantheta = tan(theta);
	}

	void EnableStoppedMaterialHandling()
	{
		stoppedMaterialHandling = true;
	}
	void DisableStoppedMaterialHandling()
	{
		stoppedMaterialHandling = false;
	}

	bool GetStoppedMaterialHandling()
	{
		return stoppedMaterialHandling;
	}

	double ZeroHeightThreshold() {return zeroHeightThreshold;}

	

protected:
	double MuIv(double Iv)
	{
		if (Iv<1e-10)
		{
			return pp.mu1;
		}
		else
		{
			return pp.mu1 + (pp.mu2 - pp.mu1)/(1+pp.I0/Iv) + Iv + 2.5*pp.phim*sqrt(Iv);
		}
	}

	// Iv of a quadratic velocity profile with depth average velocity u and thickness h
	double Iv(double u, double h)
	{
		return (3*pp.etaf)/(pp.phi*(pp.rhog-pp.rhof)*this->gcostheta)*(u/(h*h));
	}

	const double  zeroHeightThreshold, huThreshold;
	double P, rhoBulk;
	
	double thetaDeg, initThetaDeg, finalThetaDeg, theta, g, tau0, gcostheta, gsintheta, tantheta;
	MuIvParams pp;
	bool stoppedMaterialHandling;
};

template <unsigned DIM, unsigned N_UNSOLVED=0>
class SWMuIvEqnRauterBase : public SWMuIvEqnBase<DIM>
{
public:
	using SWMuIvEqnBase<DIM>::SWMuIvEqnBase;
	/// Constusing SWMuIvEqn1DFullTest::SWMuIvEqn1DFullTest;ructor, provide g, theta in degrees, and struct with MuIv friction law params

	SWMuIvEqnRauterBase(double g_, double thetaDeg_, double tau0_, double d_, double a_, double phi_rlp, double phi_rcp, double finalThetaDeg_ = -1.0, int nExtras_ = 0) :
		SWMuIvEqnBase<DIM>(g_, thetaDeg_, tau0_, finalThetaDeg_, nExtras_)
	{
		SetGTheta(g_, thetaDeg_, tau0_, d_, a_, phi_rlp, phi_rcp);
	}

	SWMuIvEqnRauterBase(double g_, double thetaDeg_, double tau0_, double d_, double a_, double phi_rlp, double phi_rcp, MuIvParams pp_, int nExtras_ = 0) 
		: SWMuIvEqnRauterBase(g_, thetaDeg_, tau0_, a, phi_rlp, phi_rcp, nExtras_)
	{
		this->SetMuIvParams(pp_);
	}

	const MuIvParams &GetMuIvParams() const
	{
		return this->pp;
	}
	
	void SteadyUniformU(double h, double &u, double &phi)
	{
		double max=1, min=1e-14, phi_eq, rho_eq, pp_eq, P_eq, Iv_eq, Iv_phi, pp_con, pp_shear;
		while (max-min > 1e-14)
		{
			phi_eq = (max+min)/2.0;
			rho_eq = phi_eq*this->pp.rhog + (1.0-phi_eq)*this->pp.rhof;
			pp_eq = (rho_eq-this->pp.rhof)*this->gcostheta*h;
			if (phi_eq > this->pp.phim)
				Iv_eq = 0;
			else {
				pp_con = a*(phi_eq-phi_rlp)/(phi_rcp-phi_eq);
				pp_shear = pp_eq-pp_con;
				Iv_phi = (this->pp.phim/phi_eq-1);
				Iv_eq = std::pow(Iv_phi,2)*pp_shear/pp_eq;
			}
			P_eq = std::max((rho_eq-this->pp.rhof)/rho_eq,1e-8);
			// std::cout <<  << std::endl;
			((this->MuIv(Iv_eq)-this->tantheta/P_eq+this->tau0/(std::max(rho_eq-this->pp.rhof,1e-8)*this->gcostheta*h)>0)?min:max)=0.5*(max+min);
		}
		phi = 0.5*(max+min);

		rho_eq = phi*this->pp.rhog + (1-phi)*this->pp.rhof;
		pp_eq = (rho_eq-this->pp.rhof)*this->gcostheta*h;
		pp_con = a*(phi_eq-phi_rlp)/(phi_rcp-phi_eq);
		pp_shear = pp_eq-pp_con;
		Iv_phi = (this->pp.phim/phi_eq-1);
		Iv_eq = std::pow(Iv_phi,2)*pp_shear/pp_eq;
		u = UFromIv(Iv_eq,h,(rho_eq-this->pp.rhof)*this->gcostheta*h);
	}

	void SteadyUniformUTheta(double alt_theta, double h, double &u, double &phi, double &pbterm)
	{
		double max=1e2, min=0, Iv_eq, phi_eq, rho_eq, P_eq, alt_gct = this->g*cos(this->alt_theta*pi/180.0), chi_eq;
		while (max-min > 1e-14)
		{
			phi_eq = this->pp.phim/(1+sqrt((max+min)/2.0));
			rho_eq = phi_eq*this->pp.rhog + (1-phi_eq)*this->pp.rhof;
			P_eq = (rho_eq-this->pp.rhof)/rho_eq;
			((this->MuIv(0.5*(max+min))-tan(alt_theta*this->pi/180.0)/P_eq+this->tau0/((rho_eq-this->pp.rhof)*this->gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		Iv_eq = 0.5*(max+min);
		phi = this->pp.phim/(1+sqrt(Iv_eq));
		rho_eq = phi*this->pp.rhof + (1-phi)*this->pp.rhof;
		chi_eq = (this->pp.rhof+3.0*rho_eq)/4.0/rho_eq;
		u = UFromIv(Iv_eq,h,(rho_eq-this->pp.rhof)*alt_gct*h);
		double pb = this->pp.rhof*alt_gct*h;
		pbterm = h*(pb-rho_eq*alt_gct*chi_eq*h);
	}
		  
	double UFromIv(double Iv, double h, double ppval)
	{
		return Iv*ppval*h/(3.0*this->pp.etaf);
	}

	void SetGTheta(double g_, double thetaDeg_, double tau0_, double d_, double a_, double phi_rlp_, double phi_rcp_)
	{
		const double pi = 3.14159265358979323846264338327950288;
		this->thetaDeg = thetaDeg_;
		this->tau0 = tau0_;
		this->d = d_;
		this->a = a_;
		this->phi_rlp = phi_rlp_;
		this->phi_rcp = phi_rcp_;

		this->theta = this->thetaDeg*pi/180.0;
		this->g = g_;
		this->gcostheta = g_*cos(this->theta);
		this->gsintheta = g_*sin(this->theta);
		this->tantheta = tan(this->theta);

		this->RegisterParameter("d", Parameter(&d));
		this->RegisterParameter("a", Parameter(&a));
		this->RegisterParameter("phi_rlp", Parameter(&phi_rlp));
		this->RegisterParameter("phi_rcp", Parameter(&phi_rcp));
	}

protected:

	// Iv of a quadratic velocity profile with depth average velocity u and thickness h
	double Iv(double u, double h, double ppval)
	{
		if (ppval < 1e-8)
			return 1e8;
		else
			return (3.0*this->pp.etaf)*(u/h)/ppval;
	}

	double d, a, phi_rlp, phi_rcp;
	// MuIvParams pp;
};

////////////////////////////////////////////////////////////////////////////////
/// 1-D SW equations with MuIv friction
class SWMuIvEqn1DRauter : public SWMuIvEqnRauterBase<3>
{
public:
	// Inherit constructors from base class
	using SWMuIvEqnRauterBase<3>::SWMuIvEqnRauterBase;
	// using SWMuIvEqnBase<3>::SWMuIvEqnBase;

	enum Variables
	{
		H = 0,
		HU = 1,
		HPHI = 2,
	};


	std::string VariableName(int d)
	{
		switch (d)
		{
		case H: return "H";
		case HU: return "HU";
		case HPHI: return "HPHI";
		default: return Equation::VariableName(d);
		}
	}

	void XConvectionFlux(double *xFlux, const double *u, const double *extras, const double *dextrasdt)
	{
		xFlux[H] = u[HU];
		// xFlux[HPHI] = u[HU];
		// xFlux[PBH] = u[HU];

		if (u[H] < this->zeroHeightThreshold) {
			xFlux[HU] = 0;
			xFlux[HPHI] = 0;
		}	
		else {
			xFlux[HU] = u[HU] * u[HU] / u[H] + 0.5 * this->gcostheta * u[H] * u[H];
			xFlux[HPHI] = u[HPHI]*u[HU]/u[H];
		}
	}

	double XWaveSpeeds(const double *u, const double *extras, const double *dextrasdt)
	{
		// No wave speed for zero height, but none either for static material
		// This stops (or reduces to machine precision) the evolution of static layers
		// that otherwise occurs even when hu = 0
		if (u[H] < this->zeroHeightThreshold || std::abs(u[HU]) < this->huThreshold)
			return 0;
		else
			return std::abs(u[HU] / u[H]) + sqrt(this->gcostheta*u[H]);
	}

	void SourceTerms(double dt, double *stvect, double *u, const double *extras, const double *dedt)
	{
		double h = u[H], hu = u[HU], hphi = u[HPHI];
		if (h < this->zeroHeightThreshold)
		{
			h = this->zeroHeightThreshold;
			hu = 0;
			hphi = h*pp.phim;
		}

		double uval = hu/h;
		double absu = std::abs(uval);
		double signU;

		double phi = hphi/h;

		double rhoBulk = pp.rhog*phi + pp.rhof*(1-phi);
		double P = (rhoBulk-this->pp.rhof)/rhoBulk;

		double phi_diff = (pp.phim/phi-1);
		double pp_shear = 3*pp.etaf*absu/h/std::pow(phi_diff,2);
		double pp_contact = std::max(a*(phi-phi_rlp)/(phi_rcp-phi),0.0);
		double ppval = pp_contact + pp_shear;
		double pb = rhoBulk*this->gcostheta*h-ppval;

		// If u is exactly equal to zero, assume (wrongly...) that the friction is upslope
		if (absu == 0)
			signU = 1;
		else 
			signU = uval/absu;

		double beta = 150.0*pp.phim*pp.phim*pp.etaf/((1.0-pp.phim)*(1.0-pp.phim)*(1.0-pp.phim)*d*d);

		double D = -2.0/beta/h*(pb-pp.rhof*this->gcostheta*h);

		// Calculate friction law
		double iv = this->Iv(absu, h, ppval);

		// mu * (rho-rho_f)/rho
		double mubf = ppval*this->MuIv(iv);
		if (ppval < 1e-8)
			mubf = 3.0*absu*pp.etaf/h;
		

		double absFriction = (1.0/rhoBulk)*(-mubf - tau0 + (rhoBulk-pp.rhof)*D*absu);

		double psi1 = D*P;
		stvect[H] += psi1;

		stvect[HPHI] += -D*phi*pp.rhof/rhoBulk;

		if (this->stoppedMaterialHandling && dt != -1)
		{
			// Add gravity terms so stvect is now gravity + pressure
			stvect[HU] += h*this->gsintheta;

			// If we are moving downstream
			if (u[HU] > this->huThreshold)
			{
				if (u[HU] + dt*(stvect[HU] + absFriction*signU) < 0)
					stvect[HU] = -u[HU]/dt;
				else
					stvect[HU] += absFriction*signU;
			}
			else if (u[HU] >= -this->huThreshold) // otherwise if static 
			{
				// apply at most enough friction to keep static, either upslope or downslope
				if (stvect[HU]>0)
				{
					stvect[HU] += std::max(absFriction, -stvect[HU]);
	 			}
				else if (stvect[HU]<0)
				{
					stvect[HU] += std::min(-absFriction, -stvect[HU]);
				}
				else
				{
					// If static and stdvect[HU] == 0 we don't apply any friction
				}
			}
			else // we are moving upstream
			{
				if (u[HU] + dt*(stvect[HU] + absFriction*signU) > 0)
					stvect[HU] = -u[HU]/dt;
				else
					stvect[HU] += absFriction*signU;
			}
		}
		else // No stopped material handling
		{
			stvect[HU] += h*this->gsintheta + absFriction*signU;
		}

		// if ((h<0.1) & (absu>2.05) & (D>0) & (h*this->gsintheta+absFriction>0)) {
		// 	double test = 1;
		// }
	}
};

class SWMuIvEqn1DRauterVis : public SWMuIvEqn1DRauter
{
public:
	using SWMuIvEqn1DRauter::SWMuIvEqn1DRauter;
	// using SWMuIvEqn1D<TRANSITION>::pp;
	// using SWMuIvEqn1D<TRANSITION>::gsintheta;
	// using SWMuIvEqn1D<TRANSITION>::gcostheta;
	// using SWMuIvEqn1D<TRANSITION>::tantheta;

	// Redefine this here, which is terrible, but otherwise two-phase lookup prevents
	// access of variable names from this class without explicit this->
	enum Variables
	{
		H = 0,
		HU = 1,
		HPHI = 2,
	};

	// Turn on diffusion terms
	constexpr const bool HasDiffusionTerms() const
	{
		return true;
	}

	void XDiffusionFlux(double *xFlux, const double *u, const double *dudx, const double *dudy)
	{
		xFlux[H] = 0;
		if (u[H] < this->zeroHeightThreshold)
		{
		xFlux[HU] = 0;
		}
		else
		{
		xFlux[HU] = this->Nu() * pow(u[H],1)*(dudx[HU]-u[HU]*dudx[H]/u[H]);
		}
		xFlux[HPHI] = 0;
	}

	constexpr const double Nu()
	{
		return 5e-3;
	}

};
#endif
