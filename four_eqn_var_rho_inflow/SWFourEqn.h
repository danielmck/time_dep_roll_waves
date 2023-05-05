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
	SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, int nExtras_ = 0) :
		Equation(DIM+1, N_UNSOLVED, nExtras_),  zeroHeightThreshold(1e-7), huThreshold(1e-14)
	{
		SetGTheta(g_, thetaDeg_, tau0_);
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
		  
	void SetGTheta(double g_, double thetaDeg_, double tau0_)
	{
		const double pi = 3.14159265358979323846264338327950288;
		thetaDeg = thetaDeg_;
		tau0 = tau0_;

		theta = thetaDeg*pi/180.0;
		g = g_;
		gcostheta = g_*cos(theta);
		gsintheta = g_*sin(theta);
		tantheta = tan(theta);

		this->RegisterParameter("theta", Parameter(&thetaDeg));
		this->RegisterParameter("g", Parameter(&g));
		this->RegisterParameter("tau0", Parameter(&tau0));
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
	
	double thetaDeg, theta, g, tau0, gcostheta, gsintheta, tantheta;
	MuIvParams pp;
	bool stoppedMaterialHandling;
};

template <unsigned DIM, unsigned N_UNSOLVED=0>
class SWMuIvEqnRhoVaryBase : public SWMuIvEqnBase<DIM>
{
public:
	using SWMuIvEqnBase<DIM>::SWMuIvEqnBase;
	/// Constusing SWMuIvEqn1DFullTest::SWMuIvEqn1DFullTest;ructor, provide g, theta in degrees, and struct with MuIv friction law params

	/// Provide g, theta in degrees.
	SWMuIvEqnRhoVaryBase(double g_, double thetaDeg_, double tau0_, double d_, double alpha_, int nExtras_ = 0) :
		SWMuIvEqnBase<DIM>(g_, thetaDeg_, tau0_, nExtras_)
	{
		SetGTheta(g_, thetaDeg_, tau0_, d_, alpha_);
	}

	SWMuIvEqnRhoVaryBase(double g_, double thetaDeg_, double tau0_, double d_, double alpha_, MuIvParams pp_, int nExtras_ = 0) 
		: SWMuIvEqnRhoVaryBase(g_, thetaDeg_, tau0_, d_, alpha_, nExtras_)
	{
		this->SetMuIvParams(pp_);
	}


	const MuIvParams &GetMuIvParams() const
	{
		return this->pp;
	}
	
	void SteadyUniformU(double h, double &u, double &phi, double &pbterm)
	{
		double max=1e2, min=0, phi_eq, rho_eq, P_eq, Iv_eq, chi_eq;
		while (max-min > 1e-14)
		{
			phi_eq = this->pp.phim/(1+sqrt(max+min));
			rho_eq = phi_eq*this->pp.rhog + (1-phi_eq)*this->pp.rhof;
			P_eq = std::max((rho_eq-this->pp.rhof)/rho_eq,1e-8);
			// std::cout <<  << std::endl;
			((this->MuIv(0.5*(max+min))-this->tantheta/P_eq+this->tau0/(std::max(rho_eq-this->pp.rhof,1e-8)*this->gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		Iv_eq = 0.5*(max+min);
		phi = this->pp.phim/(1+sqrt(Iv_eq));
		rho_eq = phi*this->pp.rhog + (1-phi)*this->pp.rhof;
		chi_eq = (this->pp.rhof+3*rho_eq)/4/rho_eq;
		u = UFromIv(Iv_eq,h,(rho_eq-this->pp.rhof)*this->gcostheta*h);
		double pb = this->pp.rhof*this->gcostheta*h;
		pbterm = h*(pb-rho_eq*this->gcostheta*chi_eq*h);
	}

	void SteadyUniformUTheta(double alt_theta, double h, double &u, double &phi, double &pbterm)
	{
		double max=1e8, min=0, Iv_eq, phi_eq, rho_eq, P_eq, alt_gct = this->g*cos(this->alt_theta*pi/180.0), chi_eq;
		while (max-min > 1e-14)
		{
			phi_eq = this->pp.phim/(1+sqrt(max+min));
			rho_eq = phi_eq*this->pp.rhof + (1-phi_eq)*this->pp.rhof;
			P_eq = (rho_eq-this->pp.rhof)/rho_eq;
			((this->MuIv(0.5*(max+min))-tan(alt_theta*this->pi/180.0)/P_eq+this->tau0/((rho_eq-this->pp.rhof)*this->gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		Iv_eq = 0.5*(max+min);
		phi = this->pp.phim/(1+sqrt(Iv_eq));
		rho_eq = phi*this->pp.rhof + (1-phi)*this->pp.rhof;
		chi_eq = (this->pp.rhof+3*rho_eq)/4/rho_eq;
		u = UFromIv(Iv_eq,h,(rho_eq-this->pp.rhof)*alt_gct*h);
		phi = this->pp.phim/(1+sqrt(Iv));
		double pb = this->pp.rhof*alt_gct*h;
		pbterm = h*(pb-rho_eq*alt_gct*chi_eq*h);
	}
		  
	double UFromIv(double Iv, double h, double ppval)
	{
		return Iv*ppval*h/(3*this->pp.etaf);
	}

	void SetGTheta(double g_, double thetaDeg_, double tau0_, double d_, double alpha_)
	{
		const double pi = 3.14159265358979323846264338327950288;
		this->thetaDeg = thetaDeg_;
		this->tau0 = tau0_;
		this->d = d_;
		this->alpha = alpha_;

		this->theta = this->thetaDeg*pi/180.0;
		this->g = g_;
		this->gcostheta = g_*cos(this->theta);
		this->gsintheta = g_*sin(this->theta);
		this->tantheta = tan(this->theta);

		this->RegisterParameter("d", Parameter(&d));
		this->RegisterParameter("alpha", Parameter(&alpha));
	}

protected:

	// Iv of a quadratic velocity profile with depth average velocity u and thickness h
	double Iv(double u, double h, double ppval)
	{
		return (3*this->pp.etaf)*(u/h)/ppval;
	}

	double d, alpha;
	// MuIvParams pp;
};

////////////////////////////////////////////////////////////////////////////////
/// 1-D SW equations with MuIv friction
class SWMuIvEqn1DRhoVary : public SWMuIvEqnRhoVaryBase<3>
{
public:
	// Inherit constructors from base class
	using SWMuIvEqnRhoVaryBase<3>::SWMuIvEqnRhoVaryBase;
	// using SWMuIvEqnBase<3>::SWMuIvEqnBase;

	enum Variables
	{
		H = 0,
		HU = 1,
		HPHI = 2,
		PBH = 3
	};


	std::string VariableName(int d)
	{
		switch (d)
		{
		case H: return "H";
		case HU: return "HU";
		case HPHI: return "HPHI";
		case PBH: return "PBH";
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
			xFlux[PBH] = 0;
		}	
		else {
			xFlux[HU] = u[HU] * u[HU] / u[H] + 0.5 * this->gcostheta * u[H] * u[H];
			xFlux[HPHI] = u[HPHI]*u[HU]/u[H];
			xFlux[PBH] = u[PBH]*u[HU]/u[H];
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
		double h = u[H], hu = u[HU], hphi = u[HPHI], pbterm = u[PBH];
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
		double chi = (this->pp.rhof+3*rhoBulk)/4/rhoBulk;
		double P = (rhoBulk-this->pp.rhof)/rhoBulk;

		double pb = pbterm/h+rhoBulk*this->gcostheta*chi*h;
		double ppval = rhoBulk*this->gcostheta*h-pb;

		// If u is exactly equal to zero, assume (wrongly...) that the friction is upslope
		if (absu == 0)
			signU = 1;
		else 
			signU = uval/absu;

		double beta = 150*pp.phim*pp.phim*pp.etaf/((1-pp.phim)*(1-pp.phim)*(1-pp.phim)*d*d);

		double D = -2/beta/h*(pb-pp.rhof*this->gcostheta*h);
		double zeta = 3/(2*alpha*h) + this->gcostheta*pp.rhof*P/4;

		// Calculate friction law
		double iv = this->Iv(absu, h, ppval);
		double tanpsi = phi - pp.phim/(1+sqrt(iv));
		if (iv<0) // in case ppval becomes negative
		{
			iv = 1;
			ppval = (rhoBulk-pp.rhof)*this->gcostheta*h;
			pb = rhoBulk*this->gcostheta*h;
			tanpsi = 0.01;
			// D = D*10;
		}
		else
		{
			tanpsi = phi - pp.phim/(1+sqrt(iv));
		}
		if (pb<0) // in case pb becomes negative
		{
			pb = pp.rhof*this->gcostheta*h/2;
			ppval = rhoBulk*this->gcostheta*h-pb;
			tanpsi = -0.01;
		}

		// mu * (rho-rho_f)/rho

		double mubf = ppval*this->MuIv(iv);

		double absFriction = (1/rhoBulk)*(-mubf - tau0 + (rhoBulk-pp.rhof)*D*uval);

		double dilatancy = 3/alpha/h*absu*tanpsi;

		double psi1 = D*P;
		double psi5 = zeta*D - dilatancy;
		stvect[H] += psi1;
		stvect[HPHI] += -D*phi*pp.rhof/rhoBulk;

		stvect[PBH] += (psi5-this->gcostheta*pp.rhof/4*psi1)*h+(pb-rhoBulk*this->gcostheta*chi*h)*psi1;

		// std::cout << absFriction << std::endl;
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

	}
};


#endif
