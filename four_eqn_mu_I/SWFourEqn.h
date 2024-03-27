#ifndef SWPOULIQUEN_EQN_H
#define SWPOULIQUEN_EQN_H

/*
  1-D Shallow water equations, MuIv drag

  dh/dt + d/dx (hu) = 0

  d/dt (hu) + d/dx (hu^2 + (g cos(theta)/2) h^2) = h g sin(theta) - mu g h cos(theta) u/|u|

*/


struct MuIParams
{
	MuIParams()
	{
		mu1 = mu2 = I0 = phim = phi = delphi = rhog = rhof = etaf = 0.0;
	}

	MuIParams(double mu1_, double mu2_, double I0_, double phim_, double phi_, double delphi_, double rhog_, double rhof_, double etaf_)
	{
		mu1 = mu1_;
		mu2 = mu2_;
		I0 = I0_;
		phi = phi_;
		phim = phim_;
		delphi = delphi_;
		rhog = rhog_;
		rhof = rhof_;
		etaf = etaf_;
	}

	double mu1, mu2, I0, phim, phi, rhog, rhof, etaf, delphi;
};

// Boyer et al mu parameters, with grain densities corresponding to rock and water, roughly, and water viscosity
MuIParams BoyerRockWater = MuIParams(0.342, 0.557, 0.069, 0.585, 0.585, 0.1, 2500, 1000, 1.0016e-3);

template <unsigned DIM, unsigned N_UNSOLVED>
class SWMuIEqnBase : public Equation
{
public:
	/// Constructor, provide g, theta in degrees, and struct with MuIv friction law params
	SWMuIEqnBase(double g_, double thetaDeg_, double tau0_, double d_, MuIParams pp_, int nExtras_ = 0) 
		: SWMuIEqnBase(g_, thetaDeg_, tau0_, d_, nExtras_)
	{
		SetMuIParams(pp_);
		solverPtr = nullptr;
	}

	/// Provide g, theta in degrees.
	SWMuIEqnBase(double g_, double thetaDeg_, double tau0_, double d_, double finalTheta_=-1, int nExtras_ = 0) :
		Equation(DIM+1, N_UNSOLVED, nExtras_),  zeroHeightThreshold(1e-7), huThreshold(1e-14)
	{
		SetGTheta(g_, thetaDeg_, tau0_, d_, finalTheta_);
		stoppedMaterialHandling = false;
		this->RegisterParameter("smh", Parameter(&stoppedMaterialHandling));
		solverPtr = nullptr;
	}

	SWMuIEqnBase(double g_, double thetaDeg_, double tau0_, double d_, double finalTheta_=-1, double change_t_=-1, int nExtras_ = 0) :
		Equation(DIM+1, N_UNSOLVED, nExtras_),  zeroHeightThreshold(1e-7), huThreshold(1e-14)
	{
		SetGTheta(g_, thetaDeg_, tau0_,d_,finalTheta_,change_t_);
		stoppedMaterialHandling = false;
		this->RegisterParameter("smh", Parameter(&stoppedMaterialHandling));
		solverPtr = nullptr;
	}

	void SetMuIParams(MuIParams pp_)
	{
		pp = pp_;
		this->RegisterParameter("mu1", Parameter(&pp.mu1));
		this->RegisterParameter("mu2", Parameter(&pp.mu2));
		this->RegisterParameter("I0", Parameter(&pp.I0));
		this->RegisterParameter("phi", Parameter(&pp.phi));
		this->RegisterParameter("phim", Parameter(&pp.phim));
		this->RegisterParameter("delphi", Parameter(&pp.delphi));
		this->RegisterParameter("rhog", Parameter(&pp.rhog));
		this->RegisterParameter("rhof", Parameter(&pp.rhof));
		this->RegisterParameter("etaf", Parameter(&pp.etaf));

		P = (rhoBulk-pp.rhof)/rhoBulk;
	}

	const MuIParams &GetMuIParams() const
	{
		return pp;
	}
	
	double SteadyUniformU(double h)
	{
		double max=1e8, min=0;
		while (max-min > 1e-14)
		{
			((MuI(I(0.5*(max+min),h))-tantheta/P+tau0/((rhoBulk-pp.rhof)*gcostheta*h)>0)?max:min)=0.5*(max+min);
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
		  
	void SetGTheta(double g_, double thetaDeg_, double tau0_, double d_,double finalThetaDeg_ = -1.0, double change_t_ = -1.0)
	{
		const double pi = 3.14159265358979323846264338327950288;
		thetaDeg = thetaDeg_;
		tau0 = tau0_;
		initThetaDeg = thetaDeg;
		if (finalThetaDeg_>0){
			finalThetaDeg = finalThetaDeg_;
		}
		else {
			finalThetaDeg = thetaDeg;
		}
		theta = thetaDeg*pi/180.0;
		change_t = change_t_;
		
		g = g_;
		gcostheta = g_*cos(theta);
		gsintheta = g_*sin(theta);
		tantheta = tan(theta);
		d=d_;

		this->RegisterParameter("theta", Parameter(&thetaDeg));
		this->RegisterParameter("g", Parameter(&g));
		this->RegisterParameter("tau0", Parameter(&tau0));
		if (finalThetaDeg_ >= 0.0)
		{
			this->RegisterParameter("initTheta", Parameter(&initThetaDeg));
			this->RegisterParameter("finalTheta", Parameter(&finalThetaDeg));
		}
		if (change_t >= 0.0)
		{
			this->RegisterParameter("change_t", Parameter(&change_t));
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

	void AlterTheta()
	{
		double t = solverPtr->Time();
		if (t>change_t){
			thetaDeg = finalThetaDeg;
		}
		else
		{
			thetaDeg = initThetaDeg + (finalThetaDeg-initThetaDeg)*t/change_t;
		}
		theta = thetaDeg*pi/180.0;
		gcostheta = g*cos(theta);
		gsintheta = g*sin(theta);
		tantheta = tan(theta);
	}

	HyperbolicSolver *solverPtr;

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
	double MuI(double I)
	{
		if (I<1e-10)
		{
			return pp.mu1;
		}
		else
		{
			return pp.mu1 + (pp.mu2 - pp.mu1)/(1+pp.I0/I);
		}
	}

	// Iv of a quadratic velocity profile with depth average velocity u and thickness h
	double I(double u, double h)
	{
		return (2.0*d)/sqrt((pp.phi*(pp.rhog-pp.rhof)*this->gcostheta))*5.0/2.0*sqrt(pp.rhog)*(u/(h*sqrt(h)));
	}

	const double  zeroHeightThreshold, huThreshold;
	double P, rhoBulk;
	
	double thetaDeg, initThetaDeg, finalThetaDeg, theta, g, tau0, gcostheta, gsintheta, tantheta, change_t, d;
	MuIParams pp;
	bool stoppedMaterialHandling;
};

template <unsigned DIM, unsigned N_UNSOLVED>
class SWMuIEqnRhoVaryBase : public SWMuIEqnBase<DIM, N_UNSOLVED>
{
public:
	using SWMuIEqnBase<DIM, N_UNSOLVED>::SWMuIEqnBase;
	/// Constusing SWMuIvEqn1DFullTest::SWMuIvEqn1DFullTest;ructor, provide g, theta in degrees, and struct with MuIv friction law params

	SWMuIEqnRhoVaryBase(double g_, double thetaDeg_, double tau0_, double d_, double alpha_, double finalThetaDeg_ = -1.0, double change_t_ = -1.0, int nExtras_ = 0) :
		SWMuIEqnBase<DIM, N_UNSOLVED>(g_, thetaDeg_, tau0_, finalThetaDeg_, change_t_, nExtras_)
	{
		SetGTheta(g_, thetaDeg_, tau0_, d_, alpha_, finalThetaDeg_, change_t_);
	}

	SWMuIEqnRhoVaryBase(double g_, double thetaDeg_, double tau0_, double d_, double alpha_, MuIParams pp_, int nExtras_ = 0) 
		: SWMuIEqnRhoVaryBase(g_, thetaDeg_, tau0_, d_, alpha_, nExtras_)
	{
		this->SetMuIParams(pp_);
	}


	const MuIParams &GetMuIParams() const
	{
		return this->pp;
	}
	
	void SteadyUniformU(double h, double &u, double &phi, double &pbterm)
	{
		double max=1e-2, min=0, phi_eq, rho_eq, P_eq, I_eq, chi_eq;
		while (max-min > 1e-14)
		{
			phi_eq = this->pp.phim-this->pp.delphi*((max+min)/2.0);
			rho_eq = phi_eq*this->pp.rhog + (1.0-phi_eq)*this->pp.rhof;
			P_eq = std::max((rho_eq-this->pp.rhof)/rho_eq,1e-8);
			// std::cout <<  << std::endl;
			((this->MuI(0.5*(max+min))-this->tantheta/P_eq+this->tau0/(std::max(rho_eq-this->pp.rhof,1e-8)*this->gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		I_eq = 0.5*(max+min);
		phi = this->pp.phim-this->pp.delphi*(I_eq);
		rho_eq = phi*this->pp.rhog + (1-phi)*this->pp.rhof;
		chi_eq = (this->pp.rhof+3.0*rho_eq)/4.0/rho_eq;
		u = UFromI(I_eq,h,(rho_eq-this->pp.rhof)*this->gcostheta*h);
		double pb = this->pp.rhof*this->gcostheta*h;
		pbterm = h*(pb-rho_eq*this->gcostheta*chi_eq*h);
	}

	void SteadyUniformUTheta(double alt_theta, double h, double &u, double &phi, double &pbterm)
	{
		double max=1e2, min=0, I_eq, phi_eq, rho_eq, P_eq, alt_gct = this->g*cos(alt_theta*pi/180.0), chi_eq;
		const double pi = 3.14159265358979323846264338327950288;
		while (max-min > 1e-14)
		{
			phi_eq = this->pp.phim-this->delphi*((max+min)/2.0);
			rho_eq = phi_eq*this->pp.rhog + (1-phi_eq)*this->pp.rhof;
			P_eq = (rho_eq-this->pp.rhof)/rho_eq;
			((this->MuI(0.5*(max+min))-tan(alt_theta*pi/180.0)/P_eq+this->tau0/((rho_eq-this->pp.rhof)*this->gcostheta*h)>0)?max:min)=0.5*(max+min);
		}
		I_eq = 0.5*(max+min);
		phi = this->pp.phim-this->delphi*(I_eq);
		rho_eq = phi*this->pp.rhog + (1-phi)*this->pp.rhof;
		chi_eq = (this->pp.rhof+3.0*rho_eq)/4.0/rho_eq;
		u = UFromI(I_eq,h,(rho_eq-this->pp.rhof)*alt_gct*h);
		double pb = this->pp.rhof*alt_gct*h;
		pbterm = h*(pb-rho_eq*alt_gct*chi_eq*h);
	}
		  
	double UFromI(double I, double h, double ppval)
	{
		return I*sqrt(ppval)*h/(2.0*5.0/2.0*this->d*sqrt(this->pp.rhog));
	}

	void SetGTheta(double g_, double thetaDeg_, double tau0_, double d_, double alpha_, double finalThetaDeg_ = -1.0, double change_t_ = -1.0)
	{
		// const double pi = 3.14159265358979323846264338327950288;
		// this->thetaDeg = thetaDeg_;
		// this->tau0 = tau0_;
		this->d = d_;
		this->alpha = alpha_;

		// this->theta = this->thetaDeg*pi/180.0;
		// this->g = g_;
		// this->gcostheta = g_*cos(this->theta);
		// this->gsintheta = g_*sin(this->theta);
		// this->tantheta = tan(this->theta);
		// this->initThetaDeg = thetaDeg_;
		// this->finalThetaDeg = finalThetaDeg_;
		// this->change_t = change_t_;
		

		this->RegisterParameter("d", Parameter(&d));
		this->RegisterParameter("alpha", Parameter(&alpha));
	}

protected:

	// Iv of a quadratic velocity profile with depth average velocity u and thickness h
	double I(double u, double h, double ppval)
	{
		if (ppval < 1e-8)
			return 1e8;
		else
			return 2.0*this->d*sqrt(this->pp.rhog)*5.0/2.0*(u/h)/sqrt(ppval);
	}

	double d, alpha;
	// MuIParams pp;
};

////////////////////////////////////////////////////////////////////////////////
/// 1-D SW equations with MuIv friction
template <unsigned N_UNSOLVED>
class SWMuIEqn1DRhoVary : public SWMuIEqnRhoVaryBase<3, N_UNSOLVED>
{
public:
	// Inherit constructors from base class
	using SWMuIEqnRhoVaryBase<3, N_UNSOLVED>::SWMuIEqnRhoVaryBase;
	// using SWMuIEqnBase<3>::SWMuIEqnBase;

	enum Variables
	{
		H = 0,
		HU = 1,
		HPHI = 2,
		PBH = 3,
		THETA = 4
	};


	std::string VariableName(int d)
	{
		switch (d)
		{
		case H: return "H";
		case HU: return "HU";
		case HPHI: return "HPHI";
		case PBH: return "PBH";
		case THETA: return "THETA";
		default: return Equation::VariableName(d);
		}
	}

	void SpatialTheta(const double *u)
	{
		const double pi = 3.14159265358979323846264338327950288;
		double spat_theta = u[THETA];
		this->thetaDeg = spat_theta;
		this->theta = this->thetaDeg*pi/180.0;
		this->gcostheta = this->g*cos(this->theta);
		this->gsintheta = this->g*sin(this->theta);
		this->tantheta = tan(this->theta);
	}

	void XConvectionFlux(double *xFlux, const double *u, const double *extras, const double *dextrasdt)
	{
		if (N_UNSOLVED == 1) {
			SpatialTheta(u);
		}
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
		if (N_UNSOLVED == 1) {
			SpatialTheta(u);
		}
		if (u[H] < this->zeroHeightThreshold || std::abs(u[HU]) < this->huThreshold)
			return 0;
		else
			return std::abs(u[HU] / u[H]) + sqrt(this->gcostheta*u[H]);
	}

	void SourceTerms(double dt, double *stvect, double *u, const double *extras, const double *dedt)
	{
		if (N_UNSOLVED == 1) {
			SpatialTheta(u);
		}
		double h = u[H], hu = u[HU], hphi = u[HPHI], pbterm = u[PBH];
		if (h < this->zeroHeightThreshold)
		{
			h = this->zeroHeightThreshold;
			hu = 0;
			hphi = h*this->pp.phim;
		}

		double uval = hu/h;
		double absu = std::abs(uval);
		double signU;

		double phi = hphi/h;

		double rhoBulk = this->pp.rhog*phi + this->pp.rhof*(1-phi);
		double chi = (this->pp.rhof+3*rhoBulk)/4/rhoBulk;
		double P = (rhoBulk-this->pp.rhof)/rhoBulk;

		double pb = pbterm/h+rhoBulk*this->gcostheta*chi*h;
		double ppval = std::max(rhoBulk*this->gcostheta*h-pb,0.0);

		// If u is exactly equal to zero, assume (wrongly...) that the friction is upslope
		if (absu == 0)
			signU = 1;
		else 
			signU = uval/absu;

		double beta = 150.0*this->pp.phim*this->pp.phim*this->pp.etaf/((1.0-this->pp.phim)*(1.0-this->pp.phim)*(1.0-this->pp.phim)*this->d*this->d);

		double D = -2.0/beta/h*(pb-this->pp.rhof*this->gcostheta*h);
		double zeta = 3.0/(2.0*this->alpha*h) + this->gcostheta*this->pp.rhof*P/4.0;


		// Calculate friction law
		double i = this->I(absu, h, ppval);
		double tanpsi = phi - (this->pp.phim-this->pp.delphi*i);

		// if (iv<0) // in case ppval becomes negative
		// {
		// 	alter=1;
		// 	iv = 1;
		// 	ppval = (rhoBulk-pp.rhof)*this->gcostheta*h;
		// 	pb = rhoBulk*this->gcostheta*h;
		// 	tanpsi = 0.01;
		// 	// D = D*10;
		// }
		// if (pb<0) // in case pb becomes negative
		// {
		// 	alter = 1;
		// 	pb = pp.rhof*this->gcostheta*h/2;
		// 	ppval = rhoBulk*this->gcostheta*h-pb;
		// 	tanpsi = -0.01;
		// }

		// mu * (rho-rho_f)/rho
		double mubf = ppval*this->MuI(i);
		if (ppval < 1e-8)
			mubf = 3.0*absu*this->pp.etaf/h; // transition to viscous fluid
		

		double absFriction = (1.0/rhoBulk)*(-mubf - this->tau0 + (rhoBulk-this->pp.rhof)*D*absu);

		double dilatancy = 1.0/this->alpha/h*absu*tanpsi*3.0;

		double psi1 = D*P;
		double psi5 = zeta*D - dilatancy;
		stvect[H] += psi1;

		stvect[HPHI] += -D*phi*this->pp.rhof/rhoBulk;

		stvect[PBH] += (psi5-this->gcostheta*this->pp.rhof*P/4.0*D)*h+(pb-rhoBulk*this->gcostheta*chi*h)*psi1;

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
					stvect[HU] += std::max(absFriction*signU, -stvect[HU]); // add *signU?
	 			}
				else if (stvect[HU]<0)
				{
					stvect[HU] += std::min(absFriction*signU, -stvect[HU]); // add *signU?
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

	void ExtrasFunction(double *dedt, const VectorArray2d *u, double t)
	{
		this->AlterTheta();
	}
};	
#endif
