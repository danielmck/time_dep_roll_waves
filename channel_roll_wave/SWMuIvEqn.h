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




template <unsigned DIM, unsigned N_UNSOLVED>
class SWMuIvEqnBase : public Equation
{
public:
	/// Constructor, provide g, theta in degrees, and struct with MuIv friction law params
	SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, MuIvParams pp_, int nExtras_ = 0) 
		: SWMuIvEqnBase(g_, thetaDeg_, tau0_, nExtras_)
	{
		SetMuIvParams(pp_);
	}

	SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, double finalThetaDeg_=-1, double change_t_=-1, int nExtras_ = 0) :
		Equation(DIM+1, N_UNSOLVED, nExtras_),  zeroHeightThreshold(1e-7), huThreshold(1e-14)
	{
		SetGTheta(g_, thetaDeg_, tau0_,finalThetaDeg_,change_t_);
		P = 0;
		rhoBulk = 0;
		stoppedMaterialHandling = false;
		this->RegisterParameter("smh", Parameter(&stoppedMaterialHandling));
		solverPtr = nullptr;
	}

	/// Provide g, theta in degrees.
	// SWMuIvEqnBase(double g_, double thetaDeg_, double tau0_, int nExtras_ = 0) :
	// 	SWMuIvEqnBase(g_, thetaDeg_, tau0_, -1.0, nExtras_)
	// {
		
	// }

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

		rhoBulk = (pp.phi*pp.rhog+(1-pp.phi)*pp.rhof);
		P = (rhoBulk-pp.rhof)/rhoBulk;

		this->RegisterParameter("rho", Parameter(&rhoBulk));
		this->RegisterParameter("P", Parameter(&P));
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
		return SteadyUniformU(h)/sqrt(gcostheta*h);
	}

	double SteadyUniformUTheta(double alt_theta, double h)
	{
		double max=1e8, min=0;
		while (max-min > 1e-14)
		{
			((MuIv(Iv(0.5*(max+min),h))-tan(alt_theta*pi/180.0)/P+tau0/((rhoBulk-pp.rhof)*g*cos(alt_theta*pi/180.0)*h)>0)?max:min)=0.5*(max+min);
		}
		return 0.5*(max+min);
	}

	double SteadyUniformFrTheta(double alt_theta, double h)
	{
		return SteadyUniformUTheta(alt_theta,h)/sqrt(g*cos(alt_theta*pi/180.0)*h);
	}
		  
	void SetGTheta(double g_, double thetaDeg_, double tau0_, double finalThetaDeg_ = -1.0, double change_t_ = 0.0)
	{
		const double pi = 3.14159265358979323846264338327950288;
		thetaDeg = thetaDeg_;
		initThetaDeg = thetaDeg;
		if (finalThetaDeg_>0){
			finalThetaDeg = finalThetaDeg_;
		}
		else {
			finalThetaDeg = thetaDeg;
		}
		tau0 = tau0_;
		change_t = change_t_;

		theta = thetaDeg*pi/180.0;
		g = g_;
		gcostheta = g_*cos(theta);
		gsintheta = g_*sin(theta);
		tantheta = tan(theta);

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

	HyperbolicSolver *solverPtr;

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
	
	double thetaDeg, initThetaDeg, finalThetaDeg, theta, g, tau0, gcostheta, gsintheta, tantheta, change_t;
	MuIvParams pp;
	bool stoppedMaterialHandling;
};

////////////////////////////////////////////////////////////////////////////////
/// 1-D SW equations with MuIv friction
template <unsigned N_UNSOLVED>
class SWMuIvEqn1D : public SWMuIvEqnBase<1, N_UNSOLVED>
{
public:
	// Inherit constructors from base class
	using SWMuIvEqnBase<1,N_UNSOLVED>::SWMuIvEqnBase;

	enum Variables
	{
		H = 0,
		HU = 1,
		THETA = 2
	};


	std::string VariableName(int d)
	{
		switch (d)
		{
		case H: return "H";
		case HU: return "HU";
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

		if (u[H] < this->zeroHeightThreshold)
			xFlux[HU] = 0;
		else
			xFlux[HU] = u[HU] * u[HU] / u[H] + 0.5 * this->gcostheta * u[H] * u[H];
	}

	double XWaveSpeeds(const double *u, const double *extras, const double *dextrasdt)
	{
		if (N_UNSOLVED == 1) {
			SpatialTheta(u);
		}
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
		if (N_UNSOLVED == 1) {
			SpatialTheta(u);
		}
		double h = u[H], hu = u[HU];

		if (h < this->zeroHeightThreshold)
		{
			h = this->zeroHeightThreshold;
			hu = 0;
		}

		double uval = hu/h;
		double absu = std::abs(uval);
		double signU;

		// If u is exactly equal to zero, assume (wrongly...) that the friction is upslope
		if (absu == 0)
			signU = 1;
		else 
			signU = uval/absu;


		// Calculate friction law
		double iv = this->Iv(absu, h);

		// mu * (rho-rho_f)/rho
		double mubf = this->P*this->MuIv(iv);

		double absFriction = -mubf * h * this->gcostheta - this->tau0/this->rhoBulk;

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

	void ExtrasFunction(double *dedt, const VectorArray2d *u, double t)
	{
		this->AlterTheta();
	}
};
#endif
