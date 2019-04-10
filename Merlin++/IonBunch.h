/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#ifndef IonBunch_h
#define IonBunch_h 1

#include "ParticleBunch.h"

namespace ParticleTracking
{
/**
 * A bunch type to hold ions.
 *
 * Charges are given in elementary charge units for example would be
 * 2 for an He++ (alpha particle). Masses in Kg (the ProtonMass
 * constant may be useful, e.g. 2*ProtonMass).
 *
 * Note that all Merlin++ bunches have a total charge value, used for collective
 * effects. This is distinct from the particle charge, which is used for single
 * particle effects, e.g. tracking.
 */

class IonBunch: public ParticleBunch
{
public:
	/**
	 * Constructs an IonBunch using the specified particle momentum,
	 * total charge, the particle array, particle charge and particle mass.
	 * Note that on exit, particles is empty.
	 */
	IonBunch(double P0, double Qtot, PSvectorArray& particles, double Qp, double Mp) :
		ParticleBunch(P0, Qtot, particles), particle_charge(Qp), particle_mass(Mp)
	{
	}

	/**
	 * Read phase space vectors from specified input stream.
	 */
	IonBunch(double P0, double Qtot, std::istream& is, double Qp, double Mp) :
		ParticleBunch(P0, Qtot, is), particle_charge(Qp), particle_mass(Mp)
	{
	}

	/**
	 * Constructs an empty ElectronBunch with the specified momentum P0,
	 * particle charge Qp, particle mass Mp, and charge per macro particle Qm
	 * (default = +1).
	 */
	IonBunch(double P0, double Qp, double Mp, double Qm = 1) :
		ParticleBunch(P0, Qm), particle_charge(Qp), particle_mass(Mp)
	{
	}

	IonBunch(size_t np, const ParticleDistributionGenerator & generator, const BeamData& beam, double Qp, double Mp,
		ParticleBunchFilter* filter = nullptr) :
		ParticleBunch(np, generator, beam, filter), particle_charge(Qp), particle_mass(Mp)
	{
	}

	virtual bool IsStable() const override;
	virtual double GetParticleMass() const override;
	virtual double GetParticleMassMeV() const override;
	virtual double GetParticleLifetime() const override;
	virtual double GetParticleCharge() const override;

private:
	double particle_charge;
	double particle_mass;
};

}

#endif
