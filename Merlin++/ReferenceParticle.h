/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef ReferenceParticle_h
#define ReferenceParticle_h 1

#include "merlin_config.h"
#include "utils.h"
#include <cassert>
#include "ParticleInfo.h"

/**
 *	A ReferenceParticle represents that particle which sits
 *	on the nominal orbit. It is responsible for maintaining
 *	the reference momentum and time (ct) for the bunch or
 *	map. ReferenceParticle cannot be instantiated, but is
 *	designed as a mixin for bunch or map-like classes.
 */
class ReferenceParticle
{
public:
	virtual ~ReferenceParticle()
	{
	}

	/**
	 *	Returns the reference momentum in GeV/c.
	 *	@return Reference momentum (GeV/c)
	 */
	double GetReferenceMomentum() const;

	/**
	 *	Returns the reference time in ct (meters).
	 *	@return Reference time in ct (m)
	 */
	double GetReferenceTime() const;

	/**
	 *	Returns either +1, 0 or -1.
	 *
	 *	@retval +1
	 *	@retval  0
	 *	@retval -1
	 */
	double GetChargeSign() const;

	/// Charge of an individual particle
	double GetParticleCharge() const;

	/// Access method: Get particle mass
	double GetParticleMass() const;

	/// Access method: Get particle mass (MeV)
	double GetParticleMassMeV() const;

	/**
	 *	Sets the reference momentum to p GeV/c. p must be
	 *	greater than zero.
	 */
	void SetReferenceMomentum(double p);

	/**
	 *	Increments the reference momentum by dp GeV/c, returning
	 *	the new value.
	 */
	double IncrReferenceMomentum(double dp);

	/**
	 *	Sets the reference time in ct (meters).
	 */
	void SetReferenceTime(double ct);

	/**
	 *	Increments the reference time by dct meters.
	 */
	double IncrReferenceTime(double dct);

protected:

	ReferenceParticle(double p, double q = 1, const ParticleInfo* ptype = nullptr);

	/**
	 *	Sets the charge sign.
	 */
	void SetChargeSign(double q);

	/**
	 * Data Members for Class Attributes
	 */

	/**
	 *	reference momentum in GeV/c
	 */
	double p0;

	/**
	 *	reference time in ct (meters)
	 */
	double ct0;

	const ParticleInfo * type = &default_particle;
};

inline ReferenceParticle::ReferenceParticle(double p, double q, const ParticleInfo* ptype) :
	p0(p), ct0(0)
{
	assert(p > 0);
	SetChargeSign(q);
	if(ptype != nullptr) // override q
	{
		type = ptype;
	}
}

inline double ReferenceParticle::GetReferenceMomentum() const
{
	return p0;
}

inline double ReferenceParticle::GetReferenceTime() const
{
	return ct0;
}

inline double ReferenceParticle::GetChargeSign() const
{
	const double &charge = type->charge;
	if(charge > 0)
		return 1;
	else if(charge < 0)
		return -1;
	return 0;
}

inline double ReferenceParticle::GetParticleCharge() const
{
	return type->charge;
}

inline double ReferenceParticle::GetParticleMass() const
{
	return type->mass;
}

inline double ReferenceParticle::GetParticleMassMeV() const
{
	return type->mass_mev;
}

inline void ReferenceParticle::SetReferenceMomentum(double p)
{
	p0 = p;
	assert(p0 > 0);
}

inline double ReferenceParticle::IncrReferenceMomentum(double dp)
{
	p0 += dp;
	assert(p0 > 0);
	return p0;
}

inline void ReferenceParticle::SetReferenceTime(double ct)
{
	ct0 = ct;
}

inline double ReferenceParticle::IncrReferenceTime(double dct)
{
	ct0 += dct;
	return ct0;
}

inline void ReferenceParticle::SetChargeSign(double q)
{
	if(q > 0)
	{
		type = &default_particle;
	}
	else if(q < 0)
	{
		type = &default_particle_negative;
	}
	else
	{
		type = &default_particle_neutral;
	}
}

#endif
