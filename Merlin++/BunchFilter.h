/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#ifndef _BUNCHFILTER_H_
#define _BUNCHFILTER_H_

#include "PSvector.h"

namespace ParticleTracking
{
/**
 * Filter to be used during bunch creation
 *
 * It is often useful to filter a particle distribution before tracking,
 * for example avoid tracking particles which are not interesting in a
 * given simulation.
 *
 * The bunch filter can be passed to ParticleBunch::ParticleBunch()
 *
 * Use ParticleBunchFilter::filter_in_realspace to choose if the filter
 * is applied in normalised or real space.
 */
class ParticleBunchFilter
{
public:

	virtual ~ParticleBunchFilter();

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	virtual bool Apply(const PSvector& v) const = 0;

	/**
	 * Filtering can either be done in normalised or real space
	 */
	bool filter_in_realspace = 1;
};

/**
 * Filter that allows particles above a given x value.
 *
 * This can be used in a loss map simulation to limit the simulation to
 * particles that will hit the collimator on the first turn.
 */
class HorizontalHaloParticleBunchFilter: public ParticleBunchFilter
{
public:

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	bool Apply(const PSvector& v) const;

	void SetHorizontalLimit(double);
	void SetHorizontalOrbit(double);

private:
	double limit;
	double orbit;
};

/**
 * Filter that allows particles above a given y value.
 *
 * This can be used in a loss map simulation to limit the simulation to
 * particles that will hit the collimator on the first turn.
 */
class VerticalHaloParticleBunchFilter: public ParticleBunchFilter
{
public:

	/**
	 *	Used by the ParticleBunch constructor object to select
	 *	vectors for inclusion in a ParticleBunch.
	 */
	bool Apply(const PSvector& v) const;

	void SetVerticalLimit(double);

private:
	double limit;
};

}   //End particle tracking namespace

#endif
