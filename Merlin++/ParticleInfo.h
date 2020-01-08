/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2020 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#ifndef ParticleInfo_h
#define ParticleInfo_h

#include "PhysicalConstants.h"

struct ParticleInfo
{
	ParticleInfo(double _mass, double _charge) :
		mass(_mass), charge(_charge)
	{
	}
	double mass;
	double charge;

	double GetMassMev() const
	{
		return mass
			   / PhysicalConstants::ElectronCharge * PhysicalConstants::SpeedOfLight * PhysicalConstants::SpeedOfLight *
			   1e-6;
	}

};

// default ultra-relativistic particles
const ParticleInfo default_particle{0, 1};
const ParticleInfo default_particle_negative{0, -1};
const ParticleInfo default_particle_neutral{0, 0};

#endif
