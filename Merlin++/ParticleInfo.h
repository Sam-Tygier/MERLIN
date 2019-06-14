/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2019 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#ifndef ParticleInfo_h
#define ParticleInfo_h

#include <map>
#include <string>

struct ParticleInfo
{
	ParticleInfo(double _mass, double _mass_mev, double _charge) :
		mass(_mass), mass_mev(_mass_mev), charge(_charge)
	{
	}
	double mass;
	double mass_mev;
	double charge;

};

// default ultra-relativistic particles
const ParticleInfo default_particle{0, 0, 1};
const ParticleInfo default_particle_negative{0, 0, -1};
const ParticleInfo default_particle_neutral{0, 0, 0};

#endif
