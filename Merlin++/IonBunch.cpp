/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "IonBunch.h"
#include "PhysicalConstants.h"
#include "PhysicalUnits.h"

using namespace PhysicalConstants;
using namespace PhysicalUnits;

namespace ParticleTracking
{
double IonBunch::GetParticleMass() const
{
	return particle_mass;
}

double IonBunch::GetParticleMassMeV() const
{
	return particle_mass / ElectronCharge * SpeedOfLight * SpeedOfLight * eV / MeV;
}

double IonBunch::GetParticleLifetime() const
{
	return 0;
}

bool IonBunch::IsStable() const
{
	return true;
}

double IonBunch::GetParticleCharge() const
{
	return particle_charge;
}

}
