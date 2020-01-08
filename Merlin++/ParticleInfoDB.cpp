/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2020 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "ParticleInfoDB.h"
#include "MerlinException.h"
#include "PhysicalConstants.h"

ParticleInfoDB::ParticleInfoDB()
{
	db.insert({"e", {PhysicalConstants::ElectronMass, -1}});
	db.insert({"e+", {PhysicalConstants::ElectronMass, 1}});
	db.insert({"p", {PhysicalConstants::ProtonMass, 1}});
	db.insert({"pbar", {PhysicalConstants::ProtonMass, -1}});
	db.insert({"muon-", {PhysicalConstants::MuonMass, -1}});
	db.insert({"muon+", {PhysicalConstants::MuonMass, 1}});
}

const ParticleInfo* ParticleInfoDB::FindParticle(std::string name)
{
	auto pos = db.find(name);
	if(pos == db.end())
	{
		throw MerlinException("Could not find particle type in ParticleInfoDB:" + name);
	}

	return &(pos->second);
}
