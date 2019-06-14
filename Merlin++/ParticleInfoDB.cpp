/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2019 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "ParticleInfoDB.h"
#include "PhysicalConstants.h"

using namespace PhysicalConstants;

const std::map<std::string, ParticleInfo> ParticleInfoDB =
{
	{"", {0, 0, 1}},  // default ultra-relativistic particle
	{"e", {ElectronMass, ElectronMassMeV, -1}},
	{"e+", {ElectronMass, ElectronMassMeV, 1}},
	{"p", {ProtonMass, ProtonMassMeV, 1}},
	{"p-", {ProtonMass, ProtonMassMeV, -1}},
	{"muon-", {MuonMass, MuonMassMeV, -1}},
	{"muon+", {MuonMass, MuonMassMeV, 1}},
};
