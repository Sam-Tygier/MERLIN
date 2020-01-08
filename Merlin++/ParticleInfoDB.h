/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2020 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#ifndef ParticleInfoDB_h
#define ParticleInfoDB_h

#include <string>
#include <map>

#include "ParticleInfo.h"

class ParticleInfoDB
{
public:
	ParticleInfoDB();

	const ParticleInfo* FindParticle(std::string name);

public:
	std::map<std::string, ParticleInfo> db;
};

#endif
