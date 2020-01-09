/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2020 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 */

#include "../tests.h"
#include <iostream>

#include "ParticleInfoDB.h"
#include "ParticleBunch.h"
#include "ProtonBunch.h"
#include "ElectronBunch.h"

using namespace std;
using namespace ParticleTracking;
int main(int argc, char* argv[])
{
	auto pdb = ParticleInfoDB();

	auto p_e = pdb.FindParticle("e");
	assert_close(p_e->mass, 9.10938215e-31, 1e-40);
	assert(p_e->charge == -1);
	assert_close(p_e->GetMassMev(), 0.510998910, 1e-9);

	auto p_epos = pdb.FindParticle("e+");
	assert_close(p_e->mass, p_epos->mass, 1e-40);

	auto p_p = pdb.FindParticle("p");
	assert_close(p_p->mass, 1.672621637e-27, 1e-40);
	assert(p_p->charge == 1);
	assert_close(p_p->GetMassMev(), 938.272013, 1e-6);
	p_p = pdb.FindParticle("pbar");
	assert_close(p_p->mass, 1.672621637e-27, 1e-40);
	assert(p_p->charge == -1);

	auto b_e = ParticleBunch(100, -10, pdb.FindParticle("e"));
	assert(b_e.GetParticleCharge() == -1);
	auto b_p = ParticleBunch(100, 10, pdb.FindParticle("p"));
	assert_close(b_p.GetParticleMassMeV(), 938.272013, 1e-6);

	auto bp_p = ProtonBunch(100, 10);
	assert_close(bp_p.GetParticleMassMeV(), 938.272013, 1e-6);

	auto be_e = ElectronBunch(100, 10);
	assert_close(be_e.GetParticleMassMeV(), 0.510998910, 1e-9);
}
