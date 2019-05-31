/*
 * Merlin++: C++ Class Library for Charged Particle Accelerator Simulations
 * Copyright (c) 2001-2018 The Merlin++ developers
 * This file is covered by the terms the GNU GPL version 2, or (at your option) any later version, see the file COPYING
 * This file is derived from software bearing the copyright notice in merlin4_copyright.txt
 */

#include "../tests.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <memory>

#include "MADInterface.h"
#include "LatticeFunctions.h"
#include "AcceleratorModel.h"
#include "AcceleratorModelConstructor.h"
#include "SectorBend.h"
#include "PhysicalConstants.h"
#include "ParticleBunch.h"
#include "ProtonBunch.h"
#include "IonBunch.h"
#include "ParticleTracker.h"
#include "SymplecticIntegrators.h"

using namespace std;
using namespace PhysicalConstants;

const string dipole_tfs = "* KEYWORD NAME S L ANGLE TILT K1L E1 E2\n"
	"$ %s %s %le %le %le %le %le %le %le\n"
	"SBEND D1 0.5 0.5 0.00174532925199432 0 0 0 0\n" // 0.5m, 0.1 deg bend
//"SBEND D2 1.3 0.8 0.01745329251994329 0 0 0 0\n" // 0.8m, 1 deg bend
;

int main(int argc, char* argv[])
{
	const auto particle_moms = {1.0, 5.0, /*100.0*/};

	const auto particles = {
		make_tuple("p", ProtonMass, 1),
		make_tuple("pbar", ProtonMass, -1),
		make_tuple("deuteron", ProtonMass * 2, 1),
		make_tuple("alpha", ProtonMass * 4, 2),
	};

	for(auto& particle : particles)
	{
		for(auto &particle_mom : particle_moms)
		{
			string pname = get<0>(particle);
			auto pmass = get<1>(particle);
			auto pcharge = get<2>(particle);
			double rigidity = particle_mom * 1e9 / SpeedOfLight / pcharge;
			cout << "\n" <<  pname << " momentum " << particle_mom << "GeV" << endl;
			cout << u8"Bρ = " << rigidity << endl;
			istringstream ins(dipole_tfs);
			MADInterface myMADinterface(&ins, particle_mom, pcharge);
			//MADInterface myMADinterface("DeveloperTools/CodeTesting/data/twiss.7.0tev.b1_new.tfs", particle_mom);

			AcceleratorModel* model = myMADinterface.ConstructModel(); // deleted by AcceleratorModelConstructor by MADInterface

			vector<SectorBend *> bends;
			model->ExtractTypedElements(bends);
			//cout << bends.size() << " bends in lattice" << endl;

			double l, b;
			for(auto& bend : bends)
			{
				l = bend->GetLength();
				b = bend->GetB0();
				double angle = b * l / rigidity;
				double angle_deg = angle / 2 / pi * 360;
				cout << "SBEND L=" << l << " m  B0=" <<  b << " T : calc angle =" << angle << " (" << angle_deg
					 << " deg)";
				cout << ((angle > 1e-3) ? "" : " BAD sign");
				cout << (abs(angle_deg - 0.1) < 1e-12 || abs(angle_deg - 1.0) < 1e-12  ? "" : " BAD val");

				cout << endl;

				unique_ptr<ParticleBunch> myBunch;

				if(pname == "p")
				{
					myBunch = make_unique<ProtonBunch>(particle_mom);
				}
				else if(pname == "pbar")
				{
					myBunch = make_unique<ProtonBunch>(particle_mom, -1);
				}
				else if(pname == "deuteron")
				{
					myBunch = make_unique<IonBunch>(particle_mom, 1, 2 * ProtonMass, 1);
				}
				else if(pname == "alpha")
				{
					myBunch = make_unique<IonBunch>(particle_mom, 2, 4 * ProtonMass, 1);
				}

				if(myBunch)
				{
					Particle p0(0), p1(0), p2(0);
					p1.dp() = 1e-3;
					p2.dp() = -1e-3;
					myBunch->AddParticle(p0);
					myBunch->AddParticle(p1);
					myBunch->AddParticle(p2);

					cout << "sign = " << myBunch->GetChargeSign() << "  tot_charge = " << myBunch->GetTotalCharge()
						 << "  particle_charge = " << myBunch->GetParticleCharge() << endl;

					auto ctor2 = make_unique<AcceleratorModelConstructor>();
					ctor2->AppendComponent(bend);
					auto model2 = ctor2->GetModel();

					auto tracker = make_unique<ParticleTracker>(model2->GetRing(), myBunch.get());

					tracker->SetIntegratorSet(new ParticleTracking::SYMPLECTIC::StdISet());
					//tracker->SetIntegratorSet(new ParticleTracking::THIN_LENS::StdISet());
					//tracker->SetIntegratorSet(new ParticleTracking::TRANSPORT::StdISet());

					tracker->Track(myBunch.get());
					p0 = myBunch->GetParticles()[0];
					p1 = myBunch->GetParticles()[1];
					p2 = myBunch->GetParticles()[2];
					cout << p0.x() << " should be 0" << ((abs(p0.x()) < 1e-12) ? "" : " BAD") << endl;
					cout << p1.x() << " should be > 0" << ((p1.x() > 1e-12) ? "" : " BAD") << endl;
					cout << p2.x() << " should be < 0" << ((p2.x() < 1e-12) ? "" : " BAD") << endl;
				}
				else
				{
					cout << "no bunch" << endl;
				}

			}
		}
	}
}
