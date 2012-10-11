#pragma once

// Anisotropic Ward BRDF, based on the original paper.
// Simply add it to the returned bsdf of a material to use it.

#include "spectrum.h"
#include "reflection.h"

class AnisoWardBrdf: public BxDF 
{
public:
	AnisoWardBrdf(const Spectrum &rs, float ax, float ay);
	Spectrum f(const Vector &wo, const Vector &wi) const;

private:
	Spectrum Rs;
	float Ax, Ay;
	float InvAx2, InvAy2;
	float FourPiAxAy;
};