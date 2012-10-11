#include "AnisoWardBrdf.h"

AnisoWardBrdf::AnisoWardBrdf(const Spectrum &rs, float ax, float ay)
	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), Rs(rs), Ax(ax), Ay(ay)
{
	FourPiAxAy = (4.0f*M_PI*Ax*Ay);
	InvAx2 = 1.0f/(Ax*Ax);
	InvAy2 = 1.0f/(Ay*Ay);
}

Spectrum AnisoWardBrdf::f(const Vector &wo, const Vector &wi) const
{
	Spectrum specular(0.0f);	
	const Vector wh = wi + wo;
	if(wh.z <= 0.0f)
		return specular;
	float cosi_coso = wi.z*wo.z;
	if(cosi_coso <= 0.0f)
		return specular;
	const float expTerm = expf(-1.0f * (wh.x*wh.x*InvAx2 + wh.y*wh.y*InvAy2) / (wh.z*wh.z));
	cosi_coso = sqrtf(cosi_coso);
	specular = Rs*expTerm/(cosi_coso*FourPiAxAy); 
	return specular;
}