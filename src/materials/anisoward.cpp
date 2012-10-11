
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
	This file is an extension to pbrt.  It implements a material that
	uses the anisotropic Ward model.  Ben Heasly wrote anisoward.cpp 
	and anisoward.h based on PBRT's plastic.c and plastic.h.  
	AnisoWardBrdf.cpp and AnisoWardBrdf.h were downloaded from the 
	University of Freiburg at 
	http://cg.informatik.uni-freiburg.de/course_notes/SS12/ACG/.

 */


// materials/anisoward.cpp*
#include "stdafx.h"
#include "materials/anisoward.h"
#include "spectrum.h"
#include "reflection.h"
#include "paramset.h"
#include "texture.h"

#include "AnisoWardBrdf.h"

// AnisoWardMaterial Method Definitions
BSDF *AnisoWardMaterial::GetBSDF(const DifferentialGeometry &dgGeom,
                               const DifferentialGeometry &dgShading,
                               MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);
    Spectrum kd = Kd->Evaluate(dgs).Clamp();
    BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
    Spectrum ks = Ks->Evaluate(dgs).Clamp();
    float au = alphaU->Evaluate(dgs);
	float av = alphaV->Evaluate(dgs);
    BxDF *spec = BSDF_ALLOC(arena, AnisoWardBrdf)(ks, au, av);
    bsdf->Add(diff);
    bsdf->Add(spec);
    return bsdf;
}


AnisoWardMaterial *CreateAnisoWardMaterial(const Transform &xform,
        const TextureParams &mp) {
    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    Reference<Texture<float> > au = mp.GetFloatTexture("alphaU", 0.1f);
	Reference<Texture<float> > av = mp.GetFloatTexture("alphaV", 0.1f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTexture("bumpmap", 0.f);
    return new AnisoWardMaterial(Kd, Ks, au, av, bumpMap);
}


