
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_ANISOWARD_H
#define PBRT_MATERIALS_ANISOWARD_H

// materials/anisoward.h*
#include "pbrt.h"
#include "material.h"

// AnisoWardMaterial Declarations
class AnisoWardMaterial : public Material {
public:
    // AnisoWardMaterial Public Methods
    AnisoWardMaterial(Reference<Texture<Spectrum> > kd,
                    Reference<Texture<Spectrum> > ks,
                    Reference<Texture<float> > au,
                    Reference<Texture<float> > av,
                    Reference<Texture<float> > bump)
        : Kd(kd), Ks(ks), alphaU(au), alphaV(av), bumpMap(bump) {
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom,
                  const DifferentialGeometry &dgShading,
                  MemoryArena &arena) const;
private:
    // AnisoWardMaterial Private Data
    Reference<Texture<Spectrum> > Kd, Ks;
    Reference<Texture<float> > alphaU, alphaV, bumpMap;
};


AnisoWardMaterial *CreateAnisoWardMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_ANISOWARD_H
