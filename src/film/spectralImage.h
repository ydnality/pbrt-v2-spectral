
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

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

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_FILM_IMAGE_H
#define PBRT_FILM_IMAGE_H

// film/image.h*
#include "pbrt.h"
#include "film.h"
#include "sampler.h"
#include "filter.h"
#include "paramset.h"

const bool debugMode = false;    //Andy added this flag for debugging

// ImageFilm Declarations
class SpectralImageFilm : public Film {
public:
    // SpectralImageFilm Public Methods
    SpectralImageFilm(int xres, int yres, Filter *filt, const float crop[4],
              const string &filename, bool openWindow);
    ~SpectralImageFilm() {
        delete pixels;
        delete filter;
        delete[] filterTable;
    }
    void AddSample(const CameraSample &sample, const Spectrum &L, const Ray &currentRay);
    void Splat(const CameraSample &sample, const Spectrum &L);
    void GetSampleExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void GetPixelExtent(int *xstart, int *xend, int *ystart, int *yend) const;
    void WriteImage(float splatScale);
    void UpdateDisplay(int x0, int y0, int x1, int y1, float splatScale);
    void ParseConversionMatrix(string filename);     //Andy Added: function that parses the conversion matrix file

private:
    //Andy added to allow for conversion matrix file parsing
    float * conversionMatrix;     //actual matrix storage
    int nCMRows;
    int nCMCols; 
    float * waveSpecify;

    // SpectralImageFilm Private Data
    Filter *filter;
    float cropWindow[4];
    string filename;
    int xPixelStart, yPixelStart, xPixelCount, yPixelCount;
    //Andy: modified this to allow for multispectral
    struct Pixel {
        Pixel() {
            for (int i = 0; i < nSpectralSamples; ++i) c[i] = splatC[i] = 0.f;
            weightSum = 0.f;
        }
        float c[nSpectralSamples];   //changed here to allow for >3 channels
        float weightSum;
        float splatC[nSpectralSamples];   //changed here
        float pad;
        float Z;   //added this
    };
    BlockedArray<Pixel> *pixels;
    float *filterTable;
};


SpectralImageFilm *CreateSpectralImageFilm(const ParamSet &params, Filter *filter);

#endif // PBRT_FILM_IMAGE_H
