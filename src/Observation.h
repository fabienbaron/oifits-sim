
/// \file Observation.h

/// \class Observation
/// Stores information about an observation.

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <complex>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
// Include headers for the unordered map.  Note, this may need to be just
// <unordered_map> if compiled in MSVS.
#include <tr1/unordered_map>

extern "C" {
#include "exchange.h"
#include "oifile.h"
#include "random.h"
}

class Array;
class Target;
class Baseline;
class NoiseModel;
class Combiner;
class SpectralMode;
class UVPoint;
class Station;
class Triplet;
class Quadruplet;

using namespace std;
typedef std::tr1::unordered_map<string, string> BLNameHash;
typedef std::tr1::unordered_map<string, string> TNameHash;

enum ObsType { HOUR_ANGLE = 0, DESCRIPTIVE = 1, OIFITS = 2 };

class Observation {
protected:
  Array *mArray;
  double mJD;
  bool mbHasVIS;
  bool mbHasV2;
  bool mbHasT3;
  bool mbHasT4;

  ObsType mObsType;

  vector<Station *> mStations;
  vector<Baseline *> mBaselines;
  vector<Triplet *> mTriplets;
  vector<Quadruplet *> mQuadruplets;

  vector<Station *> FindStations(string telescopes);
  vector<Baseline *> FindBaselines(vector<Station *> stations, string exclude_baselines);
  vector<Triplet *> FindTriplets(vector<Station *> stations, string exclude_baselines);
  vector<Quadruplet *> FindQuadruplets(vector<Station *> stations, string exclude_baselines);

public:
  Observation(void);
  ~Observation(void);

  // static vector <Observation*> ReadObservations(Array * array, string
  // filename, string comment_chars);

  static vector<Observation *> ParseCommandLine(Array *array, char *argv[], int i, int argc, string comment_chars);

private:
  static vector<Observation *> ParseCommandLineObs(Array *array, char *argv[], int i, int argc);
  static vector<Observation *> ImportFile(Array *array, string filename, string comment_chars);

public:
  virtual void WriteVis(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode,
                        Target *target, NoiseModel *noisemodel, Rand_t random_seed) = 0;
  virtual void WriteVis2(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode,
                         Target *target, NoiseModel *noisemodel, Rand_t random_seed) = 0;
  virtual void WriteT3(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode,
                       Target *target, NoiseModel *noisemodel, Rand_t random_seed) = 0;
  virtual void WriteT4(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode,
                       Target *target, NoiseModel *noisemodel, Rand_t random_seed) = 0;

  int GetNumStations(void);
  Station *GetStation(int sta_index);
  bool HasVIS(void);
  bool HasV2(void);
  bool HasT3(void);
  bool HasT4(void);

  ObsType GetObsType(void);
};

#endif // OBSERVATION_H
