/// \file Baseline.h
/// Header file for the Baseline class.

#ifndef BASELINE_H
#define BASELINE_H

// Include headers for the unordered map.  Note, this may need to be just
// <unordered_map> if compiled in MSVS.
#include <complex>
#include <string>
#include <tr1/unordered_map>
#include <vector>

// Header files for other libraries
extern "C" {
#include "exchange.h"
}

using namespace std;

typedef std::tr1::unordered_map<std::string, complex<double>> VisHash;
typedef std::tr1::unordered_map<std::string, double> Vis2Hash;

// Forward class declarations:
class UVPoint;
class Array;
class Station;
class Target;

/// \class Baseline simulator.h
/// \brief A class representing a baseline.
class Baseline {
private:
  double xyz[3]; // True XYZ coordinates

  string name;
  int indicies[2];
  VisHash mVisValues;   // Stores computed visibility values
  Vis2Hash mVis2Errors; // Stores computed/stored visibility squared error values.

public:
  Baseline(void);
  Baseline(Array *array, int station1, int station2);
  Baseline(Station *station1, Station *station2);

  /// \todo Rewrite this function to work with the new class definition.
  // double Geometric_OPD(double hour_angle, double source_declination);

private:
  complex<double> ComputeVisibility(Target &target, UVPoint uv, double wavelength, double dwavelength);
  string GetHashKey(Target &target, UVPoint uv);

public:
  string GetName(void);
  complex<double> GetVisibility(Target &target, double hour_angle, double wavelength, double dwavelength);
  complex<double> GetVisibility(Target &target, UVPoint uv_coords, double wavelength, double dwavelength);
  void SetHashKey(complex<double> vis, Target &target, UVPoint uv);
  void ComputeTargetVisibilities(int nuv, complex<double> *cvis, Target &target, UVPoint *uv_list);

  double GetVis2(Target &target, double hour_angle, double wavelength, double dwavelength);
  double GetVis2(Target &target, UVPoint uv_coords, double wavelength, double dwavelength);

  UVPoint UVcoords(double hour_angle, double source_declination);

  int GetStationID(int station_num);
};

// Hash data type
typedef std::tr1::unordered_map<std::string, Baseline *> BaselineHash;

vector<Baseline> ComputeBaselines(vector<Station> &stations);
BaselineHash ComputeBaselineHash(vector<Baseline> &baselines);

#endif // #ifndef BASELINE_H
