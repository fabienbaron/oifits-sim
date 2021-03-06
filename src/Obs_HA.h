#include "Observation.h"
#include <complex>
#ifndef OBS_HA_H
#define OBS_HA_H

class Obs_HA : public Observation {
private:
  double mHA;
  bool mComputeHA;

private:
  double GetLocalSiderealTime(double jd_high, double jd_low, double ee, double array_long);
  double GetSiderealTime(double jd_high, double jd_low, double ee);
  void ComputeTargetVisibilities(const int nuv, complex<double> *cvis, Target &target, const UVPoint *uv_list);

public:
  Obs_HA(Array *array, vector<Station *> stations, string exclude_baselines);
  Obs_HA(Array *array, double hour_angle, string telescopes, string exclude_baselines);
  Obs_HA(Array *array, double mjd, double time, string telescopes, string exclude_baselines);

  static vector<Observation *> MakeObservations(Array *array, double start, double stop, double every, string telescopes);
  static vector<Observation *> ReadObservation_HA(Array *array, vector<string> lines, int i);
  static vector<Observation *> ReadObservation_Descriptive(Array *array, vector<string> lines, int i);

  void WriteVis(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode, Target *target,
                NoiseModel *noisemodel, Rand_t random_seed);
  void WriteVis2(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode, Target *target,
                 NoiseModel *noisemodel, Rand_t random_seed);
  void WriteT3(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode, Target *target,
               NoiseModel *noisemodel, Rand_t random_seed);
  void WriteT4(fitsfile *outfile, UVPoint **puv_list, complex<double>* &cvis, Array *array, Combiner *combiner, SpectralMode *spec_mode, Target *target,
               NoiseModel *noisemodel, Rand_t random_seed);

  double GetHA(double targ_ra);
};

#endif // OBS_HA_H
