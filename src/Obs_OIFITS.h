#include "Observation.h"

#ifndef OBS_OIFITS_H
#define OBS_OIFITS_H

class Obs_OIFITS : public Observation 
{
  private:
    string mstrFilename;
    
  public:
    Obs_OIFITS(string filename);
  
    static vector <Observation*> ReadObservation_OIFITS(string filename);

    oi_vis  GetVis( UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    oi_vis2 GetVis2(UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    oi_t3   GetT3(  UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    oi_t4   GetT4(  UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
};

#endif // OBS_OIFITS_H
