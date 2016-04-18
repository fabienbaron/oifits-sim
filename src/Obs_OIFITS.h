#include "Observation.h"

#ifndef OBS_OIFITS_H
#define OBS_OIFITS_H

class Obs_OIFITS : public Observation
{
  private:
    string mstrFilename;
    int nvis_tables, nv2_tables, nt3_tables, nt4_tables;
    long nt4, nt3, nv2, nvis;
    int target_id; // the first target in the OIFITS file
  public:
    Obs_OIFITS(string filename);
    static vector <Observation*> ReadObservation_OIFITS(string filename);

    void WriteAuxTables(fitsfile* outfile, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target);
    void  WriteVis(fitsfile* outfile, UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    void WriteVis2(fitsfile* outfile,UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    void WriteT3(fitsfile* outfile,  UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
    void WriteT4(fitsfile* outfile,  UVPoint** uv_list, complex<double>** cvis,Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed);
};

#endif // OBS_OIFITS_H
