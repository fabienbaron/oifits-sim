#include "Obs_OIFITS.h"

#include "fitsio.h"
#include "UVPoint.h"
#include "Baseline.h"
#include "Array.h"
#include "Common.h"
#include "SpectralMode.h"

/// Reads an OIFITS data file and creates a series of observation based upon the data.
/// Note, presently this function only supports ONE array, combiner, spectral mode per OIFITS file.
vector <Observation*> Obs_OIFITS::ReadObservation_OIFITS(string filename)
{
  vector<Observation*> observations;

  // We don't attempt to read anything from the OIFITS file, just set the filename.
  observations.push_back(new Obs_OIFITS(filename) );

  return observations;
}

Obs_OIFITS::Obs_OIFITS(string filename)
{
  mObsType = OIFITS;
  this->mstrFilename = filename;

  // attempt to read OIFITS input file
  fitsfile *fptr;
  int status = 0;

  nvis_tables = 0; nv2_tables =0; nt3_tables =0; nt4_tables=0;
  nt4 =0; nt3 =0; nv2 =0; nvis=0;
  oi_vis vis_table;
  oi_vis2 vis2_table;
  oi_t3 t3_table;
  oi_t4 t4_table;

  // Vis tables
  status = 0;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status) throw std::runtime_error("Could not read OIFITS file.");

  read_next_oi_vis(fptr, &vis_table, &status);
  while (status == 0)
  {
    nvis += vis_table.numrec * vis_table.nwave;
    nvis_tables++;
    //printf("Vis Table\t %i\t Points %ld Channels %d\n", nvis_tables, nvis, vis_table.nwave);
    free_oi_vis(&vis_table);
    read_next_oi_vis(fptr, &vis_table, &status);
  }
  fits_close_file(fptr, &status);

  printf("OIFITS import -- OI_VIS  \tTables %i\t Entries %ld\n", nvis_tables, nvis);

  // VIS2 tables
  status = 0;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status) throw std::runtime_error("Could not read OIFITS file.");
  read_next_oi_vis2(fptr, &vis2_table, &status);
  while (status == 0)
  {
    nv2 += vis2_table.numrec * vis2_table.nwave;
    nv2_tables++;
    free_oi_vis2(&vis2_table);
    read_next_oi_vis2(fptr, &vis2_table, &status);
  }

  fits_close_file(fptr, &status);

  printf("OIFITS import -- OI_VIS2 \tTables %i\t Entries %ld\n", nv2_tables, nv2);

  // T3 tables
  status = 0;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status) throw std::runtime_error("Could not read OIFITS file.");
  read_next_oi_t3(fptr, &t3_table, &status);
  while (status == 0)
  {
    nt3 += t3_table.numrec * t3_table.nwave;
    nt3_tables++;
    //printf("T3 Table %i. Total T3 Points %ld\n",nt3_tables, nt3);
    free_oi_t3(&t3_table);
    read_next_oi_t3(fptr, &t3_table, &status);
  }
  fits_close_file(fptr, &status);

  printf("OIFITS import -- OI_T3   \tTables %i\t Entries %ld\n", nt3_tables, nt3);

  // T4 tables
  status = 0;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  if(status) throw std::runtime_error("Could not read OIFITS file.");
  read_next_oi_t4(fptr, &t4_table, &status);
  while (status == 0)
  {
    nt4 += t4_table.numrec * t4_table.nwave;
    nt4_tables++;
    //printf("T4 Table %i. Total T4 Points %ld\n",nt4_tables, nt4);
    free_oi_t4(&t4_table);
    read_next_oi_t4(fptr, &t4_table, &status);
  }
  fits_close_file(fptr, &status);

  mbHasVIS = false; //(nvis_tables > 0);
  mbHasV2 = (nv2_tables > 0);
  mbHasT3 = (nt3_tables > 0);
  mbHasT4 = false; //( nt4_tables > 0);
}


void Obs_OIFITS::WriteVis(fitsfile* outfile,UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{

}

/// \todo This is a really hacky solution.  Come up with a better method.  Perhaps use
/// the Baseline::GetOI_Vis2_record routine?
void Obs_OIFITS::WriteVis2(fitsfile* outfile, UVPoint** uv_list, complex<double>** cvis,  Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  // TODO
  // Temporary code -- needs to be moved to Obs_OIFITS to a location where input and output filenames are known
// also, we only read ONE array file
    int status = 0, status2 =0;
    fitsfile *infile;
    fits_open_file(&infile, mstrFilename.c_str(), READONLY, &status);
    oi_array oi_arr;
    char* arrname;
    read_next_oi_array(infile, &oi_arr, &status);
    write_oi_array(outfile, oi_arr, 1, &status);
    free_oi_array(&oi_arr);
    fits_close_file(infile, &status);


    UVPoint uv;
    Baseline * baseline;
    oi_wavelength wave;
    double wavelength, dwavelength;
    bool valid_v2;
    long nv2_valid = 0;

    status = 0;
    status2 = 0;
    fitsfile *fptr;
    fitsfile *fptr2;
    fits_open_file(&fptr, mstrFilename.c_str(), READONLY, &status);
    fits_open_file(&fptr2, mstrFilename.c_str(), READONLY, &status2);
    oi_vis2 vis2_table;
    for(int k=0;k<nv2_tables; k++)
    {
    read_next_oi_vis2(fptr, &vis2_table, &status);
    read_oi_wavelength(fptr2, vis2_table.insname, &wave, &status2);
    printf("V2 TABLE %d \t\tARRAY: %20s \t\t INSTRUMENT: %20s\n", k, vis2_table.arrname, vis2_table.insname);
    for (long i = 0; i < vis2_table.numrec; i++)
    {
      //baseline = array->GetBaseline((vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
      //if(baseline == NULL) printf("Station Indexes: %i %i \n", (vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
      for (long j = 0; j < vis2_table.nwave; j++)
      {
        uv.u = (vis2_table.record[i]).ucoord;
        uv.v = (vis2_table.record[i]).vcoord;
        wavelength = wave.eff_wave[j];
        dwavelength = wave.eff_band[j];
        uv.Scale(1./wavelength);

        valid_v2 =  !(((vis2_table.record[i]).flag[j] != 0)
        || isnan((vis2_table.record[i]).vis2data[j])
        || isnan((vis2_table.record[i]).vis2err[j])
        || isinf((vis2_table.record[i]).vis2data[j])
        || isinf((vis2_table.record[i]).vis2err[j])
        || ((vis2_table.record[i]).vis2data[j] == DOUBLENULLVALUE)
        || ((vis2_table.record[i]).vis2data[j] == FLOATNULLVALUE)
        || ((vis2_table.record[i]).vis2err[j] == DOUBLENULLVALUE)
        || ((vis2_table.record[i]).vis2err[j] == FLOATNULLVALUE)
        || (vis2_table.record[i]).vis2err[j] <= 0);

        if (valid_v2 == TRUE)
        {
          //baseline->GetVis2(*target, uv, wavelength, dwavelength) + v2_err * Rangauss(random_seed);
          (vis2_table.record[i]).vis2data[j] = 0.12345;
          ++nv2_valid;
        }
        else
        {
          // Force-flag it
          (vis2_table.record[i]).vis2data[j] = nan("");
          (vis2_table.record[i]).vis2err[j] = nan("");
          (vis2_table.record[i]).flag[j] = true;
        }
      }
    }
    // Table finished -- now write it to output file
    write_oi_wavelength(outfile, wave, 1, &status);
    write_oi_vis2(outfile, vis2_table, 1, &status);
    free_oi_wavelength(&wave);
    free_oi_vis2(&vis2_table);
  }
  fits_close_file(fptr, &status);
  fits_close_file(fptr2, &status2);
  printf("Total number of V2: %ld \t Valid V2: %ld\n", nv2, nv2_valid);
  fflush(stdout);
}


void Obs_OIFITS::WriteT3(fitsfile* outfile, UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  oi_t3 t3_table;

  UVPoint uv1;
  UVPoint uv2;
  UVPoint uv3;
  Triplet * triplet;
  oi_wavelength wave;
  double wavelength, dwavelength;
  bool valid_t3amp, valid_t3phi;
  long nt3_valid = 0;
  long irecord =0;
  int status = 0;
  int status2 = 0;
  fitsfile *fptr;
  fitsfile *fptr2;
  fits_open_file(&fptr, mstrFilename.c_str(), READONLY, &status);
  fits_open_file(&fptr2, mstrFilename.c_str(), READONLY, &status2);

  for(int k=0;k<nt3_tables; k++)
  {
    read_next_oi_t3(fptr, &t3_table, &status);
    read_oi_wavelength(fptr2, t3_table.insname, &wave, &status2);
    printf("T3 TABLE %d \t\tARRAY: %20s \t\tINSTRUMENT: %20s\n", k, t3_table.arrname, t3_table.insname);

    for (long i = 0; i < t3_table.numrec; i++)
    {
      // Get the triplet
      //triplet = array->GetTriplet((t3_table.record[i]).sta_index[0], (t3_table.record[i]).sta_index[1], (t3_table.record[i]).sta_index[2]);

      // Allocate memory for the vis2data, vis2error, and flag:
      for (long j = 0; j < t3_table.nwave; j++)
      {
        // Reset the UV coordinates
        uv1.u = (t3_table.record[i]).u1coord;
        uv1.v = (t3_table.record[i]).v1coord;
        uv2.u = (t3_table.record[i]).u2coord;
        uv2.v = (t3_table.record[i]).v2coord;
        uv3.u = uv1.u + uv2.u;
        uv3.v = uv1.v + uv2.v;

        wavelength = wave.eff_wave[j];
        dwavelength = wave.eff_band[j];
        uv1.Scale(1./wavelength);
        uv2.Scale(1./wavelength);
        uv3.Scale(1./wavelength);

        valid_t3amp = !(((t3_table.record[i]).t3amperr[j] <= 0)
        || isnan((t3_table.record[i]).t3amp[j])
        || isnan((t3_table.record[i]).t3amperr[j])
        || isinf((t3_table.record[i]).t3amp[j])
        || isinf((t3_table.record[i]).t3amperr[j])
        || ((t3_table.record[i]).t3amp[j] == DOUBLENULLVALUE)
        || ((t3_table.record[i]).t3amp[j] == FLOATNULLVALUE)
        || ((t3_table.record[i]).t3amperr[j] == DOUBLENULLVALUE)
        || ((t3_table.record[i]).t3amperr[j] == FLOATNULLVALUE)
        || ((t3_table.record[i]).flag[j] != 0));

        valid_t3phi =  !(((t3_table.record[i]).t3phierr[j] <= 0)
        || isnan((t3_table.record[i]).t3phi[j])
        || isnan((t3_table.record[i]).t3phierr[j])
        || isinf((t3_table.record[i]).t3phi[j])
        || isinf((t3_table.record[i]).t3phierr[j])
        || ((t3_table.record[i]).t3phi[j] == DOUBLENULLVALUE)
        || ((t3_table.record[i]).t3phi[j] == FLOATNULLVALUE)
        || ((t3_table.record[i]).t3phierr[j] == DOUBLENULLVALUE)
        || ((t3_table.record[i]).t3phierr[j] == FLOATNULLVALUE)
        || ((t3_table.record[i]).flag[j] != 0));

        // Compute the bispectra
        //  bis = triplet->GetT3(*target, uv1, uv2, uv3, wavelength, dwavelength);
        (t3_table.record[i]).flag[j] = ! ( ( (t3_table.record[i]).flag[j] == FALSE) && (valid_t3amp == TRUE) && (valid_t3amp == TRUE) );
        if((t3_table.record[i]).flag[j]  == FALSE)
          nt3_valid++;
        if (valid_t3amp == TRUE)
        {
          (t3_table.record[i]).t3amp[j] = 0.12324;
        }
        else
        {
          (t3_table.record[i]).t3amp[j] = nan("");
          (t3_table.record[i]).t3amperr[j] = nan("");
        }

        if (valid_t3phi == TRUE)
        {
          (t3_table.record[i]).t3phi[j] = 1.23456;

        }
        else
        {
          (t3_table.record[i]).t3phi[j] = nan("");
          (t3_table.record[i]).t3phierr[j] = nan("");

        }

      }
    }
    write_oi_wavelength(outfile, wave, 1, &status);
    write_oi_t3(outfile, t3_table, 1, &status);
    free_oi_wavelength(&wave);
    free_oi_t3(&t3_table);
  }
  fits_close_file(fptr, &status);
  fits_close_file(fptr2, &status2);
}

void Obs_OIFITS::WriteT4(fitsfile* outfile, UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  //  oi_t4* dummy;
  //  return *dummy;

}
