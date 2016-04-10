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


  // get the ID of the first target in the file
  oi_target targets;
  fits_open_file(&fptr, filename.c_str(), READONLY, &status);
  read_oi_target(fptr, &targets, &status);
  if(targets.ntarget > 1)
  {
    throw std::runtime_error("Can only use single target files.");
  }
  else
  {
    target_id = (targets.targ[0]).target_id;
  }
  fits_close_file(fptr, &status);
  free_oi_target(&targets);


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

  mbHasVIS =(nvis_tables > 0);
  mbHasV2 = (nv2_tables > 0);
  mbHasT3 = false; //(nt3_tables > 0);
  mbHasT4 = false; //( nt4_tables > 0);
}


oi_vis Obs_OIFITS::GetVis(UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{

}

/// \todo This is a really hacky solution.  Come up with a better method.  Perhaps use
/// the Baseline::GetOI_Vis2_record routine?
oi_vis2 Obs_OIFITS::GetVis2(UVPoint** uv_list, complex<double>** cvis,  Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  oi_vis2 vis2_table; // to store imported/input data
  oi_vis2 * outvis2 = (oi_vis2*) malloc(sizeof(oi_vis2));
  outvis2->revision = 1;
  strncpy(outvis2->date_obs, "2014-01-01", 11);
  strncpy(outvis2->arrname, "FAKE", FLEN_VALUE);
  strncpy(outvis2->insname, "SIM", FLEN_VALUE);
  outvis2->numrec = 0;
  outvis2->record = (oi_vis2_record *) malloc(nv2 * sizeof(oi_vis2_record));

  UVPoint uv;
  Baseline * baseline;
  oi_wavelength wave;
  double wavelength, dwavelength;
  bool valid_v2;
  long nv2_valid = 0;
  long irecord =0;
  int status = 0;
  int status2 = 0;
  fitsfile *fptr;
  fitsfile *fptr2;
  fits_open_file(&fptr, mstrFilename.c_str(), READONLY, &status);
  fits_open_file(&fptr2, mstrFilename.c_str(), READONLY, &status2);

  for(int k=0;k<nv2_tables; k++)
  {
    printf("V2 TABLE %d\n", k);
    read_next_oi_vis2(fptr, &vis2_table, &status);
    read_oi_wavelength(fptr2, vis2_table.insname, &wave, &status2);
    for (long i = 0; i < vis2_table.numrec; i++)
    {
      baseline = array->GetBaseline((vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
      printf("i:%ld irecord: %ld\n", i, irecord);
      if(baseline == NULL) printf("Station Indexes: %i %i \n", (vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
      outvis2->record[irecord].target_id = 1;
      outvis2->record[irecord].time = (vis2_table.record[i]).time;
      outvis2->record[irecord].mjd = (vis2_table.record[i]).mjd;
      outvis2->record[irecord].int_time = (vis2_table.record[i]).int_time;
      outvis2->record[irecord].ucoord = (vis2_table.record[i]).ucoord;
      outvis2->record[irecord].vcoord = (vis2_table.record[i]).vcoord;
      outvis2->record[irecord].sta_index[0] = (vis2_table.record[i]).sta_index[0];
      outvis2->record[irecord].sta_index[1] = (vis2_table.record[i]).sta_index[1];
      // pre-allocate memory for the vis2data, vis2error, and flag
      // because some points may be NaN/invalid, the actual final size may be smaller
      outvis2->record[irecord].vis2data = (double *) malloc(vis2_table.nwave * sizeof(double));
      outvis2->record[irecord].vis2err = (double *) malloc(vis2_table.nwave * sizeof(double));
      outvis2->record[irecord].flag = (BOOL *) malloc(vis2_table.nwave * sizeof(BOOL));
      outvis2->nwave = vis2_table.nwave;
      outvis2->numrec +=1;
      printf("Nwave: %d\n", vis2_table.nwave);
      // TODO: think about nwave
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
          outvis2->record[irecord].vis2data[j] = (vis2_table.record[i]).vis2data[j]; //baseline->GetVis2(*target, uv, wavelength, dwavelength) + v2_err * Rangauss(random_seed);
          outvis2->record[irecord].vis2err[j] = (vis2_table.record[i]).vis2err[j];
          outvis2->record[irecord].flag[j] = false;
          ++nv2_valid;
        }
        else
        {
          //printf("OIFITS import -- V2 point invalid, Record: %ld\t Wav: %ld Flag: %d V2: %f V2err: %f \n", i, j, (vis2_table.record[i]).flag[j], (vis2_table.record[i]).vis2data[j], (vis2_table.record[i]).vis2err[j]);
          outvis2->record[irecord].vis2data[j] = nan(""); //baseline->GetVis2(*target, uv, wavelength, dwavelength) + v2_err * Rangauss(random_seed);
          outvis2->record[irecord].vis2err[j] = nan("");
          outvis2->record[irecord].flag[j] = true;
        }
      }
      ++irecord;
    }
    free_oi_wavelength(&wave);
    free_oi_vis2(&vis2_table);
  }
  fits_close_file(fptr, &status);
  fits_close_file(fptr2, &status2);
  printf("Total number of V2: %ld \t Valid V2: %ld", nv2, nv2_valid);
  fflush(stdout);
  return *outvis2;
}


oi_t3 Obs_OIFITS::GetT3(UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  oi_t3 t3_table;
  oi_t3 * outt3 = (oi_t3*) malloc(sizeof(oi_t3));
  outt3->revision = 1;
  strncpy(outt3->date_obs, "2014-01-01", 11);
  strncpy(outt3->arrname, "TODO", FLEN_VALUE);
  strncpy(outt3->insname, "TODO", FLEN_VALUE);
  outt3->numrec = nt3;
  outt3->record = (oi_t3_record *) malloc(nt3 * sizeof(oi_t3_record));

  UVPoint uv1;
  UVPoint uv2;
  UVPoint uv3;
  Triplet * triplet;
  //  Baseline * baseline;
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
    printf("T3 TABLE %d\n", k);
    read_next_oi_t3(fptr, &t3_table, &status);
    read_oi_wavelength(fptr2, t3_table.insname, &wave, &status2);
    for (long i = 0; i < t3_table.numrec; i++)
    {
      // Get the triplet
      triplet = array->GetTriplet((t3_table.record[i]).sta_index[0],
      (t3_table.record[i]).sta_index[1], (t3_table.record[i]).sta_index[2]);
      outt3->record[irecord].target_id = 1;
      outt3->record[irecord].time = (t3_table.record[i]).time;
      outt3->record[irecord].mjd = (t3_table.record[i]).mjd;
      outt3->record[irecord].int_time = (t3_table.record[i]).int_time;
      outt3->record[irecord].u1coord = (t3_table.record[i]).u1coord;
      outt3->record[irecord].v1coord = (t3_table.record[i]).v1coord;
      outt3->record[irecord].u2coord = (t3_table.record[i]).u2coord;
      outt3->record[irecord].v2coord = (t3_table.record[i]).v2coord;
      outt3->record[irecord].sta_index[0] = (t3_table.record[i]).sta_index[0];
      outt3->record[irecord].sta_index[1] = (t3_table.record[i]).sta_index[1];
      outt3->record[irecord].sta_index[2] = (t3_table.record[i]).sta_index[2];

      // Allocate memory for the vis2data, vis2error, and flag:
      outt3->record[irecord].t3amp = (double *) malloc(t3_table.nwave * sizeof(double));
      outt3->record[irecord].t3amperr = (double *) malloc(t3_table.nwave * sizeof(double));
      outt3->record[irecord].t3phi = (double *) malloc(t3_table.nwave * sizeof(double));
      outt3->record[irecord].t3phierr = (double *) malloc(t3_table.nwave * sizeof(double));
      outt3->record[irecord].flag = (char *) malloc(t3_table.nwave * sizeof(char));

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
        //outt3->record[irecord].t3amp[j] = abs(bis) + amp_err * Rangauss(random_seed);
        //outt3->record[irecord].t3phi[j] = (arg(bis) * 180 / PI) + phi_err * Rangauss(random_seed);
        outt3->record[irecord].flag[j] = !( (valid_t3amp == TRUE) || (valid_t3amp == TRUE)  );
        if(outt3->record[irecord].flag[j]  == FALSE) nt3_valid++;

        // TODO: take into account original flags. Also do this for other observables

        if (valid_t3amp == TRUE)
        {
          outt3->record[irecord].t3amp[j] = (t3_table.record[i]).t3amperr[j];
          outt3->record[irecord].t3amperr[j] = (t3_table.record[i]).t3amperr[j];
        }
        else
        {
          outt3->record[irecord].t3amp[j] = nan("");
          outt3->record[irecord].t3amperr[j] = nan("");
        }

        if (valid_t3phi == TRUE)
        {
          outt3->record[irecord].t3phi[j] = (t3_table.record[i]).t3phi[j];
          outt3->record[irecord].t3phierr[j] = (t3_table.record[i]).t3phierr[j];
        }
        else
        {
          outt3->record[irecord].t3phi[j] = nan("");
          outt3->record[irecord].t3phierr[j] = nan("");
        }

      }
      ++irecord;
    }
    free_oi_wavelength(&wave);
    free_oi_t3(&t3_table);
  }
  fits_close_file(fptr, &status);
  fits_close_file(fptr2, &status2);
  // Note, this memory object contains pointers and should be freed by the OIFITSLIB free_oi_t3.
  return *outt3;
}

oi_t4 Obs_OIFITS::GetT4(UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
  //  oi_t4* dummy;
  //  return *dummy;

}
