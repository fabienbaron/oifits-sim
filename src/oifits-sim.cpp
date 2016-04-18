/*
* oifits-sim.cpp
*
*      Authors: bkloppenborg, fbaron
*/
#include "oifits-sim.h"

#include <complex>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Array.h"
#include "Combiner.h"
#include "NoiseModel.h"
#include "NoiseModel_Tatulli2006.h"
#include "Obs_HA.h"
#include "Obs_OIFITS.h"
#include "Observation.h"
#include "SpectralMode.h"
#include "Target.h"

extern "C" {
#include "random.h"
}

using namespace std;

// Prints out help describing the options on the command line
void PrintHelp() {
  string usage = "The OIFITS Simulator\n"
                 "Usage: \n"
                 " oifits-sim arguments \n"
                 " \n"
                 "Arguments: \n"
                 " -h     Prints this message \n"
                 " -t     Target definition file \n"
                 " -i     Input image \n"
                 " -o     Output OIFITS file \n"
                 " \n"
                 "There are two Simulation Options: \n"
                 "(1) from an existing OIFITS file, copying uncertainties from the real "
                 "data: \n"
                 " -a     The array used.\n"
                 " -m     Spectral Mode for the combiner (combiner must be specified "
                 "before -m).\n"
                 " -d     Input OIFITS file \n"
                 " \n"
                 "or (2) from a list of observations at an array with a specific combiner "
                 "in\n"
                 "which the noise is estimated as described in the documentation: \n"
                 " -a     The array used.\n"
                 " -c     The combiner.\n"
                 " -m     Spectral Mode for the combiner (combiner must be specified "
                 "before -m).\n"
                 " -obs   Observation Definition File \n"
                 " \n"
                 "Note: Some of the parameters also have double-dash overrides that \n"
                 "immediately follow other parameters that change default/averaged values "
                 "\n"
                 "that are found in configuration files.  See the documentation for "
                 "details.";

  cout << usage << "\n";
}

// The main routine.  Basically just used to parse out some parameters before
// handing
// things off to other functions.
int main(int argc, char *argv[]) {
  // Define characters that are used as comments:
  // TODO: Update the documentation to reflect these characters are all valid
  // comment chars.
  const string comment_chars("\\#~$&£%");

  // Define a few variables for the simulation.
  Target target;
  string output_fname;

  Array array;
  Combiner combiner;
  SpectralMode spec_mode;
  vector<Observation *> observation_list;
  bool oifits_mode = FALSE;
  bool target_info = FALSE;
  bool image_info = FALSE;
  bool spectral_info = false;
  bool output_info = false;
  bool array_info = false;
  bool combiner_info = false;
  bool obs_info = false;

  // TODO: For now we only have one noise model, so we load it by default
  NoiseModel *noisemodel = new NoiseModel_Tatulli2006();

  if (argc == 1)
    PrintHelp();

  cout << "Here" << endl;

  for (int i = 1; i < argc; i++) {
    // First see if the user is requesting help:
    if (strcmp(argv[i], "-h") == 0) {
      PrintHelp();
    }

    // We need to know some information about the target:
    if ((strcmp(argv[i], "-t") == 0) && (i < argc - 1)) {
      target.ImportFile(string(argv[i + 1]), comment_chars);
      target.ParseFileOptions(argv, i + 2, argc);
      target_info = true;
    }

    // And the image file we're going to use for the simulation:
    if ((strcmp(argv[i], "-i") == 0) && (i < argc - 1)) {
      target.SetImage(string(argv[i + 1]));
      target.ParseImageOptions(argv, i + 2, argc);
      image_info = true;
    }

    // And where we're going to write the output:
    if ((strcmp(argv[i], "-o") == 0) && (i < argc - 1)) {
      try {
        output_fname = string(argv[i + 1]);
        output_info = true;
      } catch (...) {
        throw std::runtime_error("Invalid Output File Definition");
      }
    }

    if ((strcmp(argv[i], "-d") == 0) && (i < argc - 1)) {
      observation_list = Obs_OIFITS::ReadObservation_OIFITS(string(argv[i + 1]));
      oifits_mode = true;
      printf("Command line interpreted as requiring OIFITS mode\n");
    }

    if ((strcmp(argv[i], "-a") == 0) && (i < argc - 1)) {
      array.ImportFile(string(argv[i + 1]), comment_chars);
      array.ParseOptions(argv, i, argc);
      array_info = true;
    }

    // the combiner
    if ((strcmp(argv[i], "-c") == 0) && (i < argc - 1)) {
      combiner.ImportFile(string(argv[i + 1]), comment_chars);
      combiner_info = true;
    }

    // the spectral mode.  Note, the combiner must be specified first so we can
    // check that
    // the two are indeed to be used with eachother.
    if ((strcmp(argv[i], "-m") == 0) && (i < argc - 1)) {
      if (combiner.name == "") {
        cout << "The combiner, -c, must be specified before the spectral mode, "
                "-m.\n";
        exit(0);
      }

      spec_mode.ImportFile(string(argv[i + 1]), combiner.GetName(), comment_chars);
      spectral_info = true;
    }

    // Now for observation_list
    if ((strcmp(argv[i], "-obs") == 0) && (i < argc - 1)) {
      // First ensure that the array has been defined.  If not, quit.
      if (array.GetArrayName() == "") {
        cout << "The array, -a, must be specified before the observation list, "
                "-obs.\n";
        exit(0);
      }

      // Observations are a little funny.  We permit both files and command line
      // options.
      // Just pass things off to the observation class so it can decide what to
      // do:
      observation_list = Observation::ParseCommandLine(&array, argv, i, argc, comment_chars);
      obs_info = true;
    }
  }

  if ((target_info && image_info && output_info)                                          // basic requirement for both modes
      && ((oifits_mode)                                                                   // OIFITS mode
          || (!oifits_mode && spectral_info && array_info && combiner_info && obs_info))) // HA or descriptive mode
    run_sim(&target, &array, &combiner, &spec_mode, noisemodel, observation_list, output_fname);
  else
    printf("Something is missing on the command line -- quitting !\n");

  // Clean up memory
  delete noisemodel;

  return 0;
}

void run_sim(Target *target, Array *array, Combiner *combiner, SpectralMode *spec, NoiseModel *noisemodel, vector<Observation *> observation_list,
             string output_filename) {
  // Setup openmp
  // omp_set_num_threads(omp_get_max_threads());

  // Pull up the random number generator.
  static Rand_t random_seed;

  // See if the output file aready exists:
  // TODO: this could be bad/unexpected behavior.
  if (FileExists(output_filename)) {
    cout << "Warning, output file aready exists!  Overwriting..." << endl;
    output_filename = "!" + output_filename;
  }

  // Open up the OIFITS file.
  fitsfile *outfile;
  int status = 0;
  fits_create_file(&outfile, output_filename.c_str(), &status);
  if (status) {
    fits_report_error(stderr, status);
    return;
  }

  // Now compute the vis2 records and t3s:
  int n_observations = observation_list.size();

  cout << "Simulating N Observations: " << observation_list.size() << endl;

  for (unsigned int i = 0; i < n_observations; i++) {
    Observation *observation = observation_list.back();
    ObsType type = observation->GetObsType();
    UVPoint *uv_list = NULL;      // shared between vis, v2, t3, t4...
    complex<double> *cvis = NULL; // shared between vis, v2, t3, t4...

    // Do a dymamic cast to get the subclass object back
    if (type == HOUR_ANGLE || type == DESCRIPTIVE) {
      Obs_HA *observation = dynamic_cast<Obs_HA *>(observation_list.back());
      printf("Obs HA: %d / %d - Hour Angle: %3.3f \n", i + 1, n_observations, observation->GetHA(target->right_ascension));
      // SETUP HELPER OIFITS TABLES
      if (i == 0) {
        oi_array oi_arr = array->GetOIArray();
        write_oi_array(outfile, oi_arr, 1, &status);
        free_oi_array(&oi_arr);

        oi_target oi_targ = target->GetOITarget();
        write_oi_target(outfile, oi_targ, &status);
        free_oi_target(&oi_targ);

        oi_wavelength oi_wave = spec->GetOIWavelength();
        write_oi_wavelength(outfile, oi_wave, 1, &status);
        free_oi_wavelength(&oi_wave);
        vector<double> wavenumbers = spec->GetWavenumbers();
      }
    } else //(type == OIFITS)
    {
      printf("Obs OIFITS: %d / %d \n", i + 1, n_observations);
      Obs_OIFITS *observation = dynamic_cast<Obs_OIFITS *>(observation_list.back());

      // SETUP HELPER OIFITS TABLES (COPYING THE ORIGINAL ONES)
      observation->WriteAuxTables(outfile, array, combiner, spec, target);
    }

    if (observation->HasVIS())
      observation->WriteVis(outfile, &uv_list, cvis, array, combiner, spec, target, noisemodel, random_seed);
    
    if (observation->HasV2())
      observation->WriteVis2(outfile, &uv_list, cvis, array, combiner, spec, target, noisemodel, random_seed);

    if (observation->HasT3())
      observation->WriteT3(outfile, &uv_list, cvis, array, combiner, spec, target, noisemodel, random_seed);

    //  if (observation->HasT4())
    //      observation->WriteT4(fprt, &uv_list, &cvis, array, combiner, spec,
    //      target, noisemodel, random_seed);

    // All done with this observation object.  Pop it off the vector and free
    // memory.
    observation_list.pop_back();
    delete observation;
  }

  if (status) {
    fits_delete_file(outfile, &status);
    return;
  } else {
    fits_close_file(outfile, &status);
  }
}
