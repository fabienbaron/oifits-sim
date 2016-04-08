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
        // Make a place to store the data as we are simulating it.
        // Note, the vis2_data vector doesn't know how to allocate oi_vis2_record* entries so we have
        // to do this manually below.
        vector<oi_vis2_record*> vis2_data;
        oi_vis2 vis2_table; // imported data
        //oi_vis2_record input_record;
        UVPoint uv;
        Baseline * baseline;
        oi_wavelength wave;
        double wavelength, dwavelength;
        bool valid_v2;
        long nv2_valid = 0;
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
                        oi_vis2_record * output = (oi_vis2_record*) malloc(sizeof(oi_vis2_record));
                        baseline = array->GetBaseline((vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
                        printf("%ld \n", i);
                        if(baseline == NULL) printf("Station Indexes: %i %i \n", (vis2_table.record[i]).sta_index[0], (vis2_table.record[i]).sta_index[1]);
                        output->target_id = 1;
                        output->time = (vis2_table.record[i]).time;
                        output->mjd = (vis2_table.record[i]).mjd;
                        output->int_time = (vis2_table.record[i]).int_time;
                        output->ucoord = (vis2_table.record[i]).ucoord;
                        output->vcoord = (vis2_table.record[i]).vcoord;
                        output->sta_index[0] = (vis2_table.record[i]).sta_index[0];
                        output->sta_index[1] = (vis2_table.record[i]).sta_index[1];
                        // pre-allocate memory for the vis2data, vis2error, and flag
                        // because some points may be NaN/invalid, the actual final size may be smaller
                        output->vis2data = (double *) malloc(vis2_table.nwave * sizeof(double));
                        output->vis2err = (double *) malloc(vis2_table.nwave * sizeof(double));
                        output->flag = (BOOL *) malloc(vis2_table.nwave * sizeof(BOOL));
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
                                        output->vis2data[j] = (vis2_table.record[i]).vis2data[j]; //baseline->GetVis2(*target, uv, wavelength, dwavelength) + v2_err * Rangauss(random_seed);
                                        output->vis2err[j] = (vis2_table.record[i]).vis2err[j];
                                        output->flag[j] = (vis2_table.record[i]).flag[j];
                                        ++nv2_valid;
                                }
                                else
                                {
                                        printf("OIFITS import -- V2 point invalid, Record: %ld\t Wav: %ld Flag: %d V2: %f V2err: %f \n", i, j, (vis2_table.record[i]).flag[j], (vis2_table.record[i]).vis2data[j], (vis2_table.record[i]).vis2err[j]);
                                        output->vis2data[j] = nan(""); //baseline->GetVis2(*target, uv, wavelength, dwavelength) + v2_err * Rangauss(random_seed);
                                        output->vis2err[j] = nan("");
                                        output->flag[j] = true;
                                }
                        }
                        // Now append the data to the output vector
                        vis2_data.push_back(output);
                }

                free_oi_wavelength(&wave);
                free_oi_vis2(&vis2_table);
        }
        fits_close_file(fptr, &status);
        fits_close_file(fptr2, &status2);
        printf("Total number of V2: %ld \t Valid V2: %ld", nv2, nv2_valid);

        /// \bug The observation date is set to all zeros by default.
        /// This is to ensure the user knows this is simulated data, but may not be compliant
        /// with the OIFITS format, or good "note taking"
        // TBD: copy this from original file
        string arrname = "";
        string insname = "";
        oi_vis2 * outvis2 = (oi_vis2*) malloc(sizeof(oi_vis2));
        outvis2->revision = 1;
        strncpy(outvis2->date_obs, "2014-01-01", 11);
        strncpy(outvis2->arrname, arrname.c_str(), FLEN_VALUE);
        strncpy(outvis2->insname, insname.c_str(), FLEN_VALUE);
        outvis2->numrec = nv2;
        outvis2->record = (oi_vis2_record *) malloc(nv2 * sizeof(oi_vis2_record));
        for(long i = 0; i < nv2; i++)
        {
                outvis2->record[i] = *vis2_data.back();
                // free memory and pop
                free(vis2_data.back());
                vis2_data.pop_back();
        }
        return *outvis2;
}

/// \todo This is a really hacky solution.  Try to find a better solution.
oi_t3 Obs_OIFITS::GetT3(UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
        // a place to store the data as we are simulating it.
        // Notice, the t3_data vector doesn't know how to deallocate the oi_t3_records
        // so we will need to manually deallocate memory below.
        vector<oi_t3_record*> t3_data;

        // init some local vars:
        int nwave = int(spec_mode->mean_wavenumber.size());
        fitsfile * fptr;
        int status = 0;

        double amp_err;
        double phi_err;
        complex<double> bis;

        // Open the input file as read-only data.
        fits_open_file(&fptr, this->mstrFilename.c_str(), READONLY, &status);
        if(status)
                throw std::runtime_error("Could not read OIFITS file.");

        // Now iterate through the OI_VIS2 tables, mirroring the sampling on the Source object
        // at the specified wavenumbers;
        oi_t3 t3;
        UVPoint uv1;
        UVPoint uv2;
        UVPoint uv3;
        oi_t3_record input_record;
        Triplet * triplet;
        double wavelength, dwavelength;

        do
        {
                // Read in the next vis2 table
                read_next_oi_t3(fptr, &t3, &status);
                if(status)
                        break;

                for(int record_id = 0; record_id < t3.numrec; record_id++)
                {
                        // Use a local var to store some information
                        input_record = t3.record[record_id];
                        oi_t3_record * output = (oi_t3_record*) malloc(sizeof(oi_t3_record));

                        // Get the triplet
                        triplet = array->GetTriplet(input_record.sta_index[0], input_record.sta_index[1], input_record.sta_index[2]);

                        // Copy some information over to the output record:
                        /// \bug Target is set to 1 by default
                        output->target_id = 1;
                        output->time = input_record.time;
                        output->mjd = input_record.mjd;
                        output->int_time = input_record.int_time;
                        output->u1coord = input_record.u1coord;
                        output->v1coord = input_record.v1coord;
                        output->u2coord = input_record.u2coord;
                        output->v2coord = input_record.v2coord;
                        output->sta_index[0] = input_record.sta_index[0];
                        output->sta_index[1] = input_record.sta_index[1];
                        output->sta_index[2] = input_record.sta_index[2];

                        // Allocate memory for the vis2data, vis2error, and flag:
                        output->t3amp = (double *) malloc(nwave * sizeof(double));
                        output->t3amperr = (double *) malloc(nwave * sizeof(double));
                        output->t3phi = (double *) malloc(nwave * sizeof(double));
                        output->t3phierr = (double *) malloc(nwave * sizeof(double));
                        output->flag = (char *) malloc(nwave * sizeof(char));

                        // Now iterate over the wavenumbers
                        for(int j = 0; j < nwave; j++)
                        {
                                // Reset the UV coordinates
                                uv1.u = input_record.u1coord;
                                uv1.v = input_record.v1coord;
                                uv2.u = input_record.u2coord;
                                uv2.v = input_record.v2coord;
                                uv3.u = uv1.u + uv2.u;
                                uv3.v = uv1.v + uv2.v;

                                // Scale them
                                wavelength = spec_mode->mean_wavelength[j];
                                dwavelength = spec_mode->delta_wavelength[j];

                                uv1.Scale(1./wavelength);
                                uv2.Scale(1./wavelength);
                                uv3.Scale(1./wavelength);

                                // Look up the error values:
                                amp_err = input_record.t3amperr[j];
                                phi_err = input_record.t3phierr[j];

                                // Compute the bispectra
                                bis = triplet->GetT3(*target, uv1, uv2, uv3, wavelength, dwavelength);

                                // Simulate the bispectrum's amplitude based on the source image
                                output->t3amp[j] = abs(bis) + amp_err * Rangauss(random_seed);
                                // Copy the error from the input file.
                                output->t3amperr[j] = amp_err;

                                // Simulate the bispectrum's phase based on the source image.
                                // Remember, phi_err is already in degrees
                                output->t3phi[j] = (arg(bis) * 180 / PI) + phi_err * Rangauss(random_seed);
                                // Copy the error from the input file
                                output->t3phierr[j] = phi_err;

                                // Lastly set the "ignore this data" flag to false
                                output->flag[j] = FALSE;
                        }

                        // Now append the data to the output vector
                        t3_data.push_back(output);

                }

        } while(status == 0);

        // close the file
        fits_close_file(fptr, &status);

        // Now convert the t3_data vector into a properly formatted OI_T3 table.
        oi_t3 * outt3 = (oi_t3*) malloc(sizeof(oi_t3));
        int ndata = int(t3_data.size());
        string arrname = array->GetArrayName();

        outt3->revision = 1;

        strncpy(outt3->date_obs, "2014-01-01", 11);
        strncpy(outt3->arrname, arrname.c_str(), FLEN_VALUE);
        strncpy(outt3->insname, spec_mode->spec_mode.c_str(), FLEN_VALUE);
        outt3->numrec = ndata;
        outt3->nwave = nwave;

        outt3->record = (oi_t3_record *) malloc(ndata * sizeof(oi_t3_record));
        for(int i = 0; i < ndata; i++)
        {
                outt3->record[i] = *t3_data.back();

                // Free memory and pop
                free(t3_data.back());
                t3_data.pop_back();
        }

        // Note, this memory object contains pointers and should be freed by the OIFITSLIB free_oi_t3.
        return *outt3;
}

oi_t4 Obs_OIFITS::GetT4(UVPoint** uv_list, complex<double>** cvis, Array * array, Combiner * combiner, SpectralMode * spec_mode, Target * target, NoiseModel * noisemodel, Rand_t random_seed)
{
        //  oi_t4* dummy;
        //  return *dummy;

}
