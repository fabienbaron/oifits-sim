/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include <cmath>

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Simulator.h"
#include "Station.h"
#include "Source.h"

Baseline::Baseline()
{
    this->xyz[0] = this->xyz[1] = this->xyz[2] = 0;
    this->name = "";
    this->indicies[0] = 0;
    this->indicies[1] = 0;
}

/// \todo Stop using NEU coordinates here, switch to station.x, station.y, and station.z.
Baseline::Baseline(Station * station1, Station * station2)
{
    // Now calculate the xyz coordinates:
    this->xyz[0] = station1->xyz[0] - station2->xyz[0];
    this->xyz[1] = station1->xyz[1] - station2->xyz[1];
    this->xyz[2] = station1->xyz[2] - station2->xyz[2];
        
    this->name = station1->GetName() + "-" + station2->GetName();
    this->indicies[0] = station1->GetIndex();
    this->indicies[1] = station2->GetIndex();
}

/// Computes the UV coordinates for the given information.
///     hour_angle : The hour angle (decimal hours)
///     source_dec : The declination of the source (radians)
///     wavenumber : Wavenumber (1/meter)
UVPoint Baseline::UVcoords(double hour_angle, double source_declination)
{
    // First convert all values into radians (they should be degrees or decimal hours (of time) before now.
    double h = hour_angle * PI / 12;
    double delta = source_declination;

    // Now compute the UV coordinates, again according to the APIS++ standards.
    UVPoint uv = UVPoint();
    uv.u = sin(h) * xyz[0]                    + cos(h) * xyz[1];
    uv.v = -sin(delta) * cos(h) * xyz[0]      + sin(delta) * sin(h) * xyz[1]    + cos(delta) * xyz[2];
    uv.w = cos(delta) * cos(h)*xyz[0]         - cos(delta) * sin(h) * xyz[1]    + sin(delta) * xyz[2];
    
    /// \note The coordinates (u,v) calculated here appear to be (-u, -v) with respect to CHARA+MIRC data.
    /// Given that the UV plane is symmetric, this doesn't matter.
    
    return uv;
}

string Baseline::GetName(void)
{
    return this->name;
}


/// Computes the complex visibility of the source at the specified UV point.
/// \todo Switch the hash lookup over to a multihash instead of using strings?
complex<double> Baseline::GetVisibility(Source & source, double hour_angle, double wavenumber)
{
    // First look up the UV coordinates
    UVPoint uv = UVcoords(hour_angle, source.declination);
    uv.Scale(wavenumber);    
    return GetVisibility(source, uv);
}

/// Returns the complex visibility of the given source
complex<double> Baseline::GetVisibility(Source & source, UVPoint uv)
{
    string hash_key = GetHashKey(source, uv);
    complex <double> visibility(0.0, 0.0);
    
    if(mVisValues.find(hash_key) != mVisValues.end())
    {
        visibility = mVisValues[hash_key];  
    }
    else
    {
        visibility = ComputeVisibility(source, uv);
        mVisValues[hash_key] = visibility;
    }
    
    return visibility;
}


// Computes a hash key from the source, hour angle, and wavenumber.
string  Baseline::GetHashKey(Source & source, UVPoint uv)
{
    /// \todo It may be necessary for the doubles coming into this function to be cast into some 
    /// finite floating point format.
    
    /// \todo This function is in common with the Triplet class, need to factor this code.
    
    std::ostringstream sstream;
    sstream << source.GetName() << '-' << uv.HashString();
    std::string str = sstream.str();
    return str;
}

/// Computes the error in the visibility
/// \todo Implement computations of error in visibility.  This function returns an error of zero for now.
double Baseline::ComputeVis2Error(Source & source, UVPoint uv)
{
    return 0.001;
}

/// Returns the error in the visibility assoicated with this source, hour angle, and wavenumber.
double Baseline::GetVis2Error(Source & source, double hour_angle, double wavenumber)
{
     // First look up the UV coordinates
    UVPoint uv = UVcoords(hour_angle, source.declination);
    uv.Scale(wavenumber);
    return GetVis2Error(source, uv);   
}

double Baseline::GetVis2Error(Source & source, UVPoint uv)
{
    string hash_key = GetHashKey(source, uv);
    double vis_error = 0.0;
    
    // First try looking up the value in the hash table:
    if(mVis2Errors.find(hash_key) != mVis2Errors.end())
    {
        vis_error = mVis2Errors[hash_key];
    }
    else
    {
        vis_error = ComputeVis2Error(source, uv);
        mVis2Errors[hash_key] = vis_error;
    }
    
    return vis_error;
}

/// Sets the visibility error for the given source, hour angle and wavenumber.  
/// Note: It is intended that this function is used to set the error from existing OIFITS data files.
void    Baseline::SetVis2Error(Source & source, double hour_angle, double wavenumber, double vis2error)
{
     UVPoint uv = UVcoords(hour_angle, source.declination);
     uv.Scale(wavenumber);
     SetVis2Error(source, uv, vis2error);   
}

/// Sets the visibility squared error based upon the source and UV coordiantes.
void    Baseline::SetVis2Error(Source & source, UVPoint uv, double vis2error)
{
    string hash_key = GetHashKey(source, uv);
    mVis2Errors[hash_key] = vis2error;
}

/// Computes the visibility at the specified UV point.
/// Note, uv should already be scaled by wavenumber.
complex<double> Baseline::ComputeVisibility(Source & source, UVPoint uv)
{
    complex <double> visibility(0.0, 0.0);

    // Point source
    if (source.source_image.GetCols() == 0)
    {
        visibility = 1.0;
    }
    else    // A resolved object
    {      
        int nx = source.source_image.GetRows();
        int ny = source.source_image.GetCols();    

        /// \todo This calculation could be farmed out to a GPU very easily.
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < ny; jj++)
            {
                visibility += source.source_image[ii][jj] * polar(1., -2.0 * PI * source.source_pixellation * milliarcsec * (uv.u * (double)ii + uv.v * (double)jj));
            }
        }
    }

    return visibility;
}

/// Returns an OIFITS oi_vis2_record for the given source, hour angle and wavenumber
//oi_vis2_record Baseline::GetVis2Record(Source & source, double hour_angle, <double> wavenumbers)
//{
//    // init local vars
//    complex<double> vis;
//    complex<double> vis_err;
//    double vis2;
//    double vis2_err;

//    // Finally assemble the oi_vis2_record:
//    oi_vis2_record record;
//    record.target_id = source.GetTargetID();
//    
//    /// \bug The time and MJD recorded here are nonsense values.
//    record.time = 0;
//    record.mjd = 0.0;
//    
//    /// \bug Integration time is set to 1-second by default regardless of instrument setting.
//    record.int_time = 1;
//    
//    // Get the UV point at the mid-point of the data.
//    UVPoint uv = UVcoords(hour_angle, source.declination, wavenumbers[wavenumbers.size()/2]);
//    
//    record.ucoord = uv.u;
//    record.vcoord = uv.v;
//    record.sta_index[0] = this->indicies[0];
//    record.sta_index[1] = this->indicies[1];

//    // Now write out the 
//    for(unsigned int i = 0; i < wavenumbers.size(); i++)
//    {
//        // First get the visibility and error.
//        vis = GetVisibility(source, hour_angle, wavenumbers[i]);
//        vis_err = GetVisError(source, hour_angle, wavenumbers[i]);
//        
//        // Now compute the powerspectra for this point:
//        vis2 = norm(vis);
//        vis2_err = norm(vis_err);

//        record.vis2data[i] = vis2;
//        record.vis2err[i] = vis2_err;
//        record.flag[i] = FALSE;
//    }

//    return record;
//}

// Returns the station ID of stations 0 or 1
int Baseline::GetStationID(int num)
{
    /// \todo This should probably query the telesocopes directly.
    return this->indicies[num];
}

double  Baseline::GetVis2(Source & source, double hour_angle, double wavenumber)
{
    complex<double> vis = this->GetVisibility(source, hour_angle, wavenumber);
    
    return norm(vis);
}

double Baseline::GetVis2(Source & source, UVPoint uv)
{
    complex<double> vis = this->GetVisibility(source, uv);
    return norm(vis);   
}

////////////////////////////////////////////////////////////////////
// Non Class Functions Below
////////////////////////////////////////////////////////////////////

/// Computes all possible baselines formed by the specified stations.
vector<Baseline> ComputeBaselines(vector<Station> & stations)
{
    int num_stations = stations.size();
    
    vector<Baseline> baselines;
    
    // Now compute all of the baselines and make a hash table for each baseline value
    for(int i = 0; i < num_stations; i++)
    {
        for(int j = i+1; j < num_stations; j++)
        {
            // Create a new baseline, append it to the list of baselines
            baselines.push_back(Baseline(&stations[i], &stations[j]));
        }
    }
    
    return baselines;
}

/// Computes a (baseline_name, baseline_object) hash table.
BaselineHash ComputeBaselineHash(vector<Baseline> & baselines)
{
    BaselineHash hash;
    
    for(unsigned int i = 0; i < baselines.size(); i++)
    {
        hash.insert(BaselineHash::value_type(baselines[i].GetName(), &baselines[i]) );
    }
    
    return hash;
}


/// \todo Rewrite this function to work with the new class definition.
//double Baseline::Geometric_OPD(double hour_angle, double source_declination,
//                               double array_latitude)
//{
//    double trad, drad, lrad;

//    trad = hour_angle * PI / (12.0 * 3600.);
//    drad = source_declination * PI / 180.0;
//    lrad = array_latitude * PI / 180.0;
//    Row < double >XYZ(3);

//    /// \todo This needs to be updated to use the z-coordinate of the telescopes.  See
//    /// functions above for more information.
//    XYZ[0] = -this->y * sin(lrad);
//    XYZ[1] = this->x;
//    XYZ[2] = this->y * cos(lrad);

//    double sin_alt, altitude, azimuth;

//    // Compute the altitude
//    sin_alt = cos(lrad) * cos(drad) * cos(trad) + sin(lrad) * sin(drad);
//    altitude = atan(sin_alt / sqrt(1 - sin_alt * sin_alt));

//    // Compute the azimuth
//    azimuth = atan(sin(trad) * cos(drad) / (cos(lrad) * sin(drad) - cos(trad)
//                                            * cos(drad) * sin(lrad)));

//    // Compute the unit vector
//    Row < double >s(3);

//    s[0] = cos(altitude) * cos(azimuth);
//    s[1] = cos(altitude) * sin(azimuth);
//    s[2] = sin(altitude);

//    // Return the OPD for the baseline
//    return s[0] * XYZ[0] + s[1] * XYZ[1] + s[2] * XYZ[2];
//}
