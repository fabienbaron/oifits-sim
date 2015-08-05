/// \file Baseline.cpp
/// Implements the required functions for the baseline class.

#include <cmath>

#include "Baseline.h"
#include "UVPoint.h"
#include "Array.h"
#include "Common.h"
#include "Station.h"
#include "Target.h"

Baseline::Baseline()
{
    this->xyz[0] = this->xyz[1] = this->xyz[2] = 0;
    this->name = "";
    this->indicies[0] = 0;
    this->indicies[1] = 0;
}

/// Computes the XYZ position of a baseline composed of two stations
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
///     target_dec : The declination of the target (radians)
///     wavenumber : Wavenumber (1/meter)
UVPoint Baseline::UVcoords(double hour_angle, double target_declination)
{
    // First convert all values into radians (they should be degrees or decimal hours (of time) before now.
    double h = hour_angle * PI / 12;
    double delta = target_declination;

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


// Computes a hash key from the target, hour angle, and wavenumber.
string  Baseline::GetHashKey(Target & target, UVPoint uv)
{
    /// \todo It may be necessary for the doubles coming into this function to be cast into some
    /// finite floating point format.

    /// \todo This function is in common with the Triplet class, need to factor this code.

    std::ostringstream sstream;
    sstream << target.GetName() << '-' << uv.HashString();
    std::string str = sstream.str();
    return str;
}



/// Computes the complex visibility of the target at the specified UV point.
/// \todo Switch the hash lookup over to a multihash instead of using strings?
complex<double> Baseline::GetVisibility(Target & target, double hour_angle, double wavelength, double dwavelength)
{
    // First look up the UV coordinates
    UVPoint uv = UVcoords(hour_angle, target.declination);
    uv.Scale(1./wavelength);
    return GetVisibility(target, uv, wavelength, dwavelength);
}

/// Returns the complex visibility of the given target
complex<double> Baseline::GetVisibility(Target & target, UVPoint uv, double wavelength, double dwavelength)
{
    string hash_key = GetHashKey(target, uv);
    complex <double> visibility(0.0, 0.0);

    if(mVisValues.find(hash_key) != mVisValues.end())
    {
        visibility = mVisValues[hash_key];
    }
    else
    {
        visibility = ComputeVisibility(target, uv, wavelength, dwavelength);
        mVisValues[hash_key] = visibility;
    }

    return visibility;
}

double  Baseline::GetVis2(Target & target, double hour_angle, double wavelength, double dwavelength)
{
    complex<double> vis = this->GetVisibility(target, hour_angle, wavelength, dwavelength);
    return norm(vis);
}

double Baseline::GetVis2(Target & target, UVPoint uv, double wavelength, double dwavelength)
{
  // note: getvisibility will return a precomputed complex visibility i
  complex<double> vis = this->GetVisibility(target, uv, double wavelength, double dwavelength);
  return norm(vis);
}

/// Computes the visibility at the specified UV point.
/// Note, uv should already be scaled by wavenumber.
complex<double> Baseline::ComputeVisibility(Target & target, UVPoint uv, double wavelength, double dwavelength)
{
    complex <double> visibility(0.0, 0.0);

    // Point target
    if (target.image.GetCols() == 0)
    {
        visibility = 1.0;
    }
    else    // A resolved object
    {
        int nx = target.image.GetRows();
        int ny = target.image.GetCols();
        for (int ii = 0; ii < nx; ii++)
        {
            for (int jj = 0; jj < ny; jj++)
            {
            visibility += target.image[ii][jj]
                        * sinc( dwavelength/wavelength * (uv.u * (double)(ii - nx / 2 ) - uv.v * (double)(jj - ny / 2)) * target.pixellation * milliarcsec)
                        * polar(1., 2.0 * PI * target.pixellation * milliarcsec * (uv.u * (double)(ii - nx / 2 ) - uv.v * (double)(jj - ny / 2)));
            }
        }
    }

    return visibility;
}

// Returns the station ID of stations 0 or 1
int Baseline::GetStationID(int num)
{
    /// \todo This should probably query the telesocopes directly.
    return this->indicies[num];
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
            printf("Creating baseline %i %i. \n", stations[i].GetIndex(), stations[j].GetIndex());

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
