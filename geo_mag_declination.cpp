/****************************************************************************
 *
 *   Copyright (c) 2014-2022 PX4 Development Team. All rights reserved.
 *   Copyright (c) 2025 Mohammad Nikanjam (https://github.com/GitMasterNikanjam). All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

/**
 * @file geo_mag_declination.cpp
 * @brief Implementation of Earth's magnetic field calculations using lookup tables.
 *
 * @details This file implements the geophysical calculations for Earth's magnetic
 *          field parameters (declination, inclination, and total field strength)
 *          using precomputed lookup tables with bilinear interpolation.
 *
 *          The implementation uses a fixed-resolution grid (10° spacing) from
 *          geo_magnetic_tables.hpp and performs bilinear interpolation to obtain
 *          smooth values at arbitrary geographic coordinates. This approach
 *          balances computational efficiency with reasonable accuracy for
 *          real-time embedded applications.
 *
 * @note The lookup table data was generated from the NOAA NCEI IGRF calculator
 *       on January 22, 2018. The current grid resolution is coarse (1° degree
 *       spacing in the original data, but interpolated on 10° grid).
 *
 * @todo Improve lookup table resolution to support more precise calculations.
 *       The current coarse resolution (full degrees) limits accuracy,
 *       particularly in regions with high magnetic field gradients.
 *
 * @see geo_mag_declination.h for the public API
 * @see geo_magnetic_tables.hpp for the underlying grid data
 *
 * @author PX4 Development Team
 * @author Mohammad Nikanjam
 */

// ######################################################################################

#include "geo_mag_declination.h"
#include "geo_magnetic_tables.hpp"

#include <math.h>
#include <stdint.h>
#include <cstddef>

// ###################################################################################

template<typename _Tp>
constexpr _Tp constrain(_Tp val, _Tp min_val, _Tp max_val)
{
	return (val < min_val) ? min_val : ((val > max_val) ? max_val : val);
}

// -----------------------------------------------------------------------------
// Helper functions for global grid (10° resolution, -90..90 lat, -180..180 lon)
// -----------------------------------------------------------------------------

/**
 * @brief Compute table index for a coordinate value on the global grid.
 *
 * @param val        Coordinate value (modified to the constrained value used for interpolation).
 * @param min        Minimum valid coordinate (e.g., -90 for latitude).
 * @param max        Maximum valid coordinate (e.g., +90 for latitude).
 * @param res        Grid resolution (SAMPLING_RES).
 * @return           Zero‑based index into the table.
 */
static unsigned get_table_index(float *val, float min, float max, float res)
{
    /* Clamp to ensure (index + 1) stays within table bounds */
    *val = constrain(*val, min, max - res);
    return static_cast<unsigned>((-(min) + *val) / res);
}

/**
 * @brief Generic bilinear interpolation for any regular grid.
 *
 * @tparam LAT_DIM   Number of latitude samples.
 * @tparam LON_DIM   Number of longitude samples.
 * @param latitude_deg   Input latitude (will be clamped).
 * @param longitude_deg  Input longitude (will be wrapped to [-180,180] if global).
 * @param table          The 2D lookup table (stored as int16_t).
 * @param lat_min        Minimum latitude of the grid.
 * @param lat_max        Maximum latitude of the grid.
 * @param lon_min        Minimum longitude of the grid.
 * @param lon_max        Maximum longitude of the grid.
 * @param res            Grid resolution (degrees).
 * @param wrap_lon       If true, longitudes are wrapped to [lon_min, lon_max] (for global grids).
 * @return               Interpolated value (still in raw table units, not scaled).
 */
template<size_t LAT_DIM, size_t LON_DIM>
static float interpolate_table(float latitude_deg, float longitude_deg,
                               const int16_t (&table)[LAT_DIM][LON_DIM],
                               float lat_min, float lat_max,
                               float lon_min, float lon_max,
                               float res, bool wrap_lon)
{
    // Clamp latitude to grid range
    latitude_deg = constrain(latitude_deg, lat_min, lat_max);

    // Longitude handling
    if (wrap_lon) {
        // Wrap to [-180, 180] range (assumes lon_min = -180, lon_max = 180)
        if (longitude_deg > lon_max) longitude_deg -= 360.0f;
        if (longitude_deg < lon_min) longitude_deg += 360.0f;
    }
    // For non‑wrapping grids (e.g., Iran) we simply clamp.
    longitude_deg = constrain(longitude_deg, lon_min, lon_max);

    // Find the lower grid point (floor)
    float min_lat = floorf(latitude_deg / res) * res;
    float min_lon = floorf(longitude_deg / res) * res;

    // Obtain indices using the grid‑specific function (clamps automatically)
    unsigned min_lat_idx = get_table_index(&min_lat, lat_min, lat_max, res);
    unsigned min_lon_idx = get_table_index(&min_lon, lon_min, lon_max, res);

    // Fetch the four surrounding grid values
    const float data_sw = table[min_lat_idx][min_lon_idx];
    const float data_se = table[min_lat_idx][min_lon_idx + 1];
    const float data_ne = table[min_lat_idx + 1][min_lon_idx + 1];
    const float data_nw = table[min_lat_idx + 1][min_lon_idx];

    // Interpolation factors (clamped to avoid rounding issues)
    const float lat_scale = constrain((latitude_deg - min_lat) / res, 0.0f, 1.0f);
    const float lon_scale = constrain((longitude_deg - min_lon) / res, 0.0f, 1.0f);

    // Bilinear interpolation
    const float data_min = lon_scale * (data_se - data_sw) + data_sw;
    const float data_max = lon_scale * (data_ne - data_nw) + data_nw;

    return lat_scale * (data_max - data_min) + data_min;
}

// -----------------------------------------------------------------------------
// Public API functions (global grid)
// -----------------------------------------------------------------------------

float get_mag_declination_degrees(float latitude_deg, float longitude_deg)
{
    return interpolate_table(latitude_deg, longitude_deg, declination_table,
                             SAMPLING_MIN_LAT, SAMPLING_MAX_LAT,
                             SAMPLING_MIN_LON, SAMPLING_MAX_LON,
                             SAMPLING_RES, true) * WMM_DECLINATION_SCALE_TO_DEGREES;
}

float get_mag_inclination_degrees(float latitude_deg, float longitude_deg)
{
    return interpolate_table(latitude_deg, longitude_deg, inclination_table,
                             SAMPLING_MIN_LAT, SAMPLING_MAX_LAT,
                             SAMPLING_MIN_LON, SAMPLING_MAX_LON,
                             SAMPLING_RES, true) * WMM_INCLINATION_SCALE_TO_DEGREES;
}

float get_mag_strength_tesla(float latitude_deg, float longitude_deg)
{
    // table stored as scaled nanotesla → convert to Tesla
    return interpolate_table(latitude_deg, longitude_deg, totalintensity_table,
                             SAMPLING_MIN_LAT, SAMPLING_MAX_LAT,
                             SAMPLING_MIN_LON, SAMPLING_MAX_LON,
                             SAMPLING_RES, true) * WMM_TOTALINTENSITY_SCALE_TO_NANOTESLA * 1e-9f;
}

float get_mag_strength_gauss(float latitude_deg, float longitude_deg)
{
    // 1 Gauss = 1e-4 Tesla
    return get_mag_strength_tesla(latitude_deg, longitude_deg) * 1e4f;
}

// -----------------------------------------------------------------------------
// Iran‑specific function (1° grid, latitude 25–40, longitude 43–63)
// -----------------------------------------------------------------------------

float get_mag_declination_degrees_iran(float latitude_deg, float longitude_deg)
{
    // Use the same generic interpolator with the Iran grid constants.
    // Note: longitude is NOT wrapped because the Iran grid is a small continuous block.
    return interpolate_table(latitude_deg, longitude_deg, declination_table_iran,
                             SAMPLING_MIN_LAT_IRAN, SAMPLING_MAX_LAT_IRAN,
                             SAMPLING_MIN_LON_IRAN, SAMPLING_MAX_LON_IRAN,
                             SAMPLING_RES_IRAN, false) * WMM_DECLINATION_SCALE_TO_DEGREES;
}