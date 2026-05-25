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

#if defined(_WIN32)
	#include "../mathlib/mathlib.h"
#else
	#include <mathlib/mathlib.h>
#endif

#include <math.h>
#include <stdint.h>

using math::constrain;

// ###################################################################################

/**
 * @brief Convert a geographic coordinate to a lookup table index.
 *
 * @details This helper function maps a geographic coordinate value (latitude or
 *          longitude) to the corresponding index in the lookup table grid.
 *          It ensures the value stays within valid table bounds and accounts
 *          for the sampling resolution.
 *
 * @param[in,out] val Pointer to the coordinate value. Will be constrained to
 *                    the valid range [min, max - SAMPLING_RES] to ensure
 *                    (index + 1) stays within table bounds.
 * @param[in] min Minimum valid value for the coordinate (e.g., -90° for latitude).
 * @param[in] max Maximum valid value for the coordinate (e.g., +90° for latitude).
 *
 * @return Zero-based index into the lookup table for the given coordinate.
 *
 * @note The function modifies *val to the constrained value, which is useful
 *       for subsequent calculations.
 *
 * @see get_table_data for usage in interpolation
 */
static unsigned get_lookup_table_index(float *val, float min, float max)
{
	/* for the rare case of hitting the bounds exactly
	 * the rounding logic wouldn't fit, so enforce it.
	 */

	/* limit to table bounds - required for maxima even when table spans full globe range */
	/* limit to (table bounds - 1) because bilinear interpolation requires checking (index + 1) */
	*val = constrain(*val, min, max - SAMPLING_RES);

	return static_cast<unsigned>((-(min) + *val) / SAMPLING_RES);
}

/**
 * @brief Retrieve and interpolate magnetic field data from the lookup table.
 *
 * @details This is the core interpolation function that performs bilinear
 *          interpolation on the 2D lookup table to obtain smooth values at
 *          arbitrary geographic coordinates. It handles coordinate wrapping
 *          for longitude, clamps latitude to valid ranges, and uses the four
 *          nearest grid points for interpolation.
 *
 * @param latitude_deg Geographic latitude in degrees.
 *                     Automatically clamped to [-90°, +90°].
 * @param longitude_deg Geographic longitude in degrees.
 *                      Automatically wrapped to [-180°, +180°].
 * @param table[LAT_DIM][LON_DIM] Pointer to the lookup table to use
 *                                 (declination, inclination, or total intensity).
 *
 * @return Interpolated value in the table's native units (scaled integers).
 *
 * @note The return value is in raw table units and must be scaled by the
 *       appropriate factor for physical units (degrees, nanoTesla, etc.).
 *
 * @see get_mag_declination_degrees
 * @see get_mag_inclination_degrees
 * @see get_mag_strength_tesla
 *
 * @warning The function uses floor-based rounding, which is deterministic
 *          but may have edge effects near grid boundaries.
 */
static float get_table_data(float latitude_deg, float longitude_deg, const int16_t table[LAT_DIM][LON_DIM])
{
	/* Clamp latitude to valid range */
	latitude_deg = math::constrain(latitude_deg, SAMPLING_MIN_LAT, SAMPLING_MAX_LAT);

	/* Wrap longitude to [-180°, +180°] range */
	if (longitude_deg > SAMPLING_MAX_LON) {
		longitude_deg -= 360.f;
	}

	if (longitude_deg < SAMPLING_MIN_LON) {
		longitude_deg += 360.f;
	}

	/* round down to nearest sampling resolution */
	float min_lat = floorf(latitude_deg / SAMPLING_RES) * SAMPLING_RES;
	float min_lon = floorf(longitude_deg / SAMPLING_RES) * SAMPLING_RES;

	/* find index of nearest low sampling point */
	unsigned min_lat_index = get_lookup_table_index(&min_lat, SAMPLING_MIN_LAT, SAMPLING_MAX_LAT);
	unsigned min_lon_index = get_lookup_table_index(&min_lon, SAMPLING_MIN_LON, SAMPLING_MAX_LON);

	const float data_sw = table[min_lat_index][min_lon_index];
	const float data_se = table[min_lat_index][min_lon_index + 1];
	const float data_ne = table[min_lat_index + 1][min_lon_index + 1];
	const float data_nw = table[min_lat_index + 1][min_lon_index];

	/* perform bilinear interpolation on the four grid corners */
	const float lat_scale = constrain((latitude_deg - min_lat) / SAMPLING_RES, 0.f, 1.f);
	const float lon_scale = constrain((longitude_deg - min_lon) / SAMPLING_RES, 0.f, 1.f);

	const float data_min = lon_scale * (data_se - data_sw) + data_sw;
	const float data_max = lon_scale * (data_ne - data_nw) + data_nw;

	return lat_scale * (data_max - data_min) + data_min;
}

float get_mag_declination_degrees(float latitude_deg, float longitude_deg)
{
	// table stored as scaled degrees
	return get_table_data(latitude_deg, longitude_deg, declination_table) * WMM_DECLINATION_SCALE_TO_DEGREES;
}

float get_mag_inclination_degrees(float latitude_deg, float longitude_deg)
{
	// table stored as scaled degrees
	return get_table_data(latitude_deg, longitude_deg, inclination_table) * WMM_INCLINATION_SCALE_TO_DEGREES;
}

float get_mag_strength_gauss(float latitude_deg, float longitude_deg)
{
	// 1 Gauss = 1e4 Tesla
	return get_mag_strength_tesla(latitude_deg, longitude_deg) * 1e4f;
}

float get_mag_strength_tesla(float latitude_deg, float longitude_deg)
{
	// table stored as scaled nanotesla
	return get_table_data(latitude_deg, longitude_deg, totalintensity_table)
	       * WMM_TOTALINTENSITY_SCALE_TO_NANOTESLA * 1e-9f;
}
