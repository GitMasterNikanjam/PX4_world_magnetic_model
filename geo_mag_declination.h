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
 * @file geo_mag_declination.h
 * @brief Geomagnetic field calculator for Earth's magnetic field parameters.
 *
 * @details This module provides functions to compute Earth's magnetic field
 *          characteristics including declination, inclination, and total field
 *          strength at any geographic location. The calculations are based on
 *          precomputed lookup tables generated from the IGRF (International
 *          Geomagnetic Reference Field) model.
 *
 *          The implementation uses bilinear interpolation on a fixed-resolution
 *          grid, making it suitable for real-time embedded systems where
 *          computational efficiency is critical. The data represents the
 *          geomagnetic field for a specific epoch and does not include
 *          secular variation.
 *
 * @note The current data was generated from the NOAA NCEI IGRF calculator
 *       on January 22, 2018. For long-term accuracy, the underlying tables
 *       should be periodically updated with newer IGRF models.
 *
 * @see https://www.ngdc.noaa.gov/geomag-web/#igrfgrid for data source
 * @see wmm_2020_grid.hpp for the WMM-2020 model grid implementation
 *
 * @author PX4 Development Team (original implementation)
 * @author Mohammad Nikanjam (updates and maintenance)
 *
 * @version 1.0
 * @date 2018-2025
 */

// ##############################################################################

#pragma once

// ##############################################################################

/**
 * @brief Calculate the magnetic declination at a given geographic location.
 *
 * @details Magnetic declination is the angle between geographic true north and
 *          magnetic north, measured eastward (positive) or westward (negative).
 *          This is also known as magnetic variation.
 *
 * @param latitude_deg  Geographic latitude in degrees.
 *                      Range: -90° (South) to +90° (North).
 * @param longitude_deg Geographic longitude in degrees.
 *                      Range: -180° (West) to +180° (East).
 *
 * @return Magnetic declination in degrees.
 *         - Positive values indicate eastward declination (magnetic north is
 *           east of true north).
 *         - Negative values indicate westward declination (magnetic north is
 *           west of true north).
 *
 * @note The function automatically wraps longitude values to the valid range.
 * @note Uses bilinear interpolation on a 10-degree resolution grid.
 *
 * @see get_mag_inclination_degrees
 * @see get_mag_strength_gauss
 *
 * @example
 * @code
 *   // Get declination for San Francisco, CA
 *   float lat = 37.7749f;
 *   float lon = -122.4194f;
 *   float declination = get_mag_declination_degrees(lat, lon);
 *   // declination ≈ 13.5° (eastward)
 * @endcode
 */
float get_mag_declination_degrees(float latitude_deg, float longitude_deg);

/**
 * @brief Calculate magnetic declination for Iran region with higher resolution.
 *
 * @details Uses a dedicated 1° resolution grid covering latitude 25°–40° and
 *          longitude 43°–63°. This provides more accurate results inside Iran
 *          than the global 10° grid.
 *
 * @param latitude_deg  Latitude in degrees (should be within 25–40 for best results).
 * @param longitude_deg Longitude in degrees (should be within 43–63).
 *
 * @return Magnetic declination in degrees.
 *
 * @note Input coordinates outside the Iran grid will be clamped to the grid
 *       boundaries. For locations far outside, use the global function.
 */
float get_mag_declination_degrees_iran(float latitude_deg, float longitude_deg);

/**
 * @brief Calculate the magnetic inclination (dip angle) at a given geographic location.
 *
 * @details Magnetic inclination is the angle between the Earth's magnetic field
 *          vector and the horizontal plane. Positive values indicate the field
 *          points downward into the Earth, which occurs in the Northern Hemisphere.
 *          Negative values indicate upward pointing, typical of the Southern Hemisphere.
 *          This parameter is essential for tilt compensation in magnetometers.
 *
 * @param latitude_deg  Geographic latitude in degrees.
 *                      Range: -90° (South) to +90° (North).
 * @param longitude_deg Geographic longitude in degrees.
 *                      Range: -180° (West) to +180° (East).
 *
 * @return Magnetic inclination (dip angle) in degrees.
 *         - Positive values: Field lines point downward (Northern Hemisphere)
 *         - Negative values: Field lines point upward (Southern Hemisphere)
 *         - Zero occurs near the magnetic equator
 *
 * @note The function automatically wraps longitude values to the valid range.
 * @note Accuracy is approximately ±0.5 degrees for most locations.
 *
 * @see get_mag_declination_degrees
 * @see get_mag_strength_gauss
 *
 * @example
 * @code
 *   // Get inclination for London, UK
 *   float lat = 51.5074f;
 *   float lon = -0.1278f;
 *   float inclination = get_mag_inclination_degrees(lat, lon);
 *   // inclination ≈ 67° (downward)
 * @endcode
 */
float get_mag_inclination_degrees(float latitude_deg, float longitude_deg);

/**
 * @brief Calculate the total magnetic field strength at a given geographic location.
 *
 * @details This function returns the magnitude of the Earth's magnetic field
 *          vector (total intensity) expressed in Gauss units. The value represents
 *          the combined effect of all magnetic field components.
 *
 * @param latitude_deg  Geographic latitude in degrees.
 *                      Range: -90° (South) to +90° (North).
 * @param longitude_deg Geographic longitude in degrees.
 *                      Range: -180° (West) to +180° (East).
 *
 * @return Total magnetic field strength in Gauss.
 *         Range: Typically 0.25 to 0.65 Gauss on Earth's surface.
 *         - Minimum occurs at South Atlantic Anomaly (≈ 0.25 G)
 *         - Maximum occurs near magnetic poles (≈ 0.65 G)
 *
 * @note 1 Gauss = 100,000 nanoTesla (nT)
 * @note The function automatically wraps longitude values to the valid range.
 *
 * @warning Gauss is a non-SI unit. For SI units, use @ref get_mag_strength_tesla.
 *
 * @see get_mag_strength_tesla for SI unit version
 * @see get_mag_declination_degrees
 * @see get_mag_inclination_degrees
 *
 * @example
 * @code
 *   // Get field strength for equatorial region
 *   float lat = 0.0f;
 *   float lon = 0.0f;
 *   float strength_gauss = get_mag_strength_gauss(lat, lon);
 *   // strength_gauss ≈ 0.31 G (31,000 nT)
 * @endcode
 */
float get_mag_strength_gauss(float latitude_deg, float longitude_deg);

/**
 * @brief Calculate the total magnetic field strength in SI units (Tesla).
 *
 * @details This function returns the magnitude of the Earth's magnetic field
 *          vector (total intensity) expressed in Tesla units. This is the SI
 *          version of @ref get_mag_strength_gauss.
 *
 * @param latitude_deg  Geographic latitude in degrees.
 *                      Range: -90° (South) to +90° (North).
 * @param longitude_deg Geographic longitude in degrees.
 *                      Range: -180° (West) to +180° (East).
 *
 * @return Total magnetic field strength in Tesla (T).
 *         Range: Approximately 25 to 65 microTesla (µT) on Earth's surface.
 *         - Minimum: ≈ 25 µT (South Atlantic Anomaly)
 *         - Maximum: ≈ 65 µT (magnetic poles)
 *
 * @note 1 Tesla = 10,000 Gauss
 * @note The Earth's magnetic field is typically 25-65 µT, which is very weak
 *       compared to typical permanent magnets (0.1-1 T).
 *
 * @see get_mag_strength_gauss for the non-SI unit version
 * @see get_mag_declination_degrees
 * @see get_mag_inclination_degrees
 *
 * @example
 * @code
 *   // Get field strength for use with SI units
 *   float lat = 40.7128f;  // New York
 *   float lon = -74.0060f;
 *   float strength_tesla = get_mag_strength_tesla(lat, lon);
 *   float strength_microtesla = strength_tesla * 1e6f;
 *   // strength_microtesla ≈ 52 µT
 * @endcode
 */
float get_mag_strength_tesla(float latitude_deg, float longitude_deg);
