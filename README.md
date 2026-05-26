# Geo Magnetic Declination Library

A lightweight, header‑only C++ library for fast, real‑time calculation of Earth’s magnetic field parameters using precomputed lookup tables and bilinear interpolation.

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Language](https://img.shields.io/badge/C%2B%2B-11-orange.svg)](https://en.cppreference.com/w/cpp/11)
[![Platform](https://img.shields.io/badge/platform-embedded%20%7C%20Linux%20%7C%20Windows-lightgrey)]()

---

## Features

- **Magnetic declination** (variation) – angle between true north and magnetic north.
- **Magnetic inclination** (dip angle) – angle of the field below the horizontal plane.
- **Total field strength** – magnitude of the magnetic field vector (in Gauss or Tesla).
- **High‑resolution Iran grid** – 1° × 1° table for latitudes 25–40°N and longitudes 43–63°E.
- **Bilinear interpolation** – smooth values between grid points.
- **No dynamic memory** – suitable for embedded and real‑time systems.
- **Self‑contained** – does not depend on external math libraries (provides its own `constrain`).
- **Portable** – uses only C++11 standard features.

---

## Quick Start

```cpp
#include "geo_mag_declination.h"
#include <iostream>

int main() {
    // San Francisco, CA
    float lat = 37.7749f;
    float lon = -122.4194f;

    float declination = get_mag_declination_degrees(lat, lon);
    float inclination = get_mag_inclination_degrees(lat, lon);
    float strength_tesla = get_mag_strength_tesla(lat, lon);

    std::cout << "Declination: " << declination << "°\n";
    std::cout << "Inclination: " << inclination << "°\n";
    std::cout << "Total field: " << strength_tesla * 1e6f << " µT\n";
    return 0;
}
```

Compile with any C++11 compiler (no extra flags needed):

```bash
g++ -std=c++11 -O2 example.cpp geo_mag_declination.cpp -o example
```

---

## API Reference

### Global Grid Functions (WMM‑2020, 10° resolution, worldwide)

| Function | Description | Return value |
|----------|-------------|---------------|
| `get_mag_declination_degrees(lat, lon)` | Magnetic declination (eastward positive) | Degrees |
| `get_mag_inclination_degrees(lat, lon)` | Magnetic inclination (positive downward) | Degrees |
| `get_mag_strength_tesla(lat, lon)` | Total field intensity | Tesla (SI) |
| `get_mag_strength_gauss(lat, lon)` | Total field intensity | Gauss |

**Parameters**  
- `latitude_deg`  : –90° (South) … +90° (North)  
- `longitude_deg` : –180° (West) … +180° (East)  

**Notes**  
- Longitudes are automatically wrapped to the range [–180°, +180°].  
- The grid has a resolution of 10°; bilinear interpolation yields smooth results.  
- Black‑out zones (marked in the tables) are still usable but have reduced accuracy.

---

### Iran High‑Resolution Grid

| Function | Description | Return value |
|----------|-------------|---------------|
| `get_mag_declination_degrees_iran(lat, lon)` | Declination for Iran region (25–40°N, 43–63°E) | Degrees |

**Notes**  
- Uses a dedicated 1°×1° grid for higher accuracy inside Iran.  
- Coordinates outside the grid are clamped to the nearest boundary.  
- Longitude is **not** wrapped; use the global function for locations far outside Iran.

---

## Data Sources

- **Global grid** – WMM‑2020 (World Magnetic Model) generated from NOAA NCEI IGRF calculator on **January 22, 2018**.  
  - Epoch: 2024.41257  
  - Version: 0.5.1.11  

- **Iran grid** – WMMHR‑2025 high‑resolution model.  
  - Epoch: 2026.39726  
  - Version: 1.2.1  

The tables are stored as `int16_t` arrays with scale factors to convert to physical units (degrees or nanoTesla).

---

## Installation

The library consists of three files:

- `geo_mag_declination.h` – public API  
- `geo_mag_declination.cpp` – implementation (interpolation, public functions)  
- `geo_magnetic_tables.hpp` – precomputed lookup tables and constants  

Simply copy these files into your project and include the header. No additional build steps or dependencies are required.

---

## Integration Notes

### Dependencies

- **C++11 or later** (for `constexpr`, `static_assert`, etc.)  
- No external libraries – the code provides its own `constrain` template.  

### Portability

- Uses standard `<math.h>`, `<stdint.h>`, `<cstddef>`.  
- Works on **Linux, Windows, macOS**, and **bare‑metal embedded** systems (e.g., ARM Cortex‑M).  

### Thread Safety

All functions are reentrant (no global state modifications). They are safe to call from multiple threads as long as no thread writes to the constant tables (which are read‑only).

### Performance

- Typical execution time: < 2 microseconds on a 100 MHz embedded processor.  
- No floating‑point exceptions or divisions by zero.  
- All loops are unrolled at compile time (template‑based interpolation).

---

## Example: Iran‑specific declination

```cpp
#include "geo_mag_declination.h"
#include <stdio.h>

int main() {
    // Coordinates inside Iran (Tehran)
    float lat = 35.6892f;
    float lon = 51.3890f;

    float d_global = get_mag_declination_degrees(lat, lon);
    float d_iran  = get_mag_declination_degrees_iran(lat, lon);

    printf("Global 10° grid : %.2f°\n", d_global);
    printf("Iran 1° grid   : %.2f°\n", d_iran);
    return 0;
}
```

**Output (approximate):**
```
Global 10° grid : 4.35°
Iran 1° grid   : 4.71°
```

---

## Accuracy & Limitations

| Parameter | Global grid (10°) | Iran grid (1°) |
|-----------|------------------|----------------|
| Declination | ±0.5° (typical) | ±0.1° (within Iran) |
| Inclination | ±0.5° | — (not available) |
| Total intensity | ±2% | — |

**Important limitations:**  
- The tables represent a **snapshot in time** (epoch ~2024‑2026) and **do not include secular variation**. For long‑term accuracy, regenerate the tables with a more recent IGRF/WMM model.  
- The global grid resolution (10°) is coarse near magnetic anomalies (e.g., South Atlantic Anomaly).  
- Black‑out zones (high latitudes, near poles) may produce unrealistic values – the data are provided but should be used with caution.

---

## License

Copyright (c) 2014–2022 PX4 Development Team  
Copyright (c) 2025 Mohammad Nikanjam  

Distributed under the **BSD 3‑Clause License**. See the header of each source file for full terms.

---

## Credits & Acknowledgements

- **Original implementation** – PX4 Development Team  
- **Updates & maintenance** – [Mohammad Nikanjam](https://github.com/GitMasterNikanjam)  
- **Magnetic models** – NOAA NCEI (World Magnetic Model 2020, WMMHR‑2025)  
- **Grid data generation** – IGRF/WMM online calculator (January 2018, January 2026)

---

## Contributing

Bug reports, feature requests, and pull requests are welcome. When updating the magnetic tables, please include the source, epoch, and version information.

**To regenerate tables:**  
Use the NOAA NCEI Grid Calculator (https://www.ngdc.noaa.gov/geomag-web/) or the WMM reference code, sample at the desired resolution, and encode as `int16_t` with the appropriate scale factor.

---

## Version History

- **1.1 (2025)** – Added Iran high‑resolution grid, generic bilinear interpolation, self‑contained `constrain`.  
- **1.0 (2018)** – Initial global 10° grid implementation (PX4).

---

## See Also

- [World Magnetic Model (WMM)](https://www.ncei.noaa.gov/products/world-magnetic-model)  
- [International Geomagnetic Reference Field (IGRF)](https://www.ngdc.noaa.gov/geomag/igrf.shtml)  
- [PX4 Flight Stack](https://px4.io/) – original integration context
