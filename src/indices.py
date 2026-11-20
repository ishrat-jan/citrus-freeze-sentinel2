"""Spectral indices from Sentinel-2 surface reflectance.

Band names follow Earth Search: blue=B02, green=B03, red=B04, rededge1=B05,
nir=B08, nir08=B8A, swir16=B11.

    NDVI  = (B8 - B4) / (B8 + B4)                          greenness (main damage index)
    EVI   = 2.5 (B8 - B4) / (B8 + 6 B4 - 7.5 B2 + 1)       greenness, less saturated in dense canopy
    NDRE  = (B8A - B5) / (B8A + B5)                        red-edge: chlorophyll content
    NDWI  = (B8A - B11) / (B8A + B11)   (Gao 1996)         leaf water content
    MCARI = ((B5 - B4) - 0.2 (B5 - B3)) * (B5 / B4)         chlorophyll absorption

Freeze-killed leaves lose chlorophyll and dry out, so all five are expected to drop.
"""

BANDS = {
    "NDVI": ["nir", "red"],
    "EVI": ["nir", "red", "blue"],
    "NDRE": ["nir08", "rededge1"],
    "NDWI": ["nir08", "swir16"],
    "MCARI": ["rededge1", "red", "green"],
}

FUNCTIONS = {
    "NDVI": lambda b: (b["nir"] - b["red"]) / (b["nir"] + b["red"]),
    "EVI": lambda b: 2.5 * (b["nir"] - b["red"]) / (b["nir"] + 6 * b["red"] - 7.5 * b["blue"] + 1),
    "NDRE": lambda b: (b["nir08"] - b["rededge1"]) / (b["nir08"] + b["rededge1"]),
    "NDWI": lambda b: (b["nir08"] - b["swir16"]) / (b["nir08"] + b["swir16"]),
    "MCARI": lambda b: ((b["rededge1"] - b["red"]) - 0.2 * (b["rededge1"] - b["green"]))
                       * (b["rededge1"] / b["red"]),
}


def bands_for(names):
    """All reflectance bands needed for a list of index names."""
    return sorted({b for n in names for b in BANDS[n]})


def compute(name, bands):
    """One index (float32 DataArray) from a dict of reflectance bands."""
    return FUNCTIONS[name](bands).astype("float32").rename(name)
