"""
Library Features:

Name:          reflex_tools_flooding
Author(s):     Mauro Arcorace (mauro.arcorace@cimafoundation.org)
               Alessandro Masoero (alessandro.masoero@cimafoundation.org)
               Lorenzo Alfieri (lorenzo.alfieri@cimafoundation.org)
               Valerio Basso
               Giulia Bruno (giulia.bruno@cimafoundation.org)
               Alessia Matanò
               Andrea Libertino (andrea.libertino@cimafoundation.org)
Date:          '20230330'
Version:       '2.1.0'
"""
from scipy.optimize import minimize
from numpy.random import rand
import numpy as np
from functools import partial
import rioxarray as rxr
import logging
import math
import os

# ----------------------------------

# ----------------------------------
# Function to spread a flood depth over a hand map
def draw_flood_map(hand,h):
    return np.where(hand < h, h - hand, np.nan)
# ----------------------------------

# ----------------------------------
# Function to spread a flood depth over a hand map and calculate the spreaded volume
def optimize_flood_depth(hand,V_obj,areacell,h):
    depth = draw_flood_map(hand,h)
    V_hand = np.nansum(depth*areacell)
    return np.abs(V_obj-V_hand)
# ----------------------------------

# ----------------------------------
## Function to optimize the flooded volume in 2 steps
def optimise_volume(stream_id, optimise_setting, d):
    print(" --> Optimize volume for stream " + str(stream_id))

    streams_gdf = d["streams_gdf"]
    stream_row = streams_gdf[streams_gdf["stream"] == stream_id]

    hand_type = "ext_no_ce"
    c_exp = "not active"

    if optimise_setting["coastal_expansion_active"]:
        if optimise_setting["coastal_expansion_manual"]:
            if stream_row.manual_ce.values[0] == 1:
                hand_type = "ext_ce"
                c_exp = "active"
        elif stream_row.next_strea.values[0] <= 0 or stream_row.gradient.values[0] <= optimise_setting["coastal_expansion_gradient_limit"]:
            hand_type = "ext_ce"
            c_exp = "active"

    print(" ---> Coastal expansion mode: " + c_exp)
    hand_filename = optimise_setting["rst_hand"].format(mask_type=hand_type, stream_id=str(stream_id))
    hand_map = rxr.open_rasterio(hand_filename) / 100

    print(" ---> First attempt optimization...")

    r_min, r_max = optimise_setting["min_water_depth"], optimise_setting["max_water_depth"]
    h_first_attempt = r_min + rand(1) * (r_max - r_min)

    # declare partial function to be optimized
    opt_function = partial(optimize_flood_depth, hand_map.squeeze().values, stream_row["V0"].values, optimise_setting["areacell_m"])

    # perform the l-bfgs-b algorithm search and draw flood map
    h_opt0 = minimize(opt_function, h_first_attempt[0], bounds=[(r_min,r_max)], method='L-BFGS-B')["x"]
    flood_extent = draw_flood_map(hand_map, h_opt0)

    print(" ---> Compute transit time...")

    try:
        roughness_coeff = stream_row["roughness"].values
        if not roughness_coeff > 0:
            raise ValueError
        print(" ---> Stream roughness value provided!")
    except:
        roughness_coeff = optimise_setting["roughness_coeff"]
        print(" ---> Use roughness value from config!")

    # compute transit time
    avg_depth = np.nanmean(flood_extent)

    v_i = (1 / roughness_coeff) * (avg_depth ** (2 / 3)) * (stream_row["gradient"].values  ** (0.5))
    if v_i > 0:
        t_ti = int(stream_row["str_len_km"] * (10 ** 3)  / v_i)
    else:
        t_ti = +np.Inf

    # Transit time should be minor than concentration time
    if t_ti >= stream_row["tconc"].values:
        t_ti = stream_row["tconc"].values

    # Rescale runoff volume based on Tc-Tt (considers both the case of trapezoidal and triangular shapes due to the cut)
    if stream_row["Q_p_t"].values - t_ti * math.tan(stream_row["slope_hydrograph"].values) < 0:
        # Triangular
        vol_scaled = (stream_row["Q_p_t"].values ** 2) / (2 * math.tan(stream_row["slope_hydrograph"].values))
    else:
        # Trapezoidal
        vol_scaled = ((stream_row["Q_p_t"].values + (stream_row["Q_p_t"].values - (t_ti * math.tan(stream_row["slope_hydrograph"].values)))) * t_ti) / 2

    opt_function = partial(optimize_flood_depth, hand_map.squeeze().values, vol_scaled, optimise_setting["areacell_m"])
    f_opt = minimize(opt_function, h_opt0, bounds=[(r_min,r_max)], method='L-BFGS-B')

    h_opt = f_opt["x"]
    diff_opt = f_opt["fun"]

    logging.info(" ---> Final optimization...")
    flood_extent = draw_flood_map(hand_map, h_opt)

    flood_map = hand_map.copy()
    flood_map.values[0,:,:] = flood_extent

    flood_map.rio.to_raster(os.path.join(optimise_setting["out_path"], "tmp", "flood_m_bas_" + str(stream_id) + ".tif"))

    return np.array((stream_id, vol_scaled, v_i[0], diff_opt), dtype='object')
# ----------------------------------

# ----------------------------------
