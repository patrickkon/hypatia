# MIT License
#
# Copyright (c) 2020 Debopam Bhattacherjee
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import math
import ephem
import pandas as pd
import os

try:
    from . import util
except (ImportError, SystemError):
    import util

# Visualizes paths between endpoints at specific time instances

EARTH_RADIUS = 6378135.0 # WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html

# CONSTELLATION GENERATION GENERAL CONSTANTS
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EPOCH = "2000-01-01 00:00:00" # WARNING NOTE: if modified epoch here, make sure to also modify epoch within tle_line1 of def generate_tles_from_scratch_manual. This is because the SAT IDs that we obtain from here, will be injected, as SAT IDs into the main_helper.py chain of functions, specifically in helper_dynamic_state.py. In short, we want the SAT IDs to match.

# ---------------------------------------------------------------------------------------------------------------

# SECTION: CONSTELLATION SPECIFIC PARAMETERS


# # STARLINK 550
# NAME = "starlink_550"

# ################################################################
# # The below constants are taken from Starlink's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MOD-20190830-00087
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = 940700 # From https://fcc.report/IBFS/SAT-MOD-20181108-00083/1569860.pdf (minimum angle of elevation: 25 deg)
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2)) # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# NUM_ORBS = 72
# NUM_SATS_PER_ORB = 22
# INCLINATION_DEGREE = 53

# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# path_file = "../../paper/satgenpy_analysis/data/starlink_550_isls_plus_grid_ground_stations_top_100_algorithm_free_one_only_over_isls/100ms_for_200s/manual/data/networkx_path_1608_to_1650.txt"



# # Starlink 560

# # GENERATION CONSTANTS
# NAME = "starlink_560"

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.04  # Altitude ~1015 km
# ALTITUDE_M = 560000  # Altitude ~560 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 6
# NUM_SATS_PER_ORB = 58
# INCLINATION_DEGREE = 97.6


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# path_file = "../../../../hypatia_plus/data/starlink_560/networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)

# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)



# KUIPER 630
# NAME = "kuiper_630"

# ################################################################
# # The below constants are taken from Kuiper's FCC filing as below:
# # [1]: https://www.itu.int/ITU-R/space/asreceived/Publication/DisplayPublication/8716
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 14.80  # Altitude ~630 km
# ALTITUDE_M = 630000  # Altitude ~630 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(30.0))  # Considering an elevation angle of 30 degrees; possible values [1]: 20(min)/30/35/45
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))  # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# NUM_ORBS = 34
# NUM_SATS_PER_ORB = 34
# INCLINATION_DEGREE = 51.9


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# path_file = "../../../../hypatia_plus/data/kuiper_630/networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)

# with open(path_file, 'w') as out:
#     s = '''0,1180-200-300-400-1177
# 11700000000000,1608-202-201-203-180-224-1605'''
#     out.write(s)



# #TELESAT 1015

# # GENERATION CONSTANTS
# NAME = "telesat_1015"

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 13.66  # Altitude ~1015 km
# ALTITUDE_M = 1015000  # Altitude ~1015 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(10.0))  # Considering an elevation angle of 10 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 27
# NUM_SATS_PER_ORB = 13
# INCLINATION_DEGREE = 98.98


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# path_file = "../../../../hypatia_plus/data/telesat_1015/networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)

# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)


# Custom - derived from starlink_550 and starlink_560

# # GENERATION CONSTANTS
# NAME = "starlink_550_72_72" # 550 altitude, 72 planes with 72 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 72
# NUM_SATS_PER_ORB = 72
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# City IDs are available in the city_detail_file.
# If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)

# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)


# # GENERATION CONSTANTS
# NAME = "starlink_550_48_48" # 550 altitude, 48 planes with 48 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 48
# NUM_SATS_PER_ORB = 48
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)


# # GENERATION CONSTANTS
# NAME = "starlink_550_100_100" # 550 altitude, 100 planes with 100 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 100
# NUM_SATS_PER_ORB = 100
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)

# # GENERATION CONSTANTS
# NAME = "starlink_550_130_130" # 550 altitude, 130 planes with 130 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 130
# NUM_SATS_PER_ORB = 130
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)

# # GENERATION CONSTANTS
# NAME = "starlink_550_150_150" # 550 altitude, 150 planes with 150 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 150
# NUM_SATS_PER_ORB = 150
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)

# # GENERATION CONSTANTS
# NAME = "starlink_550_200_200" # 550 altitude, 200 planes with 200 sats each

# ################################################################
# # The below constants are taken from Telesat's FCC filing as below:
# # [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
# ################################################################

# MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
# ALTITUDE_M = 550000  # Altitude ~550 km
# SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
# MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# # ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
# MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
# NUM_ORBS = 200
# NUM_SATS_PER_ORB = 200
# INCLINATION_DEGREE = 53


# # Input file; Generated during simulation
# # Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# # City IDs are available in the city_detail_file.
# # If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# # then offset ID is 1584 + 24 = 1608.
# num_sats = NUM_ORBS*NUM_SATS_PER_ORB
# gs_paris_id = num_sats + 24
# gs_moscow_id = num_sats + 21
# gs_src_id = gs_paris_id
# gs_dst_id = gs_moscow_id
# parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
# isExist = os.path.exists(parent_dir)
# if not isExist:
#     os.makedirs(parent_dir)
#     print("Creating non-existing dir: {}".format(parent_dir))
# path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
# with open(path_file, 'w') as out:
#     s = '''0,{src}-1-2-3-{dst}
# 11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
#     out.write(s)




# GENERATION CONSTANTS. PARIS MADRID
NAME = "starlink_550" # 550 altitude

################################################################
# The below constants are taken from Telesat's FCC filing as below:
# [1]: https://fcc.report/IBFS/SAT-MPL-20200526-00053/2378318.pdf
################################################################

MEAN_MOTION_REV_PER_DAY = 15.19  # Altitude ~550 km
ALTITUDE_M = 550000  # Altitude ~550 km
SATELLITE_CONE_RADIUS_M = ALTITUDE_M / math.tan(math.radians(25.0))  # Considering an elevation angle of 25 degrees;
MAX_GSL_LENGTH_M = math.sqrt(math.pow(SATELLITE_CONE_RADIUS_M, 2) + math.pow(ALTITUDE_M, 2))
# ISLs are not allowed to dip below 80 km altitude in order to avoid weather conditions
MAX_ISL_LENGTH_M = 2 * math.sqrt(math.pow(EARTH_RADIUS + ALTITUDE_M, 2) - math.pow(EARTH_RADIUS + 80000, 2))
NUM_ORBS = 72
NUM_SATS_PER_ORB = 22
INCLINATION_DEGREE = 53


# Input file; Generated during simulation
# Note the file_name consists of the 2 city IDs being offset by the size of the constellation
# City IDs are available in the city_detail_file.
# If city ID is X (for Paris X = 24) and constellation is Starlink_550 (1584 satellites),
# then offset ID is 1584 + 24 = 1608.
num_sats = NUM_ORBS*NUM_SATS_PER_ORB
gs_paris_id = num_sats + 24
gs_madrid_id = num_sats + 54
gs_sao_paolo = num_sats + 3
gs_fortaleza = num_sats + 98
gs_src_id = gs_sao_paolo
gs_dst_id = gs_fortaleza
parent_dir = "../../../../hypatia_plus/data/{}/".format(NAME)
isExist = os.path.exists(parent_dir)
if not isExist:
    os.makedirs(parent_dir)
    print("Creating non-existing dir: {}".format(parent_dir))
path_file = parent_dir + "networkx_path_{}_to_{}.txt".format(gs_src_id, gs_dst_id)
with open(path_file, 'w') as out:
    s = '''0,{src}-1-2-3-{dst}
11700000000000,1608-202-201-203-180-224-1605'''.format(src=gs_src_id, dst=gs_dst_id)
    out.write(s)


# -------------------------------------------------------------------------------------------------------------------

# General files needed to generate visualizations; Do not change for different simulations
topFile = "../static_html/top.html"
bottomFile = "../static_html/bottom.html"
city_detail_file = "../../paper/satellite_networks_state/input_data/ground_stations_cities_sorted_by_estimated_2025_pop_top_1000.basic.txt"

# Time in ms for which visualization will be generated
# GEN_TIME=167500  #ms
# GEN_TIME=11700000 #ms
GEN_TIME = 0 #ms


# Output directory for creating visualization html files
OUT_DIR = "../viz_output/"
OUT_HTML_FILE = OUT_DIR + NAME + "_path"

sat_objs = []
city_details = {}
paths_over_time = []


def generate_path_at_time():
    """
    Generates end-to-end path at specified time
    :return: HTML formatted string for visualization
    """
    viz_string = ""
    global src_GS
    global dst_GS
    global paths_over_time
    global OUT_HTML_FILE
    for line in open(path_file):
        print(line)
    lines = [line.rstrip('\n') for line in open(path_file)]
    for i in range(len(lines)):
        val = lines[i].split(",")
        nodes = val[1].split("-")
        paths_over_time.append((int(val[0]), nodes))
    paths_over_time.append((0, nodes))
    SEL_PATH_TIME = 0
    SEL_PATH = []
    # print("dudde",paths_over_time)
    for i in range(len(paths_over_time)):
        start_ms = round((paths_over_time[i][0]) / 1000000)
        print("satart mas", start_ms)
        start_next = 99999999999
        try:
            start_next = round((paths_over_time[i + 1][0]) / 1000000)
            # print("satart nex", start_next)
        except:
            None
        print(start_ms, GEN_TIME)
        if GEN_TIME >= start_ms and GEN_TIME < start_next:
            SEL_PATH_TIME = paths_over_time[i][0]
            SEL_PATH = paths_over_time[i][1]
            break
    print(SEL_PATH_TIME, SEL_PATH)

    shifted_epoch = (pd.to_datetime(EPOCH) + pd.to_timedelta(GEN_TIME, unit='ms')).strftime(format='%Y/%m/%d %H:%M:%S.%f')
    print(shifted_epoch)

    for i in range(len(sat_objs)):
        sat_objs[i]["sat_obj"].compute(shifted_epoch)
        viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", "+str(sat_objs[i]["alt_km"]*1000)+"), "\
                     + "label: new Cesium.LabelGraphics({position : Cesium.Cartesian3.fromDegrees(" \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublong)) + ", " \
                     + str(math.degrees(sat_objs[i]["sat_obj"].sublat)) + ", "+str(sat_objs[i]["alt_km"]*1000)+"), "\
                     + "text : '" + str(i) + "', font : '18px Helvetica', fillColor : Cesium.Color.BLUE, outlineColor : Cesium.Color.BLACK, outlineWidth : 4,}), "\
                     + "});\n"

    orbit_links = util.find_orbit_links(sat_objs, NUM_ORBS, NUM_SATS_PER_ORB)
    for key in orbit_links:
        sat1 = orbit_links[key]["sat1"]
        sat2 = orbit_links[key]["sat2"]
        viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                      + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                      + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                      + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                      + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                      + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                      + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                      + "width: 0.5, arcType: Cesium.ArcType.NONE, " \
                      + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                      + "color: Cesium.Color.GREY.withAlpha(0.3), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    for p in range(len(SEL_PATH)):
        if p == 0:
            GS = int(SEL_PATH[p]) - NUM_ORBS*NUM_SATS_PER_ORB
            print(city_details[GS]["name"])
            OUT_HTML_FILE += "_"+city_details[GS]["name"] + "_" +str(SEL_PATH[p])
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(city_details[GS]["long_deg"]) + ", " \
                          + str(city_details[GS]["lat_deg"]) + ", " \
                          + str(city_details[GS]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.GREEN.withAlpha(1),}});\n"
            dst = int(SEL_PATH[p + 1])
            viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                          + str(city_details[GS]["long_deg"]) + "," \
                          + str(city_details[GS]["lat_deg"]) + "," \
                          + str(city_details[GS]["alt_km"] * 1000) + "," \
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublat)) + "," \
                          + str(sat_objs[dst]["alt_km"] * 1000) + "]), " \
                          + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                          + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                          + "color: Cesium.Color.RED.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
        if p == len(SEL_PATH) - 1:
            GS = int(SEL_PATH[p]) - NUM_ORBS * NUM_SATS_PER_ORB
            print(city_details[GS]["name"])
            OUT_HTML_FILE += "_" + city_details[GS]["name"] + "_" + str(SEL_PATH[p])
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(city_details[GS]["long_deg"]) + ", " \
                          + str(city_details[GS]["lat_deg"]) + ", " \
                          + str(city_details[GS]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(50000.0, 50000.0, 50000.0), " \
                          + "material : Cesium.Color.GREEN.withAlpha(1),}});\n"
            src = int(SEL_PATH[p-1])
            viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                          + str(city_details[GS]["long_deg"]) + "," \
                          + str(city_details[GS]["lat_deg"]) + "," \
                          + str(city_details[GS]["alt_km"] * 1000) + "," \
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublat)) + "," \
                          + str(sat_objs[src]["alt_km"] * 1000) + "]), " \
                          + "width: 3.0, arcType: Cesium.ArcType.NONE, " \
                          + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                          + "color: Cesium.Color.RED.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
        if 0 < p < len(SEL_PATH) - 2:
            #print(SEL_PATH[p], SEL_PATH[p+1])
            src = int(SEL_PATH[p])
            dst = int(SEL_PATH[p+1])
            viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights(["\
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublong)) + ","\
                          + str(math.degrees(sat_objs[src]["sat_obj"].sublat)) + ","+str(sat_objs[src]["alt_km"]*1000)+","\
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublong)) + ","\
                          + str(math.degrees(sat_objs[dst]["sat_obj"].sublat)) + ","+str(sat_objs[dst]["alt_km"]*1000)+"]), "\
                          + "width: 3.0, arcType: Cesium.ArcType.NONE, "\
                          + "material: new Cesium.PolylineOutlineMaterialProperty({ "\
                          + "color: Cesium.Color.RED.withAlpha(1.0), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

    OUT_HTML_FILE += "_" + str(GEN_TIME) + ".html"
    print(OUT_HTML_FILE)
    return viz_string


city_details = util.read_city_details(city_details, city_detail_file)
sat_objs = util.generate_sat_obj_list(
    NUM_ORBS,
    NUM_SATS_PER_ORB,
    EPOCH,
    PHASE_DIFF,
    INCLINATION_DEGREE,
    ECCENTRICITY,
    ARG_OF_PERIGEE_DEGREE,
    MEAN_MOTION_REV_PER_DAY,
    ALTITUDE_M
)
viz_string = generate_path_at_time()
util.write_viz_files(viz_string, topFile, bottomFile, OUT_HTML_FILE)
