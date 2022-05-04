# Only generate orbital planes that contain simulated SAT IDs and GS IDs (as provided in the "node_to_sat_or_gs_mapping_file" field of our CLUSTER_CONFIG), with a MARGIN of extra orbital planes that are also generated.
# Also generates the ISL connectivity for SAT IDs in the above file, and for your choice of other SAT IDs. We assume and only support +grid for now. 

# from astropy import units as u
# from poliastro.bodies import Earth
# from poliastro.twobody import Orbit
# from astropy.time import Time
# from extractor import CZMLExtractor
import yaml
import math
try:
    from . import util
except (ImportError, SystemError):
    import util

def get_config_yaml_dict(filename):
    with open(filename, "r") as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise ValueError(exc)

CLUSTER_BASE_PATH = "../../../../k8s-vagrant-libvirt/temp/vagrant_libvirt_test_4"
CLUSTER_CONFIG_FILE = CLUSTER_BASE_PATH + "/cluster_config.yaml"
CLUSTER_CONFIG = get_config_yaml_dict(CLUSTER_CONFIG_FILE)

MARGIN = 5 # Represents extra orbits to simulate, that is a distance MARGIN from the orbits containing SAT IDs to be simulated. 

# Generate static visualizations for entire constellation (multiple shells).

EARTH_RADIUS = 6378135.0 # WGS72 value; taken from https://geographiclib.sourceforge.io/html/NET/NETGeographicLib_8h_source.html

# CONSTELLATION GENERATION GENERAL CONSTANTS
ECCENTRICITY = 0.0000001  # Circular orbits are zero, but pyephem does not permit 0, so lowest possible value
ARG_OF_PERIGEE_DEGREE = 0.0
PHASE_DIFF = True
EPOCH = "2000-01-01 00:00:00"

# Shell wise color codes
# COLOR = [[255, 0, 0, 200], [32, 128, 46, 200], [0, 0, 255, 200], [245, 66, 242, 200], [245, 126, 66, 200]]
COLOR = ['FORESTGREEN', 'CRIMSON', 'DODGERBLUE', 'PERU', 'BLUEVIOLET', 'DARKMAGENTA']
CUSTOM_GRID_SAT_IDS = {}
# CONSTELLATION SPECIFIC PARAMETERS


# STARLINK
NAME = "Starlink"

# Number of shells:
SHELL_CNTR = 1

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 15.19  # Altitude ~550000 km. The number of revolutions per day: field 8 of line 2 of TLE
ALTITUDE_M[0] = 550000  # Altitude ~550000 km
NUM_ORBS[0] = 72
NUM_SATS_PER_ORB[0] = 22
INCLINATION_DEGREE[0] = 53
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

# # For paris madrid (commented since it doesnt look nice since the orbits are getting squishy at that position)
# CUSTOM_GRID_SAT_IDS = { # key represents the center node. values represent the nodes connected to that center node. note that we do not want repeat connections. 
#     180: [157, 179, 181, 201],
#     201: [200, 202, 224],
#     224: [223, 225, 245]
# }

# For sao paolo and fortaleza
CUSTOM_GRID_SAT_IDS = { # key represents the center node. values represent the nodes connected to that center node. note that we do not want repeat connections. 
    285: [262, 284, 286, 306],
    306: [305, 307, 329],
    329: [328, 330, 350]
}

# TOTAL_SATS = NUM_SATS_PER_ORB[0] * NUM

# MEAN_MOTION_REV_PER_DAY[1] = 13.4  # Altitude ~1110 km
# ALTITUDE_M[1] = 1110000  # Altitude ~1110 km
# NUM_ORBS[1] = 32
# NUM_SATS_PER_ORB[1] = 50
# INCLINATION_DEGREE[1] = 53.8
# BASE_ID[1] = 1584
# ORB_WISE_IDS[1] = []

# MEAN_MOTION_REV_PER_DAY[2] = 13.35  # Altitude ~1130 km
# ALTITUDE_M[2] = 1130000  # Altitude ~1130 km
# NUM_ORBS[2] = 8
# NUM_SATS_PER_ORB[2] = 50
# INCLINATION_DEGREE[2] = 74
# BASE_ID[2] = 3184
# ORB_WISE_IDS[2] = []

# MEAN_MOTION_REV_PER_DAY[3] = 12.97  # Altitude ~1275 km
# ALTITUDE_M[3] = 1275000  # Altitude ~1275 km
# NUM_ORBS[3] = 5
# NUM_SATS_PER_ORB[3] = 75
# INCLINATION_DEGREE[3] = 81
# BASE_ID[3] = 3584
# ORB_WISE_IDS[3] = []

# MEAN_MOTION_REV_PER_DAY[4] = 12.84  # Altitude ~1325 km
# ALTITUDE_M[4] = 1325000  # Altitude ~1325 km
# NUM_ORBS[4] = 6
# NUM_SATS_PER_ORB[4] = 75
# INCLINATION_DEGREE[4] = 70
# BASE_ID[4] = 3959
# ORB_WISE_IDS[4] = []


"""
# TELESAT
NAME = "Telesat"
SHELL_CNTR = 2

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 13.66  # Altitude ~1015 km
ALTITUDE_M[0] = 1015000  # Altitude ~1015 km
NUM_ORBS[0] = 27
NUM_SATS_PER_ORB[0] = 13
INCLINATION_DEGREE[0] = 98.98
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

MEAN_MOTION_REV_PER_DAY[1] = 12.84  # Altitude ~1325 km
ALTITUDE_M[1] = 1325000  # Altitude ~1325 km
NUM_ORBS[1] = 40
NUM_SATS_PER_ORB[1] = 33
INCLINATION_DEGREE[1] = 50.88
BASE_ID[1] = 351
ORB_WISE_IDS[1] = []
"""

"""
# KUIPER
NAME = "kuiper"
################################################################
# The below constants are taken from Kuiper's FCC filing as below:
# [1]: https://www.itu.int/ITU-R/space/asreceived/Publication/DisplayPublication/8716
################################################################

SHELL_CNTR = 3

MEAN_MOTION_REV_PER_DAY = [None]*SHELL_CNTR
ALTITUDE_M = [None]*SHELL_CNTR
NUM_ORBS = [None]*SHELL_CNTR
NUM_SATS_PER_ORB = [None]*SHELL_CNTR
INCLINATION_DEGREE = [None]*SHELL_CNTR
BASE_ID = [None]*SHELL_CNTR
ORB_WISE_IDS = [None]*SHELL_CNTR

MEAN_MOTION_REV_PER_DAY[0] = 14.80  # Altitude ~630 km
ALTITUDE_M[0] = 630000  # Altitude ~630 km
NUM_ORBS[0] = 34
NUM_SATS_PER_ORB[0] = 34
INCLINATION_DEGREE[0] = 51.9
BASE_ID[0] = 0
ORB_WISE_IDS[0] = []

MEAN_MOTION_REV_PER_DAY[1] = 14.86  # Altitude ~610 km
ALTITUDE_M[1] = 610000  # Altitude ~610 km
NUM_ORBS[1] = 36
NUM_SATS_PER_ORB[1] = 36
INCLINATION_DEGREE[1] = 42
BASE_ID[1] = 1156
ORB_WISE_IDS[1] = []

MEAN_MOTION_REV_PER_DAY[2] = 14.93  # Altitude ~590 km
ALTITUDE_M[2] = 590000  # Altitude ~590 km
NUM_ORBS[2] = 28
NUM_SATS_PER_ORB[2] = 28
INCLINATION_DEGREE[2] = 33
BASE_ID[2] = 2452
ORB_WISE_IDS[2] = []
"""


# General files needed to generate visualizations; Do not change for different simulations
topFile = "../static_html/top.html"
bottomFile = "../static_html/bottom.html"
city_detail_file = "../../paper/satellite_networks_state/input_data/ground_stations_cities_sorted_by_estimated_2025_pop_top_1000.basic.txt"
city_details = {}

# Output directory for creating visualization html files
OUT_DIR = "../viz_output/"
# JSON_NAME  = NAME+"_5shell.json"
# OUT_JSON_FILE = OUT_DIR + JSON_NAME
OUT_HTML_FILE = OUT_DIR + "subset_" + NAME + ".html"

# START = Time(EPOCH, scale="tdb")
# END = START + (10*60) * u.second
# sample_points = 10
# extractor = CZMLExtractor(START, END, sample_points)

# viewer.entities.add({ polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([-100.60303790784423,-0.12458382439979618,550000,-100.60303790784423,8.896941233608787,550000]), width: 2.0, arcType: Cesium.ArcType.GEODESIC, material: new Cesium.PolylineOutlineMaterialProperty({ color: Cesium.Color.BLUE, outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});
# var redSphere = viewer.entities.add({name : 'SAT', position: Cesium.Cartesian3.fromDegrees(-100.60303790784423, -0.12458382439979618, 550000), ellipsoid : {radii : new Cesium.Cartesian3(70000.0, 70000.0, 70000.0), material : Cesium.Color.BLACK.withAlpha(1),}});
# viewer.entities.add({ polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([-95.2040329161376,7.1123809968927665,550000,-89.57043141700855,14.239714131652336,550000]), width: 8.0, arcType: Cesium.ArcType.GEODESIC, material: new Cesium.PolylineOutlineMaterialProperty({ color: Cesium.Color.TOMATO, outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});


def get_mapped_sats_and_gs():
    """Get SAT IDs and corrected GS IDs that are included in the simulation"""
    gs_ids = []
    sat_ids = []
    NODE_TO_SAT_OR_GS_MAPPING_FILE = CLUSTER_BASE_PATH + "/" + CLUSTER_CONFIG["node_to_sat_or_gs_mapping_file"]
    with open(NODE_TO_SAT_OR_GS_MAPPING_FILE, 'r') as mapping:
        for line in mapping:
            row = line.split()
            if "g" in row[0]:
                gs_ids.append(int(row[1]))
            else:
                sat_ids.append(int(row[1]))
    return sat_ids, gs_ids

def get_orbits_to_display(sat_ids, NUM_ORBS, NUM_SATS_IN_ORB):
    orbit_ids = []
    for sat in sat_ids:
        orbit_ids.append(int(sat/NUM_SATS_IN_ORB))
        if NUM_ORBS < int(sat/NUM_SATS_IN_ORB):
            raise ValueError("Orbit ID does not exist. This orbit exceeds the number of orbits available for the given constellation.")
    return orbit_ids

def generate_satellite_trajectories():
    """
    Generates and adds satellite orbits to visualization.
    :return: viz_string
    """
    viz_string = ""
    for i in range(0, SHELL_CNTR):
        sat_objs = util.generate_sat_obj_list(
            NUM_ORBS[i],
            NUM_SATS_PER_ORB[i],
            EPOCH,
            PHASE_DIFF,
            INCLINATION_DEGREE[i],
            ECCENTRICITY,
            ARG_OF_PERIGEE_DEGREE,
            MEAN_MOTION_REV_PER_DAY[i],
            ALTITUDE_M[i]
        )
        sat_ids, gs_ids = get_mapped_sats_and_gs()
        orbit_ids_to_display = get_orbits_to_display(sat_ids, NUM_ORBS[i], NUM_SATS_PER_ORB[i])
        # SECTION: generate SAT IDs as points
        for j in range(len(sat_objs)):
            # Only display SAT IDs that are simulated:
            cur_orbit = int(j/NUM_SATS_PER_ORB[i])
            orbit_found = False
            for orbit in orbit_ids_to_display:
                if orbit + MARGIN >= cur_orbit and orbit - MARGIN <= cur_orbit:
                    orbit_found = True
            if not orbit_found:
                continue
            # if cur_orbit not in orbit_ids_to_display:
            #     continue
            # print(orbit_ids_to_display)
            sat_objs[j]["sat_obj"].compute(EPOCH)
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                          + str(math.degrees(sat_objs[j]["sat_obj"].sublong)) + ", " \
                          + str(math.degrees(sat_objs[j]["sat_obj"].sublat)) + ", " + str(
                sat_objs[j]["alt_km"] * 1000) + "), " \
                          + "ellipsoid : {radii : new Cesium.Cartesian3(30000.0, 30000.0, 30000.0), " \
                          + "material : Cesium.Color.BLACK.withAlpha(1),}});\n"
        # SECTION: generate lines linking SAT IDs in the same orbital plane
        orbit_links = util.find_orbit_links(sat_objs, NUM_ORBS[i], NUM_SATS_PER_ORB[i])
        for key in orbit_links:
            sat1 = orbit_links[key]["sat1"]
            sat2 = orbit_links[key]["sat2"]

            # Only display connectivity for SAT IDs that are simulated:
            sat1_orbit = int(sat1/NUM_SATS_PER_ORB[i])
            sat2_orbit = int(sat2/NUM_SATS_PER_ORB[i])
            orbit1_found = False
            orbit2_found = False
            for orbit in orbit_ids_to_display:
                if orbit + MARGIN >= sat1_orbit and orbit - MARGIN <= sat1_orbit:
                    orbit1_found = True
                if orbit + MARGIN >= sat2_orbit and orbit - MARGIN <= sat2_orbit:
                    orbit2_found = True
            if not orbit1_found or not orbit2_found:
                continue
            # if sat1_orbit not in orbit_ids_to_display or sat2_orbit not in orbit_ids_to_display:
            #     continue

            viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                          + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                          + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                          + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                          + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                          + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                          + "width: 0.5, arcType: Cesium.ArcType.GEODESIC, " \
                          + "material: new Cesium.PolylineOutlineMaterialProperty({ " \
                          + "color: Cesium.Color."+COLOR[i]+".withAlpha(0.4), outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
        # SECTION: generate desired ISLs:
        # use a custom set of SAT IDs instead of that used in the simulation (which is used to generate the orbital planes)
        if CUSTOM_GRID_SAT_IDS:
            for sat1, values in CUSTOM_GRID_SAT_IDS.items():
                for sat2 in values:
                    viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                        + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                        + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                        + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                        + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                        + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                        + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                        + "width: 8, arcType: Cesium.ArcType.GEODESIC, " \
                        + "material: new Cesium.PolylineOutlineMaterialProperty({ color: Cesium.Color.TOMATO, outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"
        else:
            grid_links = util.find_grid_links(sat_objs, NUM_ORBS[i], NUM_SATS_PER_ORB[i])
            reverse_grid_links = util.find_reverse_grid_links(sat_objs, NUM_ORBS[i], NUM_SATS_PER_ORB[i])
            for key in grid_links:
                sat1 = grid_links[key]["sat1"]
                sat2 = grid_links[key]["sat2"]

                if sat1 not in sat_ids: 
                    continue

                viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                            + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                            + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                            + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                            + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                            + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                            + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                            + "width: 8, arcType: Cesium.ArcType.GEODESIC, " \
                            + "material: new Cesium.PolylineOutlineMaterialProperty({ color: Cesium.Color.TOMATO, outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

                sat1 = reverse_grid_links[key]["sat1"]
                sat2 = reverse_grid_links[key]["sat2"]

                viz_string += "viewer.entities.add({name : '', polyline: { positions: Cesium.Cartesian3.fromDegreesArrayHeights([" \
                    + str(math.degrees(sat_objs[sat1]["sat_obj"].sublong)) + "," \
                    + str(math.degrees(sat_objs[sat1]["sat_obj"].sublat)) + "," \
                    + str(sat_objs[sat1]["alt_km"] * 1000) + "," \
                    + str(math.degrees(sat_objs[sat2]["sat_obj"].sublong)) + "," \
                    + str(math.degrees(sat_objs[sat2]["sat_obj"].sublat)) + "," \
                    + str(sat_objs[sat2]["alt_km"] * 1000) + "]), " \
                    + "width: 8, arcType: Cesium.ArcType.GEODESIC, " \
                    + "material: new Cesium.PolylineOutlineMaterialProperty({ color: Cesium.Color.TOMATO, outlineWidth: 0, outlineColor: Cesium.Color.BLACK})}});"

        #SECTION: generate GS IDs as points:
        GS_COLOR = ["GREY", "INDIGO"]
        for index, gs in enumerate(gs_ids):
            GS = gs - NUM_ORBS[i] * NUM_SATS_PER_ORB[i]

            gs_color = GS_COLOR[int((index+1)/len(GS_COLOR))]
            
            # print(city_details[GS]["name"])
            # With city name label:
            # viz_string += "var redSphere = viewer.entities.add({name : ''," \
            #                 + "label: {text: '" + city_details[GS]["name"] + "', font : '12pt monospace', style: Cesium.LabelStyle.FILL_AND_OUTLINE, outlineWidth : 2, verticalOrigin : Cesium.VerticalOrigin.BOTTOM, pixelOffset : new Cesium.Cartesian2(0, 10)}, position: Cesium.Cartesian3.fromDegrees(" \
            #                 + str(city_details[GS]["long_deg"]) + ", " \
            #                 + str(city_details[GS]["lat_deg"]) + ", " \
            #                 + str(city_details[GS]["alt_km"] * 1000) + "), " \
            #                 + "ellipsoid : {radii : new Cesium.Cartesian3(60000.0, 60000.0, 60000.0), " \
            #                 + "material : Cesium.Color.GREY.withAlpha(1),}});\n"
            viz_string += "var redSphere = viewer.entities.add({name : '', position: Cesium.Cartesian3.fromDegrees(" \
                + str(city_details[GS]["long_deg"]) + ", " \
                + str(city_details[GS]["lat_deg"]) + ", " \
                + str(city_details[GS]["alt_km"] * 1000) + "), " \
                + "ellipsoid : {radii : new Cesium.Cartesian3(90000.0, 90000.0, 90000.0), " \
                + "material : Cesium.Color." + gs_color + ".withAlpha(1),}});\n"

    return viz_string


def write_viz_files():
    """
    Writes JSON and TML files to the output folder
    :return: None
    """
    writer_html = open(OUT_HTML_FILE, 'w')
    with open(topFile, 'r') as fi:
        writer_html.write(fi.read())
    writer_html.write(viz_string)
    with open(bottomFile, 'r') as fb:
        writer_html.write(fb.read())
    writer_html.close()

city_details = util.read_city_details(city_details, city_detail_file)
viz_string = generate_satellite_trajectories()
write_viz_files()