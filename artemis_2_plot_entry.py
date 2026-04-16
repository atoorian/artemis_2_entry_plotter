import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import CartesianRepresentation, CartesianDifferential, GCRS, ITRS, EarthLocation
from astropy.time import Time, TimeDelta
import astropy.units as u
import simplekml

# --- Configuration ---
USE_INERTIAL_VELOCITY = True 
LABEL_INTERVAL = 100 
TARGET_TEMP_MIN = 70.0
TARGET_TEMP_MAX = 5000.0
SPECIFIC_ALTITUDES = [22000, 21000, 20500] # Updated list

FEET_TO_METERS = 0.3048
MPS_TO_MPH = 2.23694

def get_std_atmo_density(alt_m):
    """US Standard Atmosphere 1976 density model."""
    Re = 6356766.0
    h = (Re * alt_m) / (Re + alt_m)
    layers = [
        [0, -0.0065, 288.15, 101325.0], [11000, 0.0, 216.65, 22632.1],
        [20000, 0.001, 216.65, 5474.89], [32000, 0.0028, 228.65, 868.019],
        [47000, 0.0, 270.65, 110.906], [51000, -0.0028, 270.65, 66.9389],
        [71000, -0.002, 214.65, 3.95642], [84852, 0.0, 186.95, 0.3734]
    ]
    R, g0 = 287.05, 9.80665
    if h > 84852: return 0.000001 * np.exp(-(h-84852) / 7500.0)
    curr = layers[0]
    for l in layers:
        if h >= l[0]: curr = l
        else: break
    h0, L, T0, P0 = curr
    T = T0 + L * (h - h0)
    if L == 0:
        P = P0 * np.exp(-g0 * (h - h0) / (R * T0))
    else:
        P = P0 * (T / T0)**(-g0 / (R * L))
    return P / (R * T)

def process_mission_data(input_file, output_kml, output_plot):
    kml = simplekml.Kml()
    ls = kml.newlinestring(name="Artemis II Path")
    folder = kml.newfolder(name="Velocity Labels")
    target_folder = kml.newfolder(name="High-Res Waypoints")
    
    base_epoch = Time("2026-01-01T00:00:00", scale='utc')
    data = np.loadtxt(input_file, skiprows=2)
    
    processed_points, intensities = [], []

    for row in data:
        t_sec, pos_ft, vel_fps, alt_ft = row[0], row[1:4], row[4:7], row[7]
        current_time = base_epoch + TimeDelta(t_sec * u.s)
        
        p_m, v_m = pos_ft * FEET_TO_METERS * u.m, vel_fps * FEET_TO_METERS * u.m / u.s
        ici = GCRS(CartesianRepresentation(p_m, differentials=CartesianDifferential(v_m)), obstime=current_time)
        itrs = ici.transform_to(ITRS(obstime=current_time))
        
        v_frame = ici.velocity if USE_INERTIAL_VELOCITY else itrs.velocity
        v_mps = np.sqrt(v_frame.d_x.value**2 + v_frame.d_y.value**2 + v_frame.d_z.value**2)
        
        loc = EarthLocation.from_geocentric(itrs.x, itrs.y, itrs.z)
        rho = get_std_atmo_density(loc.height.to(u.m).value)
        q_raw = np.sqrt(rho) * (v_mps**3)
        intensities.append(q_raw)
        
        processed_points.append({
            'lla': (loc.lon.value, loc.lat.value, loc.height.to(u.m).value), 
            'v_mph': v_mps * MPS_TO_MPH, 'alt_ft': alt_ft, 'q_raw': q_raw
        })

    # Scaling and Annotation Preparation
    q_min, q_max = min(intensities), max(intensities)
    altitudes_arr = np.array([p['alt_ft'] for p in processed_points])
    target_indices = [np.abs(altitudes_arr - target).argmin() for target in SPECIFIC_ALTITUDES]

    plot_data, path_coords = [], []
    for i, p in enumerate(processed_points):
        temp_f = TARGET_TEMP_MIN + (TARGET_TEMP_MAX - TARGET_TEMP_MIN) * ((p['q_raw'] - q_min) / (q_max - q_min))
        p['temp_f'] = temp_f
        path_coords.append(p['lla'])
        plot_data.append(p)

        # KML Point Creation
        if i in target_indices:
            name = f"ALT: {int(p['alt_ft'])}ft | {int(p['v_mph'])}mph | {int(temp_f)}F"
            tpnt = target_folder.newpoint(name=name)
            tpnt.coords = [p['lla']]
            tpnt.altitudemode = simplekml.AltitudeMode.absolute
            tpnt.style.labelstyle.color = simplekml.Color.yellow

    # --- Plotting with 3-Way Label Spread ---
    df = pd.DataFrame(plot_data)
    fig, ax1 = plt.subplots(figsize=(15, 9))
    ax1.plot(df['alt_ft'], df['v_mph'], color='blue', label='Velocity')
    ax1.set_xlim(df['alt_ft'].max(), df['alt_ft'].min())
    
    ax2 = ax1.twinx()
    ax2.plot(df['alt_ft'], df['temp_f'], color='red', alpha=0.3)
    
    # Offsets: 1: Up-Left, 2: Up-Right, 3: Down
    manual_offsets = [(-110, 80), (110, 80), (0, -110)]
    
    for i, idx in enumerate(target_indices):
        p = processed_points[idx]
        ax1.annotate(
            f"ALT: {int(p['alt_ft'])} ft\n{int(p['v_mph'])} mph\n{int(p['temp_f'])}°F", 
            xy=(p['alt_ft'], p['v_mph']), 
            xytext=manual_offsets[i],
            textcoords='offset points',
            ha='center',
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="black", alpha=0.8),
            arrowprops=dict(arrowstyle='->', connectionstyle="arc3,rad=0.1", color='black')
        )

    plt.title('Artemis II Entry: Extended Data Waypoints')
    plt.savefig(output_plot, dpi=300)
    ls.coords = path_coords
    ls.altitudemode = simplekml.AltitudeMode.absolute
    kml.save(output_kml)

process_mission_data('2026.04.10 - Post-RTC3 to EI.txt', 'Artemis_II_Entry.kml', 'Entry_Profile.png')