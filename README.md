# artemis_2_entry_plotter


My friend kept sending me conspiracy theory reels flat earther influencers.  Pretty sure he is just trolling me, but just in case I decided to make a tool to generate an Google Earth KML file from the Artemis II Entry ephemeris file on the NASA page.

https://www.nasa.gov/missions/artemis/artemis-2/track-nasas-artemis-ii-mission-in-real-time/

The script takes an ephemeris file in the M50 inertial reference frame, converts it to a geodetic reference frame, and converts the units to ft and mph.  It then generates a KML file of the position information that can be plotted in Google Earth.

It also generates a 2D plot that shows ground track velocity magnitude vs. altitude.

Finally it generates an estimated heat shield peak temperature profile assuming an initial temperature of 70F and a peak temperature of 5000F per the Artemis II website.  This part is very notional and not based on real data.

<img width="1121" height="532" alt="image" src="https://github.com/user-attachments/assets/7d666b3f-23c2-4fc1-8977-5f090dd703d7" />

https://earth.google.com/earth/d/1ZuynG8EvG_PnrAVlfprr-r49wpf0DxnL?usp=sharing
