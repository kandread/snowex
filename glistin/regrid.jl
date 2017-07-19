# Regrid from SCH to Lat/Long coordinates

"""Calculate latitude and longitude for given point."""
function latlon(slat, slon, heading, dist, ang, Rearth)
    heading += ang
    lat = asin(sin(slat)*cos(dist/Rearth)+cos(slat)*sin(dist/Rearth)*cos(heading))
    lon = slon+atan2(sin(dist/Rearth)*sin(heading), cos(slat)*cos(dist/Rearth)-sin(slat)*sin(dist/Rearth)*cos(heading))
    return (lat, lon)
end

"""
Translates a location in sc local frame to latitude/longitude given a peg point.

Inputs:
plat, plon in radians
pradius in meters
s is along-track shift in meters
c is cross-track shift in meters

Outputs:
lat is pixel latitude in radians
lon is pixel longitude in radians
"""
function schtolatlon(plat, plon, pheading, pradius, s, c)
    flat, flon = latlon(plat, plon, pheading, s, 0, pradius)
    lat, lon = latlon(flat, flon, pheading, c, -pi/2, pradius)
    return (lat, lon)
end
