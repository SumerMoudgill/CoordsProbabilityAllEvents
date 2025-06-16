
from ligo.gracedb.rest import GraceDb
import urllib.request
import sys
from astropy.time import Time
import math
g = GraceDb(service_url="https://gracedb.ligo.org/api")

triggers_250611=['S250611a', 'S250611b', 'S250611d', 'S250611e','S250611f','S250611g','S250611h', 'S250611u', 'S250611ac','S250611ag']
gpstimes = []
for i in range(len(triggers_250611)):
    #r = client.superevent(triggers_250611[i])
    #gpstime=r['gpstime']
    for event in g.events('superevent: '+triggers_250611[i]):
        gpstime=event['gpstime']
        print("gpstime =", gpstime)
    gpstimes.append(gpstime)

print("gpstimes =", gpstimes)


def gps2jd(timeGPS):
    t=Time(timeGPS, format="gps")
    t.format="jd"
    #print(t.value)
    return t.value
jdtimes = []
for i in gpstimes:
    jdtimes.append(gps2jd(i))
print("jdtimes =", jdtimes)
O4c_250611_list_0=[]
for i in range(len(triggers_250611)):
    O4c_250611_list_0.append([triggers_250611[i], 0, jdtimes[i], "H1L1"])
print("O4c_250611_list_0 =", O4c_250611_list_0)
