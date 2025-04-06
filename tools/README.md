# Horizon API Ephemeris/STK Planetary Ephemeris Converter

Queries planetary or asteroid ephemeris using the Horizon API and converts the results to an STK planetary format (.pe)

## Heliocentric Input
|         Parameter          |              Value                   |
|----------------------------|--------------------------------------|
| Dates                      |  2029-01-01  to 2030-01-01           |
| Table format               |  x,y,z,vx,vy,vz                      |
| Step Size                  | 10 minutes                           |
| CSV format                 | YES                                  |
| Coordinate Center          |  Solar Barycentric (500@0)           |
| Reference Frame            | ICRF                                 |


## Geocentric Input
|         Parameter          |              Value                   |
|----------------------------|--------------------------------------|
| Dates                      |  2029-01-01  to 2030-01-01           |
| Table format               |  x,y,z,vx,vy,vz                      |
| Step Size                  | 10 minutes                           |
| CSV format                 | YES                                  |
| Coordinate Center          | Geocentric (500)                     |
| Reference Frame            | ICRF                                 |


## Request the asteroid and the Earth position/velocity (Heliocentric)
```sh
#Request the asteroid 2024 YR4 position/velocity information between 2029-01-01 and 2030-01-01 and convert it to STK Planetary Ephemeris 
curl -s "https://ssd.jpl.nasa.gov/api/horizons.api?\
format=text\
&COMMAND='2024%20YR4'\
&OBJ_DATA=NO\
&MAKE_EPHEM=YES\
&EPHEM_TYPE=VECTORS\
&VEC_TABLE=2\
&START_TIME=2029-01-01\
&STOP_TIME=2030-01-01\
&STEP_SIZE='10m'\
&CSV_FORMAT=YES\
&CENTER='500@0'\
&REF_SYSTEM='ICRF'" | python horizons2pe.py

# Request the asteroid Apophis position/velocity information between 2029-01-01 and 2030-01-01 and convert it to STK Planetary Ephemeris 
curl -s "https://ssd.jpl.nasa.gov/api/horizons.api?\
format=text\
&COMMAND='Apophis'\
&OBJ_DATA=NO\
&MAKE_EPHEM=YES\
&EPHEM_TYPE=VECTORS\
&VEC_TABLE=2\
&START_TIME=2029-01-01\
&STOP_TIME=2030-01-01\
&STEP_SIZE='10m'\
&CSV_FORMAT=YES\
&CENTER='500@0'\
&REF_SYSTEM='ICRF'" | python horizons2pe.py

# Request the Earth's position/velocity information  between 2029-01-01 and 2030-01-01 and convert it to STK Planetary Ephemeris 
curl -s "https://ssd.jpl.nasa.gov/api/horizons.api?\
format=text\
&COMMAND='399'\
&OBJ_DATA=NO\
&MAKE_EPHEM=YES\
&EPHEM_TYPE=VECTORS\
&VEC_TABLE=2\
&START_TIME=2029-01-01\
&STOP_TIME=2030-01-01\
&STEP_SIZE='10m'\
&CSV_FORMAT=YES\
&CENTER='500@0'\
&REF_SYSTEM='ICRF'" | python horizons2pe.py

```
## Example Output
```sh
stk.v.12.0.0
BEGIN Ephemeris
NumberOfEphemerisPoints 52561
Units km/sec
EphemerisJ2000SciJedPosVel
2462137.500000000 -7.968847732904576E+07 1.374228615435446E+08 -9.232632884020425E+06 -2.429957177342700E+01 -1.061985482613310E+01 1.359831878726281E-02
2462137.506944444 -7.970305660096838E+07 1.374164888188204E+08 -9.232624670537107E+06 -2.429800126380108E+01 -1.062256088516606E+01 1.377995785100028E-02
2462137.513888889 -7.971763493051520E+07 1.374101144705260E+08 -9.232616348072268E+06 -2.429643052069193E+01 -1.062526672720390E+01 1.396159060787561E-02
2462137.520833333 -7.973221231754611E+07 1.374037384987914E+08 -9.232607916629657E+06 -2.429485954412588E+01 -1.062797235224323E+01 1.414321705709476E-02
2462137.527777778 -7.974678876192199E+07 1.373973609037464E+08 -9.232599376213051E+06 -2.429328833412910E+01 -1.063067776028081E+01 1.432483719788635E-02
2462137.534722222 -7.976136426350090E+07 1.373909816855222E+08 -9.232590726826273E+06 -2.429171689072811E+01 -1.063338295131289E+01 1.450645102943504E-02
...
END Ephemeris
```

## Request YR4 and Apophis asteroids' position/velocity (Geocentric)

```sh
#Request the asteroid 2024 YR4 position/velocity information between 2029-01-01 and 2030-01-01 and convert it to STK Planetary Ephemeris 
curl -s "https://ssd.jpl.nasa.gov/api/horizons.api?\
format=text\
&COMMAND='2024%20YR4'\
&OBJ_DATA=NO\
&MAKE_EPHEM=YES\
&EPHEM_TYPE=VECTORS\
&VEC_TABLE=2\
&START_TIME=2029-01-01\
&STOP_TIME=2030-01-01\
&STEP_SIZE='10m'\
&CSV_FORMAT=YES\
&CENTER='500'\
&REF_SYSTEM='ICRF'" | python horizons2pe.py

# Request the asteroid Apophis position/velocity information between 2029-01-01 and 2030-01-01 and convert it to STK Planetary Ephemeris 
curl -s "https://ssd.jpl.nasa.gov/api/horizons.api?\
format=text\
&COMMAND='Apophis'\
&OBJ_DATA=NO\
&MAKE_EPHEM=YES\
&EPHEM_TYPE=VECTORS\
&VEC_TABLE=2\
&START_TIME=2029-01-01\
&STOP_TIME=2030-01-01\
&STEP_SIZE='10m'\
&CSV_FORMAT=YES\
&CENTER='500'\
&REF_SYSTEM='ICRF'" | python horizons2pe.py

```
