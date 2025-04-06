import sys


class StkEphemerisWriter:
    """ Writes out the ephemeris in STK planetary ephemeris format.
    """

    def __init__(self, ephemeris):
        self._ephemeris = ephemeris

    def writeall(self):
        sys.stdout.write("stk.v.12.0.0\n")
        sys.stdout.write("BEGIN Ephemeris\n")
        sys.stdout.write(
            "NumberOfEphemerisPoints {}\n".format(len(self._ephemeris)))
        sys.stdout.write("Units km/sec\n")
        sys.stdout.write("EphemerisJ2000SciJedPosVel\n")
        for line in self._ephemeris:
            sys.stdout.write("{}\n".format(" ".join(line)))
        sys.stdout.write("END Ephemeris\n")


class HorizonOutputReader:
    """ Reads the ephemeris from the output produced by a request to the JPL Horizons Api 
    """

    def __init__(self, fs):
        self._fs = fs
        self._ephemeris = []
        self._rf = "J2000"
        self._units = "KM-S"
        self._central_body = ""
        self._target_body = ""

    def parse_posvel(self, line):
        parts = line.split(',')
        jd = parts[0]
        _ = parts[1]
        pos = [parts[2].strip(), parts[3].strip(), parts[4].strip()]
        vel = [parts[5].strip(), parts[6].strip(), parts[7].strip()]
        self._ephemeris.append([jd.strip()] + pos + vel)

    def readall(self):
        """ Read the entire contents of the input stream and parse the ephemeris.
        """
        for line in self._fs:
            if line.startswith("Target body name"):
                self._target_body = line.split(":")[1]
            elif line.startswith("Center body name"):
                self._central_body = line.split(":")[1]
            elif line.startswith("Output units"):
                self._units = line.split(":")[1]
            elif line.startswith("$$SOE"):
                for line in self._fs:
                    if line.startswith("$$EOE"):
                        break
                    self.parse_posvel(line)

    def dump(self):
        """ Dumps the parsed contents of the input stream to the standard output.
        """
        sys.stdout.write("Reference Frame: {}\n".format(self._rf))
        sys.stdout.write("Units: {}\n".format(self._units))
        sys.stdout.write("Central body: {}\n".format(self._central_body))
        sys.stdout.write("Target body: {}\n".format(self._target_body))
        for line in self._ephemeris:
            sys.stdout.write("{}\n".format(",".join(line)))


# Read the horizon API ephemeris report from piped standard input
rd = HorizonOutputReader(sys.stdin)
rd.readall()
# Convert the Horizons API ephemeris report to STK Planetary Ephemeris
wr = StkEphemerisWriter(rd._ephemeris)
# Write the converted STK ephemeris to the standard output
wr.writeall()

