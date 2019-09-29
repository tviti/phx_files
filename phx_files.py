# Python module with functions for reading and writing various Phoenix data
# files.
#
# These aren't all fully implemented yet, but the goal is for this software
# package to be able to read files with .TBL, .TSn, and .MT extensions. It
# would be nice to also be able to read .CSC and .CLB files, but it appears
# that the format for these is proprietary without open documentation, so for
# now we're stuck with using the Phoenix windows cmdline utility SYSCAL.exe for
# converting these to system transfer functions. Thankfully, SYSCAL.exe will
# run without a license for Phoenix, and runs on *NIX platforms through wine.
#
# The system for storing MT data in .TSn files described in [1] revolves around
# three central data objects: Tables, TAGs and records.  A record is a window
# of data recorded with a specific set of MTU acquisition settings, with the
# acquisition settings (e.g. window start and end times, sampling rate, number
# of samples) stored in a TAGs and TBLs.  A complete set of samples (one for
# each MTU channel) comprises a scan, with multiple scans comprising a single
# record.
#
# Table files (extension .TBL) contain survey parameters for a single MT
# sounding (e.g. dipole lengths, sensor orientations, declination...).  Each
# MTU deployment can result in a number of .TBL files describing a set of
# contiguous, uninterrupted soundings. A single sounding usually records data
# at multiple sampling rates simultaneously, the data from each of which is
# stored in the multiple .TSn files (i.e. each .TBL has associated with it
# multiple .TSn files).
#
# References:
#    [1] Rogers, S. (2005). Data Processing User Guide V3.0. Toronto, ON
#    Canada: Phoenix Geophysics Limited
#    [2] Phoenix Geophysics Limited (N.D.), MTU Time Series Format. Toronto, ON
#    Canada: Phoenix Geophysics Limited
#    [3] Phoenix Geophysics Limited (2009), SysCal User Guide V1.1. Toronto, ON
#    Canada: Phoenix Geophysics Limited
#    [4] Phoenix Geophysics Limited (N.D.). MTC-150 Operating Instructions
#    Ver2. Toronto, ON Canada: Phoenix Geophysics Limited
#    [5] Phoenix Geophysics Limited (2015). Instrument and Sensor Calibration:
#    Concepts and Utility Programs. Toronto, ON Canada: Phoenix Geophysics
#    Limited
#
# Author: Taylor Viti
# Date Created: 05/12/18
#
import struct
import numpy as np
import subprocess
import datetime
from sys import platform

from IPython.core.debugger import set_trace

# Numpy dtype object specifying the fields present in data block
TAG_dt = np.dtype([("second", np.uint8),
                   ("minute", np.uint8),
                   ("hour", np.uint8),
                   ("day", np.uint8),
                   ("month", np.uint8),
                   ("year", np.uint8),  # Last two digits
                   ("day_of_week", np.uint8),
                   ("century", np.uint8),
                   ("SN", np.uint16),  # Instrument serial number
                   ("SC", np.uint16),  # Number of scans in the record
                   ("CH", np.uint8),  # Number of channels per scan
                   ("TG", np.uint8),  # Tag length (in bytes?)
                   ("ST", np.uint8),  # Status code (see [1] p. 192)
                   ("SA", np.uint8),  # bit-wise saturation flags
                   ("reserved1", np.uint8),
                   ("BY", np.uint8),  # Sample length in bytes (always 3 [1])
                   ("SR", np.uint16),  # Sampling rate
                   ("SR_units", np.uint8),  # Sample rate units ([1] p. 192)
                   ("CK", np.uint8),  # Clock status
                   ("TE", np.uint32),  # Clock error in micro secs
                   ("reserved2", np.uint8, (6, ))])  # reserved (MUST BE ZERO)


def read_TSn(fname):
    """
    Read in all records and TAGS from a TSn file created by a Phoenix
    Geophysics MTU receiver. Note that this function assumes all TAGS follow
    the 32 byte TAG standard described in [1]. This function will NOT work for
    reading older timeseries files that follow the older 16 byte TAG standard
    (extensions .TSH and .TSL).

    The records returned will be in integer valued instrument units, and must
    be converted to field units.

    Args:
        fname (str): path to the TSn file

    Returns:
        tuple: tuple containing
            TAGs (list): list of TAG arrays as returned by read_record
            recs (list): list of record arrays as returned by read_record

    """

    # Init empty TAG and records lists
    TAGs = []
    recs = []

    with open(fname, mode="rb") as ts_file:
        data = ts_file.read()

    N_bytes = len(data)
    i = 0
    while i + 32 < N_bytes:
        TAG = np.frombuffer(data[i:(i+32)], dtype=TAG_dt)
        SC = TAG["SC"][0]
        CH = TAG["CH"][0]

        # Increment the iterator to the end of the TAG
        i = i + 32

        # Find the end of the record
        i_end = i + 3*SC*CH
        if i_end > N_bytes:
            i_end = N_bytes

        # Read in ALL the scans at once!
        rec = [data[j:(j+3)] for j in range(i, i_end, 3)]
        # Convert from 24 bit two's complement to normal integers
        rec = [struct.unpack("<i", samp + (b"\x00" if samp[2] < 128 else b"\xff"))[0] for samp in rec]

        # Reshape the rec based on the number of recs ACTUALLY read
        rec = np.array(rec)
        rec.shape = (SC, CH)

        TAGs.append(TAG)
        recs.append(rec)

        # increment the iterator to the next rec
        i = i_end

    return (TAGs, recs)


# def write_TBL(fname, TBL, grps, smphs, types):
#     """ Write a Phoenix .TBL file.

#     Note that unless all of the necessary fields are present in the .TBL file,
#     there's a good chance it will be useless for running an actual MTU. Thus,
#     this function is only intended to be useful for editing .TBL files that
#     have a few incorrect settings.

#     Args:

#     Returns:
#     """

#     with open(fname, mode="wb") as file:
#         for k in tbl.keys():

#             # Write the 5 character field code
#             bytestr = bytes(k, "utf-8")
#             struct.pack("I")


def read_TBL(fname):
    """
    Read a .TBL file created by a Phoenix Geophysics MTU receiver.

    The actual datatypes associated with each integer type code is referenced
    from [2].

    grps and smphs are only useful if you intend to write to a new .TBL

    The Phoenix .TBL files are always terminated with a ascii ETX, which this
    function excludes from the output dictionaries

    Args:
        fname (str): path to the .TBL file to read

    Returns:
        TBL : dict of table entries
        grps: dict of group codes
        smphs: dict of semaphore codes
        types : dict of integer type codes

    """

    TBL = dict()
    grps = dict()
    smphs = dict()
    types = dict()
    with open(fname, mode="rb") as file:
        data = file.read()

    stride = 25  # number of bytes per record of data
    if len(data) % stride != 0:
        raise Exception("Record length does not evenly divide data file!")

    N_rec = int(len(data)/stride)

    for i in range(N_rec):
        i0 = stride*i
        i1 = stride*(i+1)
        rec = data[i0:i1]

        Code = rec[0:5].split(b"\0")[0].decode("utf-8")
        # Last record is always signified with an ETX. Exit if we've hit it
        if Code == '\x03':
            break

        # Grp and Smph are only used by the MTU
        Grp = struct.unpack("<H", rec[5:7])[0]
        Smph = struct.unpack("<i", rec[7:11])[0]

        Type = rec[11]
        if Type in (2, 3, 4):  # String types
            V = rec[12:].split(b"\0")[0].decode("utf-8")
        elif Type == 1:  # Double prec float
            V = struct.unpack("<1d", rec[12:20])[0]
        elif Type == 5:  # AMX type code
            V = struct.unpack("<8B", rec[12:20])
            # Convert to a proper datetime obj
            # V = AMX_to_datetime(V)
        else:  # Everything else is an int
            V = struct.unpack("<1i", rec[12:16])[0]

        TBL[Code] = V
        grps[Code] = Grp
        smphs[Code] = Smph
        types[Code] = Type

    # Convert lat and long to something more sensible
    LATG = TBL["LATG"]
    NS = LATG[-1]

    LNGG = TBL["LNGG"]
    EW = LNGG[-1]

    TBL["LATG"] = ((-1 if NS == "S" else 1)*int(LATG[0:2]), float(LATG[2:-2]))
    TBL["LNGG"] = ((-1 if EW == "W" else 1)*int(LNGG[0:3]), float(LNGG[3:-2]))

    return (TBL, grps, smphs, types)

# def read_TBL_txt(fname):
#     """ Read .txt files created by the program PRNTTBL.EXE.
#     PRNTTBL will convert Phoenix .TBL files into an ascii .txt
#     file.

#     Each line in a TBL.txt is formatted as
#         code, , , type, val

#     Args:
#         fname (str): path to the .txt file to be read

#     Returns:
#         TBL (np.array): array of table entries (dtype=TBL_dt)
#     """

#     tbl_file = open(fname)

#     # Create an empty TBL file container
#     TBL = np.empty((1), dtype=TBL_dt)

#     for line in tbl_file:
#         code = line.split(",")[0]

#         # The NS/EW in lat and lon are seperated from the rest
#         # of the value field by a comma, so we have to take the
#         # entire set of trailing elements of the split to make
#         # sure we can catch it. For all other fields, the 4th
#         # element actually is the -1th element.
#         val = line.split(",")[4:-1]
#         # Recombine the elements comprising value
#         val = " ".join(val)

#         # Trim whitespace from the strings
#         code = code.strip()
#         val = val.strip()

#         if code in TBL.dtype.names:
#             TBL[code] = val

#     tbl_file.close()

#     return TBL


def read_cts(fname):
    """Read the data from a .cts file created by SYSCAL.EXE

    A .cts file contains transfer functions for each recording channel,
    including the sensor and box responses. The data is intended to be used for
    performing a deconvolution on the records from each channel, in order to
    correct out the instrument's response from the data.

    If the function returns 0 for field_type, then the real and imaginary parts
    of the responses have units of (V)^-1.

    If instead the field_type is 1, then the E channels have units of (V/m)^-1,
    and the H channels have units of T^-1.

    The return value G follows the same format as the .cts file, but with the
    channel responses stored as np.complex objects (wheras in the .cts, they're
    stored as seperate columns for the real and imaginary components).

    The first column of G is the frequencies that SYSCAL evaluated the
    responses at (this is controlled via the .PFC file used as input to
    SYSCAL. See [3]).

    The second column of G is the level that frequency corresponds to (i.e. the
    elements of column 3 are the numbers n of the .TSn files that each
    frequency corresponds to).

    The remaining columns of G are the instrument responses for each channel,
    where G[:, 2] is the complex valued response for channel 1, G[:, 3] is the
    complex valued response for channel 2, etc.

    Args:
        fname (str): path to the .cts file

    Returns:
        G (np.ndarray): complex (freqs x 2+channels) array of channel TFs
        field_type (bool): 0 for box cal, 1 for full system cal

    """

    # Read the header from the file
    header = np.genfromtxt(fname, delimiter=",", max_rows=1)
    field_type = bool(header[2])
    N_ch = int(header[3])  # number of channels on the instrument

    # Read the body from the file
    body = np.genfromtxt(fname, delimiter=",", skip_header=1)

    # Create and populate the response matrix
    G = np.empty((body.shape[0], N_ch+2), dtype=np.complex)
    G[:, 0] = body[:, 0]
    G[:, 1] = body[:, 1]
    G[:, 2:] = body[:, 2:-2:2] + 1j*body[:, 3::2]
    # for (i, re, im) in zip(range(2, G.shape[1]), body[:, 2::2]):
    #     G[:, i] = body[:, i] + 1j*body[:, i + 1]

    return (G, field_type)


def get_syscal(f, TBL_fname, CLB_dir, CLC_dir,
               SYSCAL_fname="./bin/SysCal_V7.exe",
               PFC_fname="/tmp/TMP.PFC",
               cts_fname="/tmp/TMP.cts",
               field_type=1,
               levels=None):
    """Use SYSCAL.exe to compute the system response for a specific set of
    acquisition parameters and frequencies.

    REQUIRES A WORKING WINE INSTALLATION!

    SYSCAL.exe CAN ONLY HANDLE 1144 FREQUENCIES AT ONCE! Above that, it will
    still output a response file, but with only 1144 frequencies.

    It's also probably better to supply the levels yourself (using the optional
    arg), rather than letting this function pick them for you, since you will
    likely be processing each level independently.

    NOTE: This script is currently written to be run in either OS X or Linux
    using wine to call SYSCAL.exe.  I don't think it'd be that difficult to
    make this function clever enough to know whether it's being called from a
    UNIX or Windows platform.

    This function is expecting SYSCAL.exe to reside in ./bin/.  It's probably
    not kosher to be redistributing SYSCAL.exe, so if this package is ever to
    be distributed, the user will probably have to place SYSCAL.exe in ./bin/
    (or specify the location explicitely).

    The SYSCAL decides which frequencies to calculate the response at based on
    a parameter file with extension .PFC. This function will fill a brand new
    .PFC with the frequencies contained in the argument f.  The .PFC must also
    specify the "level" corresponding to each frequency.  Since this library
    was written explicitely to read/write data from an MTU5A with MTC-150
    magnetometers, the levels are are assigned based on [4] p. 2 (which are
    specific to an MTU5A + MTC-150 combo).

    Args:
        f (np.ndarray): array of frquencies (Hz) to evaluate G at
        TBL_fname (str): path to the .TBL file for the acquisition
        CLC_dir (str): path to the directory storing the .CLC files
        CLB_dir (str): path to the directory storing the .CLB files
        SYSCAL_fname (str): path to the SYSCAL executable
        PFC_fname (str): path where a .PFC file input to SYSCAL will be created
        cts_fname (str): path where the SYSCAL result will be created
        field_type (bool): type of response to calc (0 for box, 1 for sys)
        levels (array like): int array of recording levels for each freq in f

    Returns:
        G (np.ndarray): complex (freqs x 2+channels) array of channel TFs

    """

    if platform is "win32" or platform is "cygwin":
        raise NotImplementedError("This function isn't (yet) written" +
                                  "for win32 platforms!")

    # Construct the levels array (assumes sensors are MTC-150 coils)
    if levels is None:
        TS2_ind = np.logical_and(10400 >= f, f >= 900)
        TS3_ind = np.logical_and(900 > f, f >= 39)
        TS4_ind = 39 > f
        levels = np.empty(f.shape, dtype=int)
        levels[TS2_ind] = 2
        levels[TS3_ind] = 3
        levels[TS4_ind] = 4

    # Create the requisite .PFC file
    try:
        PFC_file = open(PFC_fname, mode="w")
    except IOError:
        print("Something went wrong opening the file!")

    PFC_file.write("FldType, %i,\n" % field_type)
    for freq, lev in zip(f, levels):
        PFC_file.write("Frequency, %i, %f\n" % (lev, freq))

    PFC_file.close()

    # Make sure the cal dirs point to INSIDE the dir
    if(CLB_dir[-1] != "/"):
        CLB_dir = CLB_dir + "/"
    if(CLC_dir[-1] != "/"):
        CLC_dir = CLC_dir + "/"

    # Since we assume that we're on a *NIX box, we need to convert
    # all paths for the cmd to windows format, using winepath.
    winepath_cmd = ["winepath", "-w"]
    PFC_wfname = subprocess.check_output(winepath_cmd + [PFC_fname])
    TBL_wfname = subprocess.check_output(winepath_cmd + [TBL_fname])
    CLB_wdir = subprocess.check_output(winepath_cmd + [CLB_dir])
    CLC_wdir = subprocess.check_output(winepath_cmd + [CLC_dir])
    cts_wfname = subprocess.check_output(winepath_cmd + [cts_fname])

    # For some reason, the resultant strings have an extra newline. Kill it.
    PFC_wfname = PFC_wfname.strip(b"\n")
    TBL_wfname = TBL_wfname.strip(b"\n")
    CLB_wdir = CLB_wdir.strip(b"\n")
    CLC_wdir = CLC_wdir.strip(b"\n")
    cts_wfname = cts_wfname.strip(b"\n")

    # Construct the cmd for calling SYSCAL
    SYSCAL_cmd = ["wine", SYSCAL_fname, PFC_wfname, TBL_wfname,
                  CLB_wdir + b"*.CLB", CLC_wdir + b"*.CLC", cts_wfname]

    # Call SYSCAL
    subprocess.check_call(SYSCAL_cmd)
    # proc = subprocess.Popen(SYSCAL_cmd)
    # proc.communicate("0")
    # proc.wait()

    G, field_type_out = read_cts(cts_fname)

    return G


def rec_to_fields(rec, TBL, mt_units=False):
    """
    Convert a record in instrument units to EM field units field units: V/m for
    E channels, T for H channels, following the formula provided in [5] p 6.

    rec is assumed to follow the format [e_x, e_y, h_x, h_y, h_z], although
    arbitrary orderings of the channels could be incorporated by using the
    channel indices stored in the .TBL file.

    This routine is currently only designed to work with data created by an
    MTU5A, although other types of hardware could be incorporated as well. The
    .TBL file's HW field stores type of hardware that created the data.

    Args:
        rec (np.ndarray): a (scans x channels) record of data
        TBL (np.ndarray): the associated TAG object for rec
        mt_units (bool): Whether or not to return mt units of mV/km and nT

    Returns:
        field (np.ndarray): a record of data in field units
    """

    # Compute intermediate scale factors for E and H channels
    E_0 = TBL["FSCV"]/(2**23)/TBL["EGN"]
    H_0 = TBL["FSCV"]/(2**23)/TBL["HGN"]/TBL["HATT"]/TBL["HNOM"]*1.0e-9

    if mt_units:
        E_0 = E_0*1.0e6
        H_0 = H_0*1.0e9

    # Do the conversion
    field = np.empty(rec.shape)
    field[:, 0] = E_0/TBL["EXLN"]*rec[:, 0]
    field[:, 1] = E_0/TBL["EYLN"]*rec[:, 1]
    field[:, 2:] = H_0*rec[:, 2:]

    return field


def get_TAG_datetime(TAG):
    """
    Returns the record start time from the tag TAG as a datetime object.

    It might be useful to also have this return the clock offset from GPS time.

    Args:
        TAG (np.ndarray): a single TAG (as in those returned from read_TSn)

    returns:
        stim (datetime.datetime): datetime object containing the record start.
    """

    stim = datetime.datetime(int(TAG["century"])*100 + int(TAG["year"]),
                             TAG["month"],
                             TAG["day"],
                             TAG["hour"],
                             TAG["minute"],
                             TAG["second"])

    return stim


def AMX_to_datetime(amx):
    """
    Convert the weird Phoenix AMX datetime format to a Python datetime object

    If the input amx array contains invalid values for a datetime object
    (e.g. a negative year) then the output datetime object will be set to
    datetime.date.min

    Args:
        amx: 8 element long int array following the Amx format in [2]
    Returns:
        dt: datetime.datetime object
    """
    # Some of the TABLE entries that use the AMX datetype will end up with
    # invalid datetimes. Just set these to the min allowable datetime
    try:
        dt = datetime.datetime(second=amx[0],
                               minute=amx[1],
                               hour=amx[2],
                               day=amx[3],
                               month=amx[4],
                               year=amx[5])
    except ValueError:
        dt = datetime.date.min

    return dt
