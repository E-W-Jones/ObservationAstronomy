import math
import numpy as np


def sex_to_dec(coords):
    """
    Convert from sexegesimal coordinates to decimal degrees.

    If you want to convert just one coordinate, just put the other as 0.

    Parameters
    ----------
    coords : A list/tuple of the RA and Dec, as a string (in order to get ':')

    Returns
    -------
    A list, the first element is the decimal degree RA, the second is Declination
    """
    right_asc_inp, declination_inp = coords

    right_asc = [abs(float(x)) for x in right_asc_inp.split(":")]
    right_asc = [ang/(60**i) for (i,ang) in enumerate(right_asc)]
    right_asc = 15*sum(right_asc)

    declination = [abs(float(x)) for x in declination_inp.split(":")]
    declination = sum([ang/(60**i) for (i,ang) in enumerate(declination)])

    if right_asc_inp.startswith('-'):
        right_asc *= -1
    if declination_inp.startswith('-'):
        declination *= -1

    return [right_asc, declination]

def dec_to_sex(coords):
    """
    Convert from decimal degrees to sexegesimal coordinates.

    If you want to convert just one coordinate, just put the other as 0.

    Parameters
    ----------
    coords : A list of the decimal degree coordinates, either floats, ints, or strings should work.

    Returns
    -------
    A list of the RA and Dec, as

    """
    right_asc_inp, declination_inp = coords
    right_asc = abs(float(right_asc_inp)/15)  # Converts from degrees to decimal hours
    h, delta = divmod(right_asc, 1)  # Splits into the integer hours, and the decimal to be converted to minutes
    m, sig = divmod(60*delta, 1)  # Splits into integer minutes, and decimal to be converted into seconds
    s = 60*sig
    right_asc = "{0:n}:{1:n}:{2}".format(h, m, s)

    declination = abs(float(declination_inp))
    d, delta = divmod(declination, 1)
    m, sig = divmod(60*delta, 1)
    s = 60*sig
    declination = "{0:n}:{1:n}:{2}".format(d, m, s)

    if type(right_asc_inp) is str:
        if right_asc_inp.startswith('-'):
            right_asc *= -1
    else:
        if right_asc_inp < 0:
            right_asc *= -1

    if type(declination_inp) is str:
        if declination_inp.startswith('-'):
            declination *= -1
    else:
        if declination_inp < 0:
            declination *= -1

    return [right_asc, declination]

def resolution(wavelength, diameter):
    """
    Calculate the resolution of a telescope, in degrees and arcseconds.

    Uses the approximation resolution ~ wavelength / diameter.

    Parameters
    ----------
    wavelength : the wavelength of light, in metres, as a float or int.
    diameter : the diameter of the telescope, in metres, as a float or int.

    Returns
    -------
    A tuple of the decimal degree representation of the resolution, as a float,
    and the sexegesimal representation, as a string.

    """
    resolution_degrees = math.degrees(wavelength / diameter)
    resolution_sex = dec_to_sex([0, resolution_degrees])[1]
    return resolution_degrees, resolution_sex

def arcsec_to_rad(angle):
    """
    Convert an angle from arcseconds to radians.

    Paramters
    ---------
    angle : the angle, in arcseconds, a float or int.

    Returns
    -------
    The angle in radians.
    """
    # pi/180 rad per degree, 1/3600 degrees per arcsec
    return angle * math.pi / (180 * 3600)

def rad_to_arcsec(angle):
    """
    Convert an angle from radians to arcseconds.

    Parameters
    ----------
    angle : the angle, in radians, a float or int.

    Returns
    -------
    The angle in arcseconds.
    """
    # 180/pi degree per rad, 3600 arcsec per degree
    return angle * (180 / math.pi) * 3600

def appmag_from_dist(M, d, A=0):
    """
    Calculate the apparent magnitude from the absolute magnitude and distance.

    Uses M = m - 5log_10(d / 10pc) - A, where:
        M is abs. magnitude
        m is apparent magnitude
        d is the distance, in parsecs
        A is the extinction

    Parameters
    ----------
    M : Absolute magnitude, as an int or float.
    d : Distance, in parsecs, as an int or float.
    A : Extinction in terms of magnitude, as an int or float. If not specified,
    the default value of 0 is used to represent no extinction.

    Returns
    -------
    m, the apparent magnitude, as a float.
    """
    return M + 5*math.log10(d / 10) + A

def absmag_from_dist(m, d, A=0):
    """
    Calculate the absolute magnitude from the apparent magnitude and distance.

    Uses M = m - 5log_10(d / 10pc) - A, where:
        M is abs. magnitude
        m is apparent magnitude
        d is the distance, in parsecs
        A is the extinction

    Parameters
    ----------
    m : Apparent magnitude, as an int or float.
    d : Distance, in parsecs, as an int or float.
    A : Extinction in terms of magnitude, as an int or float. If not specified,
    the default value of 0 is used to represent no extinction.

    Returns
    -------
    M, the absolute magnitude, as a float.
    """
    return m - 5*math.log10(d / 10)

def dist_from_mag(m, M, A=0):
    """
    Calculate the distance from the absolute and apparent magnitudes.

    Uses M = m - 5log_10(d / 10pc) - A, where:
        M is abs. magnitude
        m is apparent magnitude
        d is the distance, in parsecs
        A is the extinction

    Parameters
    ----------
    m : Apparent magnitude, as an int or float.
    M : Absolute magnitude, as an int or float.
    A : Extinction in terms of magnitude, as an int or float. If not specified,
    the default value of 0 is used to represent no extinction.

    Returns
    -------
    d, the distance in parsecs, as a float.
    """
    exponent = (m - M - A) / 5 + 1
    return math.pow(10, exponent)

def main():
    """Run some examples of the functions included."""


    # --- Coordinate Conversions ---
    print(("\n\nHere are some examples of the `sex_to_dec`"
           " and `dec_to_sex` functions:"))

    print("Going from decimal degrees to sexegesimal:")
    coords = [25, 27.55]
    RA, Dec = dec_to_sex(coords)
    print(f"\t(RA, Dec) = ({coords[0]}, {coords[1]}) -> ({RA}, {Dec})")

    print("Going from sexegesimal to decimal degrees:")
    coords = ["03:32:30.14", "-27:25:10.5"]
    RA, Dec = sex_to_dec(coords)
    print(f"\t(RA, Dec) = ({coords[0]}, {coords[1]}) -> ({RA}, {Dec})")


    # --- Resolution ---
    print("\n\nHere are some examples of the `resolution` function:")
    print("For a radio telescope:")
    λ = 0.06; D = 25
    print(f"\tλ = {λ}m, D = {D}m, returns: {resolution(λ, D)}")

    print("For a visible light telescope:")
    λ = 5E-7; D = 4.2
    print(f"\tλ = {λ}m, D = {D}m, returns: {resolution(λ, D)}")


    # --- Converting between arcseconds and radians ---
    print(("\n\nHere are some examples of the `arcsec_to_rad`"
           " and `rad_to_arcsec` functions:"))
    print("Converting arcseconds to radians:")
    θ = 0.5
    print(f'\tθ = {θ}" -> {arcsec_to_rad(θ)} rad')

    print("Converting from radians to arcseconds:")
    θ = 0.0001
    print(f'\tθ = {θ} rad -> {rad_to_arcsec(θ)}"')


    # --- Magnitudes and Distances ---
    print(("\n\nHere are some examples of the `appmag_from_dist`,"
           "`absmag_from_dist`, and `dist_from_mag` functions:"))

    print("Going from a distance and an absolute magnitude:")
    M = 15; d = 300
    print(f"\tM = {M}, d = {d} -> m = {appmag_from_dist(M=15, d=300)}")

    print("Going from a distance and an apparent magnitude:")
    m = 5; d = 100
    print(f"\tm = {m}, d = {d} -> M = {absmag_from_dist(m, d)}")

    m = 15; M = 5.4; A = 3.7
    print("Calculating distance, not accounting for extinction:")
    print(f"\tm = {m}, M = {M} -> d = {dist_from_mag(m, M)}")

    print("Calculating distance, accounting for extinction:")
    print(f"\tm = {m}, M = {M}, A = {A} -> d = {dist_from_mag(m, M, A)}")


if __name__ == "__main__":
    main()