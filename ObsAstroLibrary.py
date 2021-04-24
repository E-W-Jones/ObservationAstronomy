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
    Calculate the resolution of a telescope, in degrees/arcseconds.

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

def main():
    """Run some examples of the functions included."""

    # --- Resolution ---
    print("Here are some examples of the `resolution` function:")
    print("For a radio telescope:")
    print(f"λ = {λ:=0.06}m, D = {D:=25}m, returns:\n{resolution(λ, D)}")
    print("For a visible light telescope:")
    print(f"λ = {λ:=5E-7}m, D = {D:=4.2}m, returns:\n{resolution(λ, D)}")


if __name__ == "__main__":
    main()