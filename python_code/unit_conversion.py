import numpy as np


class UnitConversion:
    in2cm = 2.54  # 1in = 2.54cm
    ft2cm = 30.48  # 1ft = 30.48cm
    cm2m = 1 / 100

    # rpm2rads = (pi / 30)  # 60 RPM = 1 RPsecond = 1 Hz = 2*pi rad/s
    @staticmethod
    def in2m(inches):
        in2m = 0.0254
        m = in2m * inches
        return m

    @staticmethod
    def ft2m(ft):
        ft2m = 0.3048
        m = ft2m * ft
        return m

    @staticmethod
    def lb2kg(lb):
        lb2kg = 0.453592
        kg = lb2kg * lb
        return kg

    # deg2rad = np.pi / 180
    # rad2deg = 180 / np.pi

    @staticmethod
    def deg2rad(deg):
        deg2rad = np.pi / 180
        rad = deg2rad * deg
        return rad

    @staticmethod
    def rad2deg(rad):
        rad2deg = 180 / np.pi
        deg = rad2deg * rad
        return deg

    @staticmethod
    def rpm2rads(rpm):
        rpm2rads = np.pi / 30
        rads = rpm2rads * rpm
        return rads

    @staticmethod
    def rads2rpm(rads):
        rads2rpm = 30 / np.pi
        rpm = rads2rpm * rads
        return rpm