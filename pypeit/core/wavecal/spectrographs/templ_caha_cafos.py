""" Generate the wavelength templates for P200/DBSP"""
import os

from pypeit.core.wavecal import templates



# ##############################

def caha_cafos_b100_450(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_b100.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [6000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_b100_450.fits')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)


def caha_cafos_g100_635(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_g100.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [8000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_g100_635.fits')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)


def caha_cafos_r100_745(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_r100.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [9000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_r100_745.fits')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

def caha_cafos_b200_510(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_b200.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [7000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_b200_510.fits')


    #    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_1_01.json')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)


def caha_cafos_g200_625(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_g200.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [8000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_g200_625.fits')


    #    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_1_01.json')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)


def caha_cafos_b400_560(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_b400.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [8000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_b400_560.fits')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

def caha_cafos_r400_785(overwrite=False):  

    binspec = 1
    outroot = 'caha_cafos_r400.fits'
    #
    ifiles = [0]
    slits = [750]  # Be careful with the order..
    lcut = [11000]
    wfile = os.path.join(templates.template_path, 'CAHA_CAFOS/caha_cafos_r400_785.fits')

    templates.build_template(wfile, slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

if __name__ == '__main__':
#    caha_cafos_g200_625(overwrite=True)
#    caha_cafos_b100_450(overwrite=True)
#    caha_cafos_g100_635(overwrite=True)
#    caha_cafos_b400_560(overwrite=True)
#    caha_cafos_r400_785(overwrite=True)
#    caha_cafos_r100_745(overwrite=True) 
    caha_cafos_b200_510(overwrite=True)