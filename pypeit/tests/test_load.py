"""
Module to run tests loading
"""
import os
import pytest
from pypeit import specobjs
from pypeit import specobj

from pypeit.tests.tstutils import data_path

# TODO: Why is this not part of cooked?
# IF THIS TEST IS FAILING BECAUSE THE SPEC1D DATAMODEL UPDATED,
#   RUN source copy_cooked_files.src in files/
#   *After* you have run the Dev Suite tests
def test_load_specobjs():
    spec_file = data_path('spec1d_r153-J0025-0312_KASTr_20150123T025323.850.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(spec_file)

    # Test
    assert isinstance(sobjs, specobjs.SpecObjs)
    assert len(sobjs[0].BOX_COUNTS) == 1200

    assert isinstance(sobjs[0], specobj.SpecObj)


