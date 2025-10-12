import numpy as np
from mafredo import Hyddb1


def test_read_ofx_conventions(data_path):
    can_we_read = Hyddb1.create_from_orcaflex_yml(
        filename=data_path / "Friesland1.yml",
        vessel_type_name="Vessel_Frieslan1_T3m4",
        iDraught=0,
    )

    can_we_read._check_dimensions()

    omegas = np.array(
        [
            0.01,
            0.13089969,
            0.14279967,
            0.15707963,
            0.17453293,
            0.19634954,
            0.22439948,
            0.26179939,
            0.31415927,
            0.34906585,
            0.39269908,
            0.40276829,
            0.41336745,
            0.42453955,
            0.43633231,
            0.44879895,
            0.46199892,
            0.47599889,
            0.49087385,
            0.50670849,
            0.52359878,
            0.54165391,
            0.56099869,
            0.58177642,
            0.60415243,
            0.62831853,
            0.64114136,
            0.65449847,
            0.66842397,
            0.68295492,
            0.6981317,
            0.71399833,
            0.73060294,
            0.74799825,
            0.76624211,
            0.78539816,
            0.80553658,
            0.82673491,
            0.8490791,
            0.87266463,
            0.8975979,
            0.92399784,
            0.95199777,
            0.9817477,
            1.01341699,
            1.04719755,
            1.08330781,
            1.12199738,
            1.16355283,
            1.20830487,
            1.25663706,
            1.30899694,
            1.36590985,
            1.42799666,
            1.4959965,
            1.57079633,
            1.65346982,
            1.74532925,
            1.84799568,
            1.96349541,
            2.0943951,
            2.24399475,
            2.41660973,
            2.61799388,
            2.85599332,
            3.14159265,
            10.0,
        ]
    )

    can_we_read.regrid_omega(omegas)

    # get a temporary folder
    import tempfile

    tempdir = tempfile.gettempdir()
    filename = tempdir + "/Friesland1.dhyd"

    can_we_read.save_as(filename)

    # now read it back
    can_we_read_again = Hyddb1.create_from(filename)

    can_we_read_again._check_dimensions()
