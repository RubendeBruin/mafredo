"""
Orcaflex / OrcaWave and Symmetry
=================================
Example of
- reading hydrodynamic data from OrcaWave.
- applying symmetry
(unrelated subjects)

Orcawave files can not be read directly. To obtain the data first import the data into an orcaflex model (vessel-type)
and then save that model as .yml

MaFreDo can now read the data from the .yml file

Symmetry is nice when a file needs to be small but may not be needed or even wanted when performing calculations. The
expand360_using_symmetry method uses the symmetry setting of the database to expand the dataset to cover the full 360
degrees of headings.

(Note: Orcaflex and OrcaWave are products by Orcina : https://www.orcina.com/orcaflex/ )
"""

from mafredo import Hyddb1

try:
    data = Hyddb1.create_from_orcaflex_yml(
        filename="./barge_100x30x4_q1.yml",
        vessel_type_name="Vessel type1",
        iDraught=0
    )

    print(data.symmetry)
    data.plot(adm=False, damp=False)

    data.expand360_using_symmetry()

    print(data.symmetry)
    data.plot(adm=False, damp=False)

except:
    pass # read-the-docs



