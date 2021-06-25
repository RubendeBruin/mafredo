from mafredo import Hyddb1

"""The orcaflex .yml example file contains two vessel types.

Both are 100x30x4 barges, but one is created using double symmetry while the other 
does not employ symmetry at all.

Applying symmetry on the first one should give the same result as the second one.

"""

def test_read_both():

    filename = r'.\files\barge_100x30x4.yml'
    vessel_type_name = 'Vessel type1'
    iDraught = 0

    data = Hyddb1.create_from_orcaflex_yml(filename, vessel_type_name, iDraught)
    data.expand360_using_symmetry()

    filename = r'.\files\barge_100x30x4.yml'
    data2 = Hyddb1.create_from_orcaflex_yml(filename, 'Full_directions')

    data.assert_allclose_to(data2)
