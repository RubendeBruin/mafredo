from mafredo import *
from numpy.testing import assert_allclose

# ----- add frequencies ----

def test_add_rao():

    # First create the object
    cylinder = Hyddb1.create_from_capytaine('./../docs/gallery/open_cylinder.nc')
    cylinder_refined = Hyddb1.create_from_capytaine('./../docs/gallery/open_cylinder_2.nc')

    rao1 = cylinder.force_rao(0)
    rao2 = cylinder_refined.force_rao(0)

    rao1.add(rao2)

    expect = [0.01      , 0.02      , 0.04      , 0.06      , 0.08      ,
       0.1       , 0.15      , 0.2       , 0.25      , 0.3       ,
       0.35      , 0.4       , 0.45      , 0.5       , 0.55      ,
       0.6       , 0.65      , 0.7       , 0.8       , 0.9       ,
       1.        , 1.1       , 1.2       , 1.4       , 1.5       ,
       1.52631579, 1.55263158, 1.57894737, 1.6       , 1.60526316,
       1.63157895, 1.65789474, 1.68421053, 1.71052632, 1.73684211,
       1.76315789, 1.78947368, 1.8       , 1.81578947, 1.84210526,
       1.86842105, 1.89473684, 1.92105263, 1.94736842, 1.97368421,
       2.        , 4.        ]

    assert_allclose(rao1['omega'].values , expect)


def test_add_database():

    # First create the object
    cylinder = Hyddb1.create_from_capytaine('./../docs/gallery/open_cylinder.nc')
    cylinder_refined = Hyddb1.create_from_capytaine('./../docs/gallery/open_cylinder_2.nc')
    cylinder_further_refined = Hyddb1.create_from_capytaine('./../docs/gallery/open_cylinder_3.nc')

    # merge
    cylinder.add(cylinder_refined)
    cylinder.add(cylinder_further_refined)

    expect = [0.01      ,0.02       ,0.04       ,0.06      , 0.08       ,0.1,
             0.15       ,0.2        ,0.25       ,0.3       , 0.35       ,0.4,
             0.45       ,0.5        ,0.55       ,0.6       , 0.65       ,0.7,
             0.8        ,0.9        ,1.         ,1.1       , 1.2        ,1.4,
             1.6        ,1.8        ,2.         ,4.        , 1.5        ,1.52631579,
             1.55263158 ,1.57894737 ,1.60526316 ,1.63157895, 1.65789474 ,1.68421053,
             1.71052632 ,1.73684211 ,1.76315789 ,1.78947368, 1.81578947 ,1.84210526,
             1.86842105 ,1.89473684 ,1.92105263 ,1.94736842, 1.97368421 ,2.,
             1.6        ,1.60512821 ,1.61025641 ,1.61538462, 1.62051282 ,1.62564103,
             1.63076923 ,1.63589744 ,1.64102564 ,1.64615385, 1.65128205 ,1.65641026,
             1.66153846 ,1.66666667 ,1.67179487 ,1.67692308, 1.68205128 ,1.68717949,
             1.69230769 ,1.6974359  ,1.7025641  ,1.70769231, 1.71282051 ,1.71794872,
             1.72307692 ,1.72820513 ,1.73333333 ,1.73846154, 1.74358974 ,1.74871795,
             1.75384615 ,1.75897436 ,1.76410256 ,1.76923077, 1.77435897 ,1.77948718,
             1.78461538 ,1.78974359 ,1.79487179 ,1.8       ]


    assert_allclose(cylinder.frequencies, expect)
