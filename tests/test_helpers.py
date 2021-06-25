from mafredo.helpers import f10

def test_f10():
    numbers = [1, 2, -1, 17.234, 1e10, 1e10 - 9, 1e10 + 1, 1325123551512511.0, -1.1641532182693481e-10,
               -1325123551512511.0,
               7.275957614183426e-12]

    for n in numbers:
        s = f10(n)
        assert len(s) == 10
