import datreant.core as dtr
import mdsynthesis as mds


def test_discover(tmpdir):
    with tmpdir.as_cwd():

        sims = ('a dir/inky',
                'something/blinky',
                'pinky',
                'something/else/clyde')

        treants = ('something/pacman',)

        sims = [mds.Sim(name) for name in sims]
        treants = [dtr.Treant(name) for name in treants]

        b = mds.discover('.')

        assert len(b) == 4

        for sim in sims:
            assert sim in b

        for treant in treants:
            assert treant not in b

        dtrb = dtr.discover('.')

        assert len(dtrb) == 5

        for treant in sims + treants:
            assert treant in dtrb
