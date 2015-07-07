from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from interface import h5nb6xxInterface
from interface import h5nb6xx

class h5nb6xxInterfaceTests(TestWithMPI):

    def test1(self):
        instance = h5nb6xxInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()

