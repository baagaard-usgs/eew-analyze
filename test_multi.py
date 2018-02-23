import multiprocessing
import numpy
import math
import time

def fn(obj, value):
    print("Starting {}.".format(value))
    time.sleep(4-0.1*value)
    print("{} {}".format(obj.name, value))
    return


class Foo(object):

    def __init__(self, nthreads):
        self.nthreads = nthreads
        self.name = "hello"
        return

    def run(self):
        values = numpy.arange(0.0, 20.01, 1.0)
        nvalues = values.shape[-1]
        
        pool = multiprocessing.Pool(self.nthreads)
        for x in values:
            pool.apply_async(fn, args=(self, x))
        pool.close()
        pool.join()
        return

if __name__ == "__main__":
    f = Foo(8)
    f.run()
    
