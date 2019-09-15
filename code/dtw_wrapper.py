import numpy as np
import dill

from os.path import exists
from time import time, sleep
from threading import Thread, Lock

class DtwWrapper:
    """Wrapper for dtw functions
    Number of series must divide by 4

    """
    num_th = 4

    def __init__(self, items, items_hash, dtw_function, distance_function,
        external,
        external_distances_path,
        ch_num=3,
        dtw_args={}):
        print(external, external_distances_path)
        self.items = items
        self.chanel_num = ch_num
        str_args = "".join(["{}{}".format(key, dtw_args[key]) for key in dtw_args])
        self.dtw_name = "{0}{1}{2}{3}{4}".format(dtw_function.__name__, distance_function.__name__, str_args, items_hash, ch_num)
        
        if external:
            self.distances = np.genfromtxt(external_distances_path)

            return

        self.series_distance = self.dtw_dist(dtw_function, distance_function, dtw_args)
        if exists("../data/distances/{0}.csv".format(self.dtw_name)):
            print("Loaded")
            self.distances = np.genfromtxt("../data/distances/{0}.csv".format(self.dtw_name))
        else:
            self.distances = np.full((len(items), len(items)), -1., dtype=float)

    def dist(self, x_index, y_index):
        x_index = int(x_index)
        y_index = int(y_index)
        if self.distances[x_index, y_index] != -1:
            return self.distances[x_index, y_index]

        self.distances[x_index, y_index] = self.series_distance(
            self.items[x_index],
            self.items[y_index]
        )

        return self.distances[x_index, y_index]

    def dtw_dist(self, dtw_function, distance_function, dtw_args):
        return lambda x, y: dtw_function(x, y, distance_function, **dtw_args)[0]

    def dump(self):
        np.savetxt("../data/distances/{0}.csv".format(self.dtw_name), X=self.distances)

