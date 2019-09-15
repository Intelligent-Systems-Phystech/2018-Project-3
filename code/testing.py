import datetime
import dill
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.stats import zscore
from data_processing import DataIterator
from models import Autoregression
from dtw_wrapper import DtwWrapper


default_params = {
    "nrow": 100000,
    "window_size": 10,
    "element_length": 200,
    "path": "../data/clustering/accelerometer2.csv",
    "overlap": 0,
    "shuffle": True,
    "sample_size": 1000,
    "chanel_num": 3,
    "repeat_num": 1
}


class TestFactory:

    def __init__(self, it=None, params=default_params, random_state=-1):
        self.data_path = params["path"]
        self.window_size = params["window_size"]
        self.element_length = params["element_length"]
        self.standart_sample_size = params["sample_size"]
        data = pd.read_csv(params["path"], header=0, nrows=params["nrow"])
        if it is not None:
            self.data_iterator = it
        else:
            self.data_iterator = DataIterator(data, self.element_length, params["shuffle"], random_state=random_state)
        self.chanel_num = params["chanel_num"]
        self.repeat_num = params["repeat_num"]
        self.results = []
        self.X = None
        self.classification_label = None
        counter = 0
        folders = os.listdir("../data/results/")
        while str(counter) in folders:
            counter += 1
        self.res_dir = "../data/results/{0}".format(counter)
        os.mkdir("../data/results/{0}".format(counter))


    def get_n(self, n):
        items = []
        labels = []
        infos = []
        for i in range(n):
            item, label, info = next(self.data_iterator)
            items.append(item)
            labels.append(label)
            infos.append(info)

        return items, np.array(labels), infos


    def test_dtw(self, dtw_function, distance_function, description=None, sample_size=-1, dump_result=False, dtw_args={}, cluster_dist="complete",
        external=False, external_distances_path=None):
        """
            work with cluster_dist from: complete, average, weighted
            the best is complete now
        """

        if sample_size < 0:
            sample_size = self.standart_sample_size

        if self.X is None:
            self.set_sample(sample_size)

        self.dtw_wrapper = DtwWrapper(self.X, hash(self.infos), dtw_function, distance_function,
            dtw_args=dtw_args,
            external=external,
            external_distances_path=external_distances_path)

        start_time = time.time()
        # Kostyl 
        # TODO: Fix this shit
        f = lambda x, y: self.dtw_wrapper.dist(x[0], y[0])

        for i in range(self.repeat_num):
            Z = linkage(
                np.linspace(0, len(self.X) - 1, len(self.X), dtype=int).reshape(-1, 1), cluster_dist, metric=f)
        end_time = time.time()
        t = "{0:.3f}".format((end_time - start_time) / self.repeat_num)
        print("Elapsed time: {0}".format(t))

        info = ClusteredInfoDTW(self.X, Z, self.element_length, self.dtw(dtw_function, distance_function, dtw_args),
                                self.chanel_num, label=self.classification_label)
        if dump_result:
            with open("{0}/{1}".format(self.res_dir, description), "wb") as f:
                dill.dump(info, f)

        self.results.append((description, datetime.datetime.now(), t))
 
        return info


    def ar_clustering(self, window_size=10, dump_result=False, description=None, normalization=False):
        start_time = time.time()
        ar_models = []
        for (i, x) in enumerate(self.X):
            ar = Autoregression(x[:, 0: self.chanel_num], window_size, normalization)
            ar.fit()
            ar_models.append(ar.coeffecients())
            print("Trained: {0}".format(i + 1))

        coeffs = []
        for model in ar_models:
            coeff = [np.array(np.concatenate([x.detach().numpy() for x in chanel], axis=0)) for chanel in model]
            coeffs.append(np.concatenate(coeff, axis=0))
        coeffs = np.array(coeffs)
        end_time = time.time()
        Z = linkage(coeffs, "ward")
        t = "{0:.3f}".format((end_time - start_time) / self.repeat_num)
        print("Elapsed time: {0}".format(t))
        info = ClusteredInfoAR(self.X, Z, self.element_length, coeffs, self.chanel_num, label=self.classification_label)
        if dump_result:
            with open("{0}/{1}".format(self.res_dir, description), "wb") as f:
                dill.dump(info, f)

        return info


    def dtw(self, dtw_function, distance_function, dtw_args):
        return lambda x, y: (dtw_function(x, y, distance_function, **dtw_args))


    def dump_result(self):
        with open("{0}/results.pkl".format(self.res_dir), "wb") as f:
            dill.dump(self.results, f)


    @staticmethod
    def load(file):
        with open(file, "rb") as f:
            return dill.load(f)


    def set_sample(self, n=-1, items=None, labels=None):
        if items is not None and labels is not None:
            self.X = items
            self.classification_label = labels
            return items, labels

        items, labels, infos = self.get_n(n)
        self.X = items
        self.classification_label = labels
        self.infos = tuple(infos)
        
        return items, labels


class ClusteredInfo:

    def __init__(self, X, Z, element_length, chanel_num, description=None, label=None):
        self.X = X
        self.Z = Z
        self.element_length = element_length
        self.description = description
        self.count = len(X)
        self.chanel_num = chanel_num
        if label is not None:
            self.label = label
    

    @staticmethod
    def load(f):
        with open(f, "rb") as f:
            return dill.load(f)


    def cluster(self, cluster_num):
        self.cluster_num = cluster_num
        self.clusters_labels = fcluster(self.Z, cluster_num, criterion="maxclust")
        unique_elements, counts_elements = np.unique(self.clusters_labels, return_counts=True)
        self.stats = pd.DataFrame(index=unique_elements, data=counts_elements, columns=["count"])
        self.stats = self.stats.sort_values(by="count", ascending=False)


    def visualize(self):
        fig = plt.figure(figsize=(25, 10))
        dn = dendrogram(self.Z)


    def show_chanel(self, ch=1, label=1, max_num=5):
        idxs = np.where(self.clusters_labels == label)[0]
        time = np.linspace(0, self.element_length - 1, self.element_length)
        if max_num is None:
            num_graphs = len(idxs)
        else:
            num_graphs = np.min([max_num, len(idxs)])
        fig, ax = plt.subplots(num_graphs, 1, sharex=True, squeeze=False, figsize=(15, 2.5 * num_graphs),
                            constrained_layout=False)
        fig.suptitle("Chanel {0} of {1} cluster".format(ch, label), y=1, fontsize = 14);
        for i in range(num_graphs):
            ax[i][0].plot(time, self.X[idxs[i]][:, ch]);
            ax[i][0].grid()
                
        plt.tight_layout();


    def clusters_compare_table(self, label, num_series=3, max_chanels=3, z_normalize=False):
        idxs = np.where(self.clusters_labels == label)[0]
        num_series = min(len(idxs), num_series)
        chanels = min(self.chanel_num, max_chanels)
        if num_series == 0:
            return

        fig, ax = plt.subplots(num_series, chanels, sharex=True, squeeze=False, figsize=(20, 2.5 * num_series), constrained_layout=False)
        fig.suptitle("Cluster {0}".format(label), y=1.01, fontsize = 14)
        
        t = np.linspace(0, self.element_length - 1, self.element_length, dtype=int)
        for df_id in range(num_series):
            for ch in range(chanels):
                x = self.X[idxs[df_id]][:, ch]
                if z_normalize:
                    x = zscore(x)
                ax[df_id][ch - 1].plot(t, x)
                if df_id == 0:
                    ax[df_id][ch].set_xlabel("ch{}".format(ch), fontsize=14)
                if ch == 1:
                    ax[df_id][ch].set_ylabel("ts{}".format(df_id + 1),fontsize=14)
                ax[df_id][ch].xaxis.set_label_position("top")
                ax[df_id][ch].xaxis.label.set_color("red")
                ax[df_id][ch].yaxis.label.set_color("red")
        
        plt.tight_layout();


    def comparing_at_one(self, label, num_series=3, max_chanel=3, z_normalize=False):
        idxs = np.where(self.clusters_labels == label)[0]
        num_series = min(len(idxs), num_series)
        if num_series == 0:
            return
        showed_chanels_num = min(self.chanel_num, max_chanel)
        if showed_chanels_num == 0:
            return

        fig, axs = plt.subplots(showed_chanels_num, 1, sharex=True, squeeze=False, figsize=(15, 2.5 * showed_chanels_num), constrained_layout=False)
        fig.suptitle("Cluster {0}".format(label), y=1.01, fontsize=14)

        for i in range(num_series):
            df_id = idxs[i]
            for (chanel_id, ax) in enumerate(axs):
                ax[0].set_title("Chanel {0}".format(chanel_id))
                y = self.X[df_id][:, chanel_id]
                if z_normalize:
                    y = zscore(y)
                ax[0].plot(y, label="y ts:{0}".format(df_id))

        plt.tight_layout();


class ClusteredInfoAR(ClusteredInfo):
    def __init__(self, X, Z, element_length, coeffs, chanel_num, description=None, label=None):
        ClusteredInfo.__init__(self, X, Z, element_length, chanel_num, description, label)
        self.coeffs = coeffs


class ClusteredInfoDTW(ClusteredInfo):

    def __init__(self, X, Z, element_length, dtw, chanel_num, description=None, label=None):
        ClusteredInfo.__init__(self, X, Z, element_length, chanel_num, description, label)
        self.paths = dict()
        self.dtw = dtw
  
    def allignment_to_random(self, label, num_series=5, max_chanel=3, show_real_y=False, z_normalize=False):
        idxs = np.where(self.clusters_labels == label)[0]
        num_series = min(len(idxs), num_series)
        if num_series == 0:
            return
        showed_chanels_num = min(self.chanel_num, max_chanel)
        if showed_chanels_num == 0:
            return
        align_to_id = np.random.choice(idxs)
        x_fixes = self.X[align_to_id]

        fig, axs = plt.subplots(showed_chanels_num, 1, sharex=True, squeeze=False, figsize=(15, 2.5 * showed_chanels_num), constrained_layout=False)
        fig.suptitle("Cluster {0}".format(label), y=1.01, fontsize = 14)

        for (chanel_id, ax) in enumerate(axs):
            x = x_fixes[:, chanel_id]
            if z_normalize:
                x = zscore(x)
            ax[0].plot(x, "black", label="X", linewidth=3)

        for i in range(num_series):
            df_id = idxs[i]
            if df_id == align_to_id:
                continue
            if align_to_id in self.paths and df_id in self.paths[align_to_id]:
                path = self.paths[align_to_id][df_id]
            elif df_id in self.paths and align_to_id in self.paths[df_id]:
                path = (self.paths[df_id][align_to_id][1], self.paths[df_id][align_to_id][0])
            else:
                if df_id not in self.paths:
                    self.paths[df_id] = dict()
                path = self.dtw(x_fixes[:, 0 : self.chanel_num], 
                    self.X[df_id][:, 0 : self.chanel_num])[1]
                self.paths[df_id][align_to_id] = path

            for (chanel_id, ax) in enumerate(axs):
                y = self.X[df_id][:, chanel_id]
                if z_normalize:
                    y = zscore(y)
                y_new = pd.DataFrame([y[i - 1] for i in path[1]], index=path[0])
                y_new = y_new.groupby(y_new.index).mean().values.reshape(-1)
                ax[0].plot(y_new, label="y_new ts:{0}".format(df_id))
                if show_real_y:
                    ax[0].plot(y, label="y ts:{0}".format(df_id))
                    for (map_x, map_y) in np.array(path).transpose():
                        plt.plot([map_x, map_y], [x[map_x], y[map_y]], "black", linewidth=0.3)
                # ax[0].legend(loc=1)

        plt.tight_layout();
