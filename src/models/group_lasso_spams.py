import numpy as np
import pandas
import spams


class SPAMSClassifier:
    def __init__(self):
        pass

    def fit(self, X, y, params, groups=None):
        self.n_samples = X.shape[0]
        self.n_features = X.shape[1]

        W0 = np.zeros((self.n_features, 1), dtype=np.float, order="F")

        if groups != None:
            W, optim_info = spams.fistaFlat(y, X, W0, True, groups=np.array(groups, dtype=np.int32), **params)
        else:
            W, optim_info = spams.fistaFlat(y, X, W0, True, **params)

        self.W = W

    def predict(self, X):
        predicted = np.dot(X, self.W)
        np.place(predicted, predicted >= 0, [1])
        np.place(predicted, predicted < 0, [-1])
        return predicted

    def get_selected_miRNAs(self, mirna_list):
        selected_miRNAs = pandas.DataFrame([(mirna_list[i], self.W[i][0]) for i in self.W.nonzero()[0]],
                                           columns=['miRNA', 'coef'])
        selected_miRNAs = selected_miRNAs.reindex(
            selected_miRNAs['coef'].abs().sort(inplace=False, ascending=False).index)

        return selected_miRNAs
