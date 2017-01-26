import numpy as np
import pandas
import spams


class SPAMSClassifier:
    def __init__(self):
        pass

    def fit(self, X_train, Y_train, params, groups=None):
        print "FISTA Flat:", params
        self.n_samples = X_train.shape[0]
        self.n_features = X_train.shape[1]
        self.groups = groups
        self.params = params

        if params['loss'] == 'multi-logistic':
            W0 = np.zeros((self.n_features, len(np.unique(Y_train))), dtype=np.float, order="F")
        else:
            W0 = np.zeros((self.n_features, 1), dtype=np.float, order="F")

        if groups is not None:
            W, optim_info = spams.fistaFlat(Y_train, X_train, W0, True, groups=np.array(groups, dtype=np.int32),
                                            **params)
        else:
            W, optim_info = spams.fistaFlat(Y_train, X_train, W0, True, **params)

        self.W = W

        return optim_info

    def predict(self, X_test):
        if self.params['loss'] == 'multi-logistic':
            predicted = np.dot(X_test, self.W)
            predicted = np.asanyarray(predicted.argmax(axis=1), order="F")
        else:
            predicted = np.dot(X_test, self.W)
            np.place(predicted, predicted >= 0, [1])
            np.place(predicted, predicted < 0, [-1])

        return predicted

    def get_selected_miRNAs(self, mirna_list):
        if self.params['loss'] == 'multi-logistic':
            selected_miRNAs = pandas.DataFrame([(mirna_list[i], max(self.W[i].max(), self.W[i].min(), key=abs)) for i in
                                                np.unique(self.W.nonzero()[0])],
                                               columns=['miRNA', 'coef'])
        else:
            selected_miRNAs = pandas.DataFrame([(mirna_list[i], self.W[i][0]) for i in self.W.nonzero()[0]],
                                               columns=['miRNA', 'coef'])

        selected_miRNAs = selected_miRNAs.reindex(
            selected_miRNAs['coef'].abs().sort(inplace=False, ascending=False).index)

        return selected_miRNAs
