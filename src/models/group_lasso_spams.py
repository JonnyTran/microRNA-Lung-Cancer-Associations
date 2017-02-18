import numpy as np
import pandas
import parsimony.algorithms as algorithms
import parsimony.estimators as estimators
import parsimony.utils.start_vectors as start_vectors
import spams
from sklearn.base import BaseEstimator, ClassifierMixin


class LRSGLWrapper(BaseEstimator, ClassifierMixin):
    def __init__(self, l1=0.1, l2=0., gl=0.5, A=None, max_iter=1000):
        self.l1 = l1
        self.l2 = l2
        self.gl = gl
        self.A = A
        self.max_iter = max_iter

    def fit(self, X, y):

        if self.A is not None:
            self.model = estimators.LogisticRegressionL1L2GL(l1=self.l1, l2=self.l2, gl=self.gl, A=self.A,
                                                         algorithm=algorithms.proximal.FISTA(),
                                                         class_weight='auto',
                                                         algorithm_params=dict(max_iter=self.max_iter),
                                                         mean=False)
        else:
            self.model = estimators.LassoLogisticRegression(l=self.l1,
                                                            algorithm=algorithms.proximal.FISTA(),
                                                            class_weight='auto',
                                                            algorithm_params=dict(max_iter=self.max_iter),
                                                            mean=False)
        beta = start_vectors.ZerosStartVector().get_vector(X.shape[1])
        self.model.fit(X, y, beta=beta)
        return self

    def predict(self, X):
        return self.model.predict(X)

    def predict_proba(self, X):
        prob = np.hstack((1 - self.model.predict_probability(X), self.model.predict_probability(X)))
        return prob

    def decision_function(self, X):
        prob = np.hstack((self.model.predict_probability(X), 1 - self.model.predict_probability(X)))
        return prob


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
