# Third-party modules
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from lazypredict.Supervised import LazyClassifier

# scikit-learn modules
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import NearestCentroid
from sklearn.naive_bayes import GaussianNB
from sklearn.decomposition import PCA
from sklearn.metrics import (
    precision_score,
    recall_score,
    accuracy_score,
    balanced_accuracy_score,
    roc_auc_score,
    f1_score,
    confusion_matrix,
    ConfusionMatrixDisplay
)
from sklearn.model_selection import (
    KFold, 
    cross_validate, 
    train_test_split,
    GridSearchCV
    )

# Custom modules
from logger import logger
from network import Network
from embedding import Embedding
from interaction import Interaction

class MachineLearning:

    def __init__(self, X: pd.DataFrame, y: pd.Series):
        self.seed = 1999
        self.X, self.index = self.shuffle(X)
        self.y, self.index = self.shuffle(y)

        #self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size = 0.25, random_state = self.seed)

    def shuffle(self, df: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
        size = df.index.get_level_values('Index')[-1] + 1
        shuffle_index = pd.DataFrame(np.arange(size)).sample(frac = 1, random_state = self.seed)[0]
        shuffle_df = df.loc[shuffle_index]
        return shuffle_df, shuffle_index
    
    def kfold(self, k: int):
        cv = KFold(n_splits = k)
        for train_index, test_index in cv.split(self.index):

            train = self.index.iloc[train_index]
            test  = self.index.iloc[test_index]

            X_train = self.X.loc[train]
            X_test  = self.X.loc[test]
            y_train = self.y.loc[train]
            y_test  = self.y.loc[test]
            
            yield X_train, X_test, y_train, y_test

    def lazy(self, kfold: int = 5):
        
        # Initialize list of prediction data frames
        predictions_list = []

        # K-fold to ensure consistent results
        for X_train, X_test, y_train, y_test in self.kfold(k = kfold):           

            # LazyPredict
            clf = LazyClassifier(random_state = self.seed)
            models, predictions = clf.fit(X_train, X_test, y_train, y_test)
            predictions_list += [predictions]
            print(predictions)

        # Unify data frames
        df = pd.concat(predictions_list)
        df = df.groupby(df.index).mean().sort_values(by = 'Balanced Accuracy', ascending = False)
        print(df)

        # Bar plot 

        df['Balanced Accuracy'].plot.bar()
        plt.ylim(0.49, 1)
        plt.tight_layout()
        plt.show()

        return df
    
    @staticmethod
    def result_metrics():
        results = {
            'Test precision' : [],
            'Train precision' : [],
            'Test specificity' : [],
            'Train specificity' : [],
            'Test accuracy': [],
            'Train accuracy' : [],
            'Test balanced accuracy' : [],
            'Train balanced accuracy' : [],
            'Test F1 score' : [],
            'Train F1 score' : [],
        }

        return results

    @staticmethod
    def param_grids():
        # Initialize
        param_grid = {}

        # Logistic Regression
        param_grid['LogisticRegression'] = {
            'penalty' : [None, 'l1', 'l2', 'elastic net'],
            'C' : np.arange(0, 1.1, 0.1),
            'dual': [True, False]
        }

        # Decision Tree
        param_grid['DecisionTree'] = {
            'criterion' : ['gini', 'entropy', 'log_loss'],
            'splitter' : ['best', 'random'],
            'max_depth': [None, 300, 150, 100, 70, 40, 20, 10, 5, 2, 1],
            #'min_samples_split': [2, 3, 4, 10, 50, 100, 1000],
            #'min_samples_leaf' : [2, 3, 4, 10, 50, 100, 1000],
            'max_features' : [None, 'sqrt', 'log2'],
            'max_leaf_nodes' : [None, 1000, 100, 5],
            'ccp_alpha' : [0, 0.2, 0.6, 1.0, 1.4, 5]
        }

        # Random Forest
        param_grid['RandomForest'] = {
            'n_estimators' : [1, 10, 20, 50, 75, 100, 200, 500],
            'criterion' : ['gini', 'entropy', 'log_loss'],
            'max_features' : ['sqrt', 'log2', None]
        }

        # Nearest Centroid
        param_grid['NearestCentroid'] = {
            'metric' : ['euclidean', 'manhattan'],
            'shrink_threshold' : [0, 0.5, 0.7, 1, 1.2, 1.5, 1.8, 2.5, 5, 10, 100, 1000]
        }

        # Gaussian Naïve Bayes
        param_grid['GaussianNB'] = {
            'var_smoothing' : [1e-15, 1e-9, 1e-5, 1e-3, 1e-2, 0.1, 0.4, 0.6, 0.8, 0.99]
        }

        return param_grid
    
    def hyperparameter_grid_search(self, algorithm):
        # Choose algorithm
        match algorithm:
            case 'LogisticRegression':
                clf = LogisticRegression(random_state = self.seed, max_iter = 10_000)
            case 'DecisionTree':
                clf = DecisionTreeClassifier(random_state = self.seed)
            case 'RandomForest':
                clf = RandomForestClassifier(random_state = self.seed)
            case 'NearestCentroid':
                clf = NearestCentroid()
            case 'GaussianNB':
                clf = GaussianNB()

        # Grid search
        grid_search = GridSearchCV(
            estimator = clf, 
            param_grid = MachineLearning.param_grids()[algorithm], 
            scoring = 'balanced_accuracy', 
            cv = 5,
            verbose = 2,
            n_jobs = -1
            )
        
        # Seach best parameters
        grid_search.fit(self.X_train, self.y_train)

        # Logger
        logger.info(f'Best parameter set: {grid_search.best_params_}')
        logger.info(f'Best metric score: {grid_search.best_score_}')

        return grid_search.best_params_, grid_search.best_score_, grid_search.best_estimator_
    
    def evaluate_model(self, clf, split):
    
        # Fit model
        X_train, X_test, y_train, y_test = split
        clf.fit(X_train.to_numpy(), y_train.to_numpy())
        y_pred_test = clf.predict(X_test)
        y_pred_train = clf.predict(X_train)

        # Interaction classes
        interactors_train = X_train.index.get_level_values('Interaction').to_series().apply(lambda x: [x.A.uniprot_id, x.B.uniprot_id]).explode()
        interactors_test  = X_test.index.get_level_values('Interaction').to_series().apply(lambda x: [x.A.uniprot_id, x.B.uniprot_id]).explode()
        print(set(interactors_train) - set(interactors_test))
    
        # Calculate performance metrics
        results = MachineLearning.result_metrics()
        results['Test precision'] = precision_score(y_test, y_pred_test)
        results['Train precision'] = precision_score(y_train, y_pred_train)
        results['Test specificity'] = recall_score(y_test, y_pred_test, pos_label = 0)
        results['Train specificity'] = recall_score(y_train, y_pred_train, pos_label = 0)
        results['Test accuracy'] = accuracy_score(y_test, y_pred_test)
        results['Train accuracy'] = accuracy_score(y_train, y_pred_train)
        results['Test balanced accuracy'] = balanced_accuracy_score(y_test, y_pred_test)
        results['Train balanced accuracy'] = balanced_accuracy_score(y_train, y_pred_train)
        results['Test F1 score'] = f1_score(y_test, y_pred_test)
        results['Train F1 score'] = f1_score(y_train, y_pred_train)

        for i,j in results.items():
            print(i, j)

        # Confusion matrix
        cm = confusion_matrix(y_test, y_pred_test)
        ConfusionMatrixDisplay(cm).plot()
        plt.show()

        return results, cm
    
    def pca(self):
        # PCA
        pca = PCA()
        scores = pca.fit_transform(self.X)
        pc1, pc2, *_= pca.explained_variance_ratio_ * 100
        columns = [f'PC1 {pc1} %', f'PC2 {pc2} %']
        pca_df = pd.DataFrame(scores[:, 0:2], columns = columns, index = self.X.index)
        pca_df['y'] = self.y.to_list()
        pca_df['Species'] = pca_df.index.to_series().apply(lambda x: Interaction(x).species)
        
        # Plot PCA
        sns.scatterplot(
            x=columns[0], 
            y=columns[1], 
            data=pca_df,
            hue = 'Species',
            style = 'y',
            legend = True,
            )
        plt.show()

        # Scree plot
        PC_values = np.arange(pca.n_components_) + 1
        fig = plt.plot(PC_values, pca.explained_variance_ratio_, 'ro-', linewidth=2)
        plt.title('Scree Plot')
        plt.xlabel('Principal Component')
        plt.ylabel('Proportion of Variance Explained')
        plt.show()


if __name__ == '__main__':

    intact_network = Network(database = 'IntAct', type = 'MADS vs MADS', species = 'all')
    e = Embedding(intact_network)
    X = e.multiindex(e.maxabs_scaler(e.length_scaler(e.kmer_sum())) * 100)
    y = e.response(X)

    ml = MachineLearning(X, y)
    ml.pca()

    #params, score, model = ml.hyperparameter_grid_search(algorithm = 'NearestCentroid')
    #results, cm = ml.evaluate_model(model)
    #print(results)
    '''
    from pulearn import ElkanotoPuClassifier
    logreg = NearestCentroid()
    pu_estimator = ElkanotoPuClassifier(estimator=logreg)
    #pu_estimator.fit(X.to_numpy(), y.to_numpy())
    arabidopsis = n.df[(n.df['Species A'] == 3702) & (n.df['Species B'] == 3702)]
    X = n.kmer_sum_embedding(arabidopsis)
    y = arabidopsis['+-']
    ml = MachineLearning(X, y)
    print(ml.evaluate_model(pu_estimator))

    ml = MachineLearning(X, y)
    ml.pca()
    df = ml.lazy()
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'LogisticRegression')
    results, cm = ml.evaluate_model(model)
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'DecisionTree')
    results, cm = ml.evaluate_model(model)
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'NearestCentroid')
    results, cm = ml.evaluate_model(model)
    

    # Arabidopsis data
    arabidopsis = n.df[(n.df['Species A'] == 3702) & (n.df['Species B'] == 3702)]
    X = n.kmer_sum_embedding(arabidopsis)
    y = arabidopsis['+-']

    ml = MachineLearning(X, y)
    ml.pca()
    df = ml.lazy()
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'LogisticRegression')
    results, cm = ml.evaluate_model(model)
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'DecisionTree')
    results, cm = ml.evaluate_model(model)
    params, score, model = ml.hyperparameter_grid_search(algorithm = 'GaussianNB')
    results, cm = ml.evaluate_model(model)
    '''

