#!/usr/bin/env python
# coding: utf-8

# ## Can we predict better fish decay if we include human data?

# In[1]:


import numpy as np
import pandas as pd
from optimalcodon.projects.optrules import utils
from optimalcodon.codons import codon_composition, generate_codons

# sklearn
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import Normalizer, StandardScaler, OneHotEncoder, PolynomialFeatures
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import r2_score, make_scorer



# ## Prepare Data

# In[2]:


decay, orfs = utils.load_endogenous_data("../../../190115-CodonOptimalityRules/data/190416-DecayProfilesANDEndogenousDataSets/")
decay.sample(10)


# In[3]:


# In[4]:


# meta is use-less
orfs = orfs.drop(['datatype', 'specie'], axis=1)

# compute the codon composition
# this functionality can be added to the sklearn pipeline
orfs = (
    pd.concat([orfs, codon_composition(orfs.cds)], axis=1)
    .drop('cds', axis=1)
)


# merge predictors and response variable
data = pd.merge(decay, orfs, on="Gene_ID").drop(['std.error', 'Gene_ID'], axis=1)


# In[5]:



# In[6]:


# drop outliers
data = data[abs(data.decay_rate) < .5 ]


# In[7]:



# In[8]:


## Scale the response variable


# In[ ]:


def my_decay_scaler(grp):
    # scale the decay rate
    grp['decay_rate'] = (grp.decay_rate - grp.decay_rate.mean()) / grp.decay_rate.std()
    return grp
    
data = data.groupby(['datatype', 'cell_type', 'specie']).apply(my_decay_scaler)


# In[ ]:



data_train, data_test = train_test_split(data, test_size=0.3)

X_train = data_train.drop(['decay_rate'], axis=1)
y_train = data_train.decay_rate

y_test = data_test.decay_rate
X_test = data_test.drop(['decay_rate'], axis=1)


# In[ ]:


X_test.shape


# ## Data Pre-processing
# 
# ### Categorical Features
# 
# + One hot encoding
# 
# 
# ### Codon Features
# 
# + L2 normalization
# + Scaling

# In[ ]:


X_test.head()


# ## The pipeline

# In[ ]:


num_features = generate_codons()
cat_features = ['specie', 'datatype', 'cell_type']

num_pipeline = Pipeline([
    ('l2_norm', Normalizer(norm='l2')),
    ('scaler', StandardScaler())
])

cat_pipeline = Pipeline([
    ('one_hot', OneHotEncoder())
])

preprocessor_1 = ColumnTransformer(
    transformers=[
        ('num', num_pipeline, num_features),
        ('cat', cat_pipeline, cat_features)
    ]
)

full_pipeline = Pipeline(steps=[
    ('preproc1', preprocessor_1),
    ('interactions', PolynomialFeatures()),
    ('zerovar', VarianceThreshold(0)),
    ('scaler', StandardScaler()),
    ('pls', PLSRegression(n_components=10))
])


# ## Grid Search

# In[ ]:


param_grid = dict(
    pls__n_components=list(range(2, 21)),
    interactions__interaction_only = [True, False]
)

grid_search = GridSearchCV(
    full_pipeline,
    param_grid=param_grid,
    cv=10,
    verbose=10,
    n_jobs=32,
    scoring=make_scorer(r2_score)
)
grid_search.fit(X_train, y_train)


# In[ ]:


print(grid_search.best_estimator_)


# In[ ]:


y_pred = grid_search.best_estimator_.predict(X_test).reshape(1, -1)


# In[ ]:




# In[ ]:




# In[ ]:


preds = data_test.drop(generate_codons(), axis=1).assign(y_pred=y_pred[0])


# In[ ]:


preds.to_csv('preds.csv', index=False)


# In[ ]:


##
cv_res = pd.DataFrame(grid_search.cv_results_)


# In[ ]:


cv_res.to_csv("cv_res.csv", index=False)


# In[ ]:




