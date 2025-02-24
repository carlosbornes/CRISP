Installation
===============

To ensure that CRISP functions optimally, you must have the Python Environment(≥3.11) with the
following software versions or higher installed on your system:

- ase: 3.23.0
- sklearn: 1.4.2
- seaborn: 0.12.2
- joblib: 1.2.0
- statsmodels: 0.14.0
- pandas: 2.0.3
- plotly: 5.9.0
- networkx: 3.1

You can install all the required packages using the following command:

.. code-block:: bash

    pip install ase>=3.23.0 scikit-learn>=1.4.2 seaborn>=0.12.2 joblib>=1.2.0 statsmodels>=0.14.0 pandas>=2.0.3 plotly>=5.9.0 networkx>=3.1


Next, clone the CRISP repository using the git clone command:

.. code-block:: bash

    git@github.com:Indranil17/CRISP.git


Now, navigate into the CRISP directory and install the package using pip:

.. code-block:: bash

    cd CRISP
    pip install .

Finally, verify that the installation was successful by launching a 
Python shell in your terminal and importing the CRISP package:

.. code-block:: bash
    
    python
    >>> import CRISP
    >>> CRISP.atom_indices()


If no errors occured and the CRISP package is successfully imported, 
your installation is complete and ready to use.









