Developer Guide
===============

Welcome to the CRISP Developer Guide. This section provides key information for 
contributing to CRISP.

Development Setup
-----------------

1. **Clone the Repository**

   .. code-block:: bash

       git clone https://github.com/Indranil17/CRISP_HOST.git
       cd CRISP

2. **Set Up a Virtual Environment**

   .. code-block:: bash

       python -m venv venv
       source venv/bin/activate  # On Windows: `venv\Scripts\activate`

3. **Install Dependencies**

   .. code-block:: bash

       pip install -r requirements-dev.txt

4. **Install CRISP**

   .. code-block:: bash

       pip install -e .

Testing Procedures
------------------

1. **Run Tests**

   .. code-block:: bash

       pytest

2. **Check Coverage**

   .. code-block:: bash

       coverage run -m pytest
       coverage report

Contribution Guidelines
-----------------------

1. **Fork and Clone**

   .. code-block:: bash

       git clone https://github.com/your-username/CRISP.git
       cd CRISP

2. **Create a Branch**

   .. code-block:: bash

       git checkout -b feature/your-feature

3. **Implement Changes**

4. **Submit a Pull Request**

   .. code-block:: bash

       git push origin feature/your-feature


Developer Contacts
------------------

For further assistance or to discuss development-related topics, you can reach out to the following contacts:

- **Developer**: Indranil Saha [indranilsaha315@gmail.com] 
- **Developer**: Daniel Willimetz [daniel.willimetz@natur.cuni.cz]

- **Project Lead**: Lukas Grajciar [lukas.grajciar@natur.cuni.cz].

Acknowledgments
---------------

The CRISP package is developed by the Nano Materials Modelling Group at Charles University. We acknowledge the
Computational Molecular Science Python Cookiecutter version 1.1 for the skeletal structure of the package documentation.