from setuptools import setup, find_packages

setup(
    name='CRISP',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'ase',
    ],
    entry_points={
        'console_scripts': [
            'CRISP=CRISP.cli:main',  # Main CRISP command
            'CRISP-test=CRISP.tests.runner:main',  # Test command
        ],
    },
)