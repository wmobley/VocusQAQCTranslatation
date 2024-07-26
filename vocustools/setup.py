from setuptools import setup, find_packages

setup(
    name="vocustools",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'h5py',
        # Add other dependencies here
    ],
)