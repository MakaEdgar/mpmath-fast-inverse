import setuptools

setuptools.setup(
    name="mpinv",
    version="0.1",
    author="Edgar Makarov",
    author_email="e.makarov@skoltech.ru",
    description="Library to fast calculation mp.matrix invert, det and logdet",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    #url="https://gitlab.com/E.Makarov/multiprecision-matrix-inversion",
    packages=setuptools.find_packages(),
    package_data={'mpinv': ['mpinv*.so']},
    install_requires=[
        "mpmath",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
)
