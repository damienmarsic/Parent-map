import setuptools

with open('README.md', 'r',encoding='utf-8') as fh:
    long_description=fh.read()

setuptools.setup(
    name="parent-map",
    version="1.1.2",
    author="Damien Marsic",
    author_email="damien.marsic@aliyun.com",
    description="Analyze parental contributions to evolved or engineered protein or DNA sequences",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/damienmarsic/parent-map",
    package_dir={'': 'parent-map'},
    packages=setuptools.find_packages(),
    py_modules=["parent-map"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.6',
)
