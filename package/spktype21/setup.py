import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="spktype21",
    license="MIT",
    version="0.1.0",
    author="whiskie14142",
    author_email="whiskie14142@gmail.com",
    description="A supporting module for jplephem to handle data type 21",
    keywords="jplephem SPK type21",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/whiskie14142/spktype21",
    py_modules=["spktype21"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)