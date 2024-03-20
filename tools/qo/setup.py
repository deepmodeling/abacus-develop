import pathlib
import setuptools


here = pathlib.Path(__file__).parent.resolve()
readme = (here / 'README.md').read_text(encoding='utf-8')

# did not include torch and pyscf here
install_requires=["scipy", "numpy"]

setuptools.setup(
    name="abacus2qo",
    author="Yike HUANG",
    author_email="huangyk@aici.ac.cn",
    description="Quasiatomic Orbital (QO) analysis",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(include=['abacus2qo', 'abacus2qo.*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        "Programming Language :: Python :: 3.7",
    ],
    install_requires=install_requires,
    python_requires=">=3.7",
    version="0.0.2"
)