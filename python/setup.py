from distutils.core import setup

setup(
    name="protomslib",
    description="Library providing the python components of the ProtoMS simulation package",
    version="3.4",
    author="Essex Research Group",
    author_email="jessexgroup@gmail.com",
    packages=["protomslib", "protomslib.prepare", "protomslib.free_energy"],
    url="https://www.protoms.org",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.6.3",
        "matplotlib>=3.4.2",
        "six>=1.16",
        "nose>=1.3.7",
    ]
)
