from distutils.core import setup

setup(
    name="protomslib",
    description="Library providing the python components of the ProtoMS simulation package",
    version="3.4",
    author="Essex Research Group",
    author_email="jessexgroup@gmail.com",
    packages=["protomslib", "protomslib.prepare", "protomslib.free_energy"],
    url="http://www.protoms.org",
)
