# This Dockerfile can be used to create a production ready protoms environment
# suitable for distribution
# It is derived from the protoms-environment image defined in Dockerfile_test
# that sets up the necessary environment. This container simply copies across
# the repository and triggers compilation of the code. It should be run from
# a freshly cloned code repository.

# A container can be built from this file by running the below command in $PROTOMSHOME:
# docker build -t protoms .

# An interactive image can then be started with
# docker run -it protoms /bin/bash

From protoms-environment

WORKDIR /protoms-dev

ADD . /protoms-dev

RUN mkdir build
RUN cd build; \
cmake ..; \
make install
