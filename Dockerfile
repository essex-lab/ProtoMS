# This Dockerfile can be used to create a production ready protoms environment
# suitable for distribution
# It is derived from the protoms-test-environment image defined in Dockerfile_test
# that sets up the necessary environment. This container simply copies across
# the repository and triggers compilation of the code. It should be run from
# a freshly cloned code repository.

# A container can be built from this file by running the below command in $PROTOMSHOME:
# docker build -t protoms:3.4 .

# An interactive container can then be started with
# docker run -it protoms:3.4 /bin/bash

From protoms-test-environment:3.4


COPY . /home/protoms/protoms-3.4
WORKDIR /home/protoms/protoms-3.4
ENV PROTOMSHOME /home/protoms/protoms-3.4/

RUN rm -rf build
RUN mkdir build
RUN cd build; \
cmake ..; \
make install
RUN chown -R protoms:protoms /home/protoms/protoms-3.4

USER protoms
