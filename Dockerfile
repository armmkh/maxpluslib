# build with: 
# docker build . -t maxpluslib

FROM ubuntu:20.04

# the following lines are to avoid the interactive time zone question during install
ENV TZ=Europe/Amsterdam
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get install --force-yes -y libxml2-dev libxml2-utils python lcov build-essential make cmake pkg-config git doxygen graphviz

WORKDIR /usr/

# copy relevant input directories
# don't copy existing cmake cache files
COPY src ./src
COPY include ./include
COPY documentation ./documentation
COPY config ./config
COPY CMakeLists.txt .

# make release executables
RUN cmake -DCMAKE_BUILD_TYPE=Release .
RUN make
RUN mkdir /maxpluslib
RUN mkdir /maxpluslib/lib
RUN mkdir /maxpluslib/include
RUN cp -R ./build/lib/* /maxpluslib/lib
RUN cp -R ./include/maxplus /maxpluslib/include/

# execute tests
RUN cmake -DCODE_COVERAGE=ON .
RUN make
# run test and continue also if test fails
RUN mkdir /maxpluslibtest
RUN ctest > /maxpluslibtest/testlog.txt || :

#documentation
RUN mkdir /maxpluslibdoc
RUN make doc
RUN cp -R ./documentation/html /maxpluslibdoc/

# get coverage information
RUN make maxpluslibcoverage
RUN mkdir /maxpluslibcoverage
RUN cp -R ./build/coverage/* /maxpluslibcoverage/

VOLUME /maxpluslib
VOLUME /maxpluslibtest
VOLUME /maxpluslibdoc
VOLUME /maxpluslibcoverage
